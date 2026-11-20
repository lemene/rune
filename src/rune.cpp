// rune.cpp - Runs Uncover Nucleotide Errors (multithreaded C++ implementation)
//
// Subcommands (the first argument; numeric aliases 1-4 are still accepted):
//   rune gen-x   -i reads.fq -t kmers.tsv -k 19 [-o X.tsv] [options...]      generate features X
//   rune gen-y   -b aln.sam|aln.bam -x X.tsv -r ref.fna -k 19 [-o Y.tsv]     generate labels Y (SAM text or BAM; BAM needs samtools in PATH)
//   rune train   -x X.tsv -y Y.tsv -m model.txt [options...]                 train random forest
//   rune predict -m model.txt -x X.tsv [-r reads.fq | -b NB] [-o prefix]     predict from X file
//   rune predict -m model.txt -i reads.fq -t kmers.tsv -k 19 [-o prefix]     predict from raw reads
// All option parsing lives in src/options.hpp.
//
// Fixes over the Python prototype:
//   * canonical lookup (aligned with jellyfish -C tables; reverse-strand windows no longer missed)
//   * Y coordinate off-by-one (physical interval = [s, s+m+k-1), half-open)
//   * half-open boundary handling for CIGAR D ops
//   * min_len default 1 (keeps m=1 short runs; old behavior: --min-len 20)
//   * automatic threshold (median window count); strong_cov scales with coverage
//   * trees built in parallel across threads; model text compatible with scripts/export_rf.py
//
// Build: g++ -O3 -std=c++17 -pthread src/rune.cpp -o rune
// Deps: none (gz via gzip pipe, SAM/BAM via samtools view pipe)

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <thread>
#include <atomic>
#include <mutex>

#include "options.hpp"
#include "sequence.hpp"
#include "kmer.hpp"
#include "dataset.hpp"

using std::string;
using std::vector;
using std::unordered_map;

#include <fstream>

// ================= Step 3: train random forest =================

static const char* FEATURE_NAMES[5] = {
    "Length", "Read_Zero_Ratio", "Boundary_Cov_Ratio", "Region_Rel_Cov", "Threshold_Rel_Pos"
};

struct Sample { float x[5]; float y; uint32_t gid; };
struct RfNode { int8_t leaf; int32_t feat, left, right; float thr, val; };

struct RfModel {
    vector<RfNode> nodes;
    vector<int32_t> roots;
    float predict(const float* x) const {
        float acc = 0;
        for (int32_t r : roots) {
            int32_t i = r;
            while (!nodes[i].leaf) {
                const RfNode& nd = nodes[i];
                i = (x[nd.feat] <= nd.thr) ? nd.left : nd.right;
            }
            acc += nodes[i].val;
        }
        return acc / (float)roots.size();
    }
};

struct ModelWriter {
    vector<RfNode> nodes;
    int32_t emit(const RfNode& nd) { nodes.push_back(nd); return (int32_t)nodes.size() - 1; }
};

static int32_t build_tree_rec(const vector<Sample>& data, vector<uint32_t>& idx,
                              int depth, int max_depth, int min_leaf,
                              vector<float>* imp, ModelWriter& W)
{
    size_t n = idx.size();
    double sy = 0, syy = 0;
    for (uint32_t i : idx) { sy += data[i].y; syy += (double)data[i].y * data[i].y; }
    double mean = sy / n;
    double sse = syy - sy * sy / n;

    bool is_leaf = (depth >= max_depth) || (n < (size_t)(2 * min_leaf)) || (sse <= 1e-12);
    int best_f = -1; float best_thr = 0; double best_gain = 0;
    if (!is_leaf) {
        thread_local vector<std::pair<float, uint32_t>> fv;
        thread_local vector<double> psum, psumsq;
        for (int f = 0; f < 5; f++) {
            fv.clear(); fv.reserve(n);
            for (uint32_t i : idx) fv.push_back({data[i].x[f], i});
            std::sort(fv.begin(), fv.end(),
                      [](const std::pair<float,uint32_t>& a, const std::pair<float,uint32_t>& b)
                      { return a.first < b.first; });
            size_t nn = fv.size();
            psum.assign(nn + 1, 0); psumsq.assign(nn + 1, 0);
            for (size_t i = 0; i < nn; i++) {
                float y = data[fv[i].second].y;
                psum[i + 1] = psum[i] + y;
                psumsq[i + 1] = psumsq[i] + (double)y * y;
            }
            for (size_t i = (size_t)min_leaf; i + (size_t)min_leaf <= nn; i++) {
                if (fv[i - 1].first == fv[i].first) continue;
                double sl = psum[i], sl2 = psumsq[i];
                double sr = psum[nn] - sl, sr2 = psumsq[nn] - sl2;
                double sseL = sl2 - sl * sl / (double)i;
                double sseR = sr2 - sr * sr / (double)(nn - i);
                double gain = sse - sseL - sseR;
                if (gain > best_gain) {
                    best_gain = gain;
                    best_f = f;
                    best_thr = 0.5f * (fv[i - 1].first + fv[i].first);
                }
            }
        }
        if (best_f < 0 || best_gain <= 1e-12) is_leaf = true;
    }
    if (is_leaf) {
        RfNode nd{}; nd.leaf = 1; nd.val = (float)mean;
        return W.emit(nd);
    }
    if (imp) (*imp)[best_f] += (float)(best_gain * (double)n);
    vector<uint32_t> li, ri;
    li.reserve(n); ri.reserve(n);
    for (uint32_t i : idx)
        (data[i].x[best_f] <= best_thr ? li : ri).push_back(i);
    if (li.empty() || ri.empty()) {
        RfNode nd{}; nd.leaf = 1; nd.val = (float)mean;
        return W.emit(nd);
    }
    RfNode nd{}; nd.leaf = 0; nd.feat = best_f; nd.thr = best_thr;
    int32_t my = W.emit(nd);
    int32_t Lc = build_tree_rec(data, li, depth + 1, max_depth, min_leaf, imp, W);
    int32_t Rc = build_tree_rec(data, ri, depth + 1, max_depth, min_leaf, imp, W);
    W.nodes[my].left = Lc;
    W.nodes[my].right = Rc;
    return my;
}

static void run_step3(const string& x_path, const string& y_path, const string& model_out,
                      const TrainOpts& o, int nthreads)
{
    // read X (keep Start for joining with Y)
    struct RawX { uint32_t start; uint32_t gid; float x[5]; };
    vector<RawX> raw;
    unordered_map<string, uint32_t> read_ids;
    {
        FILE* f = open_text(x_path);
        string line;
        bool first = true;
        while (read_line(f, line)) {
            if (first) { first = false; continue; }
            size_t t1 = line.find('\t'); if (t1 == string::npos) continue;
            size_t t2 = line.find('\t', t1 + 1); if (t2 == string::npos) continue;
            string rn = line.substr(0, t1);
            uint32_t gid;
            auto it = read_ids.find(rn);
            if (it == read_ids.end()) { gid = (uint32_t)read_ids.size(); read_ids[rn] = gid; }
            else gid = it->second;
            RawX s{};
            s.start = (uint32_t)atoll(line.c_str() + t1 + 1);
            s.gid = gid;
            size_t q = t2 + 1;
            for (int c = 0; c < 5; c++) {
                s.x[c] = strtof(line.c_str() + q, nullptr);
                size_t nt = line.find('\t', q);
                if (nt == string::npos) break;
                q = nt + 1;
            }
            raw.push_back(s);
        }
        close_text(f, x_path);
    }
    fprintf(stderr, "[rune] X: %zu regions, %zu reads\n", raw.size(), read_ids.size());

    // read Y, join on (gid, start)
    {
        unordered_map<uint64_t, uint32_t> ymap;
        size_t skipped_cov = 0;
        FILE* f = open_text(y_path);
        string line;
        bool first = true;
        while (read_line(f, line)) {
            if (first) { first = false; continue; }
            size_t t1 = line.find('\t'); if (t1 == string::npos) continue;
            string rn = line.substr(0, t1);
            auto it = read_ids.find(rn);
            if (it == read_ids.end()) continue;
            size_t q = t1 + 1;
            long long col[5] = {0,0,0,0,0};
            for (int c = 0; c < 5; c++) {
                col[c] = atoll(line.c_str() + q);
                size_t nt = line.find('\t', q);
                if (nt == string::npos) break;
                q = nt + 1;
            }
            // columns: Region_Start, Region_End, Mismatch, Insertion, Deletion [, Covered]
            // optional 7th column Covered: 1=fully inside alignment, 0=partially covered
            // (Y underestimated); Covered=0 rows are excluded from training
            {
                int tabs = 0; size_t p6 = string::npos;
                for (size_t p = line.find('\t'); p != string::npos; p = line.find('\t', p + 1)) {
                    if (++tabs == 6) { p6 = p + 1; break; }
                }
                if (p6 != string::npos && atoll(line.c_str() + p6) == 0) { skipped_cov++; continue; }
            }
            uint64_t key = ((uint64_t)it->second << 32) | (uint64_t)(uint32_t)col[0];
            ymap[key] = (uint32_t)(col[2] + col[3] + col[4]);
        }
        close_text(f, y_path);
        if (skipped_cov) fprintf(stderr, "[rune] dropped partially covered samples (Covered=0): %zu\n", skipped_cov);
        vector<Sample> data; data.reserve(raw.size());
        size_t matched = 0;
        for (auto& s : raw) {
            uint64_t key = ((uint64_t)s.gid << 32) | s.start;
            auto it = ymap.find(key);
            if (it == ymap.end()) continue;
            Sample d{};
            for (int c = 0; c < 5; c++) d.x[c] = s.x[c];
            d.gid = s.gid;
            d.y = log1p((float)it->second);
            data.push_back(d);
            matched++;
        }
        raw.clear();
        fprintf(stderr, "[rune] joined samples: %zu\n", matched);
        if (data.size() < 100) { fprintf(stderr, "[rune] too few samples, refusing to train\n"); exit(3); }

        // split grouped by read (no leakage)
        vector<uint32_t> gids; gids.reserve(read_ids.size());
        for (auto& kv : read_ids) gids.push_back(kv.second);
        std::mt19937 rng((unsigned)o.seed);
        std::shuffle(gids.begin(), gids.end(), rng);
        size_t ntest_g = (size_t)(gids.size() * o.test);
        std::unordered_set<uint32_t> test_set(gids.begin(), gids.begin() + ntest_g);
        vector<uint32_t> tri, tei;
        for (size_t i = 0; i < data.size(); i++)
            (test_set.count(data[i].gid) ? tei : tri).push_back((uint32_t)i);
        fprintf(stderr, "[rune] train %zu / test %zu (grouped by read)\n", tri.size(), tei.size());

        // tree-level parallel training
        vector<RfNode> nodes; nodes.reserve((size_t)o.trees * 511);
        vector<int32_t> roots;
        vector<float> imp(5, 0.f);
        std::mutex wmtx;
        std::atomic<int> next_tree{0};
        auto tree_worker = [&]() {
            vector<float> limp(5, 0.f);
            while (true) {
                int t = next_tree++;
                if (t >= o.trees) break;
                std::mt19937 lrng((unsigned)o.seed * 100003u + (unsigned)t + 7u);
                vector<uint32_t> bs(tri.size());
                for (auto& b : bs) b = tri[lrng() % tri.size()];
                ModelWriter LW;
                build_tree_rec(data, bs, 0, o.depth, o.leaf, &limp, LW);
                std::lock_guard<std::mutex> lk(wmtx);
                int32_t base = (int32_t)nodes.size();
                for (auto& nd : LW.nodes) {
                    RfNode nd2 = nd;
                    if (!nd2.leaf) { nd2.left += base; nd2.right += base; }
                    nodes.push_back(nd2);
                }
                roots.push_back(base);
                for (int f = 0; f < 5; f++) imp[f] += limp[f];
            }
        };
        {
            int nt = std::max(1, std::min(nthreads, o.trees));
            vector<std::thread> ths;
            for (int t = 1; t < nt; t++) ths.emplace_back(tree_worker);
            tree_worker();
            for (auto& th : ths) th.join();
        }

        // write model (export_rf.py compatible: line 1 = tree count; internal node
        // "f thr L R"; leaf "-1 0.0 0 0 val")
        {
            FILE* fo = fopen(model_out.c_str(), "w");
            if (!fo) { fprintf(stderr, "[rune] cannot write %s\n", model_out.c_str()); exit(2); }
            fprintf(fo, "%d\n", (int)roots.size());
            for (auto& nd : nodes) {
                if (nd.leaf) fprintf(fo, "-1 0.0 0 0 %.6f\n", nd.val);
                else fprintf(fo, "%d %.6f %d %d\n", nd.feat, nd.thr, nd.left, nd.right);
            }
            fclose(fo);
            FILE* ff = fopen((model_out + ".features").c_str(), "w");
            fprintf(ff, "5\n");
            for (int i = 0; i < 5; i++) fprintf(ff, "%s\n", FEATURE_NAMES[i]);
            fclose(ff);
            fprintf(stderr, "[rune] model saved: %s (+ %s.features)\n", model_out.c_str(), model_out.c_str());
        }
        double s = 0; for (double v : imp) s += v;
        fprintf(stderr, "[rune] feature importance: ");
        for (int i = 0; i < 5; i++) fprintf(stderr, "%s=%.3f ", FEATURE_NAMES[i], s > 0 ? imp[i] / s : 0.0);
        fprintf(stderr, "\n");

        // test evaluation (original scale)
        if (!tei.empty()) {
            RfModel M; M.nodes = std::move(nodes); M.roots = std::move(roots);
            double sse = 0, ss_tot = 0, mean_y = 0, sae = 0;
            size_t n = tei.size();
            for (uint32_t i : tei) mean_y += std::expm1(data[i].y);
            mean_y /= n;
            for (uint32_t i : tei) {
                float pred = std::expm1(M.predict(data[i].x)); if (pred < 0) pred = 0;
                float truth = std::expm1(data[i].y);
                sae += std::fabs(pred - truth);
                sse += (pred - truth) * (pred - truth);
                ss_tot += (truth - mean_y) * (truth - mean_y);
            }
            fprintf(stderr, "[rune] test: MAE=%.4f RMSE=%.4f R2=%.4f\n",
                    sae / n, std::sqrt(sse / n), ss_tot > 0 ? 1 - sse / ss_tot : 0.0);
        }
    }
}

// ================= Step 4: predict =================

struct XFile {
    vector<string> read;
    vector<uint32_t> start;
    vector<std::array<float, 5>> feat;
};

static XFile load_xfile(const string& path) {
    XFile X;
    FILE* f = open_text(path);
    string line;
    if (!read_line(f, line)) { fprintf(stderr, "[rune] X file is empty\n"); exit(2); }
    vector<string> cols;
    { size_t pp = 0;
      while (true) { size_t nt = line.find('\t', pp);
          cols.push_back(line.substr(pp, nt == string::npos ? string::npos : nt - pp));
          if (nt == string::npos) break;
          pp = nt + 1;
      } }
    int col_read = -1, col_start = -1, colf[5] = {-1,-1,-1,-1,-1};
    for (size_t i = 0; i < cols.size(); i++) {
        if (cols[i] == "Read_Name") col_read = (int)i;
        else if (cols[i] == "Start") col_start = (int)i;
        for (int c = 0; c < 5; c++) if (cols[i] == FEATURE_NAMES[c]) colf[c] = (int)i;
    }
    if (col_read < 0 || col_start < 0) { fprintf(stderr, "[rune] X is missing Read_Name/Start columns\n"); exit(2); }
    for (int c = 0; c < 5; c++) if (colf[c] < 0) { fprintf(stderr, "[rune] X is missing feature column %s\n", FEATURE_NAMES[c]); exit(2); }
    while (read_line(f, line)) {
        if (line.empty()) continue;
        vector<size_t> pos, nend;
        size_t q = 0;
        while (true) {
            size_t nt = line.find('\t', q);
            pos.push_back(q);
            nend.push_back(nt);
            if (nt == string::npos) break;
            q = nt + 1;
        }
        auto cptr = [&](int c) -> const char* {
            if (c < 0 || c >= (int)pos.size()) return "0";
            return line.c_str() + pos[c];
        };
        X.read.push_back(line.substr(pos[col_read],
                     nend[col_read] == string::npos ? string::npos : nend[col_read] - pos[col_read]));
        X.start.push_back((uint32_t)atoll(cptr(col_start)));
        std::array<float, 5> fr{};
        for (int c = 0; c < 5; c++) fr[c] = strtof(cptr(colf[c]), nullptr);
        X.feat.push_back(fr);
    }
    close_text(f, path);
    return X;
}

static RfModel load_model(const string& path) {
    RfModel M;
    FILE* f = fopen(path.c_str(), "r");
    if (!f) { fprintf(stderr, "[rune] cannot open model %s\n", path.c_str()); exit(2); }
    int n_trees = 0;
    if (fscanf(f, "%d", &n_trees) != 1 || n_trees <= 0) {
        fprintf(stderr, "[rune] bad model format (first line must be the tree count)\n"); exit(2);
    }
    char linebuf[256];
    for (int t = 0; t < n_trees; t++) {
        M.roots.push_back((int32_t)M.nodes.size());
        vector<int32_t> stack;
        while (true) {
            if (!fgets(linebuf, sizeof linebuf, f)) break;
            int feat = 0, L = 0, R = 0; float thr = 0, val = 0;
            int cnt = sscanf(linebuf, "%d %f %d %d %f", &feat, &thr, &L, &R, &val);
            if (cnt < 2) continue;
            RfNode nd{};
            nd.leaf = (feat == -1);
            if (nd.leaf) nd.val = (cnt >= 5 ? val : thr);
            else { nd.feat = feat; nd.thr = thr; nd.left = L; nd.right = R; }
            M.nodes.push_back(nd);
            if (!nd.leaf) { stack.push_back(nd.right); stack.push_back(nd.left); }
            if (stack.empty()) break;   // end of this tree
            stack.pop_back();           // consume one preorder slot
        }
        // sanity check: at least a root node was read for this tree
        if ((int32_t)M.nodes.size() <= M.roots.back()) { fprintf(stderr, "[rune] incomplete model file\n"); exit(2); }
    }
    fclose(f);
    fprintf(stderr, "[rune] model: %d trees, %zu nodes\n", n_trees, M.nodes.size());
    return M;
}

static unordered_map<string, long long> read_lengths(const string& reads_path, int nthreads) {
    unordered_map<string, long long> len;
    std::mutex mtx;
    auto w = [&](int tid) {
        SeqReader rd(reads_path);
        string name, seq;
        unordered_map<string, long long> local;
        size_t idx = 0;
        while (rd.next(name, seq)) {
            if (idx % (size_t)nthreads == (size_t)tid) local[name] = (long long)seq.size();
            idx++;
        }
        std::lock_guard<std::mutex> lk(mtx);
        len.insert(local.begin(), local.end());
    };
    vector<std::thread> ths;
    for (int t = 1; t < nthreads; t++) ths.emplace_back(w, t);
    w(0);
    for (auto& th : ths) th.join();
    return len;
}

static void run_step4_from_x(const RfModel& M, const XFile& X, long long total_bases,
                             const unordered_map<string, long long>* read_len,
                             int nthreads, const string& out_prefix)
{
    // parallel prediction
    vector<float> pred_log(X.feat.size());
    std::atomic<size_t> cursor{0};
    auto w = [&]() {
        while (true) {
            size_t i = cursor.fetch_add(8192);
            if (i >= X.feat.size()) break;
            size_t hi = std::min(X.feat.size(), i + 8192);
            for (size_t j = i; j < hi; j++) pred_log[j] = M.predict(X.feat[j].data());
        }
    };
    vector<std::thread> ths;
    for (int t = 1; t < nthreads; t++) ths.emplace_back(w);
    w();
    for (auto& th : ths) th.join();

    vector<float> pred(X.feat.size());
    double tot = 0;
    for (size_t i = 0; i < pred.size(); i++) {
        pred[i] = (float)std::expm1(pred_log[i]);
        if (pred[i] < 0) pred[i] = 0;
        tot += pred[i];
    }
    string p1 = out_prefix + "predicted_regions.tsv";
    {
        FILE* f = fopen(p1.c_str(), "w");
        if (!f) { fprintf(stderr, "[rune] cannot write %s\n", p1.c_str()); exit(2); }
        fprintf(f, "Read_Name\tStart\tLength\tPredicted_Errors\n");
        for (size_t i = 0; i < pred.size(); i++)
            fprintf(f, "%s\t%u\t%u\t%.4f\n", X.read[i].c_str(), X.start[i],
                    (unsigned)X.feat[i][0], pred[i]);
        fclose(f);
    }
    string p2 = out_prefix + "predicted_per_read.tsv";
    size_t nreads = 0;
    {
        unordered_map<string, double> acc; acc.reserve(X.read.size() * 2);
        for (size_t i = 0; i < pred.size(); i++) acc[X.read[i]] += pred[i];
        FILE* f = fopen(p2.c_str(), "w");
        if (!f) { fprintf(stderr, "[rune] cannot write %s\n", p2.c_str()); exit(2); }
        fprintf(f, "Read_Name\tPredicted_Errors\tRead_Bases\tRead_Error_Rate\n");
        vector<const string*> names; names.reserve(acc.size());
        for (auto& kv : acc) names.push_back(&kv.first);
        std::sort(names.begin(), names.end(),
                  [](const string* a, const string* b) { return *a < *b; });
        for (auto* np : names) {
            long long b = read_len ? (read_len->count(*np) ? read_len->at(*np) : 0) : 0;
            double rate = b > 0 ? acc[*np] / b : 0.0;
            fprintf(f, "%s\t%.3f\t%lld\t%.6f\n", np->c_str(), acc[*np], b, rate);
        }
        fclose(f);
        nreads = acc.size();
        fprintf(stderr, "[rune] per-read results: %s (%zu reads)\n", p2.c_str(), nreads);
    }
    printf("\n================================================\n");
    printf("          rune dataset-level error rate report\n");
    printf("================================================\n");
    printf("Predicted regions:   %zu\n", pred.size());
    printf("Predicted errors:    %.0f\n", tot);
    if (total_bases > 0) {
        double rate = tot / (double)total_bases;
        printf("Total bases:         %lld\n", total_bases);
        printf("Overall error rate:  %.6f  (accuracy %.4f%%)\n", rate, (1 - rate) * 100);
        if (rate > 0) printf("Phred Q:         %.2f\n", -10.0 * std::log10(rate));
    } else {
        printf("(total base count not provided, skipping overall error rate)\n");
    }
    printf("================================================\n");
}

static void run_step4(const Options& opts, int nthreads) {
    const string& model_path = opts.model;
    const string& x_path = opts.xFile;
    const string& reads_path = opts.fromReads ? opts.reads : opts.lenReads;
    const string& kmer_path = opts.kmerFile;
    const string& out_prefix = opts.output;
    long long total_bases = opts.totalBases;
    bool from_reads = opts.fromReads;
    int k = opts.xp.k;

    RfModel M = load_model(model_path);
    XFile X;
    if (from_reads) {
        if (kmer_path.empty() || reads_path.empty())
            { fprintf(stderr, "[rune] data mode requires -i reads -t kmers.tsv -k K\n"); exit(1); }
        KmerCodec codec(k);
        KmerTable table; table.load(kmer_path, codec);
        vector<XRow> rows;
        FeatureBuilder(opts).generate(reads_path, table, codec, &rows, "");
        X.read.reserve(rows.size()); X.start.reserve(rows.size()); X.feat.reserve(rows.size());
        for (auto& r : rows) {
            X.read.push_back(std::move(r.read));
            X.start.push_back(r.start);
            X.feat.push_back({ (float)r.m, r.zero_ratio, r.boundary_ratio, r.region_rel, r.thr_rel });
        }
    } else {
        if (x_path.empty()) { fprintf(stderr, "[rune] predict requires -x X.tsv or -i data mode\n"); exit(1); }
        X = load_xfile(x_path);
    }
    unordered_map<string, long long> rlen;
    const unordered_map<string, long long>* rlen_ptr = nullptr;
    if (!reads_path.empty()) {
        fprintf(stderr, "[rune] counting read lengths...\n");
        rlen = read_lengths(reads_path, nthreads);
        rlen_ptr = &rlen;
        if (total_bases < 0) {
            long long s = 0;
            for (auto& kv : rlen) s += kv.second;
            total_bases = s;
        }
    }
    run_step4_from_x(M, X, total_bases, rlen_ptr, nthreads, out_prefix);
}

// ================= main =================

int main(int argc, char** argv) {
    Options opts;
    if (!opts.parse(argc, argv)) return 1;
    if (opts.helpRequested) return 0;
    int nthreads = opts.threads;

    if (opts.command == Options::Command::GenX) {
        KmerCodec codec(opts.xp.k);
        KmerTable tb; tb.load(opts.kmerFile, codec);
        FeatureBuilder(opts).generate(opts.reads, tb, codec, nullptr, opts.output);
        return 0;
    }
    if (opts.command == Options::Command::GenY) {
        if (!opts.output.empty())
            if (!freopen(opts.output.c_str(), "w", stdout))
                { fprintf(stderr, "[rune] cannot write %s\n", opts.output.c_str()); return 2; }
        LabelBuilder(opts).generate(opts.sam, opts.xFile, opts.refFasta);
        return 0;
    }
    if (opts.command == Options::Command::Train) {
        run_step3(opts.xFile, opts.yFile, opts.model, opts.train, nthreads);
        return 0;
    }
    run_step4(opts, nthreads);
    return 0;
}
