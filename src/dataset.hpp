#pragma once

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "kmer.hpp"
#include "options.hpp"
#include "sequence.hpp"

static const char* X_HEADER =
    "Read_Name\tStart\tLength\tRead_Zero_Ratio\tBoundary_Cov_Ratio\tRegion_Rel_Cov\tThreshold_Rel_Pos";

struct XRow {
    std::string read;
    uint32_t start = 0, m = 0;
    float zero_ratio = 0, boundary_ratio = 0, region_rel = 0, thr_rel = 0;
};

// Builds the feature matrix X from reads and a k-mer count table.
class FeatureBuilder {
public:
    explicit FeatureBuilder(const Options& opts) : opts_(opts) {}

    // Fills `samples` in memory (predict step), otherwise writes to out_path or stdout.
    void generate(const std::string& reads_path, const KmerTable& table, const KmerCodec& codec,
                  std::vector<XRow>* samples, const std::string& out_path) const
    {
        const XParams& p = opts_.xp;
        int nthreads = opts_.threads;

        long long cov;
        if (p.global_cov > 0) cov = p.global_cov;
        else {
            fprintf(stderr, "[rune] gen-x: estimating coverage...\n");
            cov = estimate_median_cov(reads_path, table, codec, p.k);
            fprintf(stderr, "[rune] estimated coverage d = %lld\n", cov);
            if (cov < 5) fprintf(stderr, "[rune] warning: coverage too low, threshold unreliable\n");
        }
        int T = std::max(2, (int)std::floor(cov * p.t_ratio));
        int strong_cov = p.strong_cov > 0 ? (int)p.strong_cov
                                          : std::max(10, (int)std::floor(0.3 * std::max(cov, 10LL)));
        fprintf(stderr, "[rune] T=%d (t_ratio=%.2f) strong_cov=%d max_gap=%d min_len=%lld threads=%d\n",
                T, p.t_ratio, strong_cov, p.max_gap, p.min_len, nthreads);

        std::vector<XRow> all;
        std::mutex out_mtx;
        std::atomic<size_t> nreads{0};
        auto worker = [&](int tid) {
            SeqReader rd(reads_path);
            std::string name, seq;
            std::vector<uint8_t> bin; std::vector<uint32_t> cnts;
            std::vector<XRow> local; local.reserve(2048);
            size_t idx = 0;
            while (rd.next(name, seq)) {
                bool mine = (idx % (size_t)nthreads == (size_t)tid);
                idx++;
                if (!mine) continue;
                nreads++;
                process_read(name, seq, table, codec, p, T, strong_cov, bin, cnts, local);
                if (local.size() >= 4096) {
                    std::lock_guard<std::mutex> lk(out_mtx);
                    all.insert(all.end(), local.begin(), local.end());
                    local.clear();
                }
            }
            if (!local.empty()) {
                std::lock_guard<std::mutex> lk(out_mtx);
                all.insert(all.end(), local.begin(), local.end());
            }
        };
        std::vector<std::thread> ths;
        for (int t = 1; t < nthreads; t++) ths.emplace_back(worker, t);
        worker(0);
        for (auto& th : ths) th.join();
        fprintf(stderr, "[rune] processed %zu reads, extracted %zu regions\n", nreads.load(), all.size());

        std::sort(all.begin(), all.end(), [](const XRow& a, const XRow& b) {
            if (a.read != b.read) return a.read < b.read;
            return a.start < b.start;
        });

        if (samples) { *samples = std::move(all); return; }
        if (!out_path.empty()) {
            FILE* fo = fopen(out_path.c_str(), "w");
            if (!fo) { fprintf(stderr, "[rune] cannot write %s\n", out_path.c_str()); exit(2); }
            fprintf(fo, "%s\n", X_HEADER);
            std::string buf; buf.reserve(1 << 22);
            char tmp[256];
            for (auto& r : all) {
                buf += r.read; buf += '\t';
                snprintf(tmp, sizeof tmp, "%u\t%u\t%.4f\t%.4f\t%.4f\t%.4f\n",
                         r.start, r.m, r.zero_ratio, r.boundary_ratio, r.region_rel, r.thr_rel);
                buf += tmp;
                if (buf.size() > (1 << 22)) { fputs(buf.c_str(), fo); buf.clear(); }
            }
            fputs(buf.c_str(), fo);
            fclose(fo);
            fprintf(stderr, "[rune] gen-x done -> %s\n", out_path.c_str());
        } else {
            fprintf(stdout, "%s\n", X_HEADER);
            char tmp[256];
            for (auto& r : all) {
                fputs(r.read.c_str(), stdout); fputc('\t', stdout);
                snprintf(tmp, sizeof tmp, "%u\t%u\t%.4f\t%.4f\t%.4f\t%.4f\n",
                         r.start, r.m, r.zero_ratio, r.boundary_ratio, r.region_rel, r.thr_rel);
                fputs(tmp, stdout);
            }
        }
    }

private:
    static void process_read(const std::string& name, const std::string& seq,
                             const KmerTable& table, const KmerCodec& codec,
                             const XParams& p, int T, int strong_cov,
                             std::vector<uint8_t>& bin, std::vector<uint32_t>& cnts,
                             std::vector<XRow>& out)
    {
        size_t L = seq.size();
        if ((int)L < p.k) return;
        size_t nw = L - p.k + 1;
        bin.resize(nw); cnts.resize(nw);

        double cov_sum = 0; size_t cov_n = 0, n_zero = 0;
        for (size_t i = 0; i < nw; i++) {
            uint64_t code;
            bool ok = codec.encode(seq.c_str(), i, code);
            uint32_t c = ok ? table.get(code) : 0u;
            cnts[i] = c;
            if ((int)c >= T) { bin[i] = 1; cov_sum += c; cov_n++; }
            else { bin[i] = 0; n_zero++; }
        }
        double read_avg_cov = cov_n ? cov_sum / cov_n : 1e-5;
        double zero_ratio = (double)n_zero / nw;

        std::vector<std::pair<size_t, size_t>> raw;
        for (size_t w = 0; w < nw; w++) {
            if (bin[w] == 0) {
                if (!raw.empty() && raw.back().second == w - 1) raw.back().second = w;
                else raw.push_back({w, w});
            }
        }

        std::vector<std::pair<size_t, size_t>> merged;
        if (!raw.empty()) {
            size_t cs = raw[0].first, ce = raw[0].second;
            for (size_t j = 1; j < raw.size(); j++) {
                size_t ns = raw[j].first, ne = raw[j].second;
                size_t gap = ns - ce - 1;
                uint32_t gmax = 0;
                for (size_t w = ce + 1; w < ns; w++) gmax = std::max(gmax, cnts[w]);
                if ((long long)gap <= p.max_gap && gmax < (uint32_t)strong_cov) ce = ne;
                else { merged.push_back({cs, ce}); cs = ns; ce = ne; }
            }
            merged.push_back({cs, ce});
        }

        for (auto& rg : merged) {
            size_t s = rg.first, e = rg.second, m = e - s + 1;
            long long phys = (long long)m + p.k - 1;
            if (phys < p.min_len) continue;

            double region_abs = 0;
            for (size_t w = s; w <= e; w++) region_abs += cnts[w];
            region_abs /= (double)m;

            double bsum = 0; size_t bn = 0;
            for (size_t w = (s >= 5 ? s - 5 : 0); w < s; w++) { bsum += cnts[w]; bn++; }
            for (size_t w = e + 1; w <= e + 5 && w < nw; w++) { bsum += cnts[w]; bn++; }
            double boundary_abs = bn ? bsum / bn : 0.0;

            XRow r;
            r.read = name;
            r.start = (uint32_t)s;
            r.m = (uint32_t)m;
            r.zero_ratio = (float)zero_ratio;
            r.boundary_ratio = (float)(boundary_abs / read_avg_cov);
            r.region_rel = (float)(region_abs / read_avg_cov);
            r.thr_rel = (float)((region_abs - T) / T);
            out.push_back(std::move(r));
        }
    }

    const Options& opts_;
};

// Builds labels Y from alignment, the X regions and the reference genome.
class LabelBuilder {
public:
    explicit LabelBuilder(const Options& opts) : opts_(opts) {}

    void generate(const std::string& sam_path, const std::string& x_path,
                  const std::string& fa_path) const
    {
        int k = opts_.xp.k;
        int nthreads = opts_.threads;

        std::unordered_map<std::string, std::vector<Region2>> regions;
        {
            FILE* f = open_text(x_path);
            std::string line;
            bool first = true;
            while (read_line(f, line)) {
                if (first) { first = false; if (!line.empty() && !isdigit(line[0])) continue; }
                size_t t1 = line.find('\t'), t2 = line.find('\t', t1 + 1);
                if (t1 == std::string::npos || t2 == std::string::npos) continue;
                std::string rn = line.substr(0, t1);
                long long s = atoll(line.c_str() + t1 + 1);
                long long m = atoll(line.c_str() + t2 + 1);
                regions[rn].push_back({(uint32_t)s, (uint32_t)(s + m + k - 1)});
            }
            close_text(f, x_path);
            fprintf(stderr, "[rune] X regions: %zu reads\n", regions.size());
        }
        auto fai = load_fai(fa_path);

        fprintf(stdout, "Read_Name\tRegion_Start\tRegion_End\tMismatch\tInsertion\tDeletion\tCovered\n");
        FILE* f = open_text(sam_path);
        // Producer-consumer: the main thread reads SAM lines into a bounded queue;
        // workers pop and process. A tid-stride filter on a shared reader is WRONG
        // (lines are read by whichever thread wins the race, so most get dropped).
        std::deque<std::string> q;
        std::mutex q_mtx;
        std::condition_variable q_cv;
        bool producer_done = false;
        GenYStats st;
        std::mutex out_mtx;
        std::string out_buf; out_buf.reserve(1 << 22);

        std::vector<std::thread> ths;
        for (int t = 0; t < nthreads; t++)
            ths.emplace_back([&]() {
                sam_worker(fa_path, regions, fai, st, q_mtx, q_cv, q, producer_done, out_mtx, out_buf);
            });
        std::string line;
        while (read_line(f, line)) {
            {
                std::unique_lock<std::mutex> lk(q_mtx);
                q_cv.wait(lk, [&] { return q.size() < 4096; });
                q.push_back(std::move(line));
            }
            q_cv.notify_one();
        }
        close_text(f, sam_path);
        {
            std::lock_guard<std::mutex> lk(q_mtx);
            producer_done = true;
        }
        q_cv.notify_all();
        for (auto& th : ths) th.join();
        fputs(out_buf.c_str(), stdout);
        fprintf(stderr,
                "[rune] gen-y diagnostics: SAM lines=%lld (headers=%lld) records=%lld\n"
                "[rune]   skipped: unmapped=%lld secondary/supp=%lld no-such-read-in-X=%lld no-ref-entry=%lld malformed=%lld\n"
                "[rune]   regions: emitted=%lld skipped(qs>=qe)=%lld uncovered-by-cigar=%lld\n",
                st.lines.load(), st.headers.load(), st.recs.load(),
                st.unmapped.load(), st.secondary.load(), st.no_region.load(),
                st.no_ref.load(), st.malformed.load(),
                st.regions_ok.load(), st.regions_skipped.load(),
                st.regions_uncovered.load());
        fprintf(stderr, "[rune] gen-y done\n");
    }

private:
    struct FaidxEntry { std::string name; long long len, offset, linebases, linewidth; };
    struct Region2 { uint32_t s, e; };   // e is exclusive

    static std::unordered_map<std::string, FaidxEntry> load_fai(const std::string& fa) {
        std::unordered_map<std::string, FaidxEntry> idx;
        std::ifstream f(fa + ".fai");
        if (!f) { fprintf(stderr, "[rune] missing %s.fai (run samtools faidx %s first)\n", fa.c_str(), fa.c_str()); exit(2); }
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty()) continue;
            size_t p1 = line.find('\t');
            if (p1 == std::string::npos) continue;
            size_t p2 = line.find('\t', p1 + 1);
            size_t p3 = line.find('\t', p2 + 1);
            size_t p4 = line.find('\t', p3 + 1);
            FaidxEntry e;
            e.name = line.substr(0, p1);
            e.len = atoll(line.c_str() + p1 + 1);
            e.offset = atoll(line.c_str() + p2 + 1);
            e.linebases = atoll(line.c_str() + p3 + 1);
            e.linewidth = (p4 == std::string::npos) ? e.linebases + 1 : atoll(line.c_str() + p4 + 1);
            idx[e.name] = e;
        }
        return idx;
    }

    static std::string fa_fetch(FILE* f, const FaidxEntry& e, long long start, long long end) {
        if (start < 0) start = 0;
        if (end > e.len) end = e.len;
        if (start >= end) return "";
        std::string out; out.reserve(end - start);
        fseek64(f, e.offset + (start / e.linebases) * e.linewidth + (start % e.linebases), SEEK_SET);
        long long need = end - start;
        std::vector<char> buf(1 << 16);
        while (need > 0) {
            size_t want = (size_t)std::min<long long>(need, (long long)buf.size());
            size_t got = fread(buf.data(), 1, want, f);
            if (got == 0) break;
            for (size_t i = 0; i < got; i++) {
                char c = buf[i];
                if (c == '\n' || c == '\r') continue;
                out += (c >= 'a' && c <= 'z') ? char(c - 32) : c;
                if (--need == 0) break;
            }
        }
        return out;
    }

    static void count_errors(const std::string& seq, const std::string& refsub, const std::string& cigar,
                             long long q_start, long long q_end,
                             long long& mm, long long& ins, long long& del)
    {
        long long q_pos = 0, r_pos = 0;
        size_t ci = 0, C = cigar.size();
        while (ci < C) {
            long long len = 0;
            while (ci < C && cigar[ci] >= '0' && cigar[ci] <= '9') len = len * 10 + (cigar[ci++] - '0');
            if (ci >= C) break;
            char op = cigar[ci++];
            if (op == 'M' || op == '=' || op == 'X') {
                long long os = std::max<long long>(q_start, q_pos);
                long long oe = std::min<long long>(q_end, q_pos + len);
                if (os < oe) {
                    for (long long i = os - q_pos; i < oe - q_pos; i++) {
                        if (r_pos + i < (long long)refsub.size() && seq[q_pos + i] != refsub[r_pos + i]) mm++;
                    }
                }
                q_pos += len; r_pos += len;
            } else if (op == 'I') {
                long long os = std::max<long long>(q_start, q_pos);
                long long oe = std::min<long long>(q_end, q_pos + len);
                if (os < oe) ins += (oe - os);
                q_pos += len;
            } else if (op == 'D') {
                if (q_start <= q_pos && q_pos < q_end) del += len;
                r_pos += len;
            } else if (op == 'S') {
                q_pos += len;
            } else if (op == 'N') {
                r_pos += len;
            }
        }
    }

    // Diagnostics: where do SAM lines go during gen-y.
    struct GenYStats {
        std::atomic<long long> lines{0}, headers{0}, malformed{0}, unmapped{0},
            secondary{0}, no_region{0}, no_ref{0}, recs{0}, regions_ok{0},
            regions_skipped{0}, regions_uncovered{0};
    };

    static void sam_worker(const std::string& fa_path,
                           const std::unordered_map<std::string, std::vector<Region2>>& regions,
                           const std::unordered_map<std::string, FaidxEntry>& fai,
                           GenYStats& st,
                           std::mutex& q_mtx, std::condition_variable& q_cv,
                           std::deque<std::string>& q, bool& producer_done,
                           std::mutex& mtx, std::string& out_buf)
    {
        FILE* fa = fopen(fa_path.c_str(), "rb");   // per-thread handle avoids seek races
        if (!fa) { fprintf(stderr, "[rune] cannot open reference %s\n", fa_path.c_str()); exit(2); }
        std::string line;
        std::string local; local.reserve(1 << 20);
        while (true) {
            {
                std::unique_lock<std::mutex> lk(q_mtx);
                q_cv.wait(lk, [&] { return !q.empty() || producer_done; });
                if (q.empty()) break;
                line = std::move(q.front());
                q.pop_front();
            }
            q_cv.notify_one();
            st.lines++;
            if (line.empty() || line[0] == '@') { st.headers++; continue; }
            size_t t[10]; size_t epos[10]; int nf = 0; size_t pp = 0;
            while (nf < 10) {
                size_t nt = line.find('\t', pp);
                t[nf] = pp; epos[nf] = (nt == std::string::npos) ? std::string::npos : nt;
                if (nt == std::string::npos) { nf++; for (int z = nf; z < 10; z++) { t[z] = std::string::npos; epos[z] = std::string::npos; } break; }
                pp = nt + 1; nf++;
            }
            if (nf < 10 || t[9] == std::string::npos) { st.malformed++; continue; }
            long long flag = atoll(line.c_str() + t[1]);
            if (flag & 0x4) { st.unmapped++; continue; }
            if (flag & (0x100 | 0x800)) { st.secondary++; continue; }
            std::string qname = line.substr(t[0], epos[0] - t[0]);
            auto rit = regions.find(qname);
            if (rit == regions.end()) { st.no_region++; continue; }
            std::string rname = line.substr(t[2], epos[2] - t[2]);
            auto fit = fai.find(rname);
            if (fit == fai.end()) { st.no_ref++; continue; }
            st.recs++;
            long long pos1 = atoll(line.c_str() + t[3]);
            std::string cigar = line.substr(t[5], epos[5] - t[5]);
            std::string seq = line.substr(t[9], epos[9] - t[9]);
            for (auto& c : seq) if (c >= 'a' && c <= 'z') c = char(c - 32);
            long long refspan = 0;
            long long qpos = 0, aln_start = -1, aln_end = 0;
            for (size_t ci = 0; ci < cigar.size();) {
                long long len = 0;
                while (ci < cigar.size() && cigar[ci] >= '0' && cigar[ci] <= '9') len = len * 10 + (cigar[ci++] - '0');
                if (ci >= cigar.size()) break;
                char op = cigar[ci++];
                if (op == 'M' || op == '=' || op == 'X' || op == 'I') {
                    if (aln_start < 0) aln_start = qpos;
                    qpos += len; aln_end = qpos;
                    if (op != 'I') refspan += len;   // I consumes query but not reference
                } else if (op == 'S') {
                    qpos += len;
                } else if (op == 'D' || op == 'N') {
                    refspan += len;
                }
                // H does not consume query coordinates within SEQ; ignore
            }
            std::string refsub = fa_fetch(fa, fit->second, pos1 - 1, pos1 - 1 + refspan);
            // X regions are in ORIGINAL read coordinates; SAM SEQ is in reference
            // orientation (reverse-complemented for flag 0x10). Transform.
            long long L = (long long)seq.size();
            bool rev = (flag & 0x10) != 0;
            for (auto& rg : rit->second) {
                long long qs, qe;
                if (rev) {
                    qs = std::max<long long>(0, L - (long long)rg.e);
                    qe = std::min<long long>(L, L - (long long)rg.s);
                } else {
                    qs = std::max<long long>(0, (long long)rg.s);
                    qe = std::min<long long>(L, (long long)rg.e);
                }
                if (qs >= qe) { st.regions_skipped++; continue; }
                // CIGAR coverage test: aln_start/aln_end is the SEQ interval actually
                // consumed by the alignment (soft-clipped bases excluded). Regions
                // entirely inside the clip cannot be counted and are not emitted.
                long long ov = std::min(qe, aln_end) - std::max(qs, aln_start);
                if (ov <= 0) { st.regions_uncovered++; continue; }
                long long mm = 0, ins = 0, del = 0;
                count_errors(seq, refsub, cigar, qs, qe, mm, ins, del);
                st.regions_ok++;
                // covered: 1=region fully inside the alignment; 0=partially covered
                // (clipped part uncountable, Y underestimated; train drops covered=0)
                int covered = (qs >= aln_start && qe <= aln_end) ? 1 : 0;
                // Output ORIGINAL X coordinates so train can join on
                // (Read_Name, Start); the strand-adjusted window is internal only.
                char tmp[192];
                snprintf(tmp, sizeof tmp, "%s\t%u\t%u\t%lld\t%lld\t%lld\t%d\n",
                         qname.c_str(), rg.s, rg.e, mm, ins, del, covered);
                local += tmp;
                if (local.size() > (1 << 20)) {
                    std::lock_guard<std::mutex> lk(mtx);
                    out_buf += local; local.clear();
                }
            }
        }
        {
            std::lock_guard<std::mutex> lk(mtx);
            out_buf += local;
        }
        fclose(fa);
    }

    const Options& opts_;
};
