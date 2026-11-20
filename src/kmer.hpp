#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

#include "sequence.hpp"

// 2-bit canonical k-mer encoder. Supports k <= 32.
struct KmerCodec {
    int k;
    uint64_t mask;
    explicit KmerCodec(int kk) : k(kk), mask(kk >= 32 ? ~0ULL : ((1ULL << (2 * kk)) - 1)) {}

    bool encode(const char* s, size_t i, uint64_t& out) const {
        uint64_t fwd = 0, rc = 0;
        for (int j = 0; j < k; j++) {
            char c = s[i + j];
            int v;
            switch (c) {
                case 'A': case 'a': v = 0; break;
                case 'C': case 'c': v = 1; break;
                case 'G': case 'g': v = 2; break;
                case 'T': case 't': v = 3; break;
                default: return false;
            }
            fwd = ((fwd << 2) | (uint64_t)v) & mask;
            rc = (rc >> 2) | ((uint64_t)(3 - v) << (2 * (k - 1)));
        }
        out = fwd < rc ? fwd : rc;
        return true;
    }
};

// Canonical k-mer count table (jellyfish -C compatible).
struct KmerTable {
    std::unordered_map<uint64_t, uint32_t> m;

    void load(const std::string& path, const KmerCodec& codec) {
        FILE* f = open_text(path);
        std::string line;
        size_t n = 0;
        while (read_line(f, line)) {
            if (line.empty()) continue;
            size_t tab = line.find('\t'), sp = line.find(' ');
            size_t sep = (tab == std::string::npos) ? sp : (sp == std::string::npos ? tab : std::min(tab, sp));
            if (sep == std::string::npos) continue;
            std::string kmer = line.substr(0, sep);
            uint32_t c = (uint32_t)strtoul(line.c_str() + sep + 1, nullptr, 10);
            if ((int)kmer.size() != codec.k) continue;
            uint64_t code;
            if (codec.encode(kmer.c_str(), 0, code)) { m[code] += c; n++; }
        }
        close_text(f, path);
        fprintf(stderr, "[rune] kmer table: %zu distinct canonical k-mers\n", n);
    }

    inline uint32_t get(uint64_t code) const {
        auto it = m.find(code);
        return it == m.end() ? 0u : it->second;
    }
};

// Median of window counts, sampled from the first sample_reads reads.
inline long long estimate_median_cov(const std::string& reads_path, const KmerTable& table,
                                     const KmerCodec& codec, int k, size_t sample_reads = 4000)
{
    SeqReader rd(reads_path);
    std::string name, seq;
    std::vector<uint64_t> hist(1 << 20, 0);
    uint64_t total = 0; size_t taken = 0;
    while (taken < sample_reads && rd.next(name, seq)) {
        if ((int)seq.size() < k) continue;
        for (size_t i = 0; i + k <= seq.size(); i++) {
            uint64_t code;
            if (codec.encode(seq.c_str(), i, code)) {
                uint32_t c = table.get(code);
                if (c < hist.size()) { hist[c]++; total++; }
            }
        }
        taken++;
    }
    if (total == 0) return 0;
    uint64_t half = total / 2, acc = 0;
    for (size_t c = 0; c < hist.size(); c++) { acc += hist[c]; if (acc >= half) return (long long)c; }
    return 0;
}
