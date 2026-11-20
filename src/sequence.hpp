#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
inline void fseek64(FILE* f, long long off, int org) { _fseeki64(f, off, org); }
#else
#define POPEN popen
#define PCLOSE pclose
inline void fseek64(FILE* f, long long off, int org) { fseeko(f, (off_t)off, org); }
#endif

inline bool ends_with(const std::string& s, const char* suf) {
    size_t n = std::strlen(suf);
    return s.size() >= n && s.compare(s.size() - n, n, suf) == 0;
}

#ifdef _WIN32
#define RUNE_DEV_NULL "NUL"
#else
#define RUNE_DEV_NULL "/dev/null"
#endif

inline bool is_bam(const std::string& path) {
    return ends_with(path, ".bam") || ends_with(path, ".BAM");
}

// Unified text input: '-' is stdin; .gz through gzip; .bam through samtools view.
inline FILE* open_text(const std::string& path) {
    if (path == "-") return stdin;
    if (is_bam(path)) {
        std::string probe = std::string("samtools --version > ") + RUNE_DEV_NULL + " 2>&1";
        if (system(probe.c_str()) != 0) {
            fprintf(stderr, "[rune] reading BAM requires samtools in PATH\n");
            exit(2);
        }
        std::string cmd = "samtools view \"" + path + "\"";
        FILE* f = POPEN(cmd.c_str(), "r");
        if (!f) { fprintf(stderr, "[rune] cannot open BAM: %s\n", path.c_str()); exit(2); }
        return f;
    }
    if (ends_with(path, ".gz")) {
        std::string cmd = "gzip -dc \"" + path + "\"";
        FILE* f = POPEN(cmd.c_str(), "r");
        if (!f) { fprintf(stderr, "[rune] cannot open gz: %s\n", path.c_str()); exit(2); }
        return f;
    }
    FILE* f = fopen(path.c_str(), "r");
    if (!f) { fprintf(stderr, "[rune] cannot open: %s\n", path.c_str()); exit(2); }
    return f;
}

inline void close_text(FILE* f, const std::string& path) {
    if (path == "-") return;
    if (is_bam(path) || ends_with(path, ".gz")) PCLOSE(f);
    else fclose(f);
}

inline bool read_line(FILE* f, std::string& line) {
    line.clear();
    int c; bool any = false;
    while ((c = fgetc(f)) != EOF) {
        any = true;
        if (c == '\n') break;
        if (c != '\r') line.push_back((char)c);
    }
    return any || !line.empty();
}

// Streaming FASTA/FASTQ reader.
struct SeqReader {
    FILE* f;
    std::string path;
    std::string line;
    bool has_pending = false;
    std::string pending;

    explicit SeqReader(const std::string& p) : f(open_text(p)), path(p) {}
    ~SeqReader() { if (path != "-") close_text(f, path); }

    bool next(std::string& name, std::string& seq) {
        while (true) {
            if (has_pending) { line = pending; has_pending = false; }
            else if (!read_line(f, line)) return false;
            if (line.empty()) continue;
            char h0 = line[0];
            if (h0 == '>' || h0 == '@') {
                size_t sp = line.find_first_of(" \t");
                name = (sp == std::string::npos) ? line.substr(1) : line.substr(1, sp - 1);
                seq.clear();
                if (h0 == '>') {
                    while (read_line(f, line) && !line.empty() && line[0] != '>') seq += line;
                    if (!line.empty() && line[0] == '>') { pending = line; has_pending = true; }
                    return true;
                } else {
                    read_line(f, line); seq = line;      // SEQ
                    read_line(f, line);                  // '+'
                    read_line(f, line);                  // QUAL
                    return true;
                }
            }
            // Skip stray lines before a header.
        }
    }
};
