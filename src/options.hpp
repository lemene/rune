#pragma once

#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>

// Parameters for the gen-x step (feature extraction).
struct XParams {
    int k = 19;
    double t_ratio = 0.5;
    long long global_cov = -1;   // < 0: estimate from median window count
    int max_gap = 3;
    long long strong_cov = -1;   // < 0: max(10, 0.3 * coverage)
    long long min_len = 1;       // physical length floor; 1 keeps short segments
};

// Parameters for the train step (random forest).
struct TrainOpts {
    int trees = 100;
    int depth = 15;
    int leaf = 5;
    int seed = 42;
    double test = 0.2;
};

// All command-line parameters, centralized in one place.
// The first positional argument is a named subcommand; numeric aliases (1-4)
// are kept for backward compatibility.
struct Options {
    enum class Command { GenX, GenY, Train, Predict, Unknown };

    Command command = Command::Unknown;
    bool helpRequested = false;
    XParams xp;
    TrainOpts train;
    int threads = 0;   // 0: auto (hardware_concurrency)

    // Shared paths.
    std::string reads;      // -i
    std::string kmerFile;   // -t
    std::string output;     // -o

    // gen-y / train / predict paths.
    std::string sam;        // -b with gen-y
    std::string refFasta;   // -r with gen-y
    std::string xFile;      // -x
    std::string yFile;      // -y
    std::string model = "model.txt";  // -m
    bool modelGiven = false;

    // predict-specific.
    bool fromReads = false;          // -i with predict: generate X internally
    std::string lenReads;            // -r with predict: reads for length stats
    long long totalBases = -1;       // -b with predict

    static Command toCommand(const std::string& s) {
        if (s == "gen-x" || s == "1") return Command::GenX;
        if (s == "gen-y" || s == "2") return Command::GenY;
        if (s == "train" || s == "3") return Command::Train;
        if (s == "predict" || s == "4") return Command::Predict;
        return Command::Unknown;
    }

    bool parse(int argc, char* argv[]) {
        if (argc < 2) {
            printUsage(std::cerr);
            return false;
        }
        command = toCommand(argv[1]);
        if (command == Command::Unknown) {
            std::cerr << "[rune] unknown subcommand: " << argv[1] << "\n";
            printUsage(std::cerr);
            return false;
        }

        for (int i = 2; i < argc; ) {
            std::string key = argv[i++];
            std::string val;

            // Accept both "--key value" and "--key=value".
            if (key.rfind("--", 0) == 0) {
                size_t eq = key.find('=');
                if (eq != std::string::npos) {
                    val = key.substr(eq + 1);
                    key.resize(eq);
                }
            }
            if (key == "-h" || key == "--help") {
                helpRequested = true;
                printUsage(std::cout);
                return true;
            }
            if (val.empty() && i < argc) val = argv[i++];
            if (!apply(key, val)) {
                std::cerr << "[rune] unknown option: " << key << "\n";
                printUsage(std::cerr);
                return false;
            }
        }

        if (threads <= 0) {
            threads = std::thread::hardware_concurrency();
            if (threads <= 0) threads = 4;
        }
        if (!validate()) {
            printUsage(std::cerr);
            return false;
        }
        return true;
    }

    bool apply(const std::string& key, const std::string& val) {
        if (key == "-i" || key == "--input") {
            reads = val;
            if (command == Command::Predict) fromReads = true;
        } else if (key == "-t" || key == "--kmers") {
            kmerFile = val;
        } else if (key == "-k" || key == "--kmer-len") {
            xp.k = std::atoi(val.c_str());
        } else if (key == "-o" || key == "--output") {
            output = val;
        } else if (key == "--global-cov") {
            xp.global_cov = std::atoll(val.c_str());
        } else if (key == "--t-ratio") {
            xp.t_ratio = std::atof(val.c_str());
        } else if (key == "--max-gap") {
            xp.max_gap = std::atoi(val.c_str());
        } else if (key == "--strong-cov") {
            xp.strong_cov = std::atoll(val.c_str());
        } else if (key == "--min-len") {
            xp.min_len = std::atoll(val.c_str());
        } else if (key == "-b") {
            if (command == Command::GenY) sam = val;
            else totalBases = std::atoll(val.c_str());  // predict
        } else if (key == "--bam") {
            sam = val;
        } else if (key == "--total-bases") {
            totalBases = std::atoll(val.c_str());
        } else if (key == "-x" || key == "--x-file") {
            xFile = val;
        } else if (key == "-y" || key == "--y-file") {
            yFile = val;
        } else if (key == "-r") {
            if (command == Command::GenY) refFasta = val;
            else lenReads = val;  // predict
        } else if (key == "--ref") {
            refFasta = val;
        } else if (key == "-m" || key == "--model") {
            model = val;
            modelGiven = true;
        } else if (key == "--trees") {
            train.trees = std::atoi(val.c_str());
        } else if (key == "--depth") {
            train.depth = std::atoi(val.c_str());
        } else if (key == "--leaf") {
            train.leaf = std::atoi(val.c_str());
        } else if (key == "--seed") {
            train.seed = std::atoi(val.c_str());
        } else if (key == "--test") {
            train.test = std::atof(val.c_str());
        } else if (key == "--threads") {
            threads = std::atoi(val.c_str());
        } else {
            return false;
        }
        return true;
    }

    bool validate() const {
        switch (command) {
        case Command::GenX:
            if (reads.empty() || kmerFile.empty()) {
                std::cerr << "[rune] gen-x requires -i reads -t kmers.tsv\n";
                return false;
            }
            break;
        case Command::GenY:
            if (sam.empty() || xFile.empty() || refFasta.empty()) {
                std::cerr << "[rune] gen-y requires -b sam -x X.tsv -r ref.fna\n";
                return false;
            }
            break;
        case Command::Train:
            if (xFile.empty() || yFile.empty()) {
                std::cerr << "[rune] train requires -x X.tsv -y Y.tsv\n";
                return false;
            }
            break;
        case Command::Predict:
            if (!modelGiven) {
                std::cerr << "[rune] predict requires -m model.txt\n";
                return false;
            }
            if (!fromReads && xFile.empty()) {
                std::cerr << "[rune] predict requires -x X.tsv or -i data mode\n";
                return false;
            }
            if (fromReads && (reads.empty() || kmerFile.empty())) {
                std::cerr << "[rune] predict data mode requires -i reads -t kmers.tsv\n";
                return false;
            }
            break;
        default:
            return false;
        }
        return true;
    }

    static void printUsage(std::ostream& os) {
        os << "rune - Runs Uncover Nucleotide Errors\n"
           << "Usage: rune <command> [options]\n"
           << "Commands (numeric aliases 1-4 also accepted):\n"
           << "  gen-x    Generate feature matrix X:\n"
           << "           rune gen-x -i reads.fq -t kmers.tsv -k 19 [-o X.tsv]\n"
           << "           [--global-cov C] [--t-ratio 0.5] [--max-gap 3]\n"
           << "           [--strong-cov N] [--min-len 1] [--threads N]\n"
           << "  gen-y    Generate labels Y (SAM text or BAM; BAM needs samtools in PATH):\n"
           << "           rune gen-y -b aln.sam -x X.tsv -r ref.fna -k 19 [-o Y.tsv]\n"
           << "           rune gen-y -b aln.bam -x X.tsv -r ref.fna -k 19 [-o Y.tsv]\n"
           << "  train    Train random forest:\n"
           << "           rune train -x X.tsv -y Y.tsv -m model.txt [--trees 100]\n"
           << "           [--depth 15] [--leaf 5] [--seed 42] [--test 0.2] [--threads N]\n"
           << "  predict  Predict from X file:\n"
           << "           rune predict -m model.txt -x X.tsv [-r reads.fq | -b NB] [-o prefix]\n"
           << "           or from raw reads: rune predict -m model.txt -i reads.fq\n"
           << "                              -t kmers.tsv -k 19 [-o prefix]\n";
    }
};
