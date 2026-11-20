import argparse
import sys
import gzip
import numpy as np

def open_maybe_gzip(file_path, mode='rt'):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, mode)
    return open(file_path, mode)

def parse_fasta_fastq(file_path):
    with open_maybe_gzip(file_path, 'rt') as f:
        while True:
            header = f.readline()
            if not header: break
            if header.startswith('>'):
                read_name = header.strip().split()[0][1:]
                seq_lines = []
                while True:
                    line = f.readline()
                    if not line or line.startswith('>'):
                        if line: f.seek(f.tell() - len(line))
                        break
                    seq_lines.append(line.strip())
                yield read_name, ''.join(seq_lines)
            elif header.startswith('@'):
                read_name = header.strip().split()[0][1:]
                seq = f.readline().strip()
                f.readline()
                f.readline()
                yield read_name, seq

def merge_and_filter_regions(raw_regions, kmer_list, max_gap, strong_cov, min_len, k):
    """
    Keep the original merge strategy so extracted error intervals stay
    physically plausible.
    """
    if not raw_regions: return []
    merged = []
    curr_start, curr_end = raw_regions[0]
    
    for next_start, next_end in raw_regions[1:]:
        gap = next_start - curr_end - 1
        gap_kmers_cov = kmer_list[curr_end + 1 : next_start]
        max_cov_in_gap = max(gap_kmers_cov) if gap_kmers_cov else 0
        
        if gap <= max_gap and max_cov_in_gap < strong_cov:
            curr_end = next_end
        else:
            merged.append((curr_start, curr_end))
            curr_start, curr_end = next_start, next_end
    merged.append((curr_start, curr_end))
    
    filtered = []
    for start, end in merged:
        kmer_count = end - start + 1
        physical_len = kmer_count + k - 1
        if physical_len >= min_len:
            filtered.append((start, end))
    return filtered

def canon(kmer):
    """Canonicalize: align with jellyfish -C tables. The table stores the
    lexicographically smaller of kmer/revcomp; querying raw kmers directly misses
    ~50% of windows (reverse strand), which would all be misjudged as 0."""
    rc = kmer[::-1].translate(_RC)
    return kmer if kmer <= rc else rc

_RC = str.maketrans("ACGTN", "TGCAN")

def iter_window_counts(seq, k, kmer_counts):
    """All window counts of one read (identical lookup/canonicalization to C++ gen-x)"""
    for i in range(len(seq) - k + 1):
        yield kmer_counts.get(canon(seq[i:i+k]), 0)

def estimate_median_cov(reads_path, k, kmer_counts, sample_reads=4000):
    """Same as C++ estimate_median_cov: median window count over the first
    sample_reads reads. Error windows only pollute the low tail (f~1-exp(-eps*k)),
    so the median is nearly unbiased (~0.4% error observed)."""
    hist = {}
    total = 0
    taken = 0
    for _, seq in parse_fasta_fastq(reads_path):
        if len(seq) < k: continue
        for c in iter_window_counts(seq, k, kmer_counts):
            hist[c] = hist.get(c, 0) + 1
            total += 1
        taken += 1
        if taken >= sample_reads: break
    if total == 0: return 0.0
    half = total // 2
    acc = 0
    for c in sorted(hist):
        acc += hist[c]
        if acc >= half:
            return float(c)
    return 0.0

def process_read_features(read_name, seq, kmer_counts, global_dynamic_T, k, max_gap, strong_cov, min_len):
    """
    Compute read-level global features and extract local relative features
    for every zero run.
    """
    seq_len = len(seq)
    if seq_len < k: return

    binary_seq = []
    kmer_list = []
    
    # 1. build the 0/1 sequence
    for i in range(seq_len - k + 1):
        kmer = seq[i:i+k]
        count = kmer_counts.get(canon(kmer), 0)
        is_valid = 1 if count >= global_dynamic_T else 0
        binary_seq.append(is_valid)
        kmer_list.append(count)

    total_kmers = len(binary_seq)
    if total_kmers == 0: return

    # ================= read-level global features =================
    read_zero_ratio = binary_seq.count(0) / total_kmers
    # mean coverage over valid (>= threshold) windows of this read
    valid_counts = [c for c in kmer_list if c >= global_dynamic_T]
    read_avg_cov = np.mean(valid_counts) if valid_counts else 1e-5 
    # guard against division by zero for garbage reads
    # ========================================================

    # 2. locate runs of zeros
    raw_regions = []
    in_zero_region = False
    start_idx = 0
    for i in range(total_kmers):
        if binary_seq[i] == 0 and not in_zero_region:
            in_zero_region = True
            start_idx = i
        elif binary_seq[i] == 1 and in_zero_region:
            in_zero_region = False
            raw_regions.append((start_idx, i - 1))
    if in_zero_region:
        raw_regions.append((start_idx, total_kmers - 1))

    # 3. merge and filter
    final_regions = merge_and_filter_regions(raw_regions, kmer_list, max_gap, strong_cov, min_len, k)

    # 4. extract local relative features and output
    for start_idx, end_idx in final_regions:
        region_len = end_idx - start_idx + 1
        
        # absolute coverage inside the region
        region_counts = kmer_list[start_idx:end_idx+1]
        region_abs_cov = np.mean(region_counts)
        
        # boundary coverage (5 windows each side)
        left_win = kmer_list[max(0, start_idx - 5) : start_idx]
        right_win = kmer_list[end_idx + 1 : min(total_kmers, end_idx + 6)]
        boundary_counts = left_win + right_win
        boundary_abs_cov = np.mean(boundary_counts) if boundary_counts else 0.0
        
        # normalize by read_avg_cov
        boundary_cov_ratio = boundary_abs_cov / read_avg_cov
        region_rel_cov = region_abs_cov / read_avg_cov
        
        # position relative to the global threshold
        threshold_rel_pos = (region_abs_cov - global_dynamic_T) / global_dynamic_T

        # output: read name, start (join key), length, read-level features
        print(f"{read_name}\t{start_idx}\t{region_len}\t{read_zero_ratio:.4f}\t"
              f"{boundary_cov_ratio:.4f}\t{region_rel_cov:.4f}\t{threshold_rel_pos:.4f}")

def main():
    parser = argparse.ArgumentParser(
        description="Python reference implementation of C++ gen-x (per-segment aligned, for cross-validation)")
    parser.add_argument("-i", "--input", required=True, help="Reads FASTA/FASTQ")
    parser.add_argument("-k", "--kmer_file", required=True, help="k-mer counts file (jellyfish dump -c -t)")
    parser.add_argument("--global_cov", type=float, default=None,
                        help="coverage override (default: auto-estimated as median window count, same as C++)")
    parser.add_argument("--t_ratio", type=float, default=0.5,
                        help="threshold ratio, T = max(2, floor(d * t_ratio)) (same as C++)")
    parser.add_argument("--kmer_len", type=int, required=True, help="k-mer length")
    parser.add_argument("--max_gap", type=int, default=3, help="merge gap")
    parser.add_argument("--strong_cov", type=int, default=None,
                        help="strong signal cutoff (default max(10, floor(0.3*d)), same as C++)")
    parser.add_argument("--min_len", type=int, default=1,
                        help="min physical length (default 1: m=1 short runs contain real events, mark don't drop)")
    parser.add_argument("--sample_reads", type=int, default=4000,
                        help="reads sampled for depth estimation (same as C++)")

    args = parser.parse_args()

    print("loading k-mer table...", file=sys.stderr)
    kmer_counts = {}
    with open_maybe_gzip(args.kmer_file, 'rt') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                kmer_counts[parts[0]] = int(parts[1])

    # depth d: use --global_cov if given, else median window count (same as C++)
    if args.global_cov is not None:
        d = args.global_cov
        print(f"coverage set by --global_cov: d = {d}", file=sys.stderr)
    else:
        d = estimate_median_cov(args.input, args.kmer_len, kmer_counts, args.sample_reads)
        print(f"estimated coverage d = {d:.1f} (median window count over first {args.sample_reads} reads)", file=sys.stderr)
    if d < 5:
        print("[warn] coverage too low, threshold unreliable", file=sys.stderr)

    global_dynamic_T = max(2, int(np.floor(d * args.t_ratio)))
    strong_cov = args.strong_cov if args.strong_cov is not None else max(10, int(np.floor(0.3 * d)))
    print(f"T={global_dynamic_T} (t_ratio={args.t_ratio}) strong_cov={strong_cov} "
          f"max_gap={args.max_gap} min_len={args.min_len}", file=sys.stderr)

    print("extracting features...", file=sys.stderr)
    # X file header (absolute coverage intentionally omitted)
    print("Read_Name\tStart\tLength\tRead_Zero_Ratio\tBoundary_Cov_Ratio\tRegion_Rel_Cov\tThreshold_Rel_Pos")

    for read_name, seq in parse_fasta_fastq(args.input):
        process_read_features(read_name, seq, kmer_counts, global_dynamic_T, args.kmer_len,
                              args.max_gap, strong_cov, args.min_len)

if __name__ == "__main__":
    main()
