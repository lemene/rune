# rune — Runs Uncover Nucleotide Errors

**rune** is a reference-free, alignment-free tool for quantifying sequencing
errors in long reads. Reads are encoded as binary sequences where each position
marks whether its k-mer is low-frequency (0, putatively erroneous) or not (1).
The length of each run of zeros geometrically encodes the number of sequencing
errors inside it; per-run predictions are aggregated into per-read and
dataset-level error rates.

## Repository layout

```
rune/
├── src/                      # Main C++ implementation (single binary, four subcommands)
│   ├── rune.cpp              # main + subcommand dispatch + training/inference
│   ├── dataset.hpp           # gen-x / gen-y dataset construction
│   ├── kmer.hpp              # k-mer encoding, count table, depth estimation
│   ├── options.hpp           # CLI argument parsing
│   └── sequence.hpp          # FASTA/FASTQ/gz reading, faidx fetch
├── scripts/                  # Python reference implementations (cross-validation)
│   ├── build_feature_X.py    # Python mirror of gen-x (per-segment aligned)
│   ├── build_label_Y.py      # Python mirror of gen-y (pysam)
│   ├── train_and_save_model.py  # sklearn RF training (grouped split + Covered filter)
│   ├── export_rf.py          # sklearn model -> C++-compatible model text
│   └── predict_dataset.py    # inference with the sklearn model (rf_error_model.pkl)
└── .codebuddy/rules/rune-project.md   # coding constraints + fixed-bug list
```

## Build

```bash
g++ -O3 -std=c++17 -pthread src/rune.cpp -o rune    # or: make
```

No third-party dependencies: gzipped inputs go through a `gzip` pipe, BAM
through a `samtools view` pipe.

## Pipeline

### 0. Counting and mapping

The former `scripts/01_count.sh` and `scripts/map.sh` are merged here and no
longer exist as separate files:

```bash
# 0.1 k-mer counting (-C is mandatory: canonical form, aligned with all lookups)
zcat ../badread-r10/data/reads.fastq.gz | \
  ~/tool/kmer-counter/jellyfish-2.3.1/bin/jellyfish count \
    -s 2G -m 19 -C -t 8 -o rd_k19.jf /dev/stdin

# 0.2 Dump the count table (-c count column, -t TAB separator; -L 2 drops
#     singletons and saves ~80% disk — downstream treats missing k-mers as
#     count=0, results are unchanged)
jellyfish dump -c -t -L 2 rd_k19.jf > rd_k19.tsv

# 0.3 Alignment (needed by gen-y: SAM/BAM + reference FASTA)
minimap2 -t 8 -ax map-ont ../ref/ref.fna.gz ../badread-r10/data/reads.fastq \
  | samtools sort -@ 8 -o rd_2_ref.sorted.bam /dev/stdin \
  && samtools index rd_2_ref.sorted.bam
```

### 1. Generate features X

```bash
# C++ (recommended; depth/threshold/strong_cov are all estimated automatically)
rune gen-x -i reads.fastq -t rd_k19.tsv -k 19 -o X.tsv
#   optional: --global-cov C   --t-ratio 0.5   --min-len 1
#             --strong-cov N   --max-gap 3    --threads N
# gzipped input is supported directly (-i reads.fastq.gz)

# Python reference (output should match the C++ version byte-for-byte,
# up to 4th-decimal float rounding in ~14 lines)
python scripts/build_feature_X.py -i reads.fastq -k rd_k19.tsv --kmer_len 19 > X_py.tsv
```

### 2. Generate labels Y (requires alignment)

```bash
rune gen-y -b rd_2_ref.sorted.bam -x X.tsv -r ref.fna -k 19 -o Y.tsv
#   -b accepts .sam / .bam (BAM via a samtools pipe) / - (stdin stream)

python scripts/build_label_Y.py -b rd_2_ref.sorted.bam -x X.tsv -r ref.fna -k 19 > Y_py.tsv
```

### 3. Train

```bash
# C++ (random forest, tree-level parallelism, model stored as text)
rune train -x X.tsv -y Y.tsv -m model.txt --trees 100

# sklearn route (optional, for comparison)
python scripts/train_and_save_model.py -x X.tsv -y Y.tsv    # -> rf_error_model.pkl
python scripts/export_rf.py                                  # -> rf_model.txt (C++-readable)
```

### 4. Predict

```bash
# from an X file
rune predict -m model.txt -x X.tsv -r reads.fastq -o out_
# directly from reads (runs gen-x internally)
rune predict -m model.txt -i reads.fastq -t rd_k19.tsv -k 19 -o out_
# outputs: out_predicted_regions.tsv / out_predicted_per_read.tsv / global report
```

## Key conventions

**Coordinate join**: X's `(Start, Length)` corresponds to Y's
`(Start, Start+Length+k-1)` (the physical window extends by k-1; with k=19 the
Y interval is 18 longer than the X interval). Training joins on
`(Read_Name, X.Start) == Y.Region_Start` — do not mistake the +18 for a mismatch.

**The Covered column** (7th column of Y): 1 = region fully inside the aligned
block; 0 = region partially extends into a soft clip (Y is systematically
underestimated there). Regions entirely outside the alignment are **not
emitted**. Training on both sides discards `Covered=0` rows.

**Reverse strand**: in SAM/BAM, reverse-strand alignments (flag 0x10) store SEQ
in reference orientation. gen-y transforms regions via `[s,e) -> [L-e, L-s)`
internally for counting but **outputs original X coordinates** (so the join
works).

## Parameter alignment (Python reference ⇄ C++)

| Meaning | C++ | Python |
|---|---|---|
| Depth d | `--global-cov` (default: median window count over first 4000 reads) | `--global_cov` (same estimator) |
| Threshold T | `max(2, floor(d·t_ratio))`, t_ratio default 0.5 | same |
| Strong-signal cutoff | `--strong-cov` (default `max(10, floor(0.3·d))`) | `--strong_cov` (same) |
| Min physical length | `--min-len` default 1 (m=1 short runs contain real events; do not drop) | `--min_len` default 1 |
| Merge gap | `--max-gap` default 3 (no merge if the gap contains count ≥ strong_cov) | `--max_gap` |

## Critical fixed bugs (do not regress; see .codebuddy/rules/rune-project.md)

1. **Canonical lookup**: jellyfish `-C` tables store the lexicographically
   smaller of kmer/revcomp; every lookup must canonicalize first.
2. **Y coordinates**: physical interval is half-open `[s, s+m+k-1)`; reverse
   strands must be transformed before counting.
3. **min_len**: 68–90% of m=1 isolated short runs contain real events — mark
   them, never drop them.
4. **Multithreaded dispatch**: a shared reader plus row-index modulo filtering
   silently drops rows (the cause of "gen-y produced too few rows"); use a
   producer-consumer queue.
5. **Thread safety**: reference-sequence fetches need a per-thread file handle.
