import argparse
import sys
import pysam
import pandas as pd

def load_x_regions(x_file, kmer_len):
    """
    Read the X feature file, extract Read_Name and Start, and compute the
    physical intervals. Returns {read_name: [(phys_start, phys_end_exclusive), ...]}
    """
    regions_dict = {}
    print(f"reading X feature file: {x_file}", file=sys.stderr)
    
    # read only the needed columns to save memory
    df = pd.read_csv(x_file, sep='\t', usecols=['Read_Name', 'Start', 'Length'])
    
    for _, row in df.iterrows():
        read_name = row['Read_Name']
        kmer_start = int(row['Start'])
        kmer_count = int(row['Length'])
        
        # coordinate conversion: k-mer index -> physical base coordinates
        # window i covers bases [i, i+k-1]. A run of windows [s, s+m-1] has
        # physical interval (inclusive) [s, s+m-1+k-1]; half-open end = s+m+k-1
        phys_start = kmer_start
        phys_end_exclusive = kmer_start + kmer_count + kmer_len - 1
        
        if read_name not in regions_dict:
            regions_dict[read_name] = []
        regions_dict[read_name].append((phys_start, phys_end_exclusive))
        
    print(f"loaded {len(regions_dict)} reads containing zero regions", file=sys.stderr)
    return regions_dict

def count_errors_in_region(aln, q_start, q_end, ref_fasta):
    """
    Parse CIGAR and the reference sequence; count errors inside the physical
    interval [q_start, q_end).
    """
    mismatch = 0
    insertion = 0
    deletion = 0
    
    q_pos = 0 # current physical coordinate on the read
    r_pos = 0 # current coordinate on the reference
    
    read_seq_str = aln.query_sequence.upper()
    
    # fetch the aligned reference segment directly from the FASTA file
    ref_subseq = ""
    if ref_fasta and not aln.is_unmapped and aln.reference_name:
        # fetch returns the segment in the same orientation as query_sequence
        ref_subseq = ref_fasta.fetch(aln.reference_name, aln.reference_start, aln.reference_end).upper()
    
    for op, length in aln.cigartuples:
        # 0=M, 7==, 8=X
        if op in (0, 7, 8): 
            overlap_start = max(q_start, q_pos)
            overlap_end = min(q_end, q_pos + length)
            
            if overlap_start < overlap_end:
                if ref_fasta and ref_subseq:
                    # base-by-base comparison for mismatches
                    for i in range(overlap_start - q_pos, overlap_end - q_pos):
                        if read_seq_str[q_pos + i] != ref_subseq[r_pos + i]:
                            mismatch += 1
                else:
                    # without ref, fall back to treating M as =; only explicit X counts as mismatch
                    if op == 8: 
                        mismatch += (overlap_end - overlap_start)
            
            q_pos += length
            r_pos += length
            
        elif op == 1: # Insertion (I)
            overlap_start = max(q_start, q_pos)
            overlap_end = min(q_end, q_pos + length)
            if overlap_start < overlap_end:
                insertion += (overlap_end - overlap_start)
            q_pos += length
            
        elif op == 2: # Deletion (D)
            # deletion occurs at q_pos; q_end is half-open, use < to avoid double counting
            if q_start <= q_pos < q_end:
                deletion += length
            r_pos += length
            
        elif op == 4: # Soft clip (S)
            q_pos += length
            
        elif op == 3: # N (Skipped region)
            r_pos += length
            
    return mismatch, insertion, deletion

def main():
    parser = argparse.ArgumentParser(description="Build machine-learning Y labels (error counts) from a BAM file")
    parser.add_argument("-b", "--bam", required=True, help="input: sorted BAM file")
    parser.add_argument("-x", "--x_file", required=True, help="input: X feature TSV file (for coordinate alignment)")
    parser.add_argument("-r", "--ref", required=True, help="input: reference genome FASTA (required for exact mismatch calling)")
    parser.add_argument("-k", "--kmer_len", type=int, required=True, help="k-mer length used when building the X file")
    
    args = parser.parse_args()
    
    # 1. load X intervals
    regions_dict = load_x_regions(args.x_file, args.kmer_len)
    
    # 2. load reference genome
    print(f"loading reference genome: {args.ref}", file=sys.stderr)
    ref_fasta = pysam.FastaFile(args.ref)
    
    # 3. open BAM file
    print("processing BAM and counting errors...", file=sys.stderr)
    bam = pysam.AlignmentFile(args.bam, "rb")
    
    # Y header (Covered: 1=fully inside alignment, 0=partially covered / Y underestimated)
    print("Read_Name\tRegion_Start\tRegion_End\tMismatch\tInsertion\tDeletion\tCovered")
    
    # iterate over BAM (if sorted+indexed, bam.fetch() per chromosome is faster)
    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
            
        read_name = aln.query_name
        
        if read_name in regions_dict:
            # X regions are in original read orientation; for reverse-strand
            # alignments (flag 0x10) query_sequence is in reference orientation
            # (revcomp-ed), so the region coordinates must be transformed.
            L = aln.query_length
            # query interval actually consumed by the CIGAR (SEQ coordinates; soft clips excluded)
            a_s = aln.query_alignment_start
            a_e = aln.query_alignment_end
            for s, e in regions_dict[read_name]:
                if aln.is_reverse:
                    q_start, q_end = L - e, L - s
                else:
                    q_start, q_end = s, e
                # boundary guard: clamp to the actual read length
                q_start = max(0, q_start)
                q_end = min(aln.query_length, q_end)

                if q_start >= q_end:
                    continue

                # regions fully inside the soft clip cannot be counted; skip them
                if min(q_end, a_e) - max(q_start, a_s) <= 0:
                    continue

                covered = 1 if (q_start >= a_s and q_end <= a_e) else 0
                mm, ins, dele = count_errors_in_region(aln, q_start, q_end, ref_fasta)

                # output original X coordinates (train joins on Start);
                # the reverse-strand transform is only used internally for counting
                print(f"{read_name}\t{s}\t{e}\t{mm}\t{ins}\t{dele}\t{covered}")
                
    bam.close()
    ref_fasta.close()
    print("Y label construction finished", file=sys.stderr)

if __name__ == "__main__":
    main()
