import argparse
import sys
import gzip
import pandas as pd
import numpy as np
import joblib

def open_maybe_gzip(file_path, mode='rt'):
    """Helper: open .gz files transparently"""
    if file_path.endswith('.gz'):
        return gzip.open(file_path, mode)
    return open(file_path, mode)

def count_total_bases(fastq_or_fasta):
    """Count total sequenced bases in a FASTA/FASTQ file (supports .gz)"""
    print(f"counting total bases in: {fastq_or_fasta} ...", file=sys.stderr)
    total_bases = 0
    with open_maybe_gzip(fastq_or_fasta, 'rt') as f:
        while True:
            header = f.readline()
            if not header: break
            if header.startswith('>'):
                seq_lines = []
                while True:
                    line = f.readline()
                    if not line or line.startswith('>'):
                        if line: f.seek(f.tell() - len(line))
                        break
                    seq_lines.append(line.strip())
                total_bases += len(''.join(seq_lines))
            elif header.startswith('@'):
                seq = f.readline().strip()
                f.readline() # +
                f.readline() # quality
                total_bases += len(seq)
    
    print(f"total bases: {total_bases}", file=sys.stderr)
    return total_bases

def main():
    parser = argparse.ArgumentParser(description="Load a trained model and predict the dataset sequencing error rate")
    parser.add_argument("-m", "--model_file", required=True, help="input: trained model file (e.g. rf_error_model.pkl)")
    parser.add_argument("-f", "--feature_file", required=True, help="input: feature list matching the model (e.g. rf_model_features.pkl)")
    parser.add_argument("-x", "--x_file", required=True, help="input: X feature file to predict")
    
    # denominator for the overall error rate
    parser.add_argument("-r", "--reads_file", default=None, help="optional: FASTQ/FASTA file to compute total bases automatically")
    parser.add_argument("-b", "--total_bases", type=int, default=None, help="optional: total base count of the dataset, given directly")
    
    # output files
    parser.add_argument("-o", "--out_tsv", default="predicted_regions.tsv", help="output: detailed per-region predictions (default: predicted_regions.tsv)")
    parser.add_argument("--no_log", action="store_true", help="add this flag if the log1p transform was disabled during training")
    
    args = parser.parse_args()

    # ================= 1. load model and features =================
    print("loading model and feature list...", file=sys.stderr)
    try:
        model = joblib.load(args.model_file)
        feature_cols = joblib.load(args.feature_file)
    except Exception as e:
        print(f"failed to load model: {e}", file=sys.stderr)
        sys.exit(1)

    # ================= 2. read and predict =================
    print(f"reading feature file: {args.x_file}", file=sys.stderr)
    try:
        df = pd.read_csv(args.x_file, sep='\t')
    except Exception as e:
        print(f"failed to read data: {e}", file=sys.stderr)
        sys.exit(1)

    # check for missing features
    missing_features = [f for f in feature_cols if f not in df.columns]
    if missing_features:
        print(f"error: X file is missing feature columns: {missing_features}", file=sys.stderr)
        sys.exit(1)
        
    print(f"predicting {len(df)} regions...", file=sys.stderr)
    X = df[feature_cols].fillna(0)
    
    # predict
    pred_target = model.predict(X)
    
    # invert back to the real scale
    if not args.no_log:
        pred_errors = np.expm1(pred_target)
    else:
        pred_errors = pred_target
        
    # guard against negative predictions (rare)
    pred_errors = np.clip(pred_errors, 0, None)
    
    # add predictions back to the dataframe
    df['Predicted_Errors'] = pred_errors
    
    # save detailed predictions
    df.to_csv(args.out_tsv, sep='\t', index=False)
    print(f"detailed predictions saved to: {args.out_tsv}", file=sys.stderr)

    # ================= 3. overall error rate =================
    total_predicted_errors = df['Predicted_Errors'].sum()
    
    # try to obtain the total base count
    total_bases = 0
    if args.reads_file:
        total_bases = count_total_bases(args.reads_file)
    elif args.total_bases:
        total_bases = args.total_bases
    else:
        print("\n[warning] neither -r nor -b provided; cannot compute the overall\nerror rate (missing total base count).", file=sys.stderr)

    # ================= 4. final report =================
    print("\n" + "="*40)
    print("          dataset-level error rate report")
    print("="*40)
    print(f"predicted regions:      {len(df)}")
    print(f"total predicted errors (mismatch + ins + del): {total_predicted_errors:.0f}")
    
    if total_bases > 0:
        global_error_rate = total_predicted_errors / total_bases
        print(f"total sequenced bases:  {total_bases}")
        print(f"overall base error rate:  {global_error_rate:.6f}")
        print(f"mean base accuracy:       {(1 - global_error_rate) * 100:.4f}%")
        print(f"Phred quality (Q-score):  {-10 * np.log10(global_error_rate):.2f}")
    else:
        print("(cannot compute the error rate: total base count missing)")
    print("="*40 + "\n")

if __name__ == "__main__":
    main()
