import argparse
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from sklearn.model_selection import GroupShuffleSplit
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

def load_and_merge_data(x_file, y_file):
    """Read X and Y files and merge them"""
    print(f"reading X ({x_file}) and Y ({y_file})...")
    x_df = pd.read_csv(x_file, sep='\t')
    y_df = pd.read_csv(y_file, sep='\t')
    
    # compute the true total error count Y
    y_df['Total_Errors'] = y_df['Mismatch'] + y_df['Insertion'] + y_df['Deletion']

    # drop partially covered labels (Covered=0: region extends into soft clip, Y underestimated)
    if 'Covered' in y_df.columns:
        n0 = len(y_df)
        y_df = y_df[y_df['Covered'] == 1]
        if len(y_df) < n0:
            print(f"dropped {n0 - len(y_df)} partially covered (Covered=0) labels.", file=sys.stderr)
    
    # join: X.Start matches Y.Region_Start
    merged_df = pd.merge(
        x_df, 
        y_df[['Read_Name', 'Region_Start', 'Total_Errors']], 
        left_on=['Read_Name', 'Start'], 
        right_on=['Read_Name', 'Region_Start'], 
        how='inner'
    )
    print(f"data joined: {len(merged_df)} valid samples (regions).")
    if len(merged_df) == 0:
        print("error: joined data is empty; check X/Y column names and coordinate alignment.", file=sys.stderr)
        sys.exit(1)
    return merged_df

def main():
    # ================= CLI arguments =================
    parser = argparse.ArgumentParser(description="Train a random forest regressor to predict local sequencing error counts")
    
    # input/output path options
    parser.add_argument("-x", "--x_file", required=True, help="input: X feature TSV file")
    parser.add_argument("-y", "--y_file", required=True, help="input: Y label TSV file")
    parser.add_argument("-m", "--model_out", default="rf_error_model.pkl", help="output: path of the saved model (default: rf_error_model.pkl)")
    parser.add_argument("-f", "--features_out", default="rf_model_features.pkl", help="output: path of the feature list (default: rf_model_features.pkl)")
    parser.add_argument("-p", "--plot_out", default="rf_model_evaluation.png", help="output: path of the evaluation plot (default: rf_model_evaluation.png)")
    
    # split options
    parser.add_argument("--test_size", type=float, default=0.2, help="test set fraction (default: 0.2)")
    parser.add_argument("--seed", type=int, default=42, help="random seed (default: 42)")
    
    # random forest hyperparameters
    parser.add_argument("--n_estimators", type=int, default=100, help="number of trees (default: 100)")
    parser.add_argument("--max_depth", type=int, default=None, help="max tree depth (default: None, unlimited)")
    parser.add_argument("--min_samples_leaf", type=int, default=5, help="min samples per leaf, guards overfitting (default: 5)")
    parser.add_argument("--max_features", type=str, default="1.0", help="fraction of features considered per split (default: 1.0)")
    
    # other options
    parser.add_argument("--no_log", action="store_true", help="disable the log1p transform of y (default: enabled)")
    
    args = parser.parse_args()

    # ================= 1. load data =================
    df = load_and_merge_data(args.x_file, args.y_file)
    
    # feature columns
    feature_cols = [
        'Length', 
        'Read_Zero_Ratio', 
        'Boundary_Cov_Ratio', 
        'Region_Rel_Cov', 
        'Threshold_Rel_Pos'
    ]
    
    missing_features = [f for f in feature_cols if f not in df.columns]
    if missing_features:
        print(f"error: X file is missing feature columns: {missing_features}", file=sys.stderr)
        sys.exit(1)
    
    X = df[feature_cols].fillna(0)
    y = df['Total_Errors'].fillna(0)
    groups = df['Read_Name']

    # ================= 2. data split =================
    print(f"splitting train/test by Read_Name (test size = {args.test_size})...")
    gss = GroupShuffleSplit(n_splits=1, test_size=args.test_size, random_state=args.seed)
    train_idx, test_idx = next(gss.split(X, y, groups=groups))
    
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    
    print(f"train: {X_train.shape[0]} samples, test: {X_test.shape[0]} samples")

    # ================= 3. target transform =================
    if not args.no_log:
        y_train_target = np.log1p(y_train)
        y_test_target = np.log1p(y_test)
        print("applied log1p transform to Y labels.")
    else:
        y_train_target = y_train
        y_test_target = y_test
        print("no log1p transform; predicting raw values.")

    # ================= 4. training =================
    print(f"training random forest (trees={args.n_estimators}, max depth={args.max_depth})...")
    
    # parse max_features (numeric string -> float/int)
    max_feat = args.max_features
    if max_feat.isdigit():
        max_feat = int(max_feat)
    elif max_feat.replace('.', '', 1).isdigit():
        max_feat = float(max_feat)
    
    model = RandomForestRegressor(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        min_samples_leaf=args.min_samples_leaf,
        max_features=max_feat,
        random_state=args.seed,
        n_jobs=-1
    )
    model.fit(X_train, y_train_target)

    # ================= 5. evaluation =================
    print("evaluating model...")
    y_pred_target = model.predict(X_test)
    
    if not args.no_log:
        y_pred = np.expm1(y_pred_target)
    else:
        y_pred = y_pred_target
    
    mae = mean_absolute_error(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    r2 = r2_score(y_test, y_pred)
    
    print("\n=== random forest evaluation ===")
    print(f"MAE:  {mae:.4f}")
    print(f"RMSE: {rmse:.4f}")
    print(f"R2:   {r2:.4f}")
    print("==========================\n")

    # ================= 6. save model =================
    # persist the sklearn model with joblib
    joblib.dump(model, args.model_out)
    joblib.dump(feature_cols, args.features_out)
    print(f"model saved to: {args.model_out}")
    print(f"feature list saved to: {args.features_out}")

    # ================= 7. plot =================
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    axes[0].scatter(y_test, y_pred, alpha=0.5, color='dodgerblue', edgecolor='k', s=20)
    max_val = max(y_test.max(), y_pred.max())
    axes[0].plot([0, max_val], [0, max_val], 'r--', lw=2, label='Perfect Prediction')
    axes[0].set_xlabel('True Total Errors')
    axes[0].set_ylabel('Predicted Total Errors')
    axes[0].set_title(f'True vs Predicted (R² = {r2:.3f})')
    axes[0].legend()

    importances = model.feature_importances_
    indices = np.argsort(importances)[::-1]
    sorted_features = [feature_cols[i] for i in indices]
    
    sns.barplot(x=importances[indices], y=sorted_features, ax=axes[1], palette="viridis")
    axes[1].set_title("Feature Importances (Random Forest)")
    axes[1].set_xlabel("Relative Importance")

    plt.tight_layout()
    plt.savefig(args.plot_out, dpi=300)
    print(f"evaluation plot saved to: {args.plot_out}")

if __name__ == "__main__":
    main()
