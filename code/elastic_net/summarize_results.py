import argparse
import os

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize elastic net fold-1 metrics.")
    parser.add_argument(
        "--metrics_csv",
        type=str,
        default="/sulab/users/jqian54/enformer_fine_tuning/code/elastic_net/results/fold1_metrics.csv",
        help="Path to fold1_metrics.csv",
    )
    parser.add_argument("--top_k", type=int, default=10, help="Show top-k genes by R2.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not os.path.exists(args.metrics_csv):
        raise FileNotFoundError(f"Metrics file not found: {args.metrics_csv}")

    df = pd.read_csv(args.metrics_csv)
    ok = df[df["status"] == "ok"].copy()

    print("Elastic Net Results Summary")
    print(f"- metrics file: {args.metrics_csv}")
    print(f"- total genes: {len(df)}")
    print(f"- ok genes: {len(ok)}")
    print(f"- no_snvs genes: {(df['status'] == 'no_snvs').sum()}")
    print(f"- fit_failed genes: {(df['status'] == 'fit_failed').sum()}")

    if len(ok) == 0:
        print("- no successful genes to summarize")
        return

    median_r2 = ok["r2"].median()
    median_pcc = ok["pcc"].median()
    print(f"- median R2: {median_r2:.4f}" if not np.isnan(median_r2) else "- median R2: nan")
    print(f"- median PCC: {median_pcc:.4f}" if not np.isnan(median_pcc) else "- median PCC: nan")
    print(f"- # genes with R2 > 0.1: {(ok['r2'] > 0.1).sum()}")
    print(f"- # genes with R2 > 0.2: {(ok['r2'] > 0.2).sum()}")

    best_idx = ok["r2"].idxmax()
    best = ok.loc[best_idx]
    print(
        f"- best R2 gene: {best['gene']} (R2={best['r2']:.4f}, PCC={best['pcc'] if pd.notna(best['pcc']) else 'nan'})"
    )

    top_k = max(1, args.top_k)
    top = ok.sort_values("r2", ascending=False).head(top_k)[["gene", "r2", "pcc", "n_snvs_after_variance"]]
    print(f"\nTop {top_k} genes by R2")
    print(top.to_string(index=False))


if __name__ == "__main__":
    main()
