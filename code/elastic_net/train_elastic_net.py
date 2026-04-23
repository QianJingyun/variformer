import argparse
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import pearsonr
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline


@dataclass
class GeneResult:
    metrics: Dict[str, object]
    predictions: List[Tuple[str, str, float, float]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train per-gene elastic net models on ROSMAP.")
    parser.add_argument("--data_dir", type=str, required=True)
    parser.add_argument("--gene_file", type=str, required=True)
    parser.add_argument("--fold", type=int, default=1)
    parser.add_argument("--seq_length", type=int, default=49152)
    parser.add_argument("--n_jobs", type=int, default=-1)
    parser.add_argument("--limit", type=int, default=None)
    return parser.parse_args()


def read_lines(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def parse_gene_file(path: str) -> List[str]:
    genes: List[str] = []
    for line in read_lines(path):
        genes.append(line.split()[0])
    return genes


def define_donor_paths(data_dir: str, fold: int) -> Tuple[str, str, str]:
    donor_dir = os.path.join(data_dir, "cross_validation_folds_jq", "cv-fold-usable")
    return (
        os.path.join(donor_dir, f"person_ids-train-fold{fold}.txt"),
        os.path.join(donor_dir, f"person_ids-val-fold{fold}.txt"),
        os.path.join(donor_dir, f"person_ids-test-fold{fold}.txt"),
    )


def build_pipeline() -> Pipeline:
    return Pipeline(
        [
            ("variance_threshold", VarianceThreshold(threshold=0.0)),
            (
                "elastic_net",
                ElasticNetCV(
                    l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0],
                    cv=5,
                    max_iter=2000,
                    random_state=42,
                    n_jobs=1,
                ),
            ),
        ]
    )


def split_indices(donors: Sequence[str], train_donors: set, test_donors: set) -> Tuple[np.ndarray, np.ndarray]:
    train_idx = [i for i, d in enumerate(donors) if d in train_donors]
    test_idx = [i for i, d in enumerate(donors) if d in test_donors]
    return np.array(train_idx, dtype=int), np.array(test_idx, dtype=int)


def impute_from_train_means(X_train: np.ndarray, X_test: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    # Missing genotype calls are represented as NaN. Impute using train-only means.
    col_means = np.nanmean(X_train, axis=0)
    col_means = np.where(np.isnan(col_means), 0.0, col_means).astype(np.float32, copy=False)

    X_train_imp = X_train.copy()
    train_nan = np.isnan(X_train_imp)
    if train_nan.any():
        X_train_imp[train_nan] = np.take(col_means, np.where(train_nan)[1])

    X_test_imp = X_test.copy()
    test_nan = np.isnan(X_test_imp)
    if test_nan.any():
        X_test_imp[test_nan] = np.take(col_means, np.where(test_nan)[1])
    return X_train_imp, X_test_imp


def safe_pearsonr(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[float, float]:
    if y_true.size < 2:
        return np.nan, np.nan
    if np.std(y_true) == 0 or np.std(y_pred) == 0:
        return np.nan, np.nan
    r, p = pearsonr(y_true, y_pred)
    return float(r), float(p)


def process_gene(
    gene: str,
    geno_dir: str,
    expr_df: pd.DataFrame,
    train_donors: set,
    test_donors: set,
) -> GeneResult:
    base_metrics: Dict[str, object] = {
        "gene": gene,
        "n_snvs_raw": np.nan,
        "n_snvs_after_variance": np.nan,
        "n_train_donors": 0,
        "n_test_donors": 0,
        "r2": np.nan,
        "pcc": np.nan,
        "pcc_pval": np.nan,
        "best_alpha": np.nan,
        "best_l1_ratio": np.nan,
        "status": "fit_failed",
    }

    npz_path = os.path.join(geno_dir, f"{gene}.npz")
    if not os.path.exists(npz_path):
        return GeneResult(base_metrics, [])

    try:
        data = np.load(npz_path, allow_pickle=True)
        X = data["X"].astype(np.float32, copy=False)
        donors = [str(x) for x in data["donor_ids"].tolist()]
        base_metrics["n_snvs_raw"] = int(X.shape[1])
    except Exception:
        return GeneResult(base_metrics, [])

    if gene not in expr_df.index:
        return GeneResult(base_metrics, [])

    train_idx, test_idx = split_indices(donors, train_donors, test_donors)
    base_metrics["n_train_donors"] = int(train_idx.size)
    base_metrics["n_test_donors"] = int(test_idx.size)
    if train_idx.size == 0 or test_idx.size == 0:
        return GeneResult(base_metrics, [])

    y_row = expr_df.loc[gene, donors]
    if hasattr(y_row, "values"):
        y_all = y_row.values.astype(np.float32, copy=False)
    else:
        y_all = np.array(y_row, dtype=np.float32)

    X_train = X[train_idx, :]
    X_test = X[test_idx, :]
    y_train = y_all[train_idx]
    y_test = y_all[test_idx]

    if X_train.shape[1] == 0:
        base_metrics["status"] = "no_snvs"
        base_metrics["n_snvs_after_variance"] = 0
        return GeneResult(base_metrics, [])

    X_train_imp, X_test_imp = impute_from_train_means(X_train, X_test)
    model = build_pipeline()

    try:
        model.fit(X_train_imp, y_train)
        vt = model.named_steps["variance_threshold"]
        n_after = int(vt.get_support().sum())
        base_metrics["n_snvs_after_variance"] = n_after
        if n_after == 0:
            base_metrics["status"] = "no_snvs"
            return GeneResult(base_metrics, [])

        enet = model.named_steps["elastic_net"]
        y_pred = model.predict(X_test_imp).astype(np.float32, copy=False)
        r2 = float(r2_score(y_test, y_pred))
        pcc, pval = safe_pearsonr(y_test, y_pred)

        base_metrics.update(
            {
                "r2": r2,
                "pcc": pcc,
                "pcc_pval": pval,
                "best_alpha": float(enet.alpha_),
                "best_l1_ratio": float(enet.l1_ratio_),
                "status": "ok",
            }
        )

        preds = [(gene, donors[i], float(y_test[j]), float(y_pred[j])) for j, i in enumerate(test_idx)]
        return GeneResult(base_metrics, preds)
    except ValueError as exc:
        if "No feature in X meets the variance threshold" in str(exc):
            base_metrics["status"] = "no_snvs"
            base_metrics["n_snvs_after_variance"] = 0
            return GeneResult(base_metrics, [])
        return GeneResult(base_metrics, [])
    except Exception:
        return GeneResult(base_metrics, [])


def main() -> None:
    args = parse_args()
    if args.fold != 1:
        raise ValueError("This pipeline supports Fold-1 only (use --fold 1).")

    genes = parse_gene_file(args.gene_file)
    if args.limit is not None:
        genes = genes[: args.limit]

    expr_path = os.path.join(args.data_dir, "rosmap_wgs", "expression-rosmap.parquet")
    expr_df = pd.read_parquet(expr_path)

    train_path, valid_path, test_path = define_donor_paths(args.data_dir, args.fold)
    train_plus_valid = set(read_lines(train_path)) | set(read_lines(valid_path))
    test_donors = set(read_lines(test_path))

    # Mirror dataset donor selection behavior: keep only donors present in expression matrix columns.
    expr_donor_cols = set(expr_df.columns)
    train_plus_valid = set(sorted(d for d in train_plus_valid if d in expr_donor_cols))
    test_donors = set(sorted(d for d in test_donors if d in expr_donor_cols))

    geno_dir = os.path.join(args.data_dir, "elastic_net_genotypes")
    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)

    print(f"Genes to process: {len(genes)}")
    print(f"Train+valid donors retained: {len(train_plus_valid)}")
    print(f"Test donors retained: {len(test_donors)}")

    outputs = Parallel(n_jobs=args.n_jobs, verbose=10)(
        delayed(process_gene)(gene, geno_dir, expr_df, train_plus_valid, test_donors) for gene in genes
    )

    metrics_rows = [out.metrics for out in outputs]
    prediction_rows: List[Tuple[str, str, float, float]] = []
    for out in outputs:
        prediction_rows.extend(out.predictions)

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_path = os.path.join(results_dir, "fold1_metrics.csv")
    metrics_df.to_csv(metrics_path, index=False)

    preds_df = pd.DataFrame(prediction_rows, columns=["gene", "donor_id", "y_true", "y_pred"])
    preds_path = os.path.join(results_dir, "fold1_predictions.parquet")
    preds_df.to_parquet(preds_path, index=False)

    ok_df = metrics_df[metrics_df["status"] == "ok"]
    median_r2 = float(ok_df["r2"].median()) if len(ok_df) else np.nan
    median_pcc = float(ok_df["pcc"].median()) if len(ok_df) else np.nan
    n_r2_gt_01 = int((ok_df["r2"] > 0.1).sum()) if len(ok_df) else 0
    n_r2_gt_02 = int((ok_df["r2"] > 0.2).sum()) if len(ok_df) else 0

    print("\nTraining summary")
    print(f"- Metrics saved: {metrics_path}")
    print(f"- Predictions saved: {preds_path}")
    print(f"- Total genes: {len(metrics_df)}")
    print(f"- OK genes: {len(ok_df)}")
    print(f"- no_snvs genes: {(metrics_df['status'] == 'no_snvs').sum()}")
    print(f"- fit_failed genes: {(metrics_df['status'] == 'fit_failed').sum()}")
    print(f"- Median R2: {median_r2:.4f}" if not np.isnan(median_r2) else "- Median R2: nan")
    print(f"- Median PCC: {median_pcc:.4f}" if not np.isnan(median_pcc) else "- Median PCC: nan")
    print(f"- # genes with R2 > 0.1: {n_r2_gt_01}")
    print(f"- # genes with R2 > 0.2: {n_r2_gt_02}")


if __name__ == "__main__":
    main()
