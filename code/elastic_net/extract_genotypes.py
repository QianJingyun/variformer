import argparse
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from cyvcf2 import VCF as CyVCF2VCF  # type: ignore
except Exception:
    CyVCF2VCF = None

try:
    import pysam  # type: ignore
except Exception:
    pysam = None


@dataclass
class GeneWindow:
    gene_name: str
    chrom: str
    tss: int
    window_start: int
    window_end: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract per-gene ROSMAP genotype matrices.")
    parser.add_argument("--data_dir", type=str, required=True)
    parser.add_argument("--gene_file", type=str, required=True)
    parser.add_argument("--fold", type=int, default=1)
    parser.add_argument("--seq_length", type=int, default=49152)
    parser.add_argument("--overwrite", action="store_true")
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


def compute_gene_windows(data_dir: str, seq_length: int) -> Dict[str, GeneWindow]:
    regions_path = os.path.join(
        data_dir,
        "rosmap_wgs",
        "Enformer_genomic_regions_TSSCenteredGenes_FixedOverlapRemoval.csv",
    )
    regions = pd.read_csv(regions_path)
    windows: Dict[str, GeneWindow] = {}
    for _, row in regions.iterrows():
        gene_name = str(row["gene_name"])
        chrom = str(row["seqnames"])
        region_start = int(row["starts"])
        region_end = int(row["ends"])
        tss = region_end - (196608 // 2)
        if seq_length != 196608:
            region_start = tss - (seq_length // 2)
            region_end = tss + (seq_length // 2)
        windows[gene_name] = GeneWindow(
            gene_name=gene_name,
            chrom=chrom,
            tss=tss,
            window_start=region_start,
            window_end=region_end,
        )
    return windows


def load_expression_donors(data_dir: str) -> List[str]:
    expr_path = os.path.join(data_dir, "rosmap_wgs", "expression-rosmap.parquet")
    expr_df = pd.read_parquet(expr_path)
    return list(expr_df.columns)


def choose_vcf_backend():
    if CyVCF2VCF is not None:
        return "cyvcf2"
    if pysam is not None:
        return "pysam"
    raise RuntimeError("Neither cyvcf2 nor pysam is available.")


def resolve_chr_vcf_path(vcf_dir: str, chrom: str) -> str:
    suffix = chrom.replace("chr", "")
    path = os.path.join(vcf_dir, f"chr{suffix}.snps.vcf.gz")
    if not os.path.exists(path):
        raise FileNotFoundError(f"VCF not found for {chrom}: {path}")
    return path


def _genotype_to_dosage(gt: Sequence[int]) -> float:
    if len(gt) < 2:
        return np.nan
    a0, a1 = gt[0], gt[1]
    if a0 is None or a1 is None:
        return np.nan
    if a0 < 0 or a1 < 0:
        return np.nan
    if a0 == 0 and a1 == 0:
        return 0.0
    if (a0 == 0 and a1 == 1) or (a0 == 1 and a1 == 0):
        return 1.0
    if a0 == 1 and a1 == 1:
        return 2.0
    return np.nan


class CyVCF2Extractor:
    def __init__(self, vcf_path: str):
        self.vcf = CyVCF2VCF(vcf_path)
        self.samples = list(self.vcf.samples)
        self.contigs = set(self.vcf.seqnames)

    def _choose_query_chrom(self, chrom: str) -> str:
        if chrom in self.contigs:
            return chrom
        if chrom.startswith("chr") and chrom[3:] in self.contigs:
            return chrom[3:]
        alt = f"chr{chrom}"
        if alt in self.contigs:
            return alt
        return chrom

    def fetch(
        self, chrom: str, start_0_based: int, end_0_based: int
    ) -> Tuple[np.ndarray, List[Tuple[str, int, str, str]], int, int]:
        x_cols: List[np.ndarray] = []
        positions: List[Tuple[str, int, str, str]] = []
        skipped_non_snv = 0
        skipped_multiallelic = 0

        query_chrom = self._choose_query_chrom(chrom)
        region = f"{query_chrom}:{start_0_based + 1}-{end_0_based}"
        for var in self.vcf(region):
            if len(var.ALT) != 1:
                skipped_multiallelic += 1
                continue
            ref = str(var.REF)
            alt = str(var.ALT[0])
            if len(ref) != 1 or len(alt) != 1:
                skipped_non_snv += 1
                continue
            if not (start_0_based < int(var.POS) <= end_0_based):
                continue

            col = np.array([_genotype_to_dosage(g) for g in var.genotypes], dtype=np.float32)
            x_cols.append(col)
            positions.append((chrom, int(var.POS), ref, alt))

        if x_cols:
            X = np.stack(x_cols, axis=1).astype(np.float32, copy=False)
        else:
            X = np.zeros((len(self.samples), 0), dtype=np.float32)
        return X, positions, skipped_non_snv, skipped_multiallelic

    def close(self) -> None:
        return


class PysamExtractor:
    def __init__(self, vcf_path: str):
        self.vcf = pysam.VariantFile(vcf_path)
        self.samples = list(self.vcf.header.samples)
        self.contigs = set(self.vcf.header.contigs.keys())

    def _choose_query_chrom(self, chrom: str) -> str:
        if chrom in self.contigs:
            return chrom
        if chrom.startswith("chr") and chrom[3:] in self.contigs:
            return chrom[3:]
        alt = f"chr{chrom}"
        if alt in self.contigs:
            return alt
        return chrom

    def fetch(
        self, chrom: str, start_0_based: int, end_0_based: int
    ) -> Tuple[np.ndarray, List[Tuple[str, int, str, str]], int, int]:
        x_cols: List[np.ndarray] = []
        positions: List[Tuple[str, int, str, str]] = []
        skipped_non_snv = 0
        skipped_multiallelic = 0
        sample_names = self.samples

        query_chrom = self._choose_query_chrom(chrom)
        for rec in self.vcf.fetch(query_chrom, start_0_based, end_0_based):
            alts = rec.alts or ()
            if len(alts) != 1:
                skipped_multiallelic += 1
                continue
            ref = str(rec.ref)
            alt = str(alts[0])
            if len(ref) != 1 or len(alt) != 1:
                skipped_non_snv += 1
                continue

            dosages: List[float] = []
            for sample in sample_names:
                gt = rec.samples[sample].get("GT", None)
                if gt is None:
                    dosages.append(np.nan)
                    continue
                dosages.append(_genotype_to_dosage(gt))
            x_cols.append(np.array(dosages, dtype=np.float32))
            positions.append((chrom, int(rec.pos), ref, alt))

        if x_cols:
            X = np.stack(x_cols, axis=1).astype(np.float32, copy=False)
        else:
            X = np.zeros((len(sample_names), 0), dtype=np.float32)
        return X, positions, skipped_non_snv, skipped_multiallelic

    def close(self) -> None:
        self.vcf.close()


def build_extractor(backend: str, vcf_path: str):
    if backend == "cyvcf2":
        return CyVCF2Extractor(vcf_path)
    return PysamExtractor(vcf_path)


def subset_rows_by_donors(
    X: np.ndarray, sample_ids: Sequence[str], target_donors: Sequence[str]
) -> Tuple[np.ndarray, List[str]]:
    sample_to_idx = {s: i for i, s in enumerate(sample_ids)}
    selected = [d for d in target_donors if d in sample_to_idx]
    indices = [sample_to_idx[d] for d in selected]
    if len(indices) == 0:
        return np.zeros((0, X.shape[1]), dtype=np.float32), []
    return X[indices, :], selected


def main() -> None:
    args = parse_args()
    if args.fold != 1:
        raise ValueError("This pipeline supports Fold-1 only (use --fold 1).")

    windows = compute_gene_windows(args.data_dir, args.seq_length)
    gene_list = parse_gene_file(args.gene_file)
    if args.limit is not None:
        gene_list = gene_list[: args.limit]

    expr_donors = set(load_expression_donors(args.data_dir))
    train_path, valid_path, test_path = define_donor_paths(args.data_dir, args.fold)
    split_donors = read_lines(train_path) + read_lines(valid_path) + read_lines(test_path)
    split_donors_unique = sorted(set(split_donors))
    selected_donors = [d for d in split_donors_unique if d in expr_donors]

    print(f"Requested genes: {len(gene_list)}")
    print(f"Donors in fold files (union): {len(split_donors_unique)}")
    print(f"Donors retained (also in expression columns): {len(selected_donors)}")

    out_dir = os.path.join(args.data_dir, "elastic_net_genotypes")
    os.makedirs(out_dir, exist_ok=True)
    vcf_dir = os.path.join(args.data_dir, "rosmap_wgs", "snps_only")
    backend = choose_vcf_backend()
    print(f"VCF backend: {backend}")

    extractors: Dict[str, object] = {}
    total_snvs = 0
    fail_genes: List[str] = []
    skipped_existing = 0
    skipped_non_snv_total = 0
    skipped_multiallelic_total = 0
    processed = 0

    try:
        for i, gene in enumerate(gene_list, start=1):
            if i % 25 == 0:
                print(f"Progress: {i}/{len(gene_list)} genes")

            out_path = os.path.join(out_dir, f"{gene}.npz")
            if os.path.exists(out_path) and not args.overwrite:
                skipped_existing += 1
                continue

            if gene not in windows:
                fail_genes.append(gene)
                continue

            gw = windows[gene]
            try:
                vcf_path = resolve_chr_vcf_path(vcf_dir, gw.chrom)
                if vcf_path not in extractors:
                    extractors[vcf_path] = build_extractor(backend, vcf_path)
                extractor = extractors[vcf_path]

                X_full, positions, skipped_non_snv, skipped_multiallelic = extractor.fetch(
                    gw.chrom, gw.window_start, gw.window_end
                )
                sample_ids = extractor.samples
                X, donors = subset_rows_by_donors(X_full, sample_ids, selected_donors)

                np.savez_compressed(
                    out_path,
                    X=X.astype(np.float32, copy=False),
                    donor_ids=np.array(donors, dtype=object),
                    snv_positions=np.array(positions, dtype=object),
                    gene_name=gw.gene_name,
                    chrom=gw.chrom,
                    tss=int(gw.tss),
                    window_start=int(gw.window_start),
                    window_end=int(gw.window_end),
                )
                total_snvs += len(positions)
                skipped_non_snv_total += skipped_non_snv
                skipped_multiallelic_total += skipped_multiallelic
                processed += 1
            except Exception as exc:
                fail_genes.append(gene)
                print(f"[WARN] Failed gene {gene}: {exc}")
    finally:
        for ext in extractors.values():
            ext.close()

    mean_snvs = (total_snvs / processed) if processed > 0 else float("nan")
    print("\nExtraction summary")
    print(f"- Total genes requested: {len(gene_list)}")
    print(f"- Genes processed: {processed}")
    print(f"- Genes skipped (existing): {skipped_existing}")
    print(f"- Mean SNVs per processed gene: {mean_snvs:.2f}" if processed else "- Mean SNVs per processed gene: n/a")
    print(f"- Skipped non-SNV variants: {skipped_non_snv_total}")
    print(f"- Skipped multiallelic variants: {skipped_multiallelic_total}")
    print(f"- Genes failed: {len(fail_genes)}")
    if fail_genes:
        print(f"- Failed gene list: {','.join(fail_genes)}")


if __name__ == "__main__":
    main()
