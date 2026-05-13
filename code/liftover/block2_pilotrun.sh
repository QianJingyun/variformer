#!/bin/bash
# block2_liftover_pilot.sh
# Pilot run: lift chr21 from b37 to hg38 using Picard LiftoverVcf.
# Run on an interactive node. ~5–15 minutes.

set -euo pipefail

# --- Paths -------------------------------------------------------------------
DATA_DIR="/home/jqian54/sulab/enformer_fine_tuning/data"
LIFTOVER_DIR="${DATA_DIR}/liftover"
RAW_VCF_DIR="${DATA_DIR}/rosmap_wgs/raw_vcf"
OUT_DIR="${DATA_DIR}/rosmap_wgs/lifted_vcf"

CHROM="21"
RAW_VCF="${RAW_VCF_DIR}/DEJ_11898_B01_GRM_WGS_2017-05-15_${CHROM}.recalibrated_variants.vcf.gz"
RENAMED_VCF="${OUT_DIR}/chr${CHROM}.renamed.vcf.gz"
LIFTED_VCF="${OUT_DIR}/chr${CHROM}.hg38.vcf.gz"
REJECT_VCF="${OUT_DIR}/chr${CHROM}.hg38.rejected.vcf.gz"

CHAIN_FILE="${LIFTOVER_DIR}/hg19ToHg38.over.chain.gz"
TARGET_FA="${DATA_DIR}/hg38_genome.fa"
RENAME_MAP="${LIFTOVER_DIR}/b37_to_hg19_contigs.txt"

mkdir -p "${OUT_DIR}"

# --- 1. Rename contigs: 21 -> chr21 -----------------------------------------
# bcftools annotate --rename-chrs keeps only contigs found in the map (any others
# stay in the file with original names). We then restrict to chr21 explicitly with
# `bcftools view -r` after the rename, so anything not in the map gets dropped.
if [[ ! -f "${RENAMED_VCF}" ]]; then
    echo "[1/3] Renaming contigs in chr${CHROM}..."
    bcftools annotate \
        --rename-chrs "${RENAME_MAP}" \
        -Oz -o "${RENAMED_VCF}.tmp" \
        "${RAW_VCF}"
    # Restrict to chr21 only (drops any rogue MT/decoys silently passed through)
    bcftools view -r "chr${CHROM}" -Oz -o "${RENAMED_VCF}" "${RENAMED_VCF}.tmp"
    bcftools index -t "${RENAMED_VCF}"
    rm "${RENAMED_VCF}.tmp"
else
    echo "[1/3] Renamed VCF already present: ${RENAMED_VCF}"
fi

# Variant count before liftover
N_INPUT=$(bcftools view -H "${RENAMED_VCF}" | wc -l)
echo "    Input variants on chr${CHROM}: ${N_INPUT}"

# --- 2. LiftoverVcf ---------------------------------------------------------
# WARN_ON_MISSING_CONTIG=true: if the chain references a contig not in our target
# FASTA (e.g., some alt scaffold), warn instead of crash. Safe because our rename
# map only has primary chromosomes anyway.
#
# RECOVER_SWAPPED_REF_ALT=true: if the lifted REF allele matches what was originally
# ALT (i.e., the strand flipped), Picard will swap REF/ALT and update genotypes
# rather than rejecting the variant. Standard practice for population genetics.
#
# Heap: 16G is comfortable for one autosome. Picard's memory grows with the FASTA
# index, not the VCF size, so this is constant across chromosomes.
if [[ ! -f "${LIFTED_VCF}" ]]; then
    echo "[2/3] Running LiftoverVcf on chr${CHROM}..."
    picard -Xmx16g LiftoverVcf \
        I="${RENAMED_VCF}" \
        O="${LIFTED_VCF}" \
        CHAIN="${CHAIN_FILE}" \
        REJECT="${REJECT_VCF}" \
        R="${TARGET_FA}" \
        WARN_ON_MISSING_CONTIG=true \
        RECOVER_SWAPPED_REF_ALT=true \
        CREATE_INDEX=false
    bcftools index -t "${LIFTED_VCF}"
else
    echo "[2/3] Lifted VCF already present: ${LIFTED_VCF}"
fi

# --- 3. Sanity checks --------------------------------------------------------
echo "[3/3] Sanity checks..."
N_LIFTED=$(bcftools view -H "${LIFTED_VCF}" | wc -l)
N_REJECT=$(bcftools view -H "${REJECT_VCF}" | wc -l)
PCT=$(awk "BEGIN {printf \"%.2f\", 100*${N_LIFTED}/${N_INPUT}}")

echo
echo "  Input variants    : ${N_INPUT}"
echo "  Lifted (success)  : ${N_LIFTED}  (${PCT}%)"
echo "  Rejected          : ${N_REJECT}"
echo
echo "  Lifted contigs (should be chr21 only, possibly with a few stray alt mappings):"
bcftools view -H "${LIFTED_VCF}" | cut -f1 | sort -u
echo
echo "  Top reject reasons:"
bcftools view -H "${REJECT_VCF}" | cut -f7 | sort | uniq -c | sort -rn | head -5
echo
echo "  First 3 lifted variants (chr/pos/ref/alt):"
bcftools view -H "${LIFTED_VCF}" | head -3 | cut -f1-5
echo
echo "Pilot done. If lift rate is >97%, block 2 is a success — proceed to block 3."
