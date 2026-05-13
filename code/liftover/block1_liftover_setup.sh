#!/bin/bash
# block1_liftover_setup.sh
# Prepare chain file, target reference index, sequence dictionary, and contig-rename map
# for Picard LiftoverVcf (b37 → hg38 via UCSC hg19ToHg38 chain).
# Run this ONCE on an interactive node before any liftover jobs.

set -euo pipefail

# --- Paths -------------------------------------------------------------------
DATA_DIR="/home/jqian54/sulab/enformer_fine_tuning/data"
LIFTOVER_DIR="${DATA_DIR}/liftover"
RAW_VCF_DIR="${DATA_DIR}/rosmap_wgs/raw_vcf"
TARGET_FA="${DATA_DIR}/hg38_genome.fa"

CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_FILE="${LIFTOVER_DIR}/hg19ToHg38.over.chain.gz"
RENAME_MAP="${LIFTOVER_DIR}/b37_to_hg19_contigs.txt"

mkdir -p "${LIFTOVER_DIR}"

# --- 1. Download chain file --------------------------------------------------
if [[ ! -f "${CHAIN_FILE}" ]]; then
    echo "[1/5] Downloading hg19ToHg38 chain file..."
    wget -q -O "${CHAIN_FILE}" "${CHAIN_URL}"
    # Sanity check: should be ~220 KB, not an HTML error page
    SIZE=$(stat -c%s "${CHAIN_FILE}")
    if (( SIZE < 100000 )); then
        echo "ERROR: chain file is only ${SIZE} bytes — likely a failed download."
        exit 1
    fi
else
    echo "[1/5] Chain file already present: ${CHAIN_FILE}"
fi

# --- 2. Index target FASTA (.fai) --------------------------------------------
if [[ ! -f "${TARGET_FA}.fai" ]]; then
    echo "[2/5] Building .fai index for ${TARGET_FA}..."
    samtools faidx "${TARGET_FA}"
else
    echo "[2/5] .fai index already present."
fi

# --- 3. Build sequence dictionary (.dict) ------------------------------------
TARGET_DICT="${TARGET_FA%.fa}.dict"
if [[ ! -f "${TARGET_DICT}" ]]; then
    echo "[3/5] Building sequence dictionary..."
    picard CreateSequenceDictionary \
        R="${TARGET_FA}" \
        O="${TARGET_DICT}"
else
    echo "[3/5] Sequence dictionary already present: ${TARGET_DICT}"
fi

# --- 4. Write contig rename map (b37 → hg19/UCSC) ---------------------------
# Used by `bcftools annotate --rename-chrs` to convert "1" → "chr1", etc.
# We drop MT entirely (see block 2) since chain handling of mitochondria is messy
# and we don't need it for cis-eQTL analysis.
if [[ ! -f "${RENAME_MAP}" ]]; then
    echo "[4/5] Writing contig rename map..."
    {
        for i in {1..22}; do echo -e "${i}\tchr${i}"; done
        echo -e "X\tchrX"
        echo -e "Y\tchrY"
    } > "${RENAME_MAP}"
else
    echo "[4/5] Rename map already present: ${RENAME_MAP}"
fi

# --- 5. Sanity-check contig naming alignment --------------------------------
echo "[5/5] Sanity-checking contig naming..."
echo
echo "Input VCF contigs (from header of first VCF found):"
SAMPLE_VCF=$(ls "${RAW_VCF_DIR}"/*.vcf.gz "${RAW_VCF_DIR}"/*.vcf 2>/dev/null | head -1)
echo "  (using: ${SAMPLE_VCF})"
bcftools view -h "${SAMPLE_VCF}" | grep "^##contig" | head -5
echo
echo "Chain file source contigs (first 5):"
zcat "${CHAIN_FILE}" | grep "^chain" | awk '{print $3}' | sort -u | head -5
echo
echo "Chain file target contigs (first 5):"
zcat "${CHAIN_FILE}" | grep "^chain" | awk '{print $8}' | sort -u | head -5
echo
echo "Target FASTA contigs (first 5):"
grep "^>" "${TARGET_FA}" | head -5
echo
echo "Rename map:"
head -3 "${RENAME_MAP}"
echo "  ..."
echo
echo "Expected alignment after rename step in block 2:"
echo "  Renamed VCF       chr1, chr2, ...   ==   Chain SOURCE   chr1, chr2, ..."
echo "  Chain TARGET      chr1, chr2, ...   ==   FASTA contigs  chr1, chr2, ..."
echo
echo "Setup complete."
