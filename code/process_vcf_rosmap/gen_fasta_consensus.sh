#!/bin/bash
#SBATCH --job-name=gen_fasta
#SBATCH --output=/home/jqian54/sulab/enformer_fine_tuning/logs/rosmap_gen_fasta/%A_%a.out
#SBATCH --error=/home/jqian54/sulab/enformer_fine_tuning/logs/rosmap_gen_fasta/%A_%a.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --partition=interactive-cpu
#SBATCH --array=1-400%10
#SBATCH --mail-type=END
#SBATCH --mail-user=jingyun.qian@emory.edu


#2392%30
s
# there is 2392/2 individuals
# wc -l /home/jqian54/sulab/enformer_fine_tuning/data/rosmap_wgs/duplicated_sample_list.txt

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate variformer_env
#set -euo pipefail
#set -x
trap 'echo "Script exited at line $LINENO with status $?" >&2' ERR

BASE_DIR="/home/jqian54/sulab/enformer_fine_tuning"
DATA_DIR="$BASE_DIR/data"

BCF_DIR="$DATA_DIR/rosmap_wgs/bcf"

# TODO: ensure this FASTA matches your BCF build + contig naming (chr21 vs 21)
REF_FASTA="$DATA_DIR/GRCh37.fa"

SAMPLE_LIST="$DATA_DIR/rosmap_wgs/duplicated_sample_list.txt"

OUT_DIR="$DATA_DIR/genomes_rosmap"
mkdir -p "$OUT_DIR"

# ---- sample (MAP id) ----
RAW_SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d '"' | xargs)
SAMPLE="$RAW_SAMPLE"

if [[ -z "$SAMPLE" ]]; then
  echo "Error: No sample found on line $SLURM_ARRAY_TASK_ID of $SAMPLE_LIST" >&2
  exit 1
fi

# ---- haplotype (even/odd) ----
if (( SLURM_ARRAY_TASK_ID % 2 == 0 )); then
  HAPLOTYPE=1
else
  HAPLOTYPE=2
fi

echo "----------------------------------------"
echo "Node: $(hostname)"
echo "Task: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE}"
echo "Haplotype: H${HAPLOTYPE}"
echo "BCF_DIR: ${BCF_DIR}"
echo "REF_FASTA: ${REF_FASTA}"
echo "OUT_DIR: ${OUT_DIR}"
echo "----------------------------------------"

CHROMS=()
for i in $(seq 1 22); do CHROMS+=("chr${i}"); done
CHROMS+=("chrX" "chrY") #"chrothers")

OUTPUT_PREFIX="${SAMPLE}_consensus_H${HAPLOTYPE}"
#TMP_DIR=$(mktemp -d)
TMP_DIR="${OUT_DIR}/temp_job${SLURM_ARRAY_JOB_ID}_task${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMP_DIR"


echo "Generating per-chromosome FASTAs in: $TMP_DIR"

SUCCESS_COUNT=0
MISSING_COUNT=0
FAIL_COUNT=0

for CHR in "${CHROMS[@]}"; do
  # SNP-only: do NOT fall back
  BCF_PATH="${BCF_DIR}/${CHR}.snps.bcf"

  if [[ ! -f "$BCF_PATH" ]]; then
    echo "WARNING: Missing SNP BCF for ${CHR}: ${BCF_PATH}. Skipping." >&2
    #((MISSING_COUNT++))
    MISSING_COUNT=$((MISSING_COUNT + 1))
    continue
  fi

  OUT_CHR_FA="${TMP_DIR}/${OUTPUT_PREFIX}_${CHR}.fa"
  TEMP_CHR_FA="${TMP_DIR}/${OUTPUT_PREFIX}_${CHR}_tmp.fa" # Intermediate file

  # echo "→ ${CHR} using $(basename "$BCF_PATH")"

  # if ! /usr/bin/time -v bcftools consensus \
  #     --fasta-ref "$REF_FASTA" \
  #     --haplotype "$HAPLOTYPE" \
  #     --samples "$SAMPLE" \
  #     "$BCF_PATH" \
  #     -o "$OUT_CHR_FA" #2> "${TMP_DIR}/${OUTPUT_PREFIX}_${CHR}.consensus.stderr"
  # then
  echo "→ ${CHR} using $(basename "$BCF_PATH")"

  # 1. Extract ONLY the current chromosome from the full reference
  REF_CHR_NAME="${CHR#chr}"
  CHR_REF_FA="${TMP_DIR}/${CHR}_ref.fa"
  samtools faidx "$REF_FASTA" "$REF_CHR_NAME" > "$CHR_REF_FA"
  
  # 2. Index the temporary single-chromosome FASTA
  samtools faidx "$CHR_REF_FA"

  # 3. Run consensus using ONLY the single-chromosome reference
  if ! /usr/bin/time -v bcftools consensus \
      --fasta-ref "$CHR_REF_FA" \
      --haplotype "$HAPLOTYPE" \
      --samples "$SAMPLE" \
      "$BCF_PATH" \
      -o "$TEMP_CHR_FA"
  then
    echo "WARNING: consensus FAILED for ${CHR} (sample=$SAMPLE hap=H$HAPLOTYPE). Skipping ${CHR}." >&2
    rm -f "$OUT_CHR_FA"
    #((FAIL_COUNT++))
    FAIL_COUNT=$((FAIL_COUNT + 1))
    continue
  fi

  # CRITICAL: Rename the header from '>1' to '>chr1' so Enformer can read it
  sed "s/^>${REF_CHR_NAME}.*/>${CHR}/" "$TEMP_CHR_FA" > "$OUT_CHR_FA"
  
  # Clean up temporary single-chromosome files to save space
  rm -f "$TEMP_CHR_FA" # "$CHR_REF_FA" "$CHR_REF_FA.fai"

  if [[ ! -s "$OUT_CHR_FA" ]]; then
    echo "WARNING: Output FASTA empty for ${CHR}. Skipping." >&2
    rm -f "$OUT_CHR_FA"
    #((FAIL_COUNT++))
    FAIL_COUNT=$((FAIL_COUNT + 1))
    continue
  fi

  #((SUCCESS_COUNT++))
  SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
done

FINAL_FASTA="${TMP_DIR}/${OUTPUT_PREFIX}.fa"
FINAL_FAI="${FINAL_FASTA}.fai"

echo "----------------------------------------"
echo "Consensus summary:"
echo "  Success chrs: ${SUCCESS_COUNT}"
echo "  Missing chrs: ${MISSING_COUNT}"
echo "  Failed  chrs: ${FAIL_COUNT}"
echo "----------------------------------------"

if (( SUCCESS_COUNT == 0 )); then
  echo "ERROR: No chromosomes produced FASTA output." >&2
  mkdir -p "${OUT_DIR}/consensus_stderr"
  cp -f "${TMP_DIR}/"*.stderr "${OUT_DIR}/consensus_stderr/" 2>/dev/null || true
  # rm -rf "$TMP_DIR"
  exit 1
fi

echo "Concatenating chromosomes into ${FINAL_FASTA} ..."
: > "$FINAL_FASTA"
for CHR in "${CHROMS[@]}"; do
  CHR_FA="${TMP_DIR}/${OUTPUT_PREFIX}_${CHR}.fa"
  [[ -s "$CHR_FA" ]] && cat "$CHR_FA" >> "$FINAL_FASTA"
done

echo "Compressing and indexing final genome..."
# bgzip compresses it down to ~900MB while keeping it searchable
bgzip -f "$FINAL_FASTA" #> "${FINAL_FASTA}.gz"

# Verify the .gz was actually created before we try to index or move things
if [[ -f "${FINAL_FASTA}.gz" ]]; then
    echo "Compression successful. Removing uncompressed .fa if it still exists..."
    rm -f "$FINAL_FASTA"  # This ensures the .fa is gone
else
    echo "ERROR: bgzip failed to create ${FINAL_FASTA}.gz" >&2
    exit 1
fi

echo "Indexing..."
# samtools faidx can index bgzip files perfectly
samtools faidx "${FINAL_FASTA}.gz"

echo "Saving outputs to $OUT_DIR"
# mv -f "$FINAL_FASTA" "$OUT_DIR/"
# mv -f "$FINAL_FAI" "$OUT_DIR/"

mv -f "${FINAL_FASTA}.gz" "$OUT_DIR/"
mv -f "${FINAL_FASTA}.gz.fai" "$OUT_DIR/"
# gzipped fastas also generate a .gzi index
[[ -f "${FINAL_FASTA}.gz.gzi" ]] && mv -f "${FINAL_FASTA}.gz.gzi" "$OUT_DIR/"

# mkdir -p "${OUT_DIR}/consensus_stderr"
# cp -f "${TMP_DIR}/"*.stderr "${OUT_DIR}/consensus_stderr/" 2>/dev/null || true

rm -rf "$TMP_DIR"
echo "DONE: ${OUT_DIR}/$(basename "${FINAL_FASTA}.gz")"
