#!/bin/bash
#SBATCH --job-name=rosmap_snps
#SBATCH --output=/home/jqian54/sulab/enformer_fine_tuning/logs/rosmap_prep/rosmap_%A_%a.out
#SBATCH --error=/home/jqian54/sulab/enformer_fine_tuning/logs/rosmap_prep/rosmap_%A_%a.err
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=encore
#SBATCH --array=1-25

source ~/.bashrc
conda activate variformer_env

BASE_DIR=/home/jqian54/sulab/enformer_fine_tuning/data/rosmap_wgs
RAW_DIR=$BASE_DIR/raw_vcf
BCF_DIR=$BASE_DIR/bcf
OUT_DIR=$BASE_DIR/snps_only

mkdir -p "$BCF_DIR" "$OUT_DIR"

# ---- chromosome mapping ----
CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y others)
CHR=${CHRS[$SLURM_ARRAY_TASK_ID-1]}

infile=$(ls $RAW_DIR/*_${CHR}.recalibrated_variants.vcf.gz)

outfile_bcf=$BCF_DIR/chr${CHR}.bcf
outfile_snps_bcf=$BCF_DIR/chr${CHR}.snps.bcf
outfile_snps_vcf=$OUT_DIR/chr${CHR}.snps.vcf.gz

echo "Processing chromosome: $CHR"
echo "Input VCF: $infile"

# --- Step 1: Create binary BCF from chromosome VCF ---
if [ ! -f "$outfile_bcf" ]; then
    echo "Creating BCF for chr${CHR}"
    bcftools view \
        --output-type b \
        --output "$outfile_bcf" \
        "$infile"
    bcftools index "$outfile_bcf"
fi

# --- Step 2: Subset SNPs Only (Output as BCF) ---
if [ ! -f "$outfile_snps_bcf" ]; then
    echo "Subsetting SNPs to BCF for chr${CHR}"
    bcftools view \
        -m2 -M2 \
        --types snps \
        --output-type b \
        --output "$outfile_snps_bcf" \
        "$outfile_bcf"
    bcftools index "$outfile_snps_bcf"
fi

# --- Step 3: Subset SNPs Only (Output as VCF) ---
echo "Subsetting SNPs to VCF for chr${CHR}"
bcftools view \
    -m2 -M2 \
    --types snps \
    --output-type z \
    --output "$outfile_snps_vcf" \
    "$outfile_bcf"

bcftools index "$outfile_snps_vcf"

echo "Done chr${CHR}: $outfile_snps_vcf"
