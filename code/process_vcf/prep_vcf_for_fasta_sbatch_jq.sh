#!/bin/bash
#SBATCH --job-name=prep_vcf
#SBATCH --output=/home/jqian54/sulab/enformer_fine_tuning/logs/%j.out
#SBATCH --error=/home/jqian54/sulab/enformer_fine_tuning/logs/%j.err
#SBATCH --mem=50G
#SBATCH --time=7-00:00:00
#SBATCH --partition=su_lab


source ~/.bashrc

conda activate variformer_env

##############
# This takes the directory path
##############
# 1. Set your specific directory
vcf_dir="/home/jqian54/sulab/enformer_fine_tuning/data/GTEx_v8_SNP_array_VCF"

# 2. Define your specific Input and Output names
# infile="$vcf_dir/gtex_genotypes.vcf.gz"
# outfile_bcf="$vcf_dir/gtex_genotypes.bcf.gz"
# outfile_snps_bcf="$vcf_dir/gtex_genotypes_SNPsOnly.bcf.gz"
# outfile_snps_vcf="$vcf_dir/gtex_genotypes_SNPsOnly.vcf.gz"

infile="$vcf_dir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_from_plink.vcf.gz"
outfile_bcf="$vcf_dir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_from_plink.bcf"
outfile_snps_bcf="$vcf_dir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_from_plink_SNPsOnly.bcf"
outfile_snps_vcf="$vcf_dir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_from_plink_SNPsOnly.vcf.gz"

# --- Step 1: Create binary BCF from your VCF (Faster for processing) ---
##############
# Computers read Binary BCF much faster than text VCF.
##############
if [ ! -f "$outfile_bcf" ]; then
    echo "Creating large unphased BCF from $infile..."
    bcftools view --output-type b --output "$outfile_bcf" "$infile"
    echo "Indexing unphased BCF..."
    bcftools index "$outfile_bcf"
fi

# --- Step 2: Subset SNPs Only (Output as BCF) ---
##############
# --types snps tells the tool to throw away insertions and deletions (indels).
# It only keeps SNPs so that the length of DNA seq do not change.
##############
echo "Subsetting SNPs to BCF..."
bcftools view -m2 -M2 --types snps --output-type b --output "$outfile_snps_bcf" "$outfile_bcf"
echo "Indexing SNP-only BCF..."
bcftools index "$outfile_snps_bcf"

# --- Step 3: Subset SNPs Only (Output as VCF) ---
echo "Subsetting SNPs to VCF..."
bcftools view -m2 -M2 --types snps --output-type z --output "$outfile_snps_vcf" "$outfile_bcf"
echo "Indexing SNP-only VCF..."
bcftools index "$outfile_snps_vcf"

echo "Done! Output saved to: $outfile_snps_vcf"