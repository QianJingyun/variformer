#!/bin/bash
#SBATCH --job-name=gen_fasta
#SBATCH --output=/home/jqian54/sulab/enformer_fine_tuning/logs/%j.out
#SBATCH --error=/home/jqian54/sulab/enformer_fine_tuning/logs/%j.err
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --partition=interactive-cpu

source ~/.bashrc

conda activate variformer_env


for file in /home/jqian54/sulab/enformer_fine_tuning/data/GTEx_v8_SNP_array_VCF/*.vcf.gz; do
    echo -n "$(basename $file): "
    bcftools query -i 'PR=1' -f '%CHROM\n' "$file" | wc -l
done