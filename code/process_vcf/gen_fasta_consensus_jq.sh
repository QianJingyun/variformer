#!/bin/bash
#SBATCH --job-name=gen_fasta
#SBATCH --output=/home/jqian54/sulab/enformer_fine_tuning/logs/gen_fasta/%A_%a.out
#SBATCH --error=/home/jqian54/sulab/enformer_fine_tuning/logs/gen_fasta/%A_%a.err
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --partition=interactive-cpu
#SBATCH --array=1-1

# there are 1732 lines in the jq_duplicated_sample_list.txt

# --- 0. ENVIRONMENT SETUP ---
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate variformer_env

set -euo pipefail


# --- 1. CONFIGURATION SECTION (MODIFY THESE PATHS) ---
BASE_DIR="/home/jqian54/sulab/enformer_fine_tuning"
DATA_DIR="$BASE_DIR/data"

# Input Files
#BCF_IN="$DATA_DIR/GTEx_v8_SNP_array_VCF/gtex_genotypes_SNPsOnly.bcf.gz"
BCF_IN="$DATA_DIR/GTEx_v8_SNP_array_VCF/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_from_plink_SNPsOnly.bcf.gz"
REF_FASTA="$DATA_DIR/hg38_nochr.fa"
SAMPLE_LIST="$DATA_DIR/jq_duplicated_sample_list.txt"  

# Output Directory
OUT_DIR="$DATA_DIR/genomes"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"


#####################
# "Scratch" is a hard drive attached to a specific computer
# the logic here is a "Copy In --> Compute --> Copy Out" pattern.
# If they ran this directly on global storage (/home/jqian54/...), thousands of jobs trying to read the VCF at the same time would crash the network. By copying it to "Scratch" first, each job reads from its own private hard drive, making it much faster and safer.
#####################
# --- 2. SCRATCH SPACE SETUP ---
# # Use Slurm's local scratch if available, otherwise fallback to /scratch/user
# if [[ -z "$SLURM_TMPDIR" ]]; then
#   MY_TMP="/scratch/$USER/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID"
# else
#   MY_TMP="$SLURM_TMPDIR"
# fi
# mkdir -p "$MY_TMP"
# cd "$MY_TMP"

# echo "Running on node: $(hostname)"
# echo "Processing Array Task: $SLURM_ARRAY_TASK_ID"
# echo "Working in scratch: $MY_TMP"

# --- 3. COPY FILES TO SCRATCH ---
# Copying large files to local scratch speeds up processing and reduces network load

# if [[ -e input.bcf ]]; then
#     echo "ERROR: input.bcf already exists. Refusing to overwrite."
#     exit 1
# fi


# if [ ! -f "input.bcf" ]; then
#     echo "Copying BCF..."
#     #cp "$BCF_IN" "input.bcf"
#     ln -sf "$BCF_IN" "input.bcf"
#     #if [ -f "${BCF_IN}.csi" ]; then cp "${BCF_IN}.csi" "input.bcf.csi"; fi
#     if [ -f "${BCF_IN}.csi" ]; then ln -sf "$BCF_IN.csi" "input.bcf.csi"; fi
# fi

# bcftools view -h input.bcf > /dev/null

# if [ ! -f "ref.fa" ]; then
#     echo "Copying Reference Genome..."
#     cp "$REF_FASTA" "ref.fa"
#     # Copy index if it exists, otherwise samtools faidx will be slow generating it
#     if [ -f "${REF_FASTA}.fai" ]; then cp "${REF_FASTA}.fai" "ref.fa.fai"; fi
# fi

echo "Using BCF directly from:"
echo "  $BCF_IN"

# Fast sanity check (cheap, safe)
bcftools view -h "$BCF_IN" > /dev/null


#######################
# SGE_TASK_ID: This is the job number (e.g., 1, 2, 3, 4...).
# sed -n "Xp": This command reads the file and extracts only line number X.
#######################
# --- 4. IDENTIFY SAMPLE & HAPLOTYPE ---
# Get the sample ID from the text file using the Array Task ID as the line number
#SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
RAW_SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d '"' | xargs)
SAMPLE="${RAW_SAMPLE}_${RAW_SAMPLE}"

if [ -z "$SAMPLE" ]; then
    echo "Error: No sample found on line $SLURM_ARRAY_TASK_ID of $SAMPLE_LIST"
    exit 1
fi

#######################
# Odd Jobs : The script sets haplotype=2.
# Even Jobs : The script sets haplotype=1.

# For this to work, the duplicated_samples_list.txt must have duplicated individual_id
#######################
# Determine Haplotype based on Even/Odd logic
# Odd ID = Haplotype 2, Even ID = Haplotype 1
if [ $((SLURM_ARRAY_TASK_ID % 2)) -eq 0 ]; then
    HAPLOTYPE=1
else
    HAPLOTYPE=2
fi

echo "Target: Sample $SAMPLE | Haplotype $HAPLOTYPE"

#######################
# --- 4.5 GLOBAL REF–FASTA CONSISTENCY CHECK (RUN ONCE) ---
#######################

# # Only let ONE array task do the expensive check
# if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
#     echo "Running bcftools consensus preflight check..."

#     PREFLIGHT_LOG="consensus_mismatch.log"

#     bcftools consensus \
#         --fasta-ref ref.fa \
#         --samples "$(head -n 1 $SAMPLE_LIST)" \
#         input.bcf \
#         -o /dev/null \
#         2> "$PREFLIGHT_LOG"

#     grep -E "does not match the REF allele" "$PREFLIGHT_LOG" > consensus_mismatch_only.log

#     MISMATCH_COUNT=$(wc -l < consensus_mismatch_only.log)

#     echo "----------------------------------------"
#     echo "Consensus-detected mismatches: $MISMATCH_COUNT"
#     echo "----------------------------------------"

#     cp consensus_mismatch_only.log "$OUT_DIR/"
#     echo "$MISMATCH_COUNT" > "$OUT_DIR/consensus_mismatch_count.txt"

#     if [ "$MISMATCH_COUNT" -gt 0 ]; then
#         echo "ERROR: consensus-level REF↔FASTA mismatch detected." >&2
#         echo "Aborting FASTA generation." >&2
#         exit 0
#     fi
# fi


#######################
# bcftools consensus command generates two copies of every chromosome.
# samtools command creates a small index file (.fai) for the new FASTA file
#######################
# --- 5. RUN BCFTOOLS CONSENSUS ---
OUTPUT_FILENAME="${SAMPLE}_consensus_H${HAPLOTYPE}.fa"

echo "Running bcftools consensus..."
# bcftools consensus \
#     --fasta-ref ref.fa \
#     --haplotype "$HAPLOTYPE" \
#     --samples "$SAMPLE" \
#     input.bcf \
#     -o "$OUTPUT_FILENAME"
bcftools consensus \
    --fasta-ref "$REF_FASTA" \
    --haplotype "$HAPLOTYPE" \
    --samples "$SAMPLE" \
    "$BCF_IN" \
    -o "$OUTPUT_FILENAME"


# Check if it worked
if [ $? -ne 0 ]; then
    echo "Error: bcftools consensus failed!"
    exit 1
fi

# Index the new personal genome
echo "Indexing new genome..."
samtools faidx "$OUTPUT_FILENAME"


############### 
# This is to let the job run through

# if a REF↔FASTA mismatch happens:

# record it in the .err

# skip FASTA generation for that task

# continue to the next array task
############### 

# OUTPUT_FILENAME="${SAMPLE}_consensus_H${HAPLOTYPE}.fa"

# echo "Running bcftools consensus for sample=$SAMPLE haplotype=$HAPLOTYPE"

# bcftools consensus \
#     --fasta-ref ref.fa \
#     --haplotype "$HAPLOTYPE" \
#     --samples "$SAMPLE" \
#     input.bcf \
#     -o "$OUTPUT_FILENAME"

# BCF_STATUS=$?

# # --- HANDLE REF MISMATCH OR OTHER FAILURES ---
# if [ $BCF_STATUS -ne 0 ]; then
#     echo "ERROR: bcftools consensus failed for sample=$SAMPLE haplotype=$HAPLOTYPE" >&2
#     echo "Likely REF↔FASTA mismatch. Skipping FASTA generation." >&2
#     echo "----------------------------------------" >&2
#     exit 0   # IMPORTANT: exit cleanly so array continues
# fi

# # Index the new personal genome
# echo "Indexing new genome..."
# samtools faidx "$OUTPUT_FILENAME"

# --- 6. SAVE RESULTS / COPY OUT ---
echo "Moving results to $OUT_DIR..."
mv "$OUTPUT_FILENAME" "$OUT_DIR/"
mv "${OUTPUT_FILENAME}.fai" "$OUT_DIR/"

echo "Job Complete."