#!/bin/bash

# 1. Define your paths
DATA_DIR="/home/jqian54/sulab/enformer_fine_tuning/data"
SAMPLE_LIST="$DATA_DIR/jq_duplicated_sample_list.txt"
JOB_SCRIPT="gen_fasta_consensus_jq.sh"

# 2. Count the number of lines in your sample list automatically
# 'wc -l' counts lines, 'cut' cleans up the output
NUM_TASKS=$(wc -l < "$SAMPLE_LIST")

# Check if the file is empty
if [ "$NUM_TASKS" -eq 0 ]; then
    echo "Error: Sample list is empty or not found!"
    exit 1
fi

echo "Found $NUM_TASKS samples in list."
echo "Submitting job array for 1-$NUM_TASKS..."

# 3. Submit the job to Slurm with the calculated array size
# The --array flag here overrides whatever is in the script
sbatch --array=1-${NUM_TASKS}%30 "$JOB_SCRIPT"

