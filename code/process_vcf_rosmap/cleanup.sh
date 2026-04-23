#!/bin/bash

# Target data directory
DATA_DIR="/home/jqian54/sulab/enformer_fine_tuning/data/genomes_rosmap"
LOG_DIR="/home/jqian54/sulab/enformer_fine_tuning/logs/rosmap_gen_fasta"

echo "Starting cleanup of unfinished jobs (Task 286+)..."

# Look for logs explicitly inside LOG_DIR
for log in "${LOG_DIR}"/36131603_*.out; do
    # Extract the task ID from the log filename
    task_id=$(echo "$log" | sed 's/.*_\(.*\)\.out/\1/')
    
    # Check if the task_id is 286 or greater
    if [ "$task_id" -ge 286 ] 2>/dev/null; then
        
        # Grab the Sample and Haplotype from inside the log
        SAMPLE=$(grep "^Sample: " "$log" | awk '{print $2}')
        HAP=$(grep "^Haplotype: " "$log" | awk '{print $2}')
        
        # If both values were successfully found, print the target files
        if [ -n "$SAMPLE" ] && [ -n "$HAP" ]; then
            echo "Task $task_id unfinished. Targeting files for: ${SAMPLE}_consensus_${HAP}"
            echo " -> rm -f ${DATA_DIR}/${SAMPLE}_consensus_${HAP}.fa.gz*"
        fi
    fi
done