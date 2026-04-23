#!/bin/bash

LOG_DIR="/home/jqian54/sulab/enformer_fine_tuning/logs/rosmap_gen_fasta"
OUT_FILE="/home/jqian54/sulab/enformer_fine_tuning/data/cross_validation_folds_jq/completed_sample_list.txt"

# Ensure the output directory exists just in case
mkdir -p "$(dirname "$OUT_FILE")"

# Empty out the file if it already exists so we don't append to old data
> "$OUT_FILE"

echo "Scanning logs for completed samples..."

# Loop through all the 36131603 logs
for log in "${LOG_DIR}"/36131603_*.out; do
    
    # Check if the log contains the success message quietly (-q)
    if grep -q "DONE: /home/jqian54/sulab/enformer_fine_tuning/data/genomes_rosmap/" "$log"; then
        
        # If it does, extract the sample ID
        SAMPLE=$(grep "^Sample: " "$log" | awk '{print $2}')
        
        # Append to the output file
        if [ -n "$SAMPLE" ]; then
            echo "$SAMPLE" >> "$OUT_FILE"
        fi
    fi
done

# Sort the file alphabetically and remove duplicate sample IDs
sort -u "$OUT_FILE" -o "$OUT_FILE"

echo "Done! Extracted sample IDs saved to: $OUT_FILE"
echo "Total unique completed samples: $(wc -l < "$OUT_FILE")"