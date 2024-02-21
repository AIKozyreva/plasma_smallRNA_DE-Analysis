#!/bin/bash

# Define the destination directory
destination_dir="/mnt/projects/users/aalayeva/smallRNA/smRNA_STAR_RESULT/aligned_deduplicated"

# Iterate through each directory
for folder in {2..23}scamt_w_UMI_trimmed; do
    echo "The folder which is in processing now $folder"

    # Use find with -exec to copy files
    find ./$folder -type f \( -name "*deduplicated*Genome*SortByCoord.bam" -o -name "*deduplicated*Transcriptome*SortByCoord.bam" \) -exec cp {} $destination_dir \;

    echo "The $folder result files copying is FINISHED"
    echo -e "\n\n"
done


