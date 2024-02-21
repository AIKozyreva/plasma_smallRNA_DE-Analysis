#!/bin/bash

# List of directories
directories=("2scamt_w_UMI_trimmed" "3scamt_w_UMI_trimmed" "4scamt_w_UMI_trimmed" "5scamt_w_UMI_trimmed" "6scamt_w_UMI_trimmed" "10scamt_w_UMI_trimmed" "11scamt_w_UMI_trimmed" "12scamt_w_UMI_trimmed" "13scamt_w_UMI_trimmed" "14scamt_w_UMI_trimmed" "15scamt_w_UMI_trimmed" "16scamt_w_UMI_trimmed" "17scamt_w_UMI_trimmed" "18scamt_w_UMI_trimmed" "19scamt_w_UMI_trimmed" "20scamt_w_UMI_trimmed" "21scamt_w_UMI_trimmed" "23scamt_w_UMI_trimmed" "24scamt_w_UMI_trimmed")

# Loop through each directory
for foldername in "${directories[@]}"
do
    echo "The folder which is in processing now $foldername"

    # Step 2
    samtools sort -@ 10 -o "${foldername}/${foldername}_Aligned.toTranscriptome.sorted.out.bam" "${foldername}/${foldername}Aligned.toTranscriptome.out.bam"

    # Step 3
    echo "Finishing sorting original transcript file and starting demultiplexing it"

    # Step 4
    fumi_tools dedup --threads 10 --memory 10G -i "${foldername}/${foldername}_Aligned.toTranscriptome.sorted.out.bam" -o "${foldername}/${foldername}_deduplicated_Aligned.toTranscriptome.bam"

    # Step 5
    echo "Finishing demultiplexing transcriptome, starting genome demultiplexing"

    # Step 6
    fumi_tools dedup --threads 10 --memory 10G -i "${foldername}/${foldername}Aligned.sortedByCoord.out.bam" -o "${foldername}/${foldername}_deduplicated_Aligned.toGenome.bam"

    # Step 7
    echo "Start sorting demultiplexed result files"

    # Step 8
    samtools sort -@ 10 -o "${foldername}/${foldername}_deduplicated_Aligned.toGenome.SortByCoord.bam" "${foldername}/${foldername}_deduplicated_Aligned.toGenome.bam"

    # Step 9
    echo "Genome output is sorted"

    # Step 10
    samtools sort -@ 10 -o "${foldername}/${foldername}_deduplicated_Aligned.toTranscriptome.SortByCoord.bam" "${foldername}/${foldername}_deduplicated_Aligned.toTranscriptome.bam"

    # Step 11
    echo "Transcriptom output is sorted"
done
