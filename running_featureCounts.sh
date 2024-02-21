#!/bin/bash

#Please, run this script in the folder, where you store yout input .bam files

# Set the path to the annotation file
annotation_file="/mnt/projects/users/aalayeva/ref_human_38/gencode.v45.primary_assembly.basic.annotation.gtf"

# Set the path to the result directory
result_dir="/mnt/projects/users/aalayeva/smallRNA/smRNA_STAR_RESULT/aligned_deduplicated/featureCounts_results"

# Loop through each sample from 2 to 24
for i in {2..24}; do
    sample_name="${i}scamt"
    genome_input_file="${sample_name}_w_UMI_trimmed_deduplicated_Aligned.toGenome.SortByCoord.bam"
    transcriptome_input_file="${sample_name}_w_UMI_trimmed_deduplicated_Aligned.toTranscriptome.SortByCoord.bam"

    # Print processing message
    echo "The ${sample_name} files are in processing now"

    # Run featureCounts for genome alignment
    featureCounts -a "${annotation_file}" -o "${result_dir}/${sample_name}_dedup_Align_toGenome_FCounts.txt" "${genome_input_file}"

    # Print finishing message for genome alignment
    echo "Finishing featureCounts for ${sample_name} genome alignment"

    # Run featureCounts for transcriptome alignment
    featureCounts -a "${annotation_file}" -o "${result_dir}/${sample_name}_dedup_Align_toTranscriptome_FCounts.txt" "${transcriptome_input_file}"

    # Print finishing message for transcriptome alignment
    echo "Finishing featureCounts for ${sample_name} transcriptome alignment"
done
