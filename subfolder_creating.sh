#!/bin/bash
#this programm will for each sample create a subfolder with sample-like name with results of alignment by STAR

# input_dir-were you store files for processing
input_dir="/mnt/projects/users/aalayeva/smallRNA/trimmed"

#were you want to have subfolders with alignment output
output_base_dir="/mnt/projects/users/aalayeva/smallRNA/smRNA_STAR_RESULT"

#were star-indexed and other ref files are stored
genome_dir="~/ref_human_38/star_index_h38"
gtf_file="~/ref_human_38/gencode.v45.primary_assembly.basic.annotation.gtf"

#set sjdb_overhang parameter as *length_read*-1 
sjdb_overhang=49

#give to the program bunch of input-file names
files_to_process=("7scamt_w_UMI_trimmed.fastq" "8scamt_w_UMI_trimmed.fastq" "9scamt_w_UMI_trimmed.fastq" "scamt1_w_UMI_trimmed.fastq" "22scamt_w_UMI_trimmed.fastq")

for file in "${files_to_process[@]}"; do

    filename=$(basename "$file" .fastq)
    output_subdir="$output_base_dir/$filename"
    mkdir -p "$output_subdir"
    
    # STARcomand
    STAR --runThreadN 10 \
     --genomeDir "$genome_dir" \
     --sjdbGTFfile "$gtf_file" \
     --sjdbOverhang "$sjdb_overhang" \
     --readFilesIn "$input_dir/$file" \
     --quantMode TranscriptomeSAM \
     --quantTranscriptomeBAMcompression -1 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outFileNamePrefix "$output_subdir/$filename"
done
