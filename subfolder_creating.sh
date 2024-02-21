#!/bin/bash
#this programm will for each sample run the STAR alignmnt, create a subfolder with sample-like name with alignment output

# input_dir-were you store files for processing
input_dir="/mnt/projects/users/..."

#were you want to have subfolders with alignment output
output_base_dir="/mnt/projects/users/..."

#were star-indexed and other ref files are stored
genome_dir="mnt/projects/users/ref_human/..."
gtf_file="mnt/projects/users/ref_human/.../*your_ref_file*.gtf"

#set sjdb_overhang parameter as *length_read*-1 
sjdb_overhang=49

#give to the program bunch of input-file names. this part can be optionally automized also
files_to_process=("1s.fastq" "2s.fastq" "3s.fastq" "4s.fastq" "5s.fastq")

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
