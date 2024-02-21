#!/bin/bash
#были несколько папок, в них файлы. Из всех, нам нужны только 2 файла - результаты выравнивания на геном и транскриптом. 
#Оба файла нужно отсортировать по координатам, если он не отсортирован (транскриптомные, геномные сразу отсортированы), 
#Затем удалить дубликаты (это сделает fumi-tools, опираясь на UMI ридов), 
#Результаты fumi-tools надо снова отсортировать по координатам (sam sort это делает).
#Если у вас много папок на вход (образцов), то их тоже лучше задать через цикл или как-то ещё, как я сделала ниже.

# Get the current working directory
working_directory=$(pwd)

# Loop through subfolders with the pattern "s_rez"
for foldername in "${working_directory}"/*s_rez*/
do
    echo "The folder which is in processing now $foldername"

    # Step 2
    samtools sort -@ 10 -o "${foldername}/${foldername}_Aligned.toTranscriptome.sorted.out.bam" "${foldername}/${foldername}*Aligned.toTranscriptome.out.bam"

    # Step 3
    echo "Finishing sorting original transcript file and starting demultiplexing it"

    # Step 4
    fumi_tools dedup --threads 10 --memory 10G -i "${foldername}/${foldername}_Aligned.toTranscriptome.sorted.out.bam" -o "${foldername}/${foldername}_deduplicated_Aligned.toTranscriptome.bam"

    # Step 5
    echo "Finishing demultiplexing transcriptome, starting genome demultiplexing"

    # Step 6
    fumi_tools dedup --threads 10 --memory 10G -i "${foldername}/${foldername}*Aligned.sortedByCoord.out.bam" -o "${foldername}/${foldername}_deduplicated_Aligned.toGenome.bam"

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
