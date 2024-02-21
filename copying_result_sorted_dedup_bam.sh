#!/bin/bash
#То есть, тут у нас был ряд папочек (2-23)_trimmed, мы в каждую заходим
#берём файлы, у которых в названии паттерны из строчки 12 и копируем эти файлы в итоговую директорию, которую сами указываем.

# Define the destination directory, were you want to find your data in future
destination_dir="/mnt/projects/users/.../aligned_deduplicated"

# Iterate through each directory
for folder in {2..23}_trimmed; do
    echo "The folder which is in processing now $folder"

    # Use find with -exec to copy files
    find ./$folder -type f \( -name "*deduplicated*Genome*SortByCoord.bam" -o -name "*deduplicated*Transcriptome*SortByCoord.bam" \) -exec cp {} $destination_dir \;

    echo "The $folder result files copying is FINISHED"
    echo -e "\n\n"
done


