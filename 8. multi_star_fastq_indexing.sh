#!/bin/bash

for dir in SRR*; do
    if [ -d "$dir" ]; then
        echo -e "\033[0;32mStarting STAR alignment for $dir\033[0m"

        # Define paths
        fastq_dir="$dir/fastq"
        ref_dir="$dir/ref"
        mapped_dir="$dir/mapped"

        # Check if *.fastq files exist
        if ls "$fastq_dir"/*.fastq 1> /dev/null 2>&1; then
            echo -e "\033[0;34m*.fastq files found in $fastq_dir, starting alignment...\033[0m"

            # Loop through each *.fastq file in the fastq/ directory
            for file in "$fastq_dir"/*.fastq; do
                filename=$(basename "$file" .fastq)  # Extract file name without extension
                echo -e "\033[0;34mAligning $filename.fastq...\033[0m"

                # Run STAR alignment for *.fastq files
                STAR --runMode alignReads \
                     --genomeDir "$ref_dir" \
                     --outSAMtype BAM SortedByCoordinate \
                     --readFilesIn "$file" \
                     --runThreadN 12 \
                     --outFileNamePrefix "$mapped_dir/$filename"

                echo -e "\033[0;32mAlignment completed for $filename.fastq\033[0m"
            done

        else
            echo -e "\033[0;31mNo *.fastq files found in $fastq_dir, skipping $dir\033[0m"
        fi
    fi
done

echo -e "\033[0;32mSTAR alignment completed for *.fastq files in all SRR subfolders!\033[0m"
