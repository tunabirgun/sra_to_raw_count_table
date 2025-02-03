#!/bin/bash

for dir in SRR*; do
    if [ -d "$dir" ]; then
        echo -e "\033[0;32mRunning FastQC for trimmed files in $dir\033[0m"

        cd "$dir" || { echo -e "\033[0;31mError: Could not enter directory $dir\033[0m"; exit 1; }

        fastqc --outdir . trim_*_1.fastq.gz trim_*_2.fastq.gz

        cd ..
    fi
done

echo -e "\033[0;32mFastQC analysis completed for all trimmed files!\033[0m"
