#!/bin/bash

for dir in SRR*; do

    if [ -d "$dir" ]; then
        sra_file="$dir/$dir.sra"

    if [ -f "$sra_file" ]; then
            echo -e "\033[0;32mProcessing: $sra_file\033[0m"
			
            cd "$dir" || { echo -e "\033[0;31mError: Could not enter directory $dir\033[0m"; exit 1; }
			
            fastq-dump --split-files --gzip "$dir.sra"
			
			
            cd ..

        else
            echo -e "\033[0;31mWarning: No *.sra file found in $dir\033[0m"
        fi
    fi
done