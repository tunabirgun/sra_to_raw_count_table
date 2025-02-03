#!/bin/bash

for dir in SRR*; do

if [ -d "$dir" ]; then
		echo -e "\033[0;32mRunning MultiQC for $dir\033[0m"

		cd "$dir" || { echo -e "\033[0;31mError: Could not enter directory $dir\033[0m"; exit 1; }

	if ! 	multiqc . -o .; then
		echo -e "\033[0;31mMultiQC failed for $dir\033[0m"
            exit 1
        fi

		cd ..
fi

done

echo -e "\033[0;32mMultiQC analysis completed for all folders!\033[0m"

