#!/bin/bash

for dir in SRR*; do

	if [ -d "$dir" ]; then
		echo "Running FastQC for files in $dir"

		cd "$dir" || { echo "Error: Could not enter directory $dir"; exit 1; }

		fastqc --outdir . *.fastq.gz

	cd ..

fi
done

echo "FastQC analysis completed for all folder in the current directory."


