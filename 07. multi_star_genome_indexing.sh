#!/bin/bash

download_reference_files() {

	local folder="$1"
	echo -e "\033[0;34mChecking reference files for $folder...\033[0m"
	cd "$folder" || { echo -e "\033[0;31mError: Could not enter directory $folder\033[0m"; exit 1; }

	# Check wget
	if ! command -v wget &> /dev/null; then
		echo -e "\033[0;31mError: wget is not installed. Install it using 'sudo apt install wget' or 'conda install wget'.\033[0m"
		exit 1
	fi

	# Download reference genome FASTA file if missing
	if [[ ! -f GCF_000240135.3_ASM24013v3_genomic.fna ]]; then
		echo -e "\033[0;33mDownloading Fusarium graminearum genome...\033[0m"
		wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/135/GCF_000240135.3_ASM24013v3/GCF_000240135.3_ASM24013v3_genomic.fna.gz \
			&& gunzip GCF_000240135.3_ASM24013v3_genomic.fna.gz \
			|| { echo -e "\033[0;31mFailed to download genome file!\033[0m"; exit 1; }
	else
		echo -e "\033[0;32mGenome file already exists. Skipping download.\033[0m"
	fi

	# Download GTF annotation file if missing
	if [[ ! -f GCF_000240135.3_ASM24013v3_genomic.gtf ]]; then
		echo -e "\033[0;33mDownloading Fusarium graminearum annotation...\033[0m"
		wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/135/GCF_000240135.3_ASM24013v3/GCF_000240135.3_ASM24013v3_genomic.gtf.gz \
			&& gunzip GCF_000240135.3_ASM24013v3_genomic.gtf.gz \

	else
		echo -e "\033[0;32mGTF annotation file already exists. Skipping download.\033[0m"
	fi

	cd ..
}


for dir in SRR*; do
	if [ -d "$dir" ]; then
		echo -e "\033[0;32mSetting up directory structure for $dir\033[0m"

		# Create necessary subfolders for STAR analysis
		mkdir -p "$dir/ref" "$dir/fastq" "$dir/mapped"

		# Execute reference files downloading sub-script
		download_reference_files "$dir"

		# Find, move, check and unzip the *.fastq.gz files
		echo -e "\033[0;34mSearching for trimmed FASTQ.GZ files in $dir...\033[0m"
		find "$dir" -type f -name "trim_*.fastq.gz" -exec mv {} "$dir/fastq/" \;

		if ls "$dir/fastq/"*.fastq.gz 1> /dev/null 2>&1; then
			echo -e "\033[0;32mFASTQ.GZ files moved successfully to $dir/fastq/\033[0m"
		else
			echo -e "\033[0;31mWarning: No trimmed FASTQ.GZ files found in $dir\033[0m"
			continue
		fi
		echo -e "\033[0;34mDecompressing FASTQ.GZ files in $dir/fastq/\033[0m"
		gunzip "$dir/fastq/"*.fastq.gz

		echo -e "\033[0;34mListing decompressed FASTQ files in $dir/fastq/:\033[0m"
		ls -lh "$dir/fastq/"

		# Run STAR genome indexing
		echo -e "\033[0;34mRunning STAR genome indexing for $dir...\033[0m"

		STAR --runMode genomeGenerate \
			--genomeDir "$dir/ref/" \
			--genomeFastaFiles "$dir/GCF_000240135.3_ASM24013v3_genomic.fna" \
			--sjdbGTFfile "$dir/GCF_000240135.3_ASM24013v3_genomic.gtf" \
			--runThreadN 16

		echo -e "\033[0;32mSTAR genome indexing completed for $dir\033[0m"
	fi

done

echo -e "\033[0;32mAll SRR folders are set up and the genomes were indexed successfully!\033[0m"