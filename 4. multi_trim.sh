#!/bin/bash

for dir in SRR*; do
	if [ -d "$dir" ]; then
		 echo -e "\033[0;32mRunning fastp trimming for $dir\033[0m"

	r1="$dir/${dir}_1.fastq.gz"
	r2="$dir/${dir}_2.fastq.gz"

	trimmed_r1="$dir/trim_${dir}_1.fastq.gz"
	trimmed_r2="$dir/trim_${dir}_2.fastq.gz"

	json_report="$dir/fastp_report.json"
	html_report="$dir/fastp_report.html"

		if [[ -f "$r1" && -f "$r2" ]]; then
			echo -e "\033[0;34mTrimming: $r1 and $r2\033[0m"

            if ! fastp --in1 "$r1" --in2 "$r2" --out1 "$trimmed_r1" --out2 "$trimmed_r2" --json "$json_report" --html "$html_report"; then
                echo -e "\033[0;31mfastp failed for $dir\033[0m"
                exit 1
            fi
        else
            echo -e "\033[0;31mWarning: Paired-end files missing in $dir. Skipping...\033[0m"
        fi
    fi
done

echo -e "\033[0;32mFastp trimming completed for all folders!\033[0m"