for dir in SRR*; do
	if [ -d "$dir" ]; then
		echo -e "\033[0;32mMerging *.bam files for $dir\033[0m"


		# Define input *.bam files
		bam_r1="$dir/mapped/trim_${dir}_1Aligned.sortedByCoord.out.bam"
		bam_r2="$dir/mapped/trim_${dir}_2Aligned.sortedByCoord.out.bam"

		merged_bam="$dir/mapped/${dir}_merged.bam"

		# Check the existing *.bam files
		if [[ -f "$bam_r1" && -f "$bam_r2" ]]; then
			echo -e "\033[0;34mMerging $bam_r1 and $bam_r2 into $merged_bam...\033[0m"


		# Run samtools to merge
		samtools merge -@ 12 "$merged_bam" "$bam_r1" "$bam_r2"

		echo -e "\033[0;32mMerging completed: $merged_bam\033[0m"
	else
		echo -e "\033[0;31mWarning: One or both *.bam files are missing in $dir/mapped/. Skipping...\033[0m"
	fi
fi

done


echo -e "\033[0;32m*.bam merging completed for all SRA folders!\033[0m"
