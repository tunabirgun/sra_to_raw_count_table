#!/bin/bash

# Loop through all SRR folders
for dir in SRR*; do
    if [ -d "$dir" ]; then
        echo -e "\033[0;32mProcessing $dir...\033[0m"

        # Create required subfolders
        mkdir -p "$dir/gtf_ref" "$dir/bams"

        # Move *.gtf reference file to gtf_ref/
        gtf_file="$dir/GCF_000240135.3_ASM24013v3_genomic.gtf"
        if [ -f "$gtf_file" ]; then
            mv "$gtf_file" "$dir/gtf_ref/"
            echo -e "\033[0;34mMoved *.gtf file to $dir/gtf_ref/\033[0m"
        else
            echo -e "\033[0;31mWarning: *.gtf file not found in $dir\033[0m"
        fi

        # Move only the merged *.bam file to bams/
        merged_bam="$dir/mapped/${dir}_merged.bam"
        if [ -f "$merged_bam" ]; then
            mv "$merged_bam" "$dir/bams/"
            echo -e "\033[0;34mMoved merged *.bam file to $dir/bams/\033[0m"
        else
            echo -e "\033[0;31mWarning: Merged *.bam file not found in $dir\033[0m"
        fi

        # Count the number of non-merged *.bam files in mapped/
        bam_count=$(find "$dir/mapped/" -maxdepth 1 -type f -name "trim_*Aligned.sortedByCoord.out.bam" | wc -l)

        echo -e "\033[0;34mNon-merged *.bam file count in $dir/mapped/: $bam_count\033[0m"

        # Debugging: List *.bam files
        echo -e "\033[0;34mListing all *.bam files in $dir/mapped/: \033[0m"
        ls -lh "$dir/mapped/"*.bam 2>/dev/null || echo -e "\033[0;31mNo *.bam files found!\033[0m"

        # If merged *.bam exists and non-merged *.bam files were counted, run featureCounts
        if [[ -f "$dir/bams/${dir}_merged.bam" && "$bam_count" -gt 0 ]]; then
            echo -e "\033[0;34mRunning featureCounts with -T $bam_count on merged *.bam in $dir...\033[0m"

            featureCounts -a "$dir/gtf_ref/GCF_000240135.3_ASM24013v3_genomic.gtf" \
                          -o "$dir/bams/count.out" \
                          -T "$bam_count" \
                          "$dir/bams/${dir}_merged.bam"

            echo -e "\033[0;32mfeatureCounts completed for $dir\033[0m"
        else
            echo -e "\033[0;31mSkipping featureCounts for $dir. Either no merged *.bam or no non-merged *.bam files found.\033[0m"
        fi
    fi
done

echo -e "\033[0;32mfeatureCounts processing completed for all SRR folders!\033[0m"