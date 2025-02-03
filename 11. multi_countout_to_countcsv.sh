#!/bin/bash

# Loop through all SRR folders
for dir in SRR*; do
    if [ -d "$dir/bams" ]; then
        count_file="$dir/bams/count.out"
        csv_file="$dir/bams/count.csv"
        final_csv_file="${dir}_counts.csv"

        if [ -f "$count_file" ]; then
            echo -e "\033[0;34mConverting count.out to *.csv for $dir...\033[0m"

            # Convert tab-separated (*.out) data to *.csv format
            cat "$count_file" | tr '\t' ',' > "$csv_file"

            # Move the *.csv file to the main directory with correct naming
            mv "$csv_file" "$final_csv_file"

            echo -e "\033[0;32m*.csv file saved as $final_csv_file in the main directory\033[0m"
        else
            echo -e "\033[0;31mWarning: count.out not found in $dir\033[0m"
        fi
    fi
done

echo -e "\033[0;32mConversion to *.csv completed for all SRR folders! Files are now in the main directory.\033[0m"