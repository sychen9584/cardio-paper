#!/bin/bash

# Define base and target directories
base_directory="/home/sychen9584/projects/cardio_paper/data/raw"
target_directory="/home/sychen9584/projects/cardio_paper/data"
sorted_file="$target_directory/fragments_sorted.tsv"
compressed_file="$target_directory/fragments.tsv.gz"

# Temporary fragment files list
temp_files=()

# Iterate over scATAC* directories
for folder in "$base_directory"/scATAC*; do
    if [[ -d "$folder" ]]; then
        fragment_file="$folder/fragments.tsv.gz"
        sample_id=$(basename "$folder")  # Extract sample name as suffix (e.g., "scATAC_m12_s2")

        if [[ -f "$fragment_file" ]]; then
            temp_file="$target_directory/fragments_${sample_id}.tsv"
            temp_files+=("$temp_file")

            # Decompress and modify cell barcodes by adding a suffix (_sample_id)
            gzip -dc "$fragment_file" | awk -v suffix="_${sample_id}" 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4 suffix,$5}' - > "$temp_file"
        fi
    fi
done

# Ensure we have fragment files to merge
if [[ ${#temp_files[@]} -eq 0 ]]; then
    echo "No fragment files found in $base_directory/scATAC*"
    exit 1
fi

# Merge and sort (enforcing correct chromosome sorting)
echo "Merging and sorting fragment files..."
cat "${temp_files[@]}" | LC_ALL=C sort -k1,1V -k2,2n > "$sorted_file"

# Check if sorting worked
if [[ ! -f "$sorted_file" ]]; then
    echo "Sorting failed. Check input files."
    exit 1
fi

# Compress with multi-threaded bgzip
echo "Compressing sorted merged file..."
bgzip -@ 4 "$sorted_file"

# Ensure compression succeeded
if [[ ! -f "$compressed_file" ]]; then
    echo "Compression failed. Check bgzip installation."
    exit 1
fi

# Index with tabix
echo "Indexing sorted merged fragment file..."
tabix -p bed "$compressed_file"

# Ensure indexing succeeded
if [[ ! -f "$compressed_file.tbi" ]]; then
    echo "Indexing failed. Check tabix."
    exit 1
fi

# Clean up intermediate files
echo "Cleaning up temporary files..."
rm -f "${temp_files[@]}" "$sorted_file"

echo "Merged fragments saved to: $compressed_file"
