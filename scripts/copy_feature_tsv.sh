#!/bin/bash

# Base directory and features.tsv.gz path
base_directory="/home/sychen9584/projects/cardio_paper/data/raw"
features_file="/home/sychen9584/projects/cardio_paper/data/ref/features.tsv.gz"

# Iterate over directories starting with scRNA
for folder in "$base_directory"/scRNA*; do
    if [ -d "$folder" ]; then
        cp "$features_file" "$folder/features.tsv.gz"
        echo "Copied features.tsv.gz to $folder"
    fi
done
