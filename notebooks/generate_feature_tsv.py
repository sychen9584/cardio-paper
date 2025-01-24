# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: cardio
#     language: python
#     name: python3
# ---

# %%
import gzip 
import pandas as pd
import os

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/"
gtf_file = "ref/genes.gtf"
output_file = "ref/features.tsv.gz"


# %%
# Function to parse GTF and extract gene information
def parse_gtf(gtf_file):
    features = []
    with open(gtf_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith("#"):
                continue
            
            # Split the line into fields
            fields = line.strip().split("\t")
            
            # We only care about the "gene" feature
            if fields[2] == "gene":
                # Parse attributes (column 9)
                attributes = {key.strip(): value.strip().strip('"') for key, value in 
                              [attr.split(" ", 1) for attr in fields[8].split("; ") if attr.strip()]}
                
                # Extract gene_id and gene_name
                gene_id = attributes.get("gene_id", "unknown_gene_id")
                gene_name = attributes.get("gene_name", "unknown_gene_name")
                feature_type = "Gene Expression"  # Default feature type
                
                # Append to list
                features.append([gene_id, gene_name, feature_type])
    return features

# Parse the GTF file
features = parse_gtf(os.path.join(DATA_PATH, gtf_file))

# %%
features_df = pd.DataFrame(features, columns=["Feature ID", "Feature Name", "Feature Type"])
features_df.to_csv(os.path.join(DATA_PATH, output_file), sep='\t', header=False, index=False, compression="gzip")

print(f"features.tsv.gz created with {len(features)} entries!")
