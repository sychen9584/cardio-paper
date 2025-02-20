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

# %% [markdown]
# # Generate feature df for tss enrichment

# %%
# Define column names based on the provided BED-like structure
columns = ["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"]

# Load the file (tab-separated, without header)
bed_file = "ref/genes.bed"  # Replace with your actual file
df = pd.read_csv(os.path.join(DATA_PATH, bed_file), sep="\t", names=columns, comment="#", dtype=str)

# Extract `gene_id` and `gene_name` using regular expressions
df["gene_id"] = df["Attributes"].str.extract(r'gene_id "([^"]+)"')
df["gene_name"] = df["Attributes"].str.extract(r'gene_name "([^"]+)"')

# Filter for transcript-level entries only (to avoid multiple exon/CDS entries)
df_transcripts = df[df["Feature"] == "transcript"]
# Preserve the order by creating an index column
df_transcripts["Original_Order"] = df_transcripts.index

# Group by gene_id while keeping order based on first occurrence
df_genes = df_transcripts.groupby(["gene_id", "gene_name", "Chromosome"], sort=False).agg(
    Start=("Start", "min"),
    End=("End", "max"),
    Order=("Original_Order", "min")  # Preserve first occurrence order
).reset_index()

# Sort by the original order in the file
df_genes = df_genes.sort_values(by="Order").drop(columns=["Order"])


# %%
df_genes = df_genes[['Chromosome', "Start", "End", 'gene_id', 'gene_name']]

# %%
df_genes.to_csv(os.path.join(DATA_PATH, "ref/gene_intervals.csv"), header=True, index=False)

# %%
