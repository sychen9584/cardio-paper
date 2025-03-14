import pandas as pd
import scanpy as sc
import numpy as np
from pyjaspar import jaspardb
import pychromvar as pc

# Load scATAC-seq data stored in AnnData format
adata = sc.read_h5ad("./scATAC_all.h5ad")

# Download the reference genome (mm10) for sequence retrieval
pc.get_genome("mm10", output_dir="./")

# Define valid chromosomes (excluding non-standard contigs)
valid_chroms = [f"chr{i}" for i in range(1, 20)] + ['chrX', 'chrY']

# Filter out peaks that belong to non-standard chromosomes
adata = adata[:, adata.var['Chromosome'].isin(valid_chroms)].copy()

# Use the raw count matrix from layers
adata.X = adata.layers["counts"].copy()
adata.X = adata.X.astype(np.float32)

# Add peak sequences to AnnData object using the mm10 reference genome
pc.add_peak_seq(adata, genome_file="./mm10.fa", delimiter="[:-]")

# Compute GC bias for each peak
pc.add_gc_bias(adata)

# Compute the sum of reads per peak to identify zero-read peaks
reads_per_peak = np.sum(adata.X, axis=0)

# Identify peaks with zero reads
zero_read_peaks = np.where(reads_per_peak == 0)[1]  # Get indices

# Keep only peaks with nonzero reads
nonzero_peaks = np.where(reads_per_peak > 0)[1]  # Indices of nonzero peaks
adata = adata[:, nonzero_peaks]  # Subset AnnData object

# Generate background peaks based on GC bias and read coverage
pc.get_bg_peaks(adata)

# Fetch transcription factor (TF) motifs from JASPAR database (2020 release)
jdb_obj = jaspardb(release='JASPAR2020')
motifs = jdb_obj.fetch_motifs(
    collection='CORE',
    tax_group=['vertebrates']  # Filter motifs for vertebrates
)

# Match motifs to scATAC-seq peaks to identify TF binding sites
pc.match_motif(adata, motifs=motifs)

# Compute motif deviations (TF activity scores) per cell
dev = pc.compute_deviations(adata)

# Save the computed chromVAR deviations to an H5AD file
dev.write_h5ad("./chromvar_deviations.h5ad")
