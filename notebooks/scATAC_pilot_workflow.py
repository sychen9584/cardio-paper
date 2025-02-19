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

# %% [markdown]
# # scATAC-seq data QC and processing

# %% [markdown]
# ## Using 3 months sample 1 as pilot

# %%
# Single-cell packages
import anndata2ri
import matplotlib.pyplot as plt
import muon as mu
import numpy as np

# General helpful packages for data analysis and visualization
import pandas as pd
import scanpy as sc
import seaborn as sns
from muon import atac as ac  # the module containing function for scATAC data processing

# Packages enabling to run R code
from rpy2.robjects import pandas2ri

pandas2ri.activate()  # Automatically convert rpy2 outputs to pandas DataFrames
anndata2ri.activate()
# %load_ext rpy2.ipython


# Setting figure parameters
sc.settings.verbosity = 0
sns.set(rc={"figure.figsize": (4, 3.5), "figure.dpi": 100})
sns.set_style("whitegrid")

import os
import scipy
import anndata

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/raw"
SAMPLE_NAME = "scATAC_m3_s1"

# %% [markdown]
# # Load in Muon object

# %%
# Define the directory containing the ATAC-seq data
matrix_path = os.path.join(DATA_PATH, SAMPLE_NAME, "matrix.mtx.gz")

# Load the sparse peak-barcode matrix
mat = scipy.io.mmread(matrix_path).T.tocsc()  # Convert to CSC format for efficiency

# Load peak information (BED file)
peaks_path = os.path.join(DATA_PATH, SAMPLE_NAME, "peaks.bed.gz")
peaks = pd.read_csv(peaks_path, sep="\t", header=None, names=["chrom", "start", "end"])

# Load barcode metadata
barcodes_path = os.path.join(DATA_PATH, SAMPLE_NAME, "barcodes.tsv.gz")
barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, names=["barcode"])

# Convert barcodes and peaks into the required format
barcodes.index = barcodes["barcode"]  # Set barcodes as index
peaks["peak_id"] = ["peak_" + str(i) for i in range(len(peaks))]  # Generate unique peak IDs
peaks.index = peaks["peak_id"]  # Set peaks as index

# Create AnnData object for ATAC-seq counts
adata_atac = anndata.AnnData(
    X=mat,  # Sparse matrix
    obs=barcodes,  # Barcodes (cells)
    var=peaks  # Peaks (features)
)

# Convert AnnData to MuData (multimodal object)
mdata = mu.MuData({"atac": adata_atac})

# Print summary
print(mdata)

# %%
adata_atac.uns

# %% [markdown]
# # Quality Control

# %% [markdown]
# ## Nucleosome Signal

# %% [markdown]
# # I can't continue due to the lack of fragment files :(

# %%
# Plot fragment size distribution
#ac.pl.fragment_histogram(mdata, region="chr1:1-2000000")

# Plot nucleosome signal distribution across cells
ac.tl.nucleosome_signal(mdata, n=1e6)
mu.pl.histogram(mdata, "nucleosome_signal")

# %%
