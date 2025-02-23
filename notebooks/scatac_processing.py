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
# %load_ext autoreload
# %autoreload 2

# %%
import os
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import anndata
import subprocess
import muon as mu
from muon import atac as ac  # scATAC-seq processing

import sys
sys.path.append('../scripts')
import scatac_preprocessing

# Set plotting preferences
sc.settings.verbosity = 0
sns.set(rc={"figure.figsize": (4, 3.5), "figure.dpi": 100})
sns.set_style("whitegrid")

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/"
FIGURE_PATH = "/home/sychen9584/projects/cardio_paper/figures"
GENE_ANNOTATION_PATH = "/home/sychen9584/projects/cardio_paper/data/ref/gene_intervals.csv"
DOUBLET_SCRIPT_PATH = "/home/sychen9584/projects/cardio_paper/scripts/scDblFinder_script.R"

# %%
# Define the path to the "raw" directory
raw_data_path = os.path.join(DATA_PATH, "raw")

## Filter directories that start with "scRNA"
scATAC_samples = [
    folder for folder in os.listdir(raw_data_path) if folder.startswith("scATAC") 
    and os.path.isdir(os.path.join(raw_data_path, folder))
]

print("Folders starting with 'scATAC:", scATAC_samples)

# %% [markdown]
# # First stage: read in data and generate QC metrics

# %%
for sample in scATAC_samples:
    scatac_preprocessing.preprocess_atac_data(DATA_PATH, sample, FIGURE_PATH,
                                            gene_annot_file=GENE_ANNOTATION_PATH,
                                            script_path=DOUBLET_SCRIPT_PATH)

# %% [markdown]
# # Decide on the thresholds for filtering by drawing joint plots

# %%
filter_dict  = {}

# %%
sample = "scATAC_m24_s2"
adata = sc.read_h5ad(os.path.join(DATA_PATH, "raw", sample, "qc_metrics.h5ad"))

# %%
plot_tss_max = 20
count_cutoff_lower = 5000
count_cutoff_upper = 100000
tss_cutoff_lower = 4

g = scatac_preprocessing.fragment_tss_joint_plot(adata, count_cutoff_lower, count_cutoff_upper, tss_cutoff_lower, 20)

# %%
filter_dict[sample] = {
    "count_cutoff_lower": count_cutoff_lower,
    "count_cutoff_upper": count_cutoff_upper,
    "tss_cutoff_lower": tss_cutoff_lower
}

# %%
g.fig.savefig(os.path.join(FIGURE_PATH, sample, "preprocess/fragment_tss_joint_plot.png"), dpi=150)

# %%
filter_dict # turns out that the same thresholds work for all samples, no need to save them!

# %%
