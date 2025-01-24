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

# %% [markdown] vscode={"languageId": "plaintext"}
# # scRNA-seq data QC and processing

# %% [markdown]
# ## Using 3 months sample 1 as pilot

# %%
import os
import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/raw"
SAMPLE_NAME = "scRNA_m3_s1"

# %%
## Load in as adata object
adata = sc.read_10x_mtx(path=os.path.join(DATA_PATH, SAMPLE_NAME))

# %%
