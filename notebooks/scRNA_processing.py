# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # scRNA-seq sample preprocessing

# %%
import sys
import os
sys.path.append('../scripts')
from preprocessing import preprocess_adata

# %% [markdown]
# ## Preprocess all scRNA-seq samples

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data"
SAMPLE_NAME = "scRNA_m3_s1"
FIGURE_PATH = "/home/sychen9584/projects/cardio_paper/figures"

# %%
# Define the path to the "raw" directory
raw_data_path = os.path.join(DATA_PATH, "raw")

## Filter directories that start with "scRNA"
scRNA_samples = [
    folder for folder in os.listdir(raw_data_path) if folder.startswith("scRNA") 
    and os.path.isdir(os.path.join(raw_data_path, folder))
]

print("Folders starting with 'scRNA:", scRNA_samples)

# %%
for sample in scRNA_samples:
    preprocess_adata(data_path=DATA_PATH, sample_name=sample, figure_path=FIGURE_PATH)

# %%

# %%
