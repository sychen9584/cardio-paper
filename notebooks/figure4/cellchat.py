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
import pandas as pd
import scanpy as sc
import matplotlib.pylab as plt
import numpy as np
import liana as li
import PyComplexHeatmap as pch
import matplotlib.colors as mcolors
import sys

sys.path.append('../../scripts')
import figure_functions as fig_func

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=120,
    facecolor="white",
    frameon=True,
)

plt.rcParams["axes.grid"] = False  # Disable grids for all plots
plt.rcParams["grid.color"] = "white"  # Ensure grids are fully invisible
plt.rcParams["grid.alpha"] = 0  # Remove any transparency effects
plt.grid(False)  # Turn off grids explicitly

# %matplotlib inline

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/"
FIGURE_PATH = "/home/sychen9584/projects/cardio_paper/figures"

# %%
adata = sc.read_h5ad(os.path.join(DATA_PATH, "processed/scRNA_all.h5ad"))
adata.obs['cell_type_fine'] = adata.obs['cell_type_fine'].replace({
    'Granulocyte/Neutrophil': "Neutrophil",
    'Peri/Smooth Muscle': 'Smooth Muscle'
})

# %%
adata.obs[['month', 'sample_id']] = adata.obs['sample'].str.extract(r'(m\d+)_(s\d+)')
adata_m3 = adata[adata.obs['month'] == 'm3', :].copy()
adata_m12 = adata[adata.obs['month'] == 'm12', :].copy()
adata_m24 = adata[adata.obs['month'] == 'm24', :].copy()

# %%
adata_m3.write_h5ad(os.path.join(DATA_PATH, "processed/scRNA_m3.h5ad"))
adata_m12.write_h5ad(os.path.join(DATA_PATH, "processed/scRNA_m12.h5ad"))
adata_m24.write_h5ad(os.path.join(DATA_PATH, "processed/scRNA_m24.h5ad"))

# %% [markdown]
# # Run Cell chat

# %%
#adata_m3 = sc.read_h5ad(os.path.join(DATA_PATH, "processed/scRNA_m3.h5ad"))
#adata_m12 = sc.read_h5ad(os.path.join(DATA_PATH, "processed/scRNA_m12.h5ad"))
adata_m24 = sc.read_h5ad(os.path.join(DATA_PATH, "processed/scRNA_m24.h5ad"))

# %%
# possible to run cellchat on python but it's missing a lot of functions
#m3_cc = li.mt.cellchat(adata_m3, groupby="cell_type_fine", min_cells=10, n_jobs=4, seed=42, verbose=True, use_raw=False, inplace=False)
#m12_cc = li.mt.cellchat(adata_m12, groupby="cell_type_fine", min_cells=10, n_jobs=4, seed=42, verbose=True, use_raw=False, inplace=False)
m24_cc = li.mt.cellchat(adata_m24, groupby="cell_type_fine", min_cells=10, n_jobs=4, seed=42, verbose=True, use_raw=False, inplace=False)

# %%
#m3_cc = m3_cc.query('cellchat_pvals < 0.05')
#m12_cc = m12_cc.query('cellchat_pvals < 0.05')
m24_cc = m24_cc.query('cellchat_pvals < 0.05')

# %%
#m3_cc.to_csv(os.path.join(DATA_PATH, "cellchat/m3_cellchat.csv"))
#m12_cc.to_csv(os.path.join(DATA_PATH, "cellchat/m12_cellchat.csv"))
m24_cc.to_csv(os.path.join(DATA_PATH, "cellchat/m24_cellchat.csv"))

# %%
m24_cc.query('ligand == "S100A9" and receptor == "CD36"')

# %%
# ?li.mt.cellchat

# %%
