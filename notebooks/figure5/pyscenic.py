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
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from matplotlib import pyplot as plt
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
scenic_loom_input = os.path.join(DATA_PATH, "scenic/input.loom")
adj_output = os.path.join(DATA_PATH, "scenic/adj.tsv")

# %%
sc.pp.filter_genes(adata, min_cells=265) # expressed in at least 1% of the cells

# %%
adata.var.set_index("mouse_genesymbol", inplace=True)

# %%
# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create(scenic_loom_input, adata.X.transpose(), row_attrs, col_attrs)

# %%
# After running on AWS
