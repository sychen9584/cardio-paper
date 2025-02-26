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
# # Perturbation analysis of scRNA-seq data in the context of aging

# %%
import sys
import os
import pickle
import pertpy as pt
import scanpy as sc
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

import warnings

warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

os.environ["PYTHONWARNINGS"] = "ignore"

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=120,
    facecolor="white",
    frameon=True,
)

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data"

# %%
adata = sc.read_h5ad(os.path.join(DATA_PATH, "scRNA_all.h5ad"))
adata

# %%
del adata.layers["lognormalized"]  # Delete log-normalized data
del adata.layers["counts"]    # Delete scaled data

# %%
adata.X = adata.X.astype(np.float32)

# %%
adata.obs[['month', 'sample_num']] = adata.obs['sample'].str.extract(r'(m\d+)_(s\d+)')

# %%
adata.obs['month'] = adata.obs['month'].astype('category')

# %%
adata.obs.cell_type.value_counts()

# %%
adata.obs.month.value_counts()

# %% [markdown]
# ## Augur on cell_type

# %%
ag_rfc = pt.tl.Augur("random_forest_classifier")

# %%
loaded_data = ag_rfc.load(adata, label_col="month", cell_type_col="cell_type")
loaded_data

# %%
v_adata, v_results = ag_rfc.predict(
    loaded_data, subsample_size=20, n_threads=4, select_variance_features=True, span=1
)

v_results["summary_metrics"]

# %%
lollipop = ag_rfc.plot_lollipop(v_results)

# %%
# Define colors in hex format
colors = ["#27408B",  # royalblue4
          "#FFF8DC",  # cornsilk1
          "#FF8C00",  # DarkOrange
          "#CD4F39"]  # tomato3

# Create a custom colormap
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

# %%
sc.pl.umap(adata=v_adata, color=["augur_score"], color_map=custom_cmap)

# %%
import pickle
with open(os.path.join(DATA_PATH, "augur_cell_type_rfc.pkl"), 'wb') as f:
    pickle.dump(v_results, f)

# %% [markdown]
# ## Augur on cell_type_fine

# %%
ag_rfc = pt.tl.Augur("random_forest_classifier")

# %%
adata.obs = adata.obs.drop('cell_type', axis=1)

# %%
loaded_data = ag_rfc.load(adata, label_col="month", cell_type_col="cell_type_fine")
loaded_data

# %%
v_adata, v_results = ag_rfc.predict(
    loaded_data, subsample_size=20, n_threads=4, select_variance_features=True, span=1
)

v_results["summary_metrics"]

# %%
lollipop = ag_rfc.plot_lollipop(v_results)

# %%
# Define colors in hex format
colors = ["#27408B",  # royalblue4
          "#FFF8DC",  # cornsilk1
          "#FF8C00",  # DarkOrange
          "#CD4F39"]  # tomato3

# Create a custom colormap
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

# %%
v_adata.obs['percentile_rank'] = v_adata.obs['augur_score'].rank(pct=True)

# %%
sc.pl.umap(adata=v_adata, color=["percentile_rank"], color_map=custom_cmap)

# %%
import pickle
with open(os.path.join(DATA_PATH, "augur_cell_type_fine_rfc.pkl"), 'wb') as f:
    pickle.dump(v_results, f)

# %% [markdown]
# # Differential prioritization

# %%
ag_rfc = pt.tl.Augur("random_forest_classifier")

# %% [markdown]
# ## 12 month vs 3 month

# %%
# Default mode
adata_12 = ag_rfc.load(
    adata,
    condition_label="m3",
    treatment_label="m12",
    label_col="month"
)

# %%
_, v_results_12 = ag_rfc.predict(
    adata_12, random_state=None, n_threads=4
)
v_results_12["summary_metrics"].loc["mean_augur_score"].sort_values(
    ascending=False
)

# %%
# Permute mode
_, v_results_12_permute = ag_rfc.predict(
    adata_12,
    augur_mode="permute",
    n_subsamples=100,
    random_state=None,
    n_threads=4,
)

# %%
with open(os.path.join(DATA_PATH, "augur_cell_type_m12_v_m3_rfc.pkl"), 'wb') as f:
    pickle.dump(v_results_12, f)
    
with open(os.path.join(DATA_PATH, "augur_cell_type_m12_v_m3_permute_rfc.pkl"), 'wb') as f:
    pickle.dump(v_results_12_permute, f)

# %% [markdown]
# ## 24 months vs 12 months

# %%
# Default mode
adata_24 = ag_rfc.load(
    adata,
    condition_label="m12",
    treatment_label="m24",
    label_col="month"
)

# %%
_, v_results_24 = ag_rfc.predict(
    adata_24, random_state=None, n_threads=4
)
v_results_24["summary_metrics"].loc["mean_augur_score"].sort_values(
    ascending=False
)

# %%
# Permute mode
_, v_results_24_permute = ag_rfc.predict(
    adata_24,
    augur_mode="permute",
    n_subsamples=100,
    random_state=None,
    n_threads=4,
)

# %%
with open(os.path.join(DATA_PATH, "augur_cell_type_m24_v_m12_rfc.pkl"), 'wb') as f:
    pickle.dump(v_results_24, f)
    
with open(os.path.join(DATA_PATH, "augur_cell_type_m24_v_m12_permute_rfc.pkl"), 'wb') as f:
    pickle.dump(v_results_24_permute, f)

# %% [markdown]
# ### Visualize the results

# %%
results_12_v_3 = pickle.load(open(os.path.join(DATA_PATH, "augur_cell_type_m12_v_m3_rfc.pkl"), 'rb'))
results_24_v_12 = pickle.load(open(os.path.join(DATA_PATH, "augur_cell_type_m24_v_m12_rfc.pkl"), 'rb'))

# %%
scatter = ag_rfc.plot_scatterplot(results_24_v_12, results_12_v_3)

# %%
pvals = ag_rfc.predict_differential_prioritization(
    augur_results1=results_24_v_12,
    augur_results2=results_12_v_3,
    permuted_results1=results_24_v_12_permute,
    permuted_results2=results_12_v_3_permute,
)
pvals

# %% [markdown]
# # Clean up results files to save space

# %%
results = pickle.load(open(os.path.join(DATA_PATH, "augur_cell_type_rfc.pkl"), 'rb'))

# %%
results['summary_metrics'].to_csv(os.path.join(DATA_PATH, "augur_cell_type_rfc.csv"))

# %%
