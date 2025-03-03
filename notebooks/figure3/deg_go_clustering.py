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
import PyComplexHeatmap as pch
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

# %% [markdown]
# ## Recompute DE genes with more strigent criteria

# %%
adata = sc.read_h5ad(os.path.join(DATA_PATH, 'processed/scRNA_all.h5ad'))
adata.var.set_index('mouse_genesymbol', inplace=True)
# rename cell_type_fine names to match visualiztion in the publication
adata.obs['cell_type_fine'] = adata.obs['cell_type_fine'].replace({
    'Granulocyte/Neutrophil': "Neutrophil",
    'Peri/Smooth Muscle': 'Smooth Muscle'
})

# %%
sc.tl.rank_genes_groups(adata, groupby="cell_type_fine", method="wilcoxon", key_added="cell_type_fine_DEG", pts=True, use_raw=False)
sc.tl.rank_genes_groups(adata, groupby="cell_type", method="wilcoxon", key_added="cell_type_DEG", pts=True, use_raw=False)

# %%
sc.tl.filter_rank_genes_groups(adata, key="cell_type_fine_DEG", min_fold_change=1, min_in_group_fraction=0.5, use_raw=False, key_added="cell_type_fine_DEG_filtered")
sc.tl.filter_rank_genes_groups(adata, key="cell_type_DEG", min_fold_change=1, min_in_group_fraction=0.5, use_raw=False, key_added="cell_type_DEG_filtered")

# %%
cell_type_fine_DEGs = sc.get.rank_genes_groups_df(adata, group=None, key="cell_type_fine_DEG_filtered")
cell_type_DEGs = sc.get.rank_genes_groups_df(adata, group=None, key="cell_type_DEG_filtered")

# %%
cell_type_fine_DEGs.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_cell_type_fine.csv'))
cell_type_DEGs.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_cell_type.csv'))

# %% [markdown]
# ## Complex heatmap of DEGs

# %%
cell_type_fine_DEGs = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_cell_type_fine.csv'))
cell_type_DEGs = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_cell_type.csv'))

# %%
expr_df = sc.get.obs_df(adata, keys=['cell_type_fine', *cell_type_fine_DEGs.names.tolist()])

# %%
# sample 2000 cells for visualization
annot_df = adata.obs[['cell_type_fine']].copy().sample(n=2000, random_state=42)
annot_df.sort_values('cell_type_fine', inplace=True)
annot_df["cell_type_fine"] = annot_df["cell_type_fine"].astype(str)
downsampled_expr = expr_df.loc[annot_df.index, :].drop('cell_type_fine', axis=1)

# %%
# Get colors
cell_type_colors, cell_type_fine_colors = fig_func.get_celltype_colors()  
# identify coordinates where cell type changes; for plotting separation lines on heatmap
group_changes = np.where(annot_df["cell_type_fine"].values[:-1] != annot_df["cell_type_fine"].values[1:])[0] + 1

# %%
with plt.rc_context({"figure.dpi": 300, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10,"xtick.labelsize": 10, 
                    "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 8, "figure.titlesize": 16}):

    plt.figure(figsize=(7, 8))
    col_ha = pch.HeatmapAnnotation(
        Celltype=pch.anno_simple(annot_df.cell_type_fine, colors=cell_type_fine_colors, legend=True),  
        axis=1,  # Annotation on the x-axis (columns = cells)
        verbose=1, plot_legend=True, legend_kws={"fontsize": 8}
    )
    cm = pch.ClusterMapPlotter(
        data=downsampled_expr.T,  # Now sorted by `cell_type_fine`
        top_annotation=col_ha,
        row_cluster=False, col_cluster=False, row_dendrogram=False, col_dendrogram=False,
        col_split=annot_df["cell_type_fine"], col_split_gap=0.1, col_split_order=annot_df["cell_type_fine"].unique().tolist(),
        cmap='Reds',  vmax=2.2, label="Expression",
        xlabel="Cell Type", legend_hpad=0, rasterized=True,
        xlabel_kws=dict(color='black', fontsize=8, labelpad=0)
    )
    # 🔹 Add vertical separators at group changes
    for idx in group_changes:
        cm.ax_heatmap.axvline(idx, color='grey', linewidth=10, linestyle="-")  # Dashed separator line
    
    cm.ax.figure.savefig(os.path.join(FIGURE_PATH, "figure3a.png"), dpi=300, bbox_inches='tight', format='png')
    plt.show()

# %%
