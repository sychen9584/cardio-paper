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
from matplotlib_venn import venn3_unweighted
import matplotlib.cm as cm
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

# %% [markdown]
# ## GO analysis on cell type fine annotations

# %%
from gprofiler import GProfiler

# %%
cell_type_fine_DEGs = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_cell_type_fine.csv'))
cell_type_DEGs = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_cell_type.csv'))

# %%
deg_dict = cell_type_fine_DEGs.groupby('group')['names'].apply(list).to_dict()

# %%
gp = GProfiler(return_dataframe=True) #return pandas dataframe or plain python structures    )
go_df = gp.profile(organism='mmusculus', user_threshold=0.001, significance_threshold_method='fdr', query=deg_dict, background=cell_type_fine_DEGs['names'].unique().tolist())

# %%
go_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_cell_type_fine_GO_processed.csv'), index=False)

# %%
highlight_idx = [181, 462, 1058, 1112, 419, 761, 116, 1321, 413, 2081, 1928, 1713, 545, 656, 488, 415, 8, 804, 224]
highlight_df = go_df.loc[highlight_idx, :]

# %%
highlight_df['neg_log_p_value'] = -np.log10(highlight_df['p_value'])
highlight_df = highlight_df[['name', 'p_value', 'neg_log_p_value', 'query']]

# %%
highlight_df

# %%
highlight_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_cell_type_fine_GO_highlight.csv'), index=False)

# %% [markdown]
# ### Compare # of DEGs between cell types for venn diagram visualization

# %%
adata.obs[['month', 'sample_id']] = adata.obs['sample'].str.extract(r'(m\d+)_(s\d+)')

# %% [markdown]
# #### 3 months vs 12 months

# %%
# 🔹 Subset data for 3m vs 12m comparison
adata_3m_12m = adata[adata.obs["month"].isin(["m3", "m12"])].copy()

# %%
adata_3m_12m.obs["cell_type"].unique()

# %%
deg_results = []

# 🔹 Iterate over each unique cell type
for cell_type in ["Macrophage", "Fibroblast", "Endothelial"]:
    print(f"Running DEG analysis for {cell_type}...")

    # 🔹 Subset AnnData for the specific cell type
    adata_subset = adata_3m_12m[adata_3m_12m.obs["cell_type"] == cell_type].copy()

    # 🔹 Run DE analysis between 3m and 12m
    sc.tl.rank_genes_groups(adata_subset, groupby="month", reference="m3", method="wilcoxon", pts=True, use_raw=False)

    # 🔹 Convert results to a DataFrame
    deg_df = sc.get.rank_genes_groups_df(adata_subset, pval_cutoff=0.05, log2fc_min=1, group="m12")  # Compare 12m vs 3m
    #deg_df = deg_df.query('pct_nz_group > 0.25') # Keep genes expressed in >50% of cells in the group
    
    deg_df["cell_type"] = cell_type
    deg_results.append(deg_df)

    print(f"DEGs found for {cell_type}: {deg_df.shape[0]} genes")
    
deg_results_df = pd.concat(deg_results)

# %%
deg_results_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_3m_vs_12m_DEG.csv'), index=False)

# %% [markdown]
# #### 12 months vs 24 months

# %%
# 🔹 Subset data for 3m vs 12m comparison
adata_12m_24m = adata[adata.obs["month"].isin(["m12", "m24"])].copy()

# %%
deg_results = []

# 🔹 Iterate over each unique cell type
for cell_type in ["Macrophage", "Fibroblast", "Endothelial"]:
    print(f"Running DEG analysis for {cell_type}...")

    # 🔹 Subset AnnData for the specific cell type
    adata_subset = adata_12m_24m[adata_12m_24m.obs["cell_type"] == cell_type].copy()

    # 🔹 Run DE analysis between 12m and 24m
    sc.tl.rank_genes_groups(adata_subset, groupby="month", reference="m12", method="wilcoxon", pts=True, use_raw=False)

    # 🔹 Convert results to a DataFrame
    deg_df = sc.get.rank_genes_groups_df(adata_subset, pval_cutoff=0.05, log2fc_min=1, group="m24")  # Compare 12m vs 3m
    #deg_df = deg_df.query('pct_nz_group > 0.25') # Keep genes expressed in >50% of cells in the group
    
    deg_df["cell_type"] = cell_type
    deg_results.append(deg_df)

    print(f"DEGs found for {cell_type}: {deg_df.shape[0]} genes")
    
deg_results_df = pd.concat(deg_results)

# %%
deg_results_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_12m_vs_24m_DEG.csv'), index=False)

# %% [markdown]
# ### Venn diagram for DEGs

# %%
from matplotlib_venn import venn3_unweighted
import matplotlib.cm as cm
import figure_functions

# %%
#### 3 months vs 12 months
deg_results_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_3m_vs_12m_DEG.csv'))

# %%
macrophage = set(deg_results_df[deg_results_df['cell_type'] == 'Macrophage']['names'].values)
fibroblast = set(deg_results_df[deg_results_df['cell_type'] == 'Fibroblast']['names'].values)
endothelial = set(deg_results_df[deg_results_df['cell_type'] == 'Endothelial']['names'].values)

# %%
figure_functions.venn3_custom(macrophage, fibroblast, endothelial, labels=('Macrophage', 'Fibroblast', 'Endothelial'), title='3 months vs 12 months DEGs')
plt.show()

# %%
#### 12 months vs 24 months
deg_results_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_12m_vs_24m_DEG.csv'))

# %%
macrophage = set(deg_results_df[deg_results_df['cell_type'] == 'Macrophage']['names'].values)
fibroblast = set(deg_results_df[deg_results_df['cell_type'] == 'Fibroblast']['names'].values)
endothelial = set(deg_results_df[deg_results_df['cell_type'] == 'Endothelial']['names'].values)

# %%
figure_functions.venn3_custom(macrophage, fibroblast, endothelial, labels=('Macrophage', 'Fibroblast', 'Endothelial'), title='12 months vs 24 months DEGs')
plt.show()

# %%
