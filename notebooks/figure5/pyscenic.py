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

# %% [markdown]
# # After running SCENIC on AWS, load in resulting regulon scores

# %%
# collect SCENIC AUCell output
lf = lp.connect(os.path.join(DATA_PATH, "scenic/output.loom"), mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

# %%
auc_mtx['cell_type_fine'] = adata.obs['cell_type_fine'].copy()
regulon_celltype = auc_mtx.groupby('cell_type_fine').mean()

# %% [markdown]
# ## Figure 5A

# %%
tfs_to_plot = ['Bclaf1', "Ikzf2", "Ikzf3", "Bach2", "Irf9", "Irf7", "Zfp426", "Hey1", "Sox17", "Smad1",
               "Hltf", "Ets2", "Pparg", "Etv6", "Cux1", "Rara", "Nfkb2", "Zfp384", "Fosl1", "Pbx3", 
               "E2f1", "Maf", "Nfatc2", "Tcf7l2", "Zfp112", "Meis1", "Nr1h4", "Gli3", "Myc", "Jun", "Rela",
               "Mafg", "Hey2", "Zfp467", "Etv1", "Sox6", "Etv4", "Arntl", "Klf2", "Klf4"]

regulon_celltype.columns = regulon_celltype.columns.str.replace(r"\(\+\)", "", regex=True)
# Filter columns where the TF name appears before '('
filtered_columns = [col for col in regulon_celltype.columns if any(tf in col for tf in tfs_to_plot)]
regulon_celltype = regulon_celltype[filtered_columns]

# %%
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()

# %%
regulon_celltype_scaled = pd.DataFrame(scaler.fit_transform(regulon_celltype), columns=regulon_celltype.columns, index=regulon_celltype.index)

# %%
annot_df = pd.DataFrame(regulon_celltype_scaled.index)
annot_df["cell_type_fine"] = annot_df["cell_type_fine"].astype(str)
annot_df.index = annot_df["cell_type_fine"]
# Get colors
cell_type_colors, cell_type_fine_colors = fig_func.get_celltype_colors()  

# %%
with plt.rc_context({"figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10,"xtick.labelsize": 10, 
                    "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 8, "figure.titlesize": 16}):

    plt.figure(figsize=(6, 7))
    col_ha = pch.HeatmapAnnotation(
        Celltype=pch.anno_simple(annot_df.cell_type_fine, colors=cell_type_fine_colors, legend=True),  
        axis=1,  # Annotation on the x-axis (columns = cells)
        verbose=1, plot_legend=True, legend_kws={"fontsize": 8}
    )
    cm = pch.ClusterMapPlotter(
        data=regulon_celltype_scaled.T, 
        top_annotation=col_ha,
        row_cluster=True, col_cluster=False, row_dendrogram=False, col_dendrogram=False,
        vmin=-2, vmax=4, show_rownames=True, cmap="coolwarm", label="Regulon \n Activity")
    cm.ax.figure.savefig(os.path.join(FIGURE_PATH, "figure5/figure5a.png"), dpi=150, bbox_inches='tight', format='png')
    plt.show()

# %% [markdown]
# ## Figure 5B

# %%
adata.obs[['month', 'sample_num']] = adata.obs['sample'].str.extract(r'(m\d+)_(s\d+)')
auc_mtx['month'] = adata.obs['month'].copy()
regulon_celltype_month = auc_mtx.groupby(['cell_type_fine','month']).mean()

# %%
tfs_to_plot = ["Bclaf1", "Dbp", "Thra", "Hlf", "Rarb", "Klf13", "Stat1", "Irf9", "Irf7", "Elf2", "Fosl2", "Atf6"]
regulon_celltype_month.columns = regulon_celltype_month.columns.str.replace(r"\(\+\)", "", regex=True)
# Filter columns where the TF name appears before '('
filtered_columns = [col for col in regulon_celltype_month.columns if any(tf in col for tf in tfs_to_plot)]
regulon_celltype_month = regulon_celltype_month[filtered_columns]

# %%
regulon_celltype_month_scaled = pd.DataFrame(scaler.fit_transform(regulon_celltype_month), columns=regulon_celltype_month.columns, index=regulon_celltype_month.index)

# %%
# reorder months
month_order = ["m3", "m12", "m24"]
regulon_celltype_month_scaled = regulon_celltype_month_scaled.reindex(month_order, level=1)

# %%
annot_df = regulon_celltype_month_scaled.index.to_frame()
annot_df["cell_type_fine"] = annot_df["cell_type_fine"].astype(str)

# %%
with plt.rc_context({"figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10,"xtick.labelsize": 10, 
                    "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 8, "figure.titlesize": 16}):

    plt.figure(figsize=(5, 6))
    col_ha = pch.HeatmapAnnotation(
        Celltype=pch.anno_simple(annot_df.cell_type_fine, colors=cell_type_fine_colors, legend=True), 
        Timepoint=pch.anno_simple(annot_df.month, colors={"m3": "#f69697", "m12": "#9dd189", "m24": "#c6b0d6"}, legend = True),        
        axis=1,  # Annotation on the x-axis (columns = cells)
        verbose=1, plot_legend=True, legend_kws={"fontsize": 8}
    )
    cm = pch.ClusterMapPlotter(
        data=regulon_celltype_month_scaled.T, 
        top_annotation=col_ha,
        row_cluster=True, col_cluster=False, row_dendrogram=False, col_dendrogram=False,
        vmin=-2, vmax=4, show_rownames=True, cmap="coolwarm", label="Regulon \n Activity")
    cm.ax.figure.savefig(os.path.join(FIGURE_PATH, "figure5/figure5b.png"), dpi=150, bbox_inches='tight', format='png')
    plt.show()

# %%
