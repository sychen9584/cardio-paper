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
def compute_adata_pariwise_DEGs(adata, celltypes, groupby, group, reference, pval_cutoff=0.05, log2fc_cutoff=1, min_pts=0):
    
    deg_results = []
    
    for celltype in celltypes:
        
        print(f"Running DEG analysis for {celltype}...")
        adata_subset = adata[adata.obs["cell_type"] == celltype].copy()
        sc.tl.rank_genes_groups(adata_subset, groupby=groupby, reference=reference, method="wilcoxon", pts=True, use_raw=False)
        degs = sc.get.rank_genes_groups_df(adata_subset, group=group)
        degs = degs.query('pct_nz_group > @min_pts and pvals_adj < @pval_cutoff')
        degs['abs_lfc'] = np.abs(degs['logfoldchanges'])
        degs = degs.query('abs_lfc > @log2fc_cutoff').drop('abs_lfc', axis=1)
        
        degs["cell_type"] = celltype
        deg_results.append(degs)
        
        print(f"DEGs found for {celltype}: {degs.shape[0]} genes")
    
    deg_results = pd.concat(deg_results)
    
    return deg_results


# %%
deg_results_df = compute_adata_pariwise_DEGs(adata_3m_12m, 
                                             celltypes=["Macrophage", "Fibroblast", "Endothelial"], 
                                             groupby="month", group="m12", reference="m3", 
                                             pval_cutoff=0.05, log2fc_cutoff=0.25, min_pts=0.25)

# %%
deg_results_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_3m_vs_12m_DEG.csv'), index=False)

# %% [markdown]
# #### 12 months vs 24 months

# %%
# 🔹 Subset data for 3m vs 12m comparison
adata_12m_24m = adata[adata.obs["month"].isin(["m12", "m24"])].copy()

# %%
deg_results_df = compute_adata_pariwise_DEGs(adata_12m_24m, 
                                             celltypes=["Macrophage", "Fibroblast", "Endothelial"], 
                                             groupby="month", group="m24", reference="m12", 
                                             pval_cutoff=0.05, log2fc_cutoff=0.25, min_pts=0.25)

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
figure_functions.venn3_custom(macrophage, fibroblast, endothelial, labels=('Macrophage', 'Fibroblast', 'Endothelial'), 
                              normalize_range=(0, 2500), title='3 months vs 12 months DEGs')
plt.show()

# %%
#### 12 months vs 24 months
deg_results_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_12m_vs_24m_DEG.csv'))

# %%
macrophage = set(deg_results_df[deg_results_df['cell_type'] == 'Macrophage']['names'].values)
fibroblast = set(deg_results_df[deg_results_df['cell_type'] == 'Fibroblast']['names'].values)
endothelial = set(deg_results_df[deg_results_df['cell_type'] == 'Endothelial']['names'].values)

# %%
figure_functions.venn3_custom(macrophage, fibroblast, endothelial, labels=('Macrophage', 'Fibroblast', 'Endothelial'), 
                              normalize_range=(0, 2500), title='12 months vs 24 months DEGs')
plt.show()

# %% [markdown]
# ### Figure 3D - scRNA-seq DEG heatmap for the three cell types

# %%
deg_results = []
    
for celltype in ["Macrophage", "Fibroblast", "Endothelial"]:
    
    print(f"Running DEG analysis for {celltype}...")
    adata_subset = adata[adata.obs["cell_type"] == celltype].copy()
    sc.tl.rank_genes_groups(adata_subset, groupby='month', reference="m3", method="wilcoxon", pts=True, use_raw=False)
    degs = sc.get.rank_genes_groups_df(adata_subset, group=None)
    
    degs["cell_type"] = celltype
    deg_results.append(degs)
    
deg_results = pd.concat(deg_results)

# %%
deg_results['abs_lfc'] = np.abs(deg_results['logfoldchanges'])

# %%
macrophage_degs = deg_results.query('cell_type == "Macrophage" and pvals_adj < 0.05 and abs_lfc > 0.25 and pct_nz_group > 0.25').names.unique().tolist()
fibroblast_degs = deg_results.query('cell_type == "Fibroblast" and pvals_adj < 0.05 and abs_lfc > 0.25 and pct_nz_group > 0.25').names.unique().tolist()
endothelial_degs = deg_results.query('cell_type == "Endothelial" and pvals_adj < 0.05 and abs_lfc > 0.25 and pct_nz_group > 0.25').names.unique().tolist()

# %%
macrophage_deg_df = deg_results[['names', 'group', 'cell_type', 'logfoldchanges']].query('cell_type == "Macrophage" and names in @macrophage_degs').\
    drop(columns='cell_type').pivot(index='names', columns='group', values='logfoldchanges')
fibroblast_deg_df = deg_results[['names', 'group', 'cell_type', 'logfoldchanges']].query('cell_type == "Fibroblast" and names in @fibroblast_degs').\
    drop(columns='cell_type').pivot(index='names', columns='group', values='logfoldchanges')
endothelial_deg_df = deg_results[['names', 'group', 'cell_type', 'logfoldchanges']].query('cell_type == "Endothelial" and names in @endothelial_degs').\
    drop(columns='cell_type').pivot(index='names', columns='group', values='logfoldchanges')

# %%
macrophage_deg_df['m3'] = 0
fibroblast_deg_df['m3'] = 0
endothelial_deg_df['m3'] = 0

macrophage_deg_df = macrophage_deg_df[['m3', 'm12', 'm24']]
fibroblast_deg_df = fibroblast_deg_df[['m3', 'm12', 'm24']]
endothelial_deg_df = endothelial_deg_df[['m3', 'm12', 'm24']]

# %%
macrophage_deg_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_macrophage_logfc.csv'))
fibroblast_deg_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_fibroblast_logfc.csv'))
endothelial_deg_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_endothelial_logfc.csv'))

# %% [markdown]
# ### Kmeans clustering

# %%
from sklearn.cluster import KMeans

# %%
macrophage_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_macrophage_logfc.csv'), index_col=0)
fibroblast_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_fibroblast_logfc.csv'), index_col=0)
endothelial_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_endothelial_logfc.csv'), index_col=0)


# %%
def kmeans_elbow_plot(range_n_clusters, data, title, random_seed=42):
    # 🔹 Test different cluster numbers (2 to 10)
    range_n_clusters = range(*range_n_clusters)
    wcss = []
    for n_clusters in range_n_clusters:
        kmeans = KMeans(n_clusters=n_clusters, random_state=random_seed).fit(data)
        wcss.append(kmeans.inertia_)
        
   # 🔹 Plot the Elbow Curve
    plt.figure(figsize=(6,4))
    plt.plot(range(2, 11), wcss, marker='o', linestyle='--')
    plt.xlabel("Number of Clusters")
    plt.ylabel("WCSS (Inertia)")
    plt.title("Elbow Plot - " + title)
    plt.show()


# %%
kmeans_elbow_plot((2, 11), macrophage_deg_df, "Macrophage") # 7 clusters

# %%
kmeans_elbow_plot((2, 11), fibroblast_deg_df, "Fibroblast") # 6 clusters

# %%
kmeans_elbow_plot((2, 11), endothelial_deg_df, "Endothelial") #5 clusters


# %%
def kmeans_cluster(data, n_clusters, random_seed=42):
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_seed).fit(data)
    cluster_labels = kmeans.labels_
    return cluster_labels


# %%
macrophage_deg_df['cluster'] = kmeans_cluster(macrophage_deg_df, 7)
#fibroblast_deg_df['cluster'] = kmeans_cluster(fibroblast_deg_df, 6)
#endothelial_deg_df['cluster'] = kmeans_cluster(endothelial_deg_df, 4)

# %%
macrophage_deg_df.sort_values('cluster', inplace=True)
#fibroblast_deg_df.sort_values('cluster', inplace=True)
#endothelial_deg_df.sort_values('cluster', inplace=True)

# %%
macrophage_deg_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_macrophage_logfc.csv'))
#fibroblast_deg_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_fibroblast_logfc.csv'))
#endothelial_deg_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_endothelial_logfc.csv'))

# %% [markdown]
# ### Heatmap time

# %%
macrophage_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_macrophage_logfc.csv'), index_col=0)
fibroblast_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_fibroblast_logfc.csv'), index_col=0)
endothelial_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_endothelial_logfc.csv'), index_col=0)


# %%
def plot_log2fc_heatmap(df, title, column_names, cluster_col, cluster_order, 
                        cmap='coolwarm', vmin=-2.5, vmax=2.5, title_xpad=0.6, title_ypad=1.1, 
                        plot_legend=True, ax=None):
    """
    Plot a log2 fold-change heatmap using PyComplexHeatmap.

    Parameters:
        df (pd.DataFrame): Dataframe containing log2 fold changes.
        title (str): Title for the heatmap.
        cluster_col (str): Column containing cluster annotations.
        cluster_order (list): Custom order for cluster splitting.
        cmap (str): Colormap for the heatmap.
        vmin (float): Minimum value for colormap scaling.
        vmax (float): Maximum value for colormap scaling.
        plot_legend (bool): Whether to display the heatmap legend.
        ax (matplotlib.axes._subplots.AxesSubplot, optional): Existing axis to plot on.

    Returns:
        matplotlib.axes._subplots.AxesSubplot: The axis containing the heatmap.
    """
    
    annot_df = pd.DataFrame(df[cluster_col])
    unique_clusters = annot_df[cluster_col].unique()
    color_dict = {cluster: mcolors.to_hex(cm.get_cmap("Set3")(i % 10)) for i, cluster in enumerate(unique_clusters)}
    
    plot_df = df.drop(columns=cluster_col)
    plot_df.columns = column_names

    row_ha = pch.HeatmapAnnotation(
        cluster=pch.anno_simple(annot_df[cluster_col], colors=color_dict, add_text=True, legend=False), 
        axis=0, verbose=1, legend=False, label_side='top',
    )

    # Create figure and axis if none provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(3, 7))
        
    hm = pch.ClusterMapPlotter(
        data=plot_df,  # Now sorted by `cell_type_fine`
        left_annotation=row_ha,
        row_cluster=False, col_cluster=False, row_dendrogram=False, col_dendrogram=False,
        cmap=cmap, label="", vmin=vmin, vmax=vmax, row_split=annot_df[cluster_col], row_split_gap=0.75, row_split_order=cluster_order,
        legend_hpad=0, legend_vpad=144, rasterized=True, show_colnames=True, col_names_side='top', xticklabels_kws={'labelrotation':45},
        plot_legend=plot_legend
    )

    ax.set_title(title, fontsize=14 , x=title_xpad, y=title_ypad)
    
    return ax


# %%
plot_log2fc_heatmap(macrophage_deg_df, "Macrophage", ['3 months', '12 months', '24 months'], 'cluster', [4, 6, 0, 3, 5, 1, 2], title_ypad=1.15)
plt.show()

# %%
plot_log2fc_heatmap(fibroblast_deg_df, "Fibroblast", ['3 months', '12 months', '24 months'], 'cluster', [3, 4, 0, 5, 1, 2], title_ypad=1.15)
plt.show()

# %%
plot_log2fc_heatmap(endothelial_deg_df, "Endothelial", ['3 months', '12 months', '24 months'], 'cluster', [3, 2, 0 ,1], title_ypad=1.15)
plt.show()

# %% [markdown]
# ## GO analysis on specific clusters

# %%
from gprofiler import GProfiler

# %%
macrophage_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_macrophage_logfc.csv'), index_col=0)
fibroblast_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_fibroblast_logfc.csv'), index_col=0)
endothelial_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_endothelial_logfc.csv'), index_col=0)

# %%
cluster_dict = {
    'endo_up': endothelial_deg_df.query('cluster == 3').index.tolist(),
    'endo_down': endothelial_deg_df.query('cluster == 1').index.tolist(),
    
    'fib_up': fibroblast_deg_df.query('cluster == 3').index.tolist(),
    'fib_down': fibroblast_deg_df.query('cluster == 2').index.tolist(),
    
    'mac_up': macrophage_deg_df.query('cluster == 4').index.tolist(),
    'mac_down': macrophage_deg_df.query('cluster == 1').index.tolist()
}

# %%
gp = GProfiler(return_dataframe=True) #return pandas dataframe or plain python structures  )
go_df = gp.profile(organism='mmusculus', user_threshold=0.05, significance_threshold_method='fdr', sources=["GO:BP", "GO:CC", "KEGG"], query=cluster_dict, background=adata.var_names.tolist())

# %%
go_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_cluster_GO_processed.csv'), index=False)

# %%
highlight_idx = [1537, 2355, 2568, 230, 2629, 3448, 1156, 1570, 1779, 693, 195, 233, 213, 305, 455, 461, 189, 205]
highlight_df = go_df.loc[highlight_idx, :]

# %%
highlight_df['neg_log_p_value'] = -np.log10(highlight_df['p_value'])
highlight_df = highlight_df[['name', 'p_value', 'neg_log_p_value', 'query']]

# %%
highlight_df['query'] = pd.Categorical(highlight_df['query'], categories=['endo_up', 'endo_down', 'fib_up', 'fib_down', 'mac_up', 'mac_down'], ordered=True)
highlight_df.sort_values(['query', 'neg_log_p_value'], ascending=[True, True], inplace=True)

# %%
highlight_df

# %%
highlight_df.to_csv(os.path.join(DATA_PATH, 'deg/scRNA_cluster_GO_highlight.csv'), index=False)

# %%
