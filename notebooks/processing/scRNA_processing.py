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

# %% [markdown]
# # scRNA-seq sample preprocessing

# %%
import sys
import os
import scanpy as sc
import anndata as ad
sys.path.append('../../scripts')
from preprocessing import preprocess_adata
from scipy.sparse import csr_matrix

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# %% [markdown]
# ## Preprocess all scRNA-seq samples

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data"
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

# %% [markdown]
# ## Combine all samples into one

# %%
adata_paths = {}
for sample in scRNA_samples:
    adata_paths[sample] = os.path.join(DATA_PATH, "preprocessed", sample, f"{sample}_preprocessed.h5ad")

# %%
ad.experimental.concat_on_disk(in_files=adata_paths, out_file=os.path.join(DATA_PATH, "scRNA_all.h5ad"), label="sample")

# %% [markdown]
# ## Dimensional Reduction and Clustering the combined anndata object

# %%
adata = sc.read_h5ad('../data/scRNA_all.h5ad')

# %%
adata.obs_names_make_unique()
adata

# %%
adata.X = adata.layers["counts"].copy()

# %%
# identify variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")

# %%
sc.pl.highly_variable_genes(adata)

# %%
adata.X = adata.layers["lognormalized"].copy()

# %%
# scale highly variable genes
sc.pp.scale(adata, max_value=10)

# %%
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# %%
sc.pp.neighbors(adata, n_pcs=20)
sc.tl.umap(adata)

# %%
# check batch effects, it's minimal
sc.tl.embedding_density(adata, groupby="sample")
sc.pl.embedding_density(adata, basis="umap", key="umap_density_sample")

# %%
# try different resolution for clustering
sc.tl.leiden(adata, resolution=1, flavor="igraph", key_added="res_1")
sc.tl.leiden(adata, resolution=0.8, flavor="igraph", key_added="res_0_8")
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", key_added="res_0_5")

# %%
sc.pl.umap(adata, color=["res_1", "res_0_8", "res_0_5"], wspace=0.5)
# use res = 0.8 for easier annotation

# %%
adata.write_h5ad(os.path.join(DATA_PATH, f"scRNA_all.h5ad"))

# %% [markdown]
# ## Annotating clusters

# %%
import decoupler as dc

# %%
ref_markers = pd.read_excel("../../data/ref/cardio_geneset_ref.xlsx", sheet_name="markers")
#adata = sc.read_h5ad('../data/scRNA_all.h5ad')

# %%
# keep top 20 marker genes for each cell type
ref_makers_top20 = ref_markers.groupby('cluster').head(20)
ref_makers_top20

# %%
adata.X = adata.layers["lognormalized"].copy()

# %%
# run Multivariate Linear Model with decoupler
dc.run_mlm(mat=adata, net=ref_makers_top20, weight=None, source="cluster", target="gene", verbose=True, use_raw=False)

# %%
adata.obsm["mlm_estimate"].head()

# %%
# visualize MLM enrichment results
annotations = adata.obsm["mlm_estimate"].columns.tolist()
acts = dc.get_acts(adata=adata, obsm_key="mlm_estimate")
sc.pl.umap(
    acts,
    color=[
        "res_0_8",
        *annotations
    ],
    wspace=0.5,
    ncols=3,
)

# %%
# manually annotate res_0_8 clusters based on decoupler gene set enrichment results
# this is not going to be perfect as the a lot of the annotations seemed very arbitrary in the paper
adata.obsm["mlm_estimate"]['cluster'] = adata.obs['res_0_8']
enrichment_scores = adata.obsm["mlm_estimate"].groupby('cluster').mean()

# %%
# automatically assign res_0_8 clusters with maximum MLM score cell types
annot_map_fine = enrichment_scores.idxmax(axis=1).to_dict()

# %%
# manually adjust some clusters based on enrichment plots above
annot_map_fine['7'] = "Endo.3"
annot_map_fine['8'] = "Endo.1"
annot_map_fine['9'] = "Endo.4"
annot_map_fine['15'] = "MC.3"

annot_map_fine['0'] = "Fib.1"
annot_map_fine['1'] = "Fib.1"

# %%
adata.obs['cell_type_fine'] = adata.obs['res_0_8'].map(annot_map_fine)

# %%
annot_map = {
    "Endo": "Endothelial",
    "Fib": "Fibroblast",
    'MC': "Macrophage",
    'Granulocyte': "Neutrophil",
    "T-Cell": "T Cell",
    "B-Cell": "B Cell",
    "Peri": "Smooth Muscle"
}

adata.obs['cell_type'] = adata.obs['cell_type_fine'].str.extract(r'(' + '|'.join(annot_map.keys()) + ')')[0].map(annot_map)

# %%
sc.pl.umap(adata, color=["cell_type","cell_type_fine"], wspace=0.5)

# %%
sc.pl.violin(adata, keys=["Ctss", "Lyz2", "C1qb"], groupby="cell_type", rotation=90)

# %%
adata

# %%
adata.write_h5ad(os.path.join(DATA_PATH, f"scRNA_all.h5ad"), compression="gzip")

# %% [markdown]
# # Differentially Expressed Genes

# %%
# Convert the sparse matrix to a dense array
adata.X = adata.X.toarray()
sc.tl.rank_genes_groups(adata, groupby="cell_type", method="wilcoxon", key_added="cell_type_DEG")
sc.tl.filter_rank_genes_groups(adata, key="cell_type_DEG", min_fold_change=1)

# %%
sc.pl.rank_genes_groups_dotplot(adata,  key="cell_type_DEG", groupby="cell_type", standard_scale="var", n_genes=5)

# %%
sc.tl.rank_genes_groups(adata, groupby="cell_type_fine", method="wilcoxon", key_added="cell_type_fine_DEG")
sc.tl.filter_rank_genes_groups(adata, key="cell_type_fine_DEG", min_fold_change=1)

# %%
sc.pl.rank_genes_groups_dotplot(adata, key="cell_type_fine_DEG", groupby="cell_type_fine", standard_scale="var", n_genes=5)

# %%
import pandas as pd

def extract_deg_df(key: str, adjpval_cutoff: float = 1e-5, log2fc_cutoff: float = 1):
    """Extracts differentially expressed genes (DEGs) from Scanpy `adata.uns` results.

    Args:
        key (str): The key in `adata.uns` storing DEG results.
        adjpval_cutoff (float, optional): Adjusted p-value threshold. Defaults to 1e-5.
        log2fc_cutoff (float, optional): Absolute log2 fold-change threshold. Defaults to 1.

    Returns:
        pd.DataFrame: Filtered DataFrame of DEGs.
    """
    
    # Extract DEG results
    deg_results = adata.uns[key]

    # Get the group names (e.g., cell types)
    groups = deg_results["names"].dtype.names

    # Store DEG results for all groups
    deg_dfs = []

    for group in groups:
        df = pd.DataFrame({
            "gene": deg_results["names"][group],
            "logfoldchange": deg_results["logfoldchanges"][group].astype(float),
            "pvals": deg_results["pvals"][group].astype(float),
            "pvals_adj": deg_results["pvals_adj"][group].astype(float),
        })
        df["group"] = group  # Add the cluster/cell type label
        deg_dfs.append(df)

    # Combine all groups into a single DataFrame
    deg_df = pd.concat(deg_dfs, ignore_index=True)

    # Apply filtering (include both upregulated and downregulated genes)
    deg_df = deg_df.query('pvals_adj <= @adjpval_cutoff and abs(logfoldchange) >= @log2fc_cutoff')
    
    return deg_df



# %%
celltype_degs = extract_deg_df(key="cell_type_DEG")
celltype_fine_degs = extract_deg_df(key="cell_type_fine_DEG")

# %%
celltype_degs.to_csv(os.path.join(DATA_PATH, "deg", "scRNA_cell_type.csv"), index=False)
celltype_fine_degs.to_csv(os.path.join(DATA_PATH, "deg", "scRNA_cell_type_fine.csv"), index=False)

# %%
# save final scRNA adata object
adata.X = csr_matrix(adata.X)
adata.write_h5ad(os.path.join(DATA_PATH, f"scRNA_all.h5ad"), compression="gzip")

# %% [markdown]
# # That weird cluster from m12_s2
# Lgals3+ Cotl1+

# %%
adata = sc.read_h5ad(os.path.join(DATA_PATH, f"scRNA_all.h5ad"))

# %%
adata.X = adata.layers["lognormalized"].copy()

# %%
#adata.X = adata.X.toarray()
sc.tl.rank_genes_groups(adata, groupby="res_0_5", method="wilcoxon", key_added="test_DEG")
sc.tl.filter_rank_genes_groups(adata, key="test_DEG", min_fold_change=1)

# %%
sc.get.rank_genes_groups_df(adata, group="18", key="test_DEG").head(25)

# %%
sc.pl.umap(adata, color=["Lgals3", "Cotl1", "Ctss", "doublet_score", "res_0_5"], wspace=0.5)

# %%
