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
    frameon=True,
)
# %matplotlib inline

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/raw"
SAMPLE_NAME = "scRNA_m3_s1"

# %% [markdown]
# ## Load in  adata object

# %%

adata = sc.read_10x_mtx(path=os.path.join(DATA_PATH, SAMPLE_NAME))
adata

# %% [markdown]
# ## Quality control

# %%
# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("mt-")  # "MT-" for human, "mt-" for mouse
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))

# %%
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)

# %%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# %%
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# %%
# filter out cells with <1500 and >7500 detected genes
sc.pp.filter_cells(adata, min_genes = 1500)
sc.pp.filter_cells(adata, max_genes = 7500)
# filter out cells with >10% mitochondrial gene content
adata = adata[adata.obs['pct_counts_mt'] < 10]

# %% [markdown]
# ## Doublet removal using Scrublet

# %%
sc.pp.scrublet(adata)

# %%
sc.pl.scrublet_score_distribution(adata)

# %%
adata = adata[adata.obs['predicted_doublet'] == False]
adata

# %% [markdown]
# ## Normalize data

# %%
# Saving count data
adata.layers["counts"] = adata.X.copy()

# %%
# Normalizing to 10K counts
sc.pp.normalize_total(adata, target_sum=10000)
# Logarithmize the data
sc.pp.log1p(adata)

# %%
adata.layers['lognormalized'] = adata.X.copy()
adata.raw = adata.copy()

# %% [markdown]
# ## Identify highly variable genes

# %%
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")

# %%
sc.pl.highly_variable_genes(adata)

# %% [markdown]
# ## Scale data

# %%
adata = adata[:, adata.var['highly_variable']]
sc.pp.scale(adata, max_value=10)

# %% [markdown]
# ## Dimensional reduction and clustering

# %%
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, flavor="igraph")

# %%
sc.pl.umap(adata, color=["leiden"], legend_loc="on data")

# %%
adata.layers['scaled'] = adata.X.copy()

# %% [markdown]
# ## Re-assess quality control

# %%
sc.pl.umap(
    adata, color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"], wspace=0.5, ncols=2
)

# %% [markdown]
# ## Diffrerential expression analysis

# %%
# Restore the full dataset from adata.raw into a new AnnData object
adata_raw = adata.raw.to_adata()

# Convert the sparse matrix to a dense array
adata_raw.X = adata_raw.X.toarray()

# Use the restored dataset (adata_raw) for filtering or further analysis
sc.tl.rank_genes_groups(adata_raw, groupby="leiden", method="wilcoxon")
sc.tl.filter_rank_genes_groups(adata_raw, min_fold_change=1)

# %%
sc.pl.rank_genes_groups_dotplot(adata_raw, groupby="leiden", standard_scale="var", n_genes=5)

# %%
# visualize macrophage cluster (leiden 6)
cluster6_genes = ["Ctsc", "Lyz2", "Ctss", "C1qb", "Csf1r"]
sc.pl.umap(adata_raw, color=[*cluster6_genes, "leiden"], legend_loc="on data", frameon=False, ncols=3)


# %%
sc.pl.violin(adata_raw, keys=cluster6_genes[0:3], groupby="leiden", multi_panel=True)

# %% [markdown]
# ## Cell type annotation

# %% [markdown]
# Just trying out automated annotation methods here, will recreate manual annotations reported in the publication after merging all samplkes

# %%
import celltypist as ct
import decoupler as dc

# %% [markdown]
# ### Decoupler enrichment cell annotation

# %%
# Query Omnipath and get PanglaoDB cell type annotations'
markers = dc.get_resource(name="PanglaoDB", organism="human")
# Keep canonical markers only
markers = markers[markers["canonical_marker"]]

# Remove duplicated entires
markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]
markers.head()

# %%
# convert to human genes symbols by capitalizing the gene names
adata_raw.var['mouse_genesymbol'] = adata_raw.var_names
adata_raw.var_names = adata_raw.var_names.str.upper()

# %%
dc.run_mlm(mat=adata_raw, net=markers, weight=None, source="cell_type", target="genesymbol", verbose=True, use_raw=False)

# %%
adata_raw.obsm["mlm_estimate"].head()

# %%
# visualize MLM enrichment results
acts = dc.get_acts(adata=adata_raw, obsm_key="mlm_estimate")
sc.pl.umap(
    acts,
    color=[
        "leiden",
        "Fibroblasts",
        "Macrophages",
        "B cells",
        "Endothelial cells",
        "Smooth muscle cells"
    ],
    wspace=0.5,
    ncols=3,
)

# %%
# transfer max overrepresentation score estimates to each cluster
mean_enr = dc.summarize_acts(acts, groupby="leiden", min_std=1)
mean_enr = mean_enr.reindex(sorted(mean_enr.columns), axis=1)
annotation_dict = dc.assign_groups(mean_enr)
adata_raw.obs["dc_anno"] = [annotation_dict[clust] for clust in adata_raw.obs["leiden"]]

# %%
summary = mean_enr

# %%
sc.pl.umap(adata_raw, color=["dc_anno"])

# %% [markdown]
# ### Celltypist automatic annotation

# %%
ct.models.download_models(model=["Healthy_Adult_Heart.pkl"])

# %%
ct_model = ct.models.Model.load(model="Healthy_Adult_Heart.pkl")
predictions = ct.annotate(adata_raw, model=ct_model, majority_voting=True, over_clustering="leiden")
# convert back to anndata
adata_raw = predictions.to_adata()

# %%
sc.pl.umap(adata_raw, color=["majority_voting", "dc_anno"], ncols=1)

# %%
adata_raw.var_names = adata_raw.var["mouse_genesymbol"]
