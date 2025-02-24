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
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import anndata
import scipy
import pyranges as pr
from muon import atac as ac  # scATAC-seq processing

import sys
sys.path.append('../scripts')
import scatac_preprocessing

# Set plotting preferences
sc.settings.verbosity = 0
sns.set(rc={"figure.figsize": (4, 3.5), "figure.dpi": 100})
sns.set_style("whitegrid")

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/"
FIGURE_PATH = "/home/sychen9584/projects/cardio_paper/figures"
GENE_ANNOTATION_PATH = "/home/sychen9584/projects/cardio_paper/data/ref/gene_intervals.csv"
DOUBLET_SCRIPT_PATH = "/home/sychen9584/projects/cardio_paper/scripts/scDblFinder_script.R"

# %%
# Define the path to the "raw" directory
raw_data_path = os.path.join(DATA_PATH, "raw")

## Filter directories that start with "scRNA"
scATAC_samples = [
    folder for folder in os.listdir(raw_data_path) if folder.startswith("scATAC") 
    and os.path.isdir(os.path.join(raw_data_path, folder))
]

print("Folders starting with 'scATAC:", scATAC_samples)

# %% [markdown]
# # Read in data for each sample and generate QC metrics

# %%
for sample in scATAC_samples:
    scatac_preprocessing.preprocess_atac_data(DATA_PATH, sample, FIGURE_PATH,
                                            gene_annot_file=GENE_ANNOTATION_PATH,
                                            script_path=DOUBLET_SCRIPT_PATH)

# %% [markdown]
# # Decide on the thresholds for filtering by drawing joint plots

# %%
filter_dict  = {}

# %%
sample = "scATAC_m24_s2"
adata = sc.read_h5ad(os.path.join(DATA_PATH, "raw", sample, "qc_metrics.h5ad"))

# %%
plot_tss_max = 20
count_cutoff_lower = 5000
count_cutoff_upper = 100000
tss_cutoff_lower = 4

g = scatac_preprocessing.fragment_tss_joint_plot(adata, count_cutoff_lower, count_cutoff_upper, tss_cutoff_lower, 20)

# %%
filter_dict[sample] = {
    "count_cutoff_lower": count_cutoff_lower,
    "count_cutoff_upper": count_cutoff_upper,
    "tss_cutoff_lower": tss_cutoff_lower
}

# %%
g.fig.savefig(os.path.join(FIGURE_PATH, sample, "preprocess/fragment_tss_joint_plot.png"), dpi=150)

# %%
filter_dict # turns out that the same thresholds work for all samples, no need to save them!

# %% [markdown]
# # Filter each sample based on the decided thresholds

# %%
for sample in scATAC_samples:
    adata = sc.read_h5ad(os.path.join(DATA_PATH, "raw", sample, "qc_metrics.h5ad"))
    sample_preprocessed_path = os.path.join(DATA_PATH, "preprocessed", sample)
    os.makedirs(sample_preprocessed_path, exist_ok=True)
    
    logging.info(f"Filtering {sample}")
    scatac_preprocessing.filter_atac_adata(adata,
                                           count_cutoff_lower=5000, 
                                           count_cutoff_upper=100000,
                                           nuclesome_threshold= 4,
                                           tss_threshold=3)
    adata.X = adata.X.tocsr()
    adata.raw = adata.copy()
    adata.write(os.path.join(DATA_PATH, "preprocessed", sample, "preprocessed.h5ad"))

# %% [markdown]
# # Create a unified peak set for all samples

# %%
adata_paths = {}
for sample in scATAC_samples:
    adata_paths[sample] = os.path.join(DATA_PATH, "preprocessed", sample, "preprocessed.h5ad")

# %%
adatas = [sc.read_h5ad(adata_paths[sample]) for sample in scATAC_samples]

# %%
original_peak_list = []
for adata in adatas:
    peaks = adata.var[["chrom", "start", "end"]].copy()
    original_peak_list.append(peaks)
    
original_peaks = pd.concat(original_peak_list, ignore_index=True)

# %%
# Convert to PyRanges object for efficient merging
gr = pr.PyRanges(original_peaks.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
# Merge overlapping peaks
merged_gr = gr.merge()
print(f"Total unified peaks after merging: {len(merged_gr)}")

# %%
# Find where each original peak maps to the merged peak set
peak_mapping = gr.join(merged_gr, how="left").df
# Create a mapping dictionary {original_peak: merged_peak}
peak_mapping["original_peak"] = peak_mapping["Chromosome"].astype(str) + ":" + peak_mapping["Start"].astype(str) + "-" + peak_mapping["End"].astype(str)
peak_mapping["merged_peak"] = peak_mapping["Chromosome"].astype(str) + ":" + peak_mapping["Start_b"].astype(str) + "-" + peak_mapping["End_b"].astype(str)
# Remove unmapped peaks and keep only necessary columns
peak_mapping = peak_mapping.dropna(subset=["merged_peak"])[["original_peak", "merged_peak"]].drop_duplicates()

# %%
peak_mapping = pd.read_csv(os.path.join(DATA_PATH, "unified_peaks.csv"))

# %%
for sample in scATAC_samples:
    adata = sc.read_h5ad(os.path.join(DATA_PATH, "preprocessed", sample, "preprocessed.h5ad"))
    logging.info(f"Mapping peaks for {sample} ---------------------------------")
    adata_new = scatac_preprocessing.map_unified_peaks(adata, peak_mapping)
    adata_new.write(os.path.join(DATA_PATH, "preprocessed", sample, "unified.h5ad"))

# %% [markdown]
# # Combine all samples into one

# %%
adata_paths = {}
for sample in scATAC_samples:
    adata_paths[sample] = os.path.join(DATA_PATH, "preprocessed", sample, "unified.h5ad")

# %%
anndata.experimental.concat_on_disk(in_files=adata_paths, out_file=os.path.join(DATA_PATH, "scATAC_all.h5ad"), label="sample_name", join="outer", index_unique="_")

# %% [markdown]
# # Dimensional Reduction and Clustering the combined anndata object

# %%
adata = sc.read_h5ad('../data/scATAC_all.h5ad')

# %%
ac.pp.tfidf(adata, scale_factor=1e4)
ac.tl.lsi(adata) # latent semantic indexing

# %%
# Remove the first LSI component (captures sequencing depth variation)
adata.obsm['X_lsi'] = adata.obsm['X_lsi'][:,1:]
adata.varm["LSI"] = adata.varm["LSI"][:,1:]
adata.uns["lsi"]["stdev"] = adata.uns["lsi"]["stdev"][1:]

# %%
adata.obs[['month', 'sample']] = adata.obs['sample_name'].str.extract(r'(m\d+)_(s\d+)')

# %%
# Apply Harmony on the LSI embedding
sc.external.pp.harmony_integrate(adata, key='month', basis="X_lsi")

# %%
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=30)  # Use harmony embeddings
sc.tl.umap(adata)

# %%
# check batch effects after harmony integration
sc.tl.embedding_density(adata, groupby="sample_name")
sc.pl.embedding_density(adata, basis="umap", key="umap_density_sample_name")

# %%
sc.tl.leiden(adata, resolution=0.3)

# %%
sc.pl.umap(adata, color=["leiden", "n_features_per_cell"], legend_loc="on data")

# %%
adata.write_h5ad(os.path.join(DATA_PATH, "scATAC_all.h5ad"))

# %% [markdown]
# # Gene activity matrix

# %%
adata_atac = sc.read_h5ad(os.path.join(DATA_PATH, "scATAC_all.h5ad"))

# %%
gene_intervals = pd.read_csv('/home/sychen9584/projects/cardio_paper/data/ref/gene_intervals.csv')
gene_intervals_filtered = gene_intervals[gene_intervals['Chromosome'].str.startswith('chr')]
gene_intervals_filtered = gene_intervals_filtered[gene_intervals_filtered['Chromosome'] != "chrM"]

# %%
for sample in scATAC_samples:
    logging.info(f"Counting fragments for {sample}----------------------------------")
    adata = sc.read_h5ad(os.path.join(DATA_PATH, "preprocessed", sample, "unified.h5ad"))
    ac.tl.locate_fragments(adata, fragments=os.path.join(DATA_PATH, 'raw', sample, 'fragments.tsv.gz'))
    adata_gene = ac.tl.count_fragments_features(adata, features=gene_intervals_filtered, stranded=True)
    adata_gene.write_h5ad(os.path.join(DATA_PATH, "preprocessed", sample, "gene_activity.h5ad"))

# %%
adatas = {}
for sample in scATAC_samples:
    adata = sc.read_h5ad(os.path.join(DATA_PATH, "preprocessed", sample, "gene_activity.h5ad"))
    adata.X = adata.X.astype(np.int32)
    adatas.update({sample: adata})

# %%
adata_gene = anndata.concat(adatas, label="sample_name", index_unique="_", merge="same")

# %%
# # copy some data over from the original adata object
adata_gene.var = adata_gene.var.set_index('gene_name')
adata_gene.var.index = adata_gene.var.index.astype(str)
adata_gene.var_names_make_unique()

# %%
adata_gene.uns['leiden_colors'] = adata_atac.uns['leiden_colors'].copy()
adata_gene.obsm['X_pca_harmony'] = adata_atac.obsm['X_pca_harmony'].copy()
adata_gene.obsm['X_umap'] = adata_atac.obsm['X_umap'].copy()
adata_gene.obs['leiden'] = adata_atac.obs['leiden'].copy()

# %%
adata_gene.obs['barcode'] = adata_gene.obs.index.astype(str)

# %%
adata_gene.write_h5ad(os.path.join(DATA_PATH, "scATAC_all_gene_activity.h5ad"), compression="gzip")

# %%
adata_gene.layers['count'] = adata_gene.X.copy()
# Normalizing to 10K counts
sc.pp.normalize_total(adata_gene, target_sum=10000)
# Logarithmize the data
sc.pp.log1p(adata_gene)

# %%
# we can use the functionality of the ATAC module in muon to color plots by cut values in peaks correspoonding to a certain gene
sc.pl.umap(adata_gene, color=["Gsn", "Lyz2", "Egfl7", 'leiden'])

# %%
sc.pl.violin(adata_gene, keys=["Gsn", "Ctss", "Egfl7"], groupby="leiden")

# %%
adata_gene.write_h5ad("scATAC_all_gene_activity.h5ad", compression="gzip")

# %%
