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

# %% [markdown]
# # scATAC-seq data QC and processing

# %% [markdown]
# ## Using 3 months sample 1 as pilot

# %%
import os
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import anndata
import subprocess
import muon as mu
from muon import atac as ac  # scATAC-seq processing

import sys
sys.path.append('../scripts')
import scatac_preprocessing # replacement for muon implementation that accounts for strands

# Set plotting preferences
sc.settings.verbosity = 0
sns.set(rc={"figure.figsize": (4, 3.5), "figure.dpi": 100})
sns.set_style("whitegrid")

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/raw"
SAMPLE_NAME = "scATAC_m12_s2"

# %% [markdown]
# # Load in adata object

# %%
# Define the directory containing the ATAC-seq data
matrix_path = os.path.join(DATA_PATH, SAMPLE_NAME, "matrix.mtx.gz")

# Load the sparse peak-barcode matrix
mat = scipy.io.mmread(matrix_path).T.tocsc()  # Convert to CSC format for efficiency

# Load peak information (BED file)
peaks_path = os.path.join(DATA_PATH, SAMPLE_NAME, "peaks.bed.gz")
peaks = pd.read_csv(peaks_path, sep="\t", header=None, names=["chrom", "start", "end"])

# Load barcode metadata
barcodes_path = os.path.join(DATA_PATH, SAMPLE_NAME, "barcodes.tsv.gz")
barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, names=["barcode"])

# Convert barcodes and peaks into the required format
barcodes.index = barcodes["barcode"]  # Set barcodes as index
peaks["peak_name"] = peaks['chrom'] + ":" + peaks['start'].astype(str) + "-" + peaks['end'].astype(str)
peaks.index = peaks["peak_name"]  # Set peaks as index

# Create AnnData object for ATAC-seq counts
adata_atac = anndata.AnnData(
    X=mat,  # Sparse matrix
    obs=barcodes,  # Barcodes (cells)
    var=peaks  # Peaks (features)
)

# %% [markdown]
# ## Add fragment files to adata object

# %%
ac.tools.locate_fragments(adata_atac, fragments=os.path.join(DATA_PATH, SAMPLE_NAME, "fragments.tsv.gz"))

# %% [markdown]
# ## Peak annotations with HOMER

# %%
# Define paths
dir_path = os.path.join(DATA_PATH, SAMPLE_NAME)
bed_path = os.path.join(dir_path, "peaks.bed.gz")
bed_unzipped_path = os.path.join(dir_path, "peaks.bed")
annot_path = os.path.join(dir_path, "annotated_peaks.txt")

# %%
# Bash script with multiple commands
bash_script = f"""
gunzip -c {bed_path} > {bed_unzipped_path} && \
annotatePeaks.pl {bed_unzipped_path} mm10 > {annot_path} && \
echo "Annotation complete!"
"""

# Run the script
result = subprocess.run(bash_script, shell=True, text=True, capture_output=True)

# Check for errors
if result.returncode == 0:
    print("Annotation completed successfully!")
else:
    print("Error:", result.stderr)


# %%
# Load HOMER annotation output and wrangle into 10X tsv format
homer_df = pd.read_csv(annot_path, sep="\t")
homer_df = homer_df.iloc[:, [1, 2, 3, 14, 15, 9, 7]]
homer_df.columns = ['chrom', 'start', 'end', 'gene_id', 'gene', 'distance', 'peak_type']
homer_df['peak_type'] = homer_df['peak_type'].str.extractall(r'^(.*?)(?=\s\()').unstack()
homer_df['peak_type'] = homer_df['peak_type'].fillna('intergenic')
homer_df['peak_type'] = homer_df['peak_type'].replace({
    'intron': 'distal',
    'exon': 'distal', 
    'promoter-TSS': 'promoter',
    "5' UTR": 'distal', 'non-coding': 'distal',
    'TTS': 'distal', "3' UTR": 'distal'
})
homer_df['start'] = homer_df['start'] - 1

# %%
homer_df.to_csv(os.path.join(dir_path, "peak_annotations.tsv"), sep="\t", index=False)

# %%
ac.tools.add_peak_annotation(adata_atac, annotation=os.path.join(dir_path, "peak_annotations.tsv"))

# %%
adata_atac.write_h5ad(os.path.join(DATA_PATH, SAMPLE_NAME, "raw.h5ad"))

# %% [markdown]
# # Doublet Removal

# %%
adata_atac = sc.read_h5ad(os.path.join(DATA_PATH, SAMPLE_NAME, "raw.h5ad"))

# %%
# used a R script to run scDblFinder
dbl_scores = pd.read_csv(os.path.join(DATA_PATH, SAMPLE_NAME, "doublet_scores.csv")).set_index('barcode')
adata_atac.obs['dbl_score'] = dbl_scores['doublet_score']
adata_atac.obs['dbl_class'] = dbl_scores['doublet_class']

# %% [markdown]
# # Quality Control

# %%
# Calculate general qc metrics using scanpy
sc.pp.calculate_qc_metrics(adata_atac, percent_top=None, log1p=False, inplace=True)

# Rename columns
adata_atac.obs.rename(
    columns={
        "n_genes_by_counts": "n_features_per_cell",
        "total_counts": "total_fragment_counts",
    },
    inplace=True,
)

# %%
sc.pl.violin(adata_atac, ['total_fragment_counts', 'n_features_per_cell'], jitter=0.4, multi_panel=True)

# %% [markdown]
# ## Nucleosome Signal

# %%
# Plot nucleosome signal distribution across cells
ac.tl.nucleosome_signal(adata_atac, n=10e3* adata_atac.n_obs)
mu.pl.histogram(adata_atac, "nucleosome_signal")

# %%
# Add group labels for above and below the nucleosome signal threshold
nuc_signal_threshold = 4
adata_atac.obs["nuc_signal_filter"] = [
    "NS_FAIL" if ns > nuc_signal_threshold else "NS_PASS"
    for ns in adata_atac.obs["nucleosome_signal"]
]

# Print number cells not passing nucleosome signal threshold
adata_atac.obs["nuc_signal_filter"].value_counts()

# %% [markdown]
# ## TSS enrichment

# %%
gene_intervals = pd.read_csv('/home/sychen9584/projects/cardio_paper/data/ref/gene_intervals.csv')

# %%
tss = scatac_preprocessing.tss_enrichment(adata_atac, features=gene_intervals, n_tss=3000, random_state=666, extend_upstream=2000, extend_downstream=2000)

# %%
fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

p1 = sns.histplot(adata_atac.obs, x="tss_score", ax=axs[0])
p1.set_title("Full range")

p2 = sns.histplot(
    adata_atac.obs,
    x="tss_score",
    binrange=(0, adata_atac.obs["tss_score"].quantile(0.95)),
    ax=axs[1],
)
p2.axvline(x=3)
p2.set_title("Up to 95% percentile")

plt.suptitle("Distribution of the TSS score")

plt.tight_layout()
plt.show()

# %%
tss_threshold = 3
tss.obs["tss_filter"] = [
    "TSS_FAIL" if score < tss_threshold else "TSS_PASS"
    for score in adata_atac.obs["tss_score"]
]

# Print number cells not passing nucleosome signal threshold
tss.obs["tss_filter"].value_counts()

# %%
# Temporarily set different color palette
sns.set_palette(palette="Set1")
ac.pl.tss_enrichment(tss, color="tss_filter")
# reset color palette
sns.set_palette(palette="tab10")

# %% [markdown]
# ## Filtering cells

# %%
# log-transform total counts and add as column
adata_atac.obs["log_total_fragment_counts"] = np.log10(adata_atac.obs["total_fragment_counts"])

# %%
plot_tss_max = 20
count_cutoff_lower = 5000
lcount_cutoff_upper = 100000
tss_cutoff_lower = 4

# Scatter plot & histograms
g = sns.jointplot(
    data=adata_atac[(adata_atac.obs["tss_score"] < plot_tss_max)].obs,
    x="log_total_fragment_counts",
    y="tss_score",
    color="black",
    marker=".",
)
# Density plot including lines
g.plot_joint(sns.kdeplot, fill=True, cmap="Blues", zorder=1, alpha=0.75)
g.plot_joint(sns.kdeplot, color="black", zorder=2, alpha=0.75)

# Lines thresholds
plt.axvline(x=np.log10(count_cutoff_lower), c="red")
plt.axvline(x=np.log10(lcount_cutoff_upper), c="red")
plt.axhline(y=tss_cutoff_lower, c="red")

plt.show()

# %%
adata_atac.write_h5ad(os.path.join(DATA_PATH, SAMPLE_NAME, "qc_metrics.h5ad"))

# %%
adata_atac = sc.read_h5ad(os.path.join(DATA_PATH, SAMPLE_NAME, "qc_metrics.h5ad"))

# %% [markdown]
# ### Actual filtering

# %%
print(f"Total number of cells: {adata_atac.n_obs}")
mu.pp.filter_obs(
    adata_atac,
    "total_fragment_counts",
    lambda x: (x >= 5000) & (x <= 100000),
)
print(f"Number of cells after filtering on fragment_counts: {adata_atac.n_obs}")

# %%
mu.pp.filter_obs(adata_atac, "nucleosome_signal", lambda x: x <= 4)
print(f"Number of cells after filtering on nucleosome_signal: {adata_atac.n_obs}")

# %%
mu.pp.filter_obs(adata_atac, "dbl_class", lambda x: x == "singlet")
print(f"Number of cells after filtering doublets: {adata_atac.n_obs}")

# %%
mu.pp.filter_obs(adata_atac, "tss_score", lambda x: x >= 3)
print(f"Number of cells after on TSS enrichment score: {adata_atac.n_obs}")

# %% [markdown]
# # TF-IDF Normalization

# %%
print("Performing TF-IDF normalization...")
ac.pp.tfidf(adata_atac, scale_factor=1e4)

# %% [markdown]
# # Dimensionality Reduction using SVD (LSI components)

# %%
ac.tl.lsi(adata_atac) 

# %%
# We find the first component is typically associated with number of peaks or counts per cell so it is reasonable to remove it:
adata_atac.obsm['X_lsi'] = adata_atac.obsm['X_lsi'][:,1:]
adata_atac.varm["LSI"] = adata_atac.varm["LSI"][:,1:]
adata_atac.uns["lsi"]["stdev"] = adata_atac.uns["lsi"]["stdev"][1:]

# %%
svd_variance = adata_atac.uns["lsi"]["stdev"]

# %% [markdown]
# ## UMAP and clustering

# %%
sc.pp.neighbors(adata_atac, use_rep="X_lsi", n_neighbors=10, n_pcs=30)  # Use LSI embeddings for nearest-neighbor graph
sc.tl.umap(adata_atac)

# %%
sc.tl.leiden(adata_atac, resolution=0.8)

# %%
sc.pl.umap(adata_atac, color=["leiden", "n_features_per_cell"], legend_loc="on data")

# %%
ac.pl.umap(adata_atac, color=["Gsn", "Lyz2", "Egfl7"], average="total", use_raw=False)

# %% [markdown]
# # Gene activity matrix

# %%
gene_intervals = pd.read_csv('/home/sychen9584/projects/cardio_paper/data/ref/gene_intervals.csv')

# %%
gene_intervals_filtered = gene_intervals[gene_intervals['Chromosome'].str.startswith('chr')]
gene_intervals_filtered = gene_intervals_filtered[gene_intervals_filtered['Chromosome'] != "chrM"]

# %%
adata_atac_gene = ac.tl.count_fragments_features(adata_atac, features=gene_intervals_filtered, stranded=True)

# %%
# # copy some data over from the original adata object
adata_atac_gene.var = adata_atac_gene.var.set_index('gene_name')
adata_atac_gene.uns['leiden_colors'] = adata_atac.uns['leiden_colors'].copy()
adata_atac_gene.obsm['X_lsi'] = adata_atac.obsm['X_lsi'] .copy()
adata_atac_gene.obsm['X_umap'] = adata_atac.obsm['X_umap'] .copy()

# %%
adata_atac_gene.layers['count'] = adata_atac_gene.X.copy()
# Normalizing to 10K counts
sc.pp.normalize_total(adata_atac_gene, target_sum=10000)
# Logarithmize the data
sc.pp.log1p(adata_atac_gene)

# %%
# we can use the functionality of the ATAC module in muon to color plots by cut values in peaks correspoonding to a certain gene
sc.pl.umap(adata_atac_gene, color=["Gsn", "Lyz2", "Egfl7", 'leiden'])

# %%
adata_atac_gene.var.index = adata_atac_gene.var.index.astype(str)
adata_atac_gene.var_names_make_unique()

# %%
sc.pl.violin(adata_atac_gene, keys=["Gsn", "Lyz2", "Egfl7"], groupby="leiden")

# %%
adata_atac_gene.write_h5ad("gene_activity.h5ad")
