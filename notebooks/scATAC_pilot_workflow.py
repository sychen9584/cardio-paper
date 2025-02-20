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
# Load R-Python bridge FIRST
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
import anndata2ri

# Activate conversion BEFORE loading other libraries
pandas2ri.activate()
#anndata2ri.activate()

# Load Jupyter R magic extension
# %load_ext rpy2.ipython

# Now import other packages
import os
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

# Set plotting preferences
sc.settings.verbosity = 0
sns.set(rc={"figure.figsize": (4, 3.5), "figure.dpi": 100})
sns.set_style("whitegrid")

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/raw"
SAMPLE_NAME = "scATAC_m3_s1"

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
peaks["peak_id"] = ["peak_" + str(i) for i in range(len(peaks))]  # Generate unique peak IDs
peaks.index = peaks["peak_id"]  # Set peaks as index

# Create AnnData object for ATAC-seq counts
adata_atac = anndata.AnnData(
    X=mat,  # Sparse matrix
    obs=barcodes,  # Barcodes (cells)
    var=peaks  # Peaks (features)
)

# %%
adata_atac

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
homer_df = homer_df.iloc[:, [1, 2, 3, 15, 9, 7]]
homer_df.columns = ['chrom', 'start', 'end', 'gene', 'distance', 'peak_type']
homer_df['peak_type'] = homer_df['peak_type'].str.extractall(r'^(.*?)(?=\s\()').unstack()
homer_df['peak_type'] = homer_df['peak_type'].fillna('intergenic')
homer_df['peak_type'] = homer_df['peak_type'].replace({
    'intron': 'distal',
    'exon': 'distal', 
    'promoter-TSS': 'promoter',
    "5' UTR": 'distal', 'non-coding': 'distal',
    'TTS': 'distal', "3' UTR": 'distal'
})

# %%
homer_df.to_csv(os.path.join(dir_path, "peak_annotations.tsv"), sep="\t", index=False)

# %%
ac.tools.add_peak_annotation(adata_atac, annotation=os.path.join(dir_path, "peak_annotations.tsv"))

# %%
adata_atac.write_h5ad(os.path.join(DATA_PATH, SAMPLE_NAME, "raw.h5ad"))

# %% [markdown]
# # Doublet Removal

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
gene_intervals

# %%
gene_intervals_positive = gene_intervals.query('Strand == "-"')

# %%
tss = ac.tl.tss_enrichment(adata_atac, features=gene_intervals_positive, n_tss=3000, random_state=666, extend_upstream=2000, extend_downstream=2000)

# %%
fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

p1 = sns.histplot(
    adata_atac.obs,
    x="tss_score",
    binrange=(0, adata_atac.obs["tss_score"].quantile(0.95)),
    bins=100, 
    ax=axs[0],
)
p1.set_title("Up to 90% percentile")

p2 = sns.histplot(
    adata_atac.obs,
    x="tss_score",
    binrange=(0, adata_atac.obs["tss_score"].quantile(0.75)),
    bins=100,
    ax=axs[1],
)
p2.set_title("Up to 50% percentile")

plt.suptitle("Distribution of the TSS score")

plt.tight_layout()
plt.show()

# %%
tss_threshold = 2
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
plot_tss_max = 1
count_cutoff_lower = 10000
lcount_cutoff_upper = 60000
tss_cutoff_lower = 0.03

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
    lambda x: (x >= 10000) & (x <= 60000),
)
print(f"Number of cells after filtering on fragment_counts: {adata_atac.n_obs}")

# %%
mu.pp.filter_obs(adata_atac, "nucleosome_signal", lambda x: x <= 2)
print(f"Number of cells after filtering on nucleosome_signal: {adata_atac.n_obs}")

# %%
mu.pp.filter_obs(adata_atac, "dbl_class", lambda x: x == "singlet")
print(f"Number of cells after filtering doublets: {adata_atac.n_obs}")

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

# %%
# we can use the functionality of the ATAC module in muon to color plots by cut values in peaks correspoonding to a certain gene
ac.pl.umap(adata_atac, color=["KLF4"], average="peak_type")

# %%
adata_atac.var

# %%
