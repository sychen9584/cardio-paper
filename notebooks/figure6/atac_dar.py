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
sys.path.append('../../scripts')
import figure_functions as fig_func

# Set plotting preferences
sc.settings.verbosity = 0
sns.set(rc={"figure.figsize": (4, 3.5), "figure.dpi": 100})
sns.set_style("whitegrid")

# %% [markdown]
# # load in scATAC adata object with peaks

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/"
FIGURE_PATH = "/home/sychen9584/projects/cardio_paper/figures"

# %%
adata_annot = sc.read_h5ad(os.path.join(DATA_PATH, "processed/scATAC_scanvi_annot.h5ad"))
adata = sc.read_h5ad(os.path.join(DATA_PATH, "processed/scATAC_all.h5ad"))

# %%
adata.obs['cell_type_fine'] = adata_annot.obs['cell_type_fine'].copy()

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
adata.write_h5ad(os.path.join(DATA_PATH, "processed/scATAC_all.h5ad"))

# %%
sc.tl.rank_genes_groups(adata, groupby="cell_type", method="wilcoxon")

# %%
dar_df = sc.get.rank_genes_groups_df(adata, group=None)
dar_df = dar_df.query("logfoldchanges > 1 and pvals_adj < 0.05").copy()

# %%
dar_df.group.value_counts()

# %%
dar_df.to_csv(os.path.join(DATA_PATH, "deg/atac_dar.csv"), index=False)

# %% [markdown]
# ## Peak annotation bar plots

# %%
dar_df = pd.read_csv(os.path.join(DATA_PATH, "deg/atac_dar.csv"))

# %%
# Homer annotation on unified peaks
bed_unzipped_path = os.path.join(DATA_PATH, "processed/unified_peaks.bed")
annot_path = os.path.join(DATA_PATH, "processed/unified_raw_annotated_peaks.txt")

# %%
import subprocess

# Bash script with multiple commands
bash_script = f"""
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
homer_df = pd.read_csv(annot_path, sep="\t")
homer_df = homer_df.iloc[:, [1, 2, 3, 14, 15, 9, 7]]
homer_df.columns = ['chrom', 'start', 'end', 'gene_id', 'gene', 'distance', 'peak_type']
homer_df['peak_type'] = homer_df['peak_type'].str.extractall(r'^(.*?)(?=\s\()').unstack()
homer_df['peak_type'] = homer_df['peak_type'].fillna('intergenic')
homer_df['start'] = (homer_df['start'] - 1).astype(str)
homer_df['end'] = homer_df['end'].astype(str)
homer_df.dropna(inplace=True)

# %%
dar_df[['chrom', 'start_end']] = dar_df['names'].str.split(':', expand=True)
dar_df[['start', 'end']] = dar_df['start_end'].str.split('-', expand=True)

# Dropping the temporary column
dar_df.drop(columns=['start_end'], inplace=True)

# %%
dar_df = dar_df.merge(homer_df, on=['chrom', 'start', 'end'])

# %%
dar_df.to_csv(os.path.join(DATA_PATH, "deg/atac_dar_annotated.csv"), index=False)

# %%
fibroblast = dar_df.query("group == 'Fibroblast'").copy()
endothelial = dar_df.query("group == 'Endothelial'").copy()
macrophage = dar_df.query("group == 'Macrophage'").copy()

# %%
fibroblast = pd.DataFrame(100*fibroblast.value_counts("peak_type") / fibroblast.shape[0]).reset_index().sort_values("peak_type", ascending=False)
endothelial = pd.DataFrame(100*endothelial.value_counts("peak_type") / endothelial.shape[0]).reset_index().sort_values("peak_type", ascending=False)
macrophage = pd.DataFrame(100*macrophage.value_counts("peak_type") / macrophage.shape[0]).reset_index().sort_values("peak_type", ascending=False)

# %%
celltypes = ['Fibroblast', 'Endothelial', 'Macrophage']
palette = sns.color_palette('Set2')

# %%
fig, axs = plt.subplots(3, 1, figsize=(4, 12))

for i, (data, celltype) in enumerate(zip([fibroblast, endothelial, macrophage], celltypes)):
    wedges, _ = axs[i].pie(data['count'], startangle=140, colors=palette)  # No labels
    legend_labels = data['peak_type'] + " (" + data['count'].round(2).astype(str) + "%)"
    axs[i].legend(wedges, legend_labels, title=celltype, loc="upper right", bbox_to_anchor=(1.6, 0.95), frameon=False, fontsize=10)


# %% [markdown]
# ### venn diagram of DARs

# %%
dar_df

# %%
macrophage_dar = set(dar_df[dar_df['group'] == 'Macrophage']['names'].values)
fibroblast_dar = set(dar_df[dar_df['group'] == 'Fibroblast']['names'].values)
endothelial_dar = set(dar_df[dar_df['group'] == 'Endothelial']['names'].values)

# %%
fig_func.venn3_custom(fibroblast_dar, endothelial_dar, macrophage_dar, labels=('Fibroblast', 'Endothelial', 'Macrophage'), 
                              normalize_range=(0, 30000), title='')
plt.show()

# %% [markdown]
# ## Figure 6c

# %%
# Create bed files for each cell type and background to input into HOMER

# Background peaks
background_peaks = adata.var[['Chromosome', 'Start', 'End', 'peak_id']].reset_index(drop=True).copy()
background_bed = os.path.join(DATA_PATH, "homer_motifs/input/background.bed")
background_peaks.to_csv(background_bed, sep="\t", header=False, index=False)

# %%
# input peaks
input_peaks = dar_df.query('group in ["Endothelial", "Fibroblast", "Macrophage"]')[['chrom', 'start', 'end', 'names']].drop_duplicates().copy()
input_bed = os.path.join(DATA_PATH, "homer_motifs/input/input.bed")
input_peaks.to_csv(input_bed, sep="\t", header=False, index=False)

# %%
# visualize svg motif logos
Rfx6 = os.path.join(DATA_PATH, "homer_motifs/knownResults/known114.logo.svg")
Irf8 = os.path.join(DATA_PATH, "homer_motifs/knownResults/known36.logo.svg")
Stat4 = os.path.join(DATA_PATH, "homer_motifs/knownResults/known70.logo.svg")
Atf3 = os.path.join(DATA_PATH, "homer_motifs/knownResults/known1.logo.svg")

# %%
import cairosvg
cairosvg.svg2png(url=Rfx6, write_to=os.path.join(DATA_PATH, "homer_motifs/rfx6.png"))
cairosvg.svg2png(url=Irf8, write_to=os.path.join(DATA_PATH, "homer_motifs/irf8.png"))
cairosvg.svg2png(url=Stat4, write_to=os.path.join(DATA_PATH, "homer_motifs/stat4.png"))
cairosvg.svg2png(url=Atf3, write_to=os.path.join(DATA_PATH, "homer_motifs/atf3.png"))

# %% [markdown]
# ## Try pychromVAR

# %%
from pyjaspar import jaspardb
import pychromvar as pc

# %%
pc.get_genome("mm10", output_dir=os.path.join(DATA_PATH, "ref"))

# %%
valid_chroms = [f"chr{i}" for i in range(1, 20)] + ['chrX', 'chrY']
adata = adata[:, adata.var['Chromosome'].isin(valid_chroms)].copy()

# %%
pc.add_peak_seq(adata, genome_file=os.path.join(DATA_PATH, "ref/mm10.fa"), delimiter="[:-]")

# %%
#  estimate GC bias per peak and get the backgrounds
pc.add_gc_bias(adata)

# %%
adata.X = adata.layers["counts"].copy()

# %%
# Sum of reads per peak
reads_per_peak = np.sum(adata.X, axis=0)
# Find peaks with zero reads
zero_read_peaks = np.where(reads_per_peak == 0)[1]  # Get indices
print(f"Number of peaks with zero reads: {len(zero_read_peaks)}")

# Keep only peaks with nonzero reads
nonzero_peaks = np.where(reads_per_peak > 0)[1]  # Indices of nonzero peaks
adata = adata[:, nonzero_peaks]  # Subset AnnData object

# %%
pc.get_bg_peaks(adata)

# %%
adata.varm['bg_peaks'].shape

# %%
# extract TF motifs and perform motif matching to identify TF binding sites
jdb_obj = jaspardb(release='JASPAR2020')
motifs = jdb_obj.fetch_motifs(
    collection = 'CORE',
    tax_group = ['vertebrates'])

pc.match_motif(adata, motifs=motifs)

# %%
dev = pc.compute_deviations(adata)
dev

# %%
adata = sc.read_h5ad(os.path.join(DATA_PATH, "processed/scATAC_all.h5ad"))
chrom_deviations = sc.read_h5ad(os.path.join(DATA_PATH, "processed/chromvar_deviations.h5ad"))

# %%
chrom_deviations

# %%
df_dev = pd.DataFrame(chrom_deviations.X,
                             columns=chrom_deviations.var_names,
                             index=chrom_deviations.obs_names)

# %%
df_dev

# %%
df_dev['cell_type'] = adata.obs['cell_type'].values
df_dev['month'] = adata.obs['month'].values

# %%
df_dev_long = df_dev.melt(id_vars=['cell_type', 'month'], var_name='motif', value_name='deviation')
df_dev_long = df_dev_long.groupby(['cell_type', 'month', 'motif']).mean().reset_index()

# %%
df_dev_long['motif'] = df_dev_long['motif'].str.split('.').str[2]

# %%
motifs_to_plot = ['ATF3', "ATF6", "CEBPB", "CREB3L2", "CREM", "DBP", "E2F3", "ELF2", "ELK1", "ELK4", "ETS2", "ETV4", "EVX1",
                  "FOXO4", "GMEB2", "HES1", "HMBOX1", "HLF", "IRF1", "IRF2", "IRF3", "IRF7", "IRF9", "KLF13", "MAFG:NFE2L1",
                  "KLF4", "MXI1", "NFIA", "NFIB", "NFKB2", "NR1D1", "NRF1", "POU2F1", "RARB", "STAT1:STAT2", "STAT1", "STAT2",
                  "TCF12", "TCF7L1", "TEAD1", "THAP1", "XBP1", "YY1", "ZBTB7A"]
df_dev_long = df_dev_long.query("motif in @motifs_to_plot and cell_type in ['Fibroblast', 'Endothelial', 'Macrophage']").sort_values('motif').copy()

# %%
df_dev_long['cell_type'] = df_dev_long['cell_type'].astype(str)
df_dev_long['cell_type'] = pd.Categorical(df_dev_long['cell_type'], categories=['Endothelial', "Fibroblast", 'Macrophage'], ordered=True)

# %%
month_colors = {"m3": "#98df8a", "m12":"#FFED6F", "m24":"#ff9896"}

# %%
g = sns.catplot(
    data=df_dev_long,
    x="motif", 
    y="deviation", 
    hue="month",
    palette=month_colors,
    row="cell_type", 
    kind="bar",
    height=6, 
    aspect=3.2, 
    dodge=True  # Enables dodging for hue
)

g.set_titles("{row_name}", fontweight='bold', fontsize=14)
g.set_axis_labels("Motif", "")
g.fig.text(0.01, 0.5, "Average Motif Activity", va='center', rotation=90, fontsize=12)

for ax in g.axes.flat:
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

# %%
df_dev_long.to_csv(os.path.join(DATA_PATH, "processed/motif_activity.csv"), index=False)

# %%
