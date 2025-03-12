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
adata_subset = adata[adata.obs['cell_type'].isin(['Endothelial', 'Fibroblast', 'Macrophage'])].copy()

# %%
del adata, adata_annot

# %%
sc.tl.rank_genes_groups(adata_subset, groupby="cell_type", method="wilcoxon")

# %%
dar_df = sc.get.rank_genes_groups_df(adata_subset, group=None)
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


# %%
### venn diagram of DARs

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

# %%

# %%
