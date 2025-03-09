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
import pandas as pd

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data/"
FIGURE_PATH = "/home/sychen9584/projects/cardio_paper/figures"

# %%
motifs_mgi_fname = os.path.join(DATA_PATH, "ref/motifs-v9-nr.mgi-m0.001-o0.0.tbl")
out_tfs_mgi_fname = os.path.join(DATA_PATH, "ref/mm_mgi_tfs.txt")

# %%
df_motifs_mgi = pd.read_csv(motifs_mgi_fname, sep='\t')
mm_tfs = df_motifs_mgi.gene_name.unique()
with open(out_tfs_mgi_fname, 'wt') as f:
    f.write('\n'.join(mm_tfs) + '\n')
len(mm_tfs)

# %%
