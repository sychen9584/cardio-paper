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
# # Figure 1

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import os
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import plotly.graph_objects as go
import plotly.io as pio

from PIL import Image
from io import BytesIO
from adjustText import adjust_text
from typing import Dict, List, Optional, Set
from mpl_toolkits.axes_grid1 import make_axes_locatable

import sys
sys.path.append('../../scripts')
import figure1_functions as fig1

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=120,
    facecolor="white",
    frameon=True,
)
# %matplotlib inline

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data"
FIGURE_PATH = "/home/sychen9584/projects/cardio_paper/figures"

# %% [markdown]
# ## Figure 1B/D

# %%
adata_rna = sc.read_h5ad(os.path.join(DATA_PATH, 'processed/scRNA_all.h5ad'))
adata_atac = sc.read_h5ad(os.path.join(DATA_PATH, 'processed/scATAC_scanvi_annot.h5ad'))

# rename cell_type_fine names to match visualiztion in the publication
adata_rna.obs['cell_type_fine'] = adata_rna.obs['cell_type_fine'].replace({
    'Granulocyte/Neutrophil': "Neutrophil",
    'Peri/Smooth Muscle': 'Smooth Muscle'
})

adata_atac.obs['cell_type_fine'] = adata_atac.obs['cell_type_fine'].replace({
    'Granulocyte/Neutrophil': "Neutrophil",
    'Peri/Smooth Muscle': 'Smooth Muscle'
})

# %%
# colors for umap plots
cell_type_colors, cell_type_fine_colors = fig1.get_celltype_colors()

# ordering of cell types and samples
cell_type_order = ['Fibroblast', 'Endothelial', 'Smooth Muscle', "Macrophage", "Neutrophil", "T Cell", "B Cell"]
cell_type_fine_order = ['Fib.1', 'Fib.2', 'Fib.3', 'Fib.4', 'Fib.5', 'Fib.6',
                        'Endo.1', "Endo.2", "Endo.3", "Endo.4", "Endo.6", "Endo.7",
                        'Smooth Muscle', 'MC.1', "MC.2", 'MC.3', 'MC/B-Cell', 'Neutrophil', 'T-Cell', "B-Cell"]

# %% [markdown]
# ## Figure 1C

# %%
river_df = adata_rna.obs[['cell_type', 'cell_type_fine']].groupby(['cell_type', 'cell_type_fine'], as_index=False, observed=True).size()
river_df['cell_type'] = pd.Categorical(river_df['cell_type'], categories=reversed(cell_type_order), ordered=True)
river_df = river_df.sort_values(['cell_type', 'size'], ascending=False)

# colors for river plot
unique_celltypes = adata_rna.obs['cell_type'].unique().tolist()
cmap = plt.get_cmap("Pastel1")
hex_colors = [mcolors.to_hex(cmap(i)) for i in range(cmap.N)][0:7]
river_colors = dict(zip(unique_celltypes, hex_colors))

# Create a combined list of unique labels
label_celltype = river_df['cell_type'].unique().tolist() 
label_celltype_fine = river_df['cell_type_fine'].unique().tolist()
labels = label_celltype + label_celltype_fine

# Map labels to indices
label_to_index1 = {label: i for i, label in enumerate(label_celltype)}
river_df['source_index'] = river_df['cell_type'].map(label_to_index1)
label_to_index2 = {label: i+7 for i, label in enumerate(label_celltype_fine)}
river_df['target_index'] = river_df['cell_type_fine'].map(label_to_index2)

# colors for nodes and links
river_df['color'] = river_df['cell_type'].map(river_colors)
sankey_color_dict = {**dict(zip(river_df['cell_type'], river_df['color'])), **dict(zip(river_df['cell_type_fine'], river_df['color']))}

# hide node labels for cell types with small populaiton
hidden_labels = {'Endo.3', 'Endo.7', 'MC.2', 'MC.3', 'T Cell', 'T-Cell', 'MC/B-Cell', 'B-Cell'}
node_labels = [
    "" if label in hidden_labels else label
    for label in labels
]

# %%
fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 0,
      thickness = 100,
      line = dict(color = "black", width = 2),
      label = node_labels,
      color = [sankey_color_dict.get(label) for label in labels]
    ),
    link = dict(
      source = river_df['source_index'],
      target = river_df['target_index'],
      value = river_df['size'],
      color = [sankey_color_dict.get(label) for label in label_celltype_fine]
  ))])

fig.update_layout(
    margin=dict(l=5, r=5, t=5, b=5),
    autosize=False,
    width=450,
    height=400,
    font=dict(size=11, weight="bold")
    
)

# Convert Plotly figure to an image (PNG format)
img_bytes = pio.to_image(fig, format="png")

# Load image as a PIL object
img = Image.open(BytesIO(img_bytes))

# %% [markdown]
# ## Figure 1E

# %%
month_order = ['m3', 'm12', 'm24']

cellnum_df = adata_rna.obs[['sample', 'cell_type']]
cellnum_df[['month', 'sample']] = cellnum_df['sample'].str.split("_", expand=True).iloc[:, -2:]
cellnum_df = cellnum_df[['month', 'cell_type']].value_counts().reset_index().sort_values(['month', 'cell_type'])

cellnum_df['cell_type'] = pd.Categorical(cellnum_df['cell_type'], categories=reversed(cell_type_order), ordered=True)
cellnum_df['month'] = pd.Categorical(cellnum_df['month'], categories=month_order, ordered=True)
cellnum_df = cellnum_df.sort_values(['month', 'cell_type'])

cellnum_df_pivot = cellnum_df.pivot(index="month", columns="cell_type", values="count")
coarse_colors = [cell_type_colors[cell_type] for cell_type in cellnum_df_pivot.columns]

# %% [markdown]
# ## Figure 1F

# %%

cellnum_fine_df = adata_rna.obs[['sample', 'cell_type_fine']]
cellnum_fine_df[['month', 'sample']] = cellnum_fine_df['sample'].str.split("_", expand=True).iloc[:, -2:]
cellnum_fine_df = cellnum_fine_df[['month', 'cell_type_fine']].groupby('month').value_counts(normalize=True).reset_index().sort_values(['month', 'cell_type_fine'])

cellnum_fine_df['cell_type_fine'] = pd.Categorical(cellnum_fine_df['cell_type_fine'], categories=reversed(cell_type_fine_order), ordered=True)
cellnum_fine_df['month'] = pd.Categorical(cellnum_fine_df['month'], categories=month_order, ordered=True)
cellnum_fine_df = cellnum_fine_df.sort_values(['month', 'cell_type_fine'])

cellnum_fine_df_pivot = cellnum_fine_df.pivot(index="month", columns="cell_type_fine", values="proportion")
fine_colors = [cell_type_fine_colors[cell_type] for cell_type in cellnum_fine_df_pivot.columns]

# %% [markdown]
# ## Figure 1G

# %%
augur_cmap = fig1.get_augur_colors()

augur_celltype = pd.read_csv(os.path.join(DATA_PATH, "augur/augur_cell_type_rfc.csv"), index_col=0)
augur_celltype_map = dict(augur_celltype.loc['mean_auc'])

augur_celltype_fine = pd.read_csv(os.path.join(DATA_PATH, "augur/augur_cell_type_fine_rfc.csv"), index_col=0)
augur_celltype_fine_map = dict(augur_celltype_fine.loc['mean_auc'])

adata_rna.obs['augur_celltype'] = adata_rna.obs['cell_type'].map(augur_celltype_map)
adata_rna.obs['augur_celltype_fine'] = adata_rna.obs['cell_type_fine'].map(augur_celltype_fine_map).rank(pct=True)

# %% [markdown]
# ## Figure 1H

# %%
augur_3_v_12 = pd.read_csv(os.path.join(DATA_PATH, "augur/augur_cell_type_m3_v_m12_rfc.csv"), index_col=0)
augur_12_v_24 = pd.read_csv(os.path.join(DATA_PATH, "augur/augur_cell_type_m12_v_m24_rfc.csv"), index_col=0)

augur_df = pd.DataFrame({'m3_v_m12': augur_3_v_12.loc['mean_auc'],
                         'm12_v_m24': augur_12_v_24.loc['mean_auc']})
augur_df['delta'] = augur_df['m3_v_m12'] - augur_df['m12_v_m24']
augur_df.index.name = "cell_type"
augur_df.reset_index(inplace=True)

# %% [markdown]
# ## Main figure

# %%
with plt.rc_context({"figure.figsize": (16, 14), "figure.dpi": 150, "figure.frameon": True}):
    
    fig = plt.figure()
    gs0 = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[1, 1, 1])
    
    # first row
    gs00 = gs0[0].subgridspec(1, 3)  # 1 row, 3 columns
    
    # ax1: scRNA-seq UMAP
    ax1 = fig.add_subplot(gs00[0, 0])
    sc.pl.umap(adata_rna, color='cell_type_fine', show=False, title="", legend_loc=None, palette=cell_type_fine_colors, ax=ax1)
    fig1.repel_umap_labels(
        adata_rna,
        "cell_type_fine",
        ax=ax1,
        adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
        text_kwargs=dict(fontsize=11, weight='bold')
    )
    
    # ax2: river plot of cell types and their subtypes
    ax2 = fig.add_subplot(gs00[0, 1])
    ax2.axis("off")
    ax2.set_adjustable("box")  # Makes sure the image scales properly
    ax2.imshow(img, aspect="auto", extent=ax2.get_position().extents)  # Forces it to fit tighter
    
    # ax3: scATAC-seq UMAP
    ax3 = fig.add_subplot(gs00[0, 2])
    sc.pl.umap(adata_atac, color='cell_type_fine', show=False, title="", legend_loc=None, palette=cell_type_fine_colors, ax=ax3)
    fig1.repel_umap_labels(
        adata_atac,
        "cell_type_fine",
        ax=ax3,
        adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
        text_kwargs=dict(fontsize=11, weight='bold')
    )
    
    # second row
    gs01 = gs0[1].subgridspec(1, 4, width_ratios =[1.5, 0.5, 1.5, 0.5])
    
    # ax4: bar plot of cell type proportions
    ax4 = fig.add_subplot(gs01[0, 0])
    cellnum_df_pivot.plot(kind="barh", stacked=True, color=coarse_colors, ax=ax4)
    ax4.tick_params(axis="x", labelsize=10)
    ax4.set_yticklabels(['3 month', '12 month', '24 month'], fontsize=10)
    ax4.set_xlabel("number of cells", fontsize=10)
    ax4.set_ylabel("")
    
    handles, labels = ax4.get_legend_handles_labels()
    ax4.legend(handles[::-1], labels[::-1], title="Cell Type", loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, title_fontsize=8, frameon=False)
    fig1.label_bars(ax4, cellnum_df_pivot, celltypes={'Fibroblast', "Endothelial", "Macrophage"})
    
    # ax5: bar plot of cell type fine proportions
    ax5 = fig.add_subplot(gs01[0, 2])
    cellnum_fine_df_pivot.plot(kind="barh", stacked=True, color=fine_colors, ax=ax5)
    ax5.set_xticklabels([])
    ax5.set_yticklabels(['3 month', '12 month', '24 month'], fontsize=10)
    ax5.set_xlabel("percentage of cells", fontsize=10)
    ax5.set_ylabel("")
    
    handles, labels = ax5.get_legend_handles_labels()
    ax5.legend(handles[::-1], labels[::-1], title="Cell Type", loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, title_fontsize=8, frameon=False)
    
    # third row
    gs02 = gs0[2].subgridspec(1, 4, width_ratios=[1, 1, 0, 1]) 
    
    # ax6: UMAP of cell types colored by Augur scores
    ax6 = fig.add_subplot(gs02[0, 0])
    sc.pl.umap(adata_rna, color='augur_celltype', show=False, title="", legend_loc=None, color_map=augur_cmap, colorbar_loc=None, ax=ax6)
    fig1.repel_umap_labels(
        adata_rna,
        "cell_type",
        include=["Macrophage", "Neutrophil", "Fibroblast"],
        ax=ax6,
        adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
        text_kwargs=dict(fontsize=11, weight='bold')
    )
    fig1.augur_colorbar(ax6, "AUC", label_fontsize=10, tick_fontsize=8, pad_size=0.1, size="1.5%")
    
    # ax7: UMAP of cell type fine colored by Augur scores
    ax7 = fig.add_subplot(gs02[0, 1])
    sc.pl.umap(adata_rna, color='augur_celltype_fine', show=False, title="", legend_loc=None, color_map=augur_cmap, colorbar_loc=None, ax=ax7)
    fig1.repel_umap_labels(
        adata_rna,
        "cell_type_fine",
        include=['Fib.1', 'MC.1', "MC.2", "MC/B-Cell", 'Neutrophil'],
        ax=ax7,
        adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
        text_kwargs=dict(fontsize=11, weight='bold')
    )
    fig1.augur_colorbar(ax7, "Rank (%)", label_fontsize=10, tick_fontsize=8, pad_size=0.1, size="1.5%")
    
    # ax8: scatter plot of AUGUR perturbation scores between 3-12 months and 12-24 months
    ax8 = fig.add_subplot(gs02[0, 3])
    sns.scatterplot(
        data=augur_df,
        x="m12_v_m24",
        y="m3_v_m12",
        hue="delta",
        hue_norm=(augur_df["delta"].min(), augur_df["delta"].max()),  # Normalize color scale
        palette=augur_cmap, legend=False, s=75, ax=ax8
    )
    texts = []
    for i, row in augur_df.iterrows():
        texts.append(
            ax8.text(row["m12_v_m24"], row["m3_v_m12"], row["cell_type"], fontsize=8, color="black", ha="right", va="bottom")
        )
    adjust_text(texts, ax=ax8)
    
    x_vals = np.linspace(augur_df["m12_v_m24"].min(), augur_df["m12_v_m24"].max(), 100)  # Get min-max range
    ax8.plot(x_vals, x_vals, linestyle="--", color="grey", linewidth=1)# diagonal line
    
    ax8.tick_params(labelsize=10)
    ax8.set_xlabel("AUC 12 months vs 24 months", fontsize=10)
    ax8.set_ylabel("AUC 3 months vs 12 months", fontsize=10)
    
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax8]
    for i, ax in enumerate(axes, 1):
        ax.text(-0.1, 1.1, f"{chr(65+i)}", transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='left') # subfigure labels
        ax.grid(False)

    plt.suptitle("Figure 1. Integrated single-cell transcriptomic and epigenomic analysis of non-cardiac (non-CM) cells during aging")
    fig.subplots_adjust(wspace=0.3)
    fig.tight_layout()
    plt.show()

# %%
fig.savefig(os.path.join(FIGURE_PATH, "figure1.png"), bbox_inches='tight', dpi=300)

# %%
