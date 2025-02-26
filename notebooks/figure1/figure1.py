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

# %%
adata_rna = sc.read_h5ad(os.path.join(DATA_PATH, 'processed/scRNA_all.h5ad'))
adata_atac = sc.read_h5ad(os.path.join(DATA_PATH, 'processed/scATAC_scanvi_annot.h5ad'))

# %%
# rename cell_type_fine names to match visualiztion in the publication
adata_rna.obs['cell_type_fine'] = adata_rna.obs['cell_type_fine'].replace({
    'Granulocyte/Neutrophil': "Neutrophil",
    'Peri/Smooth Muscle': 'Smooth Muscle'
})

adata_atac.obs['cell_type_fine'] = adata_atac.obs['cell_type_fine'].replace({
    'Granulocyte/Neutrophil': "Neutrophil",
    'Peri/Smooth Muscle': 'Smooth Muscle'
})

# colors of cell type on UMAP
cell_type_fine_colors = {
# B-Cell
"B-Cell" : "#1f77b4",

# Endo types in distinct blue shades
"Endo.1" : "#2171b5", 
"Endo.2" : "#4292c6", 
"Endo.3" : "#6baed6", 
"Endo.4" : "#9ecae1", 
"Endo.5" : "#c6dbef", 
"Endo.6" : "#deebf7", 
"Endo.7" : "#08306b", 

# Fib types in distinct orange shades
"Fib.1" : "#e6550d", 
"Fib.2" : "#fd8d3c", 
"Fib.3" : "#fdae6b", 
"Fib.4" : "#fdd0a2", 
"Fib.5" : "#feedde", 
"Fib.6" : "#a63603", 

# MC types in similar colors
"MC.1" : "#31a354", 
"MC.2" : "#74c476", 
"MC.3" : "#a1d99b", 
"MC/B-Cell" : "#c7e9c0", 

# Others
"Neutrophil" : "#9467bd",
"Smooth Muscle" : "#17becf",
"T-Cell" : "#aec7e8"
}

cell_type_colors = {
"Fibroblast": "#fdae6b",
"Endothelial": "#9ecae1",
"Smooth Muscle": "#17becf",
"Macrophage": "#a1d99b",
"Neutrophil": "#9467bd",
"T Cell": "#aec7e8",
"B Cell": "#1f77b4"
}

unique_celltypes = adata_rna.obs['cell_type'].unique().tolist()
cmap = plt.get_cmap("Pastel1")
hex_colors = [mcolors.to_hex(cmap(i)) for i in range(cmap.N)][0:7]
river_colors = dict(zip(unique_celltypes, hex_colors))

adata_rna.uns['cell_type_fine_colors'] = [cell_type_fine_colors[c] for c in adata_rna.obs["cell_type_fine"].cat.categories]
adata_rna.uns['cell_type_colors'] = [cell_type_colors[c] for c in adata_rna.obs["cell_type"].cat.categories]

adata_atac.uns['cell_type_fine_colors'] = [cell_type_fine_colors[c] for c in adata_rna.obs["cell_type_fine"].cat.categories]
adata_atac.uns['cell_type_colors'] = [cell_type_colors[c] for c in adata_rna.obs["cell_type"].cat.categories]

cell_type_order = ['Fibroblast', 'Endothelial', 'Smooth Muscle', "Macrophage", "Neutrophil", "T Cell", "B Cell"]
cell_type_fine_order = ['Fib.1', 'Fib.2', 'Fib.3', 'Fib.4', 'Fib.5', 'Fib.6',
                        'Endo.1', "Endo.2", "Endo.3", "Endo.4", "Endo.6", "Endo.7",
                        'Smooth Muscle', 'MC.1', "MC.2", 'MC.3', 'MC/B-Cell',
                        'Neutrophil', 'T-Cell', "B-Cell"]
month_order = ['m3', 'm12', 'm24']


# %%
# figure 1b
def repel_umap_labels(
    adata, groupby, include=None, ax=None, adjust_kwargs=None, text_kwargs=None
):
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}
    
    if include is None:
        include = adata.obs[groupby].unique()

    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in include:
            medians[g] = np.median(adata[g_idx].obsm["X_umap"], axis=0)

    if ax is None:
        texts = [
            plt.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()
        ]
    else:
        texts = [ax.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()]

    adjust_text(texts, **adjust_kwargs)


# %%
# figure 1c
river_df = adata_rna.obs[['cell_type', 'cell_type_fine']].groupby(['cell_type', 'cell_type_fine'], as_index=False, observed=True).size()
river_df['cell_type'] = pd.Categorical(river_df['cell_type'], categories=reversed(cell_type_order), ordered=True)
river_df = river_df.sort_values(['cell_type', 'size'], ascending=False)

# %%
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

# %%
# figure 1e
cellnum_df = adata_rna.obs[['sample', 'cell_type']]
cellnum_df[['month', 'sample']] = cellnum_df['sample'].str.split("_", expand=True).iloc[:, -2:]
cellnum_df = cellnum_df[['month', 'cell_type']].value_counts().reset_index().sort_values(['month', 'cell_type'])

cellnum_df['cell_type'] = pd.Categorical(cellnum_df['cell_type'], categories=reversed(cell_type_order), ordered=True)
cellnum_df['month'] = pd.Categorical(cellnum_df['month'], categories=month_order, ordered=True)
cellnum_df = cellnum_df.sort_values(['month', 'cell_type'])

cellnum_df_pivot = cellnum_df.pivot(index="month", columns="cell_type", values="count")
coarse_colors = [cell_type_colors[cell_type] for cell_type in cellnum_df_pivot.columns]


# %%
def label_bars(ax: plt.Axes, df: pd.DataFrame, celltypes: Set[str], text_kwargs: Optional[Dict[str, str]] = None):
    """
    Annotates selected bars in a stacked bar chart.

    Parameters:
    - ax: Matplotlib axis object
    - df: DataFrame used for plotting (stacked bar format)
    - celltypes: List or set of cell types to annotate
    - text_kwargs: Dictionary for text formatting (default: bold, center-aligned)
    """

    # Default text styling if none provided
    if text_kwargs is None:
        text_kwargs = {
            'ha': 'center', 'va': 'center', 'fontsize': 10,
            'fontweight': 'bold', "color": 'black'
        }

    y_labels = list(df.index)  # Ensure proper indexing for categories
    group_labels = list(df.columns)  # Stacked categories (groups)

    # Dictionary to track bars' y-coordinates
    bar_positions = {}  
    for bar in ax.patches:
        y_pos = round(bar.get_y(), 2)  # Track bar positions
        if y_pos not in bar_positions:
            bar_positions[y_pos] = []
        bar_positions[y_pos].append(bar)

    # Correctly map bars to their respective categories
    for y_index, (y_pos, bars) in enumerate(bar_positions.items()):
        for bar, group in zip(bars, group_labels):
            if group in celltypes:  # Only annotate selected groups
                value = int(bar.get_width())  # Get bar width (cell count)

                height = bar.get_y() + bar.get_height() / 2  # Center vertically
                width = bar.get_x() + bar.get_width() / 2  # Center horizontally
                
                ax.text(width, height, f"{value}", **text_kwargs)


# %%
# figure 1f
cellnum_fine_df = adata_rna.obs[['sample', 'cell_type_fine']]
cellnum_fine_df[['month', 'sample']] = cellnum_fine_df['sample'].str.split("_", expand=True).iloc[:, -2:]
cellnum_fine_df = cellnum_fine_df[['month', 'cell_type_fine']].groupby('month').value_counts(normalize=True).reset_index().sort_values(['month', 'cell_type_fine'])

cellnum_fine_df['cell_type_fine'] = pd.Categorical(cellnum_fine_df['cell_type_fine'], categories=reversed(cell_type_fine_order), ordered=True)
cellnum_fine_df['month'] = pd.Categorical(cellnum_fine_df['month'], categories=month_order, ordered=True)
cellnum_fine_df = cellnum_fine_df.sort_values(['month', 'cell_type_fine'])

cellnum_fine_df_pivot = cellnum_fine_df.pivot(index="month", columns="cell_type_fine", values="proportion")
fine_colors = [cell_type_fine_colors[cell_type] for cell_type in cellnum_fine_df_pivot.columns]

# %%
# figure 1g
augur_colors = ["#27408B",  # royalblue4
          "#FFF8DC",  # cornsilk1
          "#FF8C00",  # DarkOrange
          "#CD4F39"]  # tomato3

# Create a custom colormap
augur_cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", augur_colors)

# %%
augur_celltype = pd.read_csv(os.path.join(DATA_PATH, "augur/augur_cell_type_rfc.csv"), index_col=0)
augur_celltype_map = dict(augur_celltype.loc['mean_auc'])

augur_celltype_fine = pd.read_csv(os.path.join(DATA_PATH, "augur/augur_cell_type_fine_rfc.csv"), index_col=0)
augur_celltype_fine_map = dict(augur_celltype_fine.loc['mean_auc'])

adata_rna.obs['augur_celltype'] = adata_rna.obs['cell_type'].map(augur_celltype_map)
adata_rna.obs['augur_celltype_fine'] = adata_rna.obs['cell_type_fine'].map(augur_celltype_fine_map).rank(pct=True)


# %%
def augur_colorbar(ax, label, label_fontsize=10, tick_fontsize=8, pad_size=0.1, size="3%"):
    """
    Adds a colorbar to a given axis, with customizable label, font sizes, and padding.

    Parameters:
    - ax (matplotlib.axes.Axes): The axis to which the colorbar is attached.
    - label (str): The text label for the colorbar.
    - label_fontsize (int): Font size for the colorbar label.
    - tick_fontsize (int): Font size for the colorbar ticks.
    - pad_size (float): Padding between the plot and the colorbar.
    - size (str): Width of the colorbar (default "3%").

    Returns:
    - cbar (matplotlib.colorbar.Colorbar): The created colorbar.
    """

    # Extract the scatter plot (mappable object)
    mappable = ax.collections[0]

    # Create a divider for better placement control
    divider = make_axes_locatable(ax)

    # Append colorbar to the right with specified size & padding
    cax = divider.append_axes("right", size=size, pad=pad_size)

    # Add colorbar
    cbar = plt.colorbar(mappable, cax=cax)

    # Set colorbar label and font size
    cbar.set_label(label, fontsize=label_fontsize)

    # Adjust tick font size
    cbar.ax.tick_params(labelsize=tick_fontsize)

    return cbar  # Return colorbar for further modifications if needed


# %%
# figure 1h
augur_3_v_12 = pd.read_csv(os.path.join(DATA_PATH, "augur/augur_cell_type_m3_v_m12_rfc.csv"), index_col=0)
augur_12_v_24 = pd.read_csv(os.path.join(DATA_PATH, "augur/augur_cell_type_m12_v_m24_rfc.csv"), index_col=0)

augur_df = pd.DataFrame({'m3_v_m12': augur_3_v_12.loc['mean_auc'],
                         'm12_v_m24': augur_12_v_24.loc['mean_auc']})
augur_df['delta'] = augur_df['m3_v_m12'] - augur_df['m12_v_m24']
augur_df.index.name = "cell_type"
augur_df.reset_index(inplace=True)

# %%
with plt.rc_context({"figure.figsize": (16, 14), "figure.dpi": 150, "figure.frameon": True}):
    fig = plt.figure()
    gs0 = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[1, 1, 1])
    
    # first row
    gs00 = gs0[0].subgridspec(1, 3)  # 1 row, 3 columns
    
    # ax1: scRNA-seq UMAP
    ax1 = fig.add_subplot(gs00[0, 0])
    sc.pl.umap(adata_rna, color='cell_type_fine', show=False, title="", legend_loc=None, palette=cell_type_fine_colors, ax=ax1)
    repel_umap_labels(
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
    repel_umap_labels(
        adata_atac,
        "cell_type_fine",
        ax=ax3,
        adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
        text_kwargs=dict(fontsize=11, weight='bold')
    )
    
    # second row
    gs01 = gs0[1].subgridspec(1, 4, width_ratios =[1.5, 0.5, 1.5, 0.5])
    
    ax4 = fig.add_subplot(gs01[0, 0])
    cellnum_df_pivot.plot(kind="barh", stacked=True, color=coarse_colors, ax=ax4)
    ax4.tick_params(axis="x", labelsize=10)
    ax4.set_yticklabels(['3 month', '12 month', '24 month'], fontsize=10)
    ax4.set_xlabel("number of cells", fontsize=10)
    ax4.set_ylabel("")
    
    handles, labels = ax4.get_legend_handles_labels()
    ax4.legend(handles[::-1], labels[::-1], title="Cell Type", loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, title_fontsize=8, frameon=False)
    label_bars(ax4, cellnum_df_pivot, celltypes={'Fibroblast', "Endothelial", "Macrophage"})
    
    ax5 = fig.add_subplot(gs01[0, 2])
    cellnum_fine_df_pivot.plot(kind="barh", stacked=True, color=fine_colors, ax=ax5)
    ax5.set_xticklabels([])
    ax5.set_yticklabels(['3 month', '12 month', '24 month'], fontsize=10)
    ax5.set_xlabel("percentage of cells", fontsize=10)
    ax5.set_ylabel("")
    
    handles, labels = ax5.get_legend_handles_labels()
    ax5.legend(handles[::-1], labels[::-1], title="Cell Type", loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, title_fontsize=8, frameon=False)
    
    # third row
    gs02 = gs0[2].subgridspec(1, 4, width_ratios=[1, 1, 0, 1])  # 1 row, 3 columns
    ax6 = fig.add_subplot(gs02[0, 0])
    sc.pl.umap(adata_rna, color='augur_celltype', show=False, title="", legend_loc=None, color_map=augur_cmap, colorbar_loc=None, ax=ax6)
    repel_umap_labels(
        adata_rna,
        "cell_type",
        include=["Macrophage", "Neutrophil", "Fibroblast"],
        ax=ax6,
        adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
        text_kwargs=dict(fontsize=11, weight='bold')
    )
    augur_colorbar(ax6, "AUC", label_fontsize=10, tick_fontsize=8, pad_size=0.1, size="1.5%")
    
    ax7 = fig.add_subplot(gs02[0, 1])
    sc.pl.umap(adata_rna, color='augur_celltype_fine', show=False, title="", legend_loc=None, color_map=augur_cmap, colorbar_loc=None, ax=ax7)
    repel_umap_labels(
        adata_rna,
        "cell_type_fine",
        include=['Fib.1', 'MC.1', "MC.2", "MC/B-Cell", 'Neutrophil'],
        ax=ax7,
        adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
        text_kwargs=dict(fontsize=11, weight='bold')
    )
    augur_colorbar(ax7, "Rank (%)", label_fontsize=10, tick_fontsize=8, pad_size=0.1, size="1.5%")
    
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
    # diagonal line
    x_vals = np.linspace(augur_df["m12_v_m24"].min(), augur_df["m12_v_m24"].max(), 100)  # Get min-max range
    ax8.plot(x_vals, x_vals, linestyle="--", color="grey", linewidth=1)
    
    ax8.tick_params(labelsize=10)
    ax8.set_xlabel("AUC 12 months vs 24 months", fontsize=10)
    ax8.set_ylabel("AUC 3 months vs 12 months", fontsize=10)
    
    
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax8]
    for i, ax in enumerate(axes, 1):
        ax.text(-0.1, 1.1, f"{chr(65+i)}", transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='left')
        ax.grid(False)

    #fig.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1, wspace=0.4, hspace=0.4)
    plt.suptitle("Figure 1. Integrated single-cell transcriptomic and epigenomic analysis of non-cardiac (non-CM) cells during aging")
    fig.subplots_adjust(wspace=0.3)
    fig.tight_layout()
    plt.show()

# %%
fig.savefig(os.path.join(FIGURE_PATH, "figure1.png"), bbox_inches='tight', dpi=300)

# %%
