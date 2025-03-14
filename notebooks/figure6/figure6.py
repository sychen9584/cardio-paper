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
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
sys.path.append('../../scripts')
import figure_functions as fig_func

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
# ### Figure 6A

# %%
dar_df = pd.read_csv(os.path.join(DATA_PATH, "deg/atac_dar_annotated.csv"))

fibroblast = dar_df.query("group == 'Fibroblast'").copy()
endothelial = dar_df.query("group == 'Endothelial'").copy()
macrophage = dar_df.query("group == 'Macrophage'").copy()

fibroblast = pd.DataFrame(100*fibroblast.value_counts("peak_type") / fibroblast.shape[0]).reset_index().sort_values("peak_type", ascending=False)
endothelial = pd.DataFrame(100*endothelial.value_counts("peak_type") / endothelial.shape[0]).reset_index().sort_values("peak_type", ascending=False)
macrophage = pd.DataFrame(100*macrophage.value_counts("peak_type") / macrophage.shape[0]).reset_index().sort_values("peak_type", ascending=False)

celltypes = ['Fibroblast', 'Endothelial', 'Macrophage']
palette = sns.color_palette('Set2')

# %% [markdown]
# ### Figure 6B

# %%
macrophage_dar = set(dar_df[dar_df['group'] == 'Macrophage']['names'].values)
fibroblast_dar = set(dar_df[dar_df['group'] == 'Fibroblast']['names'].values)
endothelial_dar = set(dar_df[dar_df['group'] == 'Endothelial']['names'].values)

# %% [markdown]
# ### Figure 6C

# %%
rfx6 = plt.imread(f"{DATA_PATH}/homer_motifs/rfx6.png")
irf8 = plt.imread(f"{DATA_PATH}/homer_motifs/irf8.png")
stat4 = plt.imread(f"{DATA_PATH}/homer_motifs/stat4.png")
atf3 = plt.imread(f"{DATA_PATH}/homer_motifs/atf3.png")

# %% [markdown]
# ### Figure 6D   

# %%
df_dev_long = pd.read_csv(os.path.join(DATA_PATH, "processed/motif_activity.csv"))

# %%
month_colors = {"m3": "#98df8a", "m12":"#FFED6F", "m24":"#ff9896"}
df_dev_long['cell_type'] = pd.Categorical(df_dev_long['cell_type'], categories=['Endothelial', "Fibroblast", 'Macrophage'], ordered=True)
g = sns.catplot(
    data=df_dev_long,
    x="motif", 
    y="deviation", 
    hue="month",
    palette=month_colors,
    row="cell_type",
    sharey=False,
    kind="bar",
    height=4, 
    aspect=5.5, 
    dodge=True  # Enables dodging for hue
)

g.set_titles("{row_name}")
g.set_axis_labels("Motif", "")
g.fig.text(0.01, 0.5, "Average Motif Activity", va='center', rotation=90, fontsize=12)

for ax in g.axes.flat:
    ax.set_title(ax.get_title(), fontsize=14, fontweight="bold", bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"))
    ax.grid(False)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

# %%
plt.show()

# %%
with plt.rc_context({"figure.figsize": (12, 14), "figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10,"xtick.labelsize": 10, 
                    "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 10, "figure.titlesize": 16}):
    
    fig = plt.figure(layout="constrained")
    ax_dict = fig.subplot_mosaic(
        """
        ABB
        ABB
        ACC
        DDD
        DDD
        DDD
        """
    )
    
    # A: pie charts
    for i, (data, celltype) in enumerate(zip([fibroblast, endothelial, macrophage], celltypes)):
        inset_ax = ax_dict['A'].inset_axes([-0.1, 0.7 - i * 0.35, 0.7, 0.3])  # Adjust y position for stacking
        wedges, _ = inset_ax.pie(data['count'], startangle=140, colors=palette)

        # Create formatted legend labels
        legend_labels = data['peak_type'] + " (" + data['count'].round(2).astype(str) + "%)"

        # Add legend inside A (outside each pie chart)
        inset_ax.legend(wedges, legend_labels, title=celltype, loc="upper right", 
                        bbox_to_anchor=(1.7, 0.95), frameon=False, fontsize=8)
        
    # Strip all axis labels, ticks, and spines from subplot A
    ax_dict['A'].set_xticks([])
    ax_dict['A'].set_yticks([])
    ax_dict['A'].grid(False)
    ax_dict['A'].set_frame_on(False)
    
    # B: Venn diagram
    fig_func.venn3_custom(fibroblast_dar, endothelial_dar, macrophage_dar, labels=('Fibroblast', 'Endothelial', 'Macrophage'), 
                              normalize_range=(0, 30000), title='', ax=ax_dict['B'])
    
    # C: Motif logos
    motifs = {"Rfx6": rfx6, "Irf8": irf8, "Stat4": stat4, "Atf3": atf3}
    grid_positions = [(0, 0), (0, 1), (1, 0), (1, 1)]  # (row, col) format
    for i, (motif_name, motif_image) in enumerate(motifs.items()):
        row, col = grid_positions[i]  # Get grid position
        x_offset = 0.05 + col * 0.45   # Adjust x based on column index
        y_offset = 0.55 - row * 0.45   # Adjust y based on row index

        # Create inset axis
        inset_ax = ax_dict['C'].inset_axes([x_offset, y_offset, 0.4, 0.4])  # [x, y, width, height]
        inset_ax.imshow(motif_image)
        inset_ax.set_title(motif_name, fontsize=14, fontweight='bold')
        inset_ax.axis("off")
        
    ax_dict['C'].set_xticks([])
    ax_dict['C'].set_yticks([])
    ax_dict['C'].grid(False)
    ax_dict['C'].set_frame_on(False)
    
    ax_dict["D"].imshow(g.fig.canvas.buffer_rgba(), aspect='auto')  # Embed catplot into the subplot
    # Hide the default axes of ax_dict['D']
    ax_dict["D"].axis("off")
    
    label_positions = {
        "A": (-0.1, 1.05),
        "B": (-0.05, 1.1),
        "C": (-0.05, 1.1),
        "D": (-0.05, 1.05)  # Slightly lower for the large subplot
    }
    
    for k, ax in ax_dict.items():
        ax.annotate(k, xy=label_positions[k], xycoords='axes fraction', 
                    fontsize=20, fontweight='bold', va='top', ha='left')

    

    plt.suptitle("Figure 6. Integrated single-cellATAC seq epigenomic analysis of cardiac of non-CM cells during aging")
    plt.show()

# %%
fig.savefig(os.path.join(FIGURE_PATH, "figure6.png"), bbox_inches='tight', dpi=300)

# %%
