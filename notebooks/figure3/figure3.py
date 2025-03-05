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
# # Figure 3

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import os
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.image as mpimg  # Load images
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

plt.rcParams["axes.grid"] = False  # Disable grids for all plots
plt.rcParams["grid.color"] = "white"  # Ensure grids are fully invisible
plt.rcParams["grid.alpha"] = 0  # Remove any transparency effects
plt.grid(False)  # Turn off grids explicitly

# %%
DATA_PATH = "/home/sychen9584/projects/cardio_paper/data"
FIGURE_PATH = "/home/sychen9584/projects/cardio_paper/figures"

# %% [markdown]
# ### Figure 3A

# %%
# Load image
image_path = os.path.join(FIGURE_PATH, "figure3a.png")  # Adjust path as needed
fig3a = mpimg.imread(image_path)

# %% [markdown]
# ### Figure 3B

# %%
go_df = pd.read_csv(os.path.join(DATA_PATH, "deg/scRNA_cell_type_fine_GO_highlight.csv"))

# %%
go_df['query'].replace({
    'Granulocyte/Neutrophil': "Neutrophil",
    'Peri/Smooth Muscle': 'Smooth Muscle'
}, inplace=True)

# rename some terms so they fit in plot
go_df.at[14, 'name'] = 'MHC class II antigen presentation'
go_df.at[13, 'name'] = 'myeloid cell activation'

# %%
cell_type_colors, cell_type_fine_colors = fig_func.get_celltype_colors()

# %%
ax = sns.barplot(data=go_df, x="neg_log_p_value", y="query", palette=cell_type_fine_colors)
ax.set_xlim(0, 20)

# 🔹 Annotate each bar correctly
for index, p in enumerate(ax.patches):
    go_term = go_df.iloc[index]['name']
    ax.annotate(go_term,
                (0.5, p.get_y() + p.get_height() / 2),  # Place at end of bar
                ha='left', va='center', fontsize=8, color="black", xytext=(5, 0), textcoords="offset points")

plt.show()


# %% [markdown]
# ### Figure 3C

# %%
#### 3 months vs 12 months
deg_results_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_3m_vs_12m_DEG.csv'))
macrophage3v12 = set(deg_results_df[deg_results_df['cell_type'] == 'Macrophage']['names'].values)
fibroblast3v12 = set(deg_results_df[deg_results_df['cell_type'] == 'Fibroblast']['names'].values)
endothelial3v12 = set(deg_results_df[deg_results_df['cell_type'] == 'Endothelial']['names'].values)

# %%
#### 12 months vs 24 months
deg_results_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_12m_vs_24m_DEG.csv'))
macrophage12v24 = set(deg_results_df[deg_results_df['cell_type'] == 'Macrophage']['names'].values)
fibroblast12v24 = set(deg_results_df[deg_results_df['cell_type'] == 'Fibroblast']['names'].values)
endothelial12v24 = set(deg_results_df[deg_results_df['cell_type'] == 'Endothelial']['names'].values)

# %% [markdown]
# ### Figure 3D

# %%
macrophage_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_macrophage_logfc.csv'), index_col=0)
fibroblast_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_fibroblast_logfc.csv'), index_col=0)
endothelial_deg_df = pd.read_csv(os.path.join(DATA_PATH, 'deg/scRNA_endothelial_logfc.csv'), index_col=0)

# %%
endo_hm = fig_func.plot_log2fc_heatmap(endothelial_deg_df, "Endothelial", ['3 months', '12 months', '24 months'], 'cluster', [3, 2, 0 ,1], title_ypad=1.1, plot_legend=False)
fib_hm = fig_func.plot_log2fc_heatmap(fibroblast_deg_df, "Fibroblast", ['3 months', '12 months', '24 months'], 'cluster', [3, 4, 0, 5, 1, 2], title_ypad=1.1, plot_legend=False)
mac_hm = fig_func.plot_log2fc_heatmap(macrophage_deg_df, "Macrophage", ['3 months', '12 months', '24 months'], 'cluster', [4, 6, 0, 3, 5, 1, 2], title_ypad=1.1, legend_vpad=260)
plt.show()

# %%
endo_img = fig_func.ax_to_image(endo_hm)
fib_img = fig_func.ax_to_image(fib_hm)
mac_img = fig_func.ax_to_image(mac_hm)

# %% [markdown]
# ## Main Figure

# %%
with plt.rc_context({"figure.figsize": (12, 18), "figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10,"xtick.labelsize": 10, 
                    "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 10, "figure.titlesize": 16}):
    fig = plt.figure()
    gs0 = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[2, 1, 2])
    
    # first row
    gs00 = gs0[0].subgridspec(1, 2, width_ratios=[1, 1])
    ax1 = fig.add_subplot(gs00[0, 0])
    ax1.imshow(fig3a, aspect="auto")
    ax1.axis('off')
    
    ax2 = fig.add_subplot(gs00[0, 1])
    ax2 = sns.barplot(data=go_df, x="neg_log_p_value", y="query", palette=cell_type_fine_colors)
    ax2.set_xlim(0, 20)
    ax2.set_xlabel("-Log10 (P-value Adjust)")
    ax2.set_ylabel("")

    # 🔹 Annotate each bar correctly
    for index, p in enumerate(ax2.patches):
        go_term = go_df.iloc[index]['name']
        ax2.annotate(go_term,
                    (0.5, p.get_y() + p.get_height() / 2),  # Place at end of bar
                    ha='left', va='center', fontsize=8, weight="bold", color="black", xytext=(5, 0), textcoords="offset points")
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    
    # second row
    gs01 = gs0[1].subgridspec(1, 2, width_ratios=[1, 1])
    
    ax3 = fig.add_subplot(gs01[0, 0])
    fig_func.venn3_custom(macrophage3v12, fibroblast3v12, endothelial3v12,
                          labels=('Macrophage', 'Fibroblast', 'Endothelial'), 
                          normalize_range=(0, 2500), title='3 months vs 12 months DEGs', ax=ax3)
    
    ax4 = fig.add_subplot(gs01[0, 1])
    fig_func.venn3_custom(macrophage12v24, fibroblast12v24, endothelial12v24,
                          labels=('Macrophage', 'Fibroblast', 'Endothelial'), 
                          normalize_range=(0, 2500), title='12 months vs 24 months DEGs', ax=ax4)
    
    # third row
    gs02 = gs0[2].subgridspec(1, 4, width_ratios=[1, 1, 1.1, 3]) 
    
    ax5 = fig.add_subplot(gs02[0, 0])
    ax5.imshow(endo_img, aspect="auto")
    ax5.axis('off')
    
    ax6 = fig.add_subplot(gs02[0, 1])
    ax6.imshow(fib_img, aspect="auto")
    ax6.axis('off')
    
    ax7 = fig.add_subplot(gs02[0, 2])
    ax7.imshow(mac_img, aspect="auto")
    ax7.axis('off')
    
    ax8 = fig.add_subplot(gs02[0, 3])
    
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    for i, ax in enumerate(axes, 1):
        ax.grid(False)
    
    axes_heading = [ax for ax in axes if ax not in [ax4, ax6, ax7]]
    for i, ax in enumerate(axes_heading, 1):
        ax.text(-0.1, 1.1, f"{chr(64+i)}", transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='left') # subfigure labels
    

    plt.suptitle("Figure 3. Differentially expressed genes and cluster analysis during aging")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.tight_layout()
    plt.show()

# %%

# %%

# %%
