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
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import json
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
# ### Figure 5a-b

# %%
fig5a = plt.imread(f"{FIGURE_PATH}/figure5/figure5a.png")
fig5b = plt.imread(f"{FIGURE_PATH}/figure5/figure5b.png")

# %% [markdown]
# ### Figure 5c-d

# %%
auc_mtx = pd.read_csv(os.path.join(DATA_PATH, "scenic/tf_auc_scores.csv"), index_col=0)
fib_scores = auc_mtx[auc_mtx['cell_type_fine'].str.startswith("Fib")]
endo_scores = auc_mtx[auc_mtx['cell_type_fine'].str.startswith("Endo")]

fib_scores = fib_scores[['Stat1', "Irf7", "Atf4", "Yy1", "cell_type_fine", "month"]]
fib_scores_long = fib_scores.melt(id_vars=['cell_type_fine', 'month'], var_name='TF', value_name='TF_activity')

endo_scores = endo_scores[["Hlf", "Irf2", "Dbp", "Elk4", "cell_type_fine", "month"]]
endo_scores_long = endo_scores.melt(id_vars=['cell_type_fine', 'month'], var_name='TF', value_name='TF_activity')

# %%
cell_type_colors, cell_type_fine_colors = fig_func.get_celltype_colors()
regulon_dict = json.load(open(os.path.join(DATA_PATH, "scenic/regulons.json"), "r"))
highlight_df = pd.read_csv(os.path.join(DATA_PATH, "scenic/go_analysis_highlight.csv"))

# %%
fig5c_l = plt.imread(os.path.join(FIGURE_PATH, "figure5/figure5c_l.png"))
fig5d_l = plt.imread(os.path.join(FIGURE_PATH, "figure5/figure5d_l.png"))

# %% [markdown]
# # Main Figure

# %%
with plt.rc_context({"figure.figsize": (12, 16), "figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10, "xtick.labelsize": 10, 
                     "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 10, "figure.titlesize": 16}):
    
    fig = plt.figure()
    
    # Define the main grid with reduced spacing
    gs0 = gridspec.GridSpec(5, 1, figure=fig, height_ratios=[1.85, 0.6, 1, 0.6, 1])

    # First row
    gs00 = gs0[0].subgridspec(1, 2, width_ratios=[1, 1])
    ax1_1 = fig.add_subplot(gs00[0])
    ax1_1.imshow(fig5a, aspect="auto")
    ax1_1.axis('off')
    
    ax1_2 = fig.add_subplot(gs00[1])
    ax1_2.imshow(fig5b, aspect="auto")
    ax1_2.axis('off')

    # Second row 
    gs01 = gs0[1].subgridspec(1, 4, width_ratios=[1, 1, 1, 1])
    ax2_1 = fig.add_subplot(gs01[0])
    ax2_1.set_box_aspect(0.8) 
    fig_func.tf_activity_lineplot(fib_scores_long, regulon_dict, 'Stat1', significance=[('m3', 'm12', 0.001), ('m3', 'm24', 0.001)],
                     legend=False, hue='cell_type_fine', palette=cell_type_fine_colors, ax=ax2_1)
    
    ax2_2 = fig.add_subplot(gs01[1])
    ax2_2.set_box_aspect(0.8) 
    fig_func.tf_activity_lineplot(fib_scores_long, regulon_dict, 'Irf7', significance=[('m3', 'm12', 0.001), ('m3', 'm24', 0.001)],
                     legend=False, hue='cell_type_fine', palette=cell_type_fine_colors, ax=ax2_2)
    
    ax2_3 = fig.add_subplot(gs01[2])
    ax2_3.set_box_aspect(0.8) 
    fig_func.tf_activity_lineplot(fib_scores_long, regulon_dict, 'Atf4', significance=[('m3', 'm12', 0.001), ('m3', 'm24', 0.001)],
                     legend=False, hue='cell_type_fine', palette=cell_type_fine_colors, ax=ax2_3)
    
    ax2_4 = fig.add_subplot(gs01[3])
    ax2_4.set_box_aspect(0.8) 
    fig_func.tf_activity_lineplot(fib_scores_long, regulon_dict, 'Yy1', significance=[('m3', 'm12', 0.001), ('m3', 'm24', 0.001)],
                     legend=False, hue='cell_type_fine', palette=cell_type_fine_colors, ax=ax2_4)
   
    # Third row
    gs02 = gs0[2].subgridspec(1, 1, width_ratios=[1])
    ax3 = fig.add_subplot(gs02[0])
    ax3.imshow(fig5c_l, aspect="auto")
    ax3.axis('off')

    # Fourth row
    gs03 = gs0[3].subgridspec(1, 4, width_ratios=[1, 1, 1, 1])
    ax4_1 = fig.add_subplot(gs03[0])
    ax4_1.set_box_aspect(0.8) 
    fig_func.tf_activity_lineplot(endo_scores_long, regulon_dict, 'Hlf', significance=[('m3', 'm12', 0.1), ('m3', 'm24', 0.001)],
                     legend=False, hue='cell_type_fine', palette=cell_type_fine_colors, ax=ax4_1)
    
    ax4_2 = fig.add_subplot(gs03[1])
    ax4_2.set_box_aspect(0.8) 
    fig_func.tf_activity_lineplot(endo_scores_long, regulon_dict, 'Irf2', significance=[('m3', 'm12', 0.1), ('m3', 'm24', 0.1)],
                     legend=False, hue='cell_type_fine', palette=cell_type_fine_colors, ax=ax4_2)
    
    ax4_3 = fig.add_subplot(gs03[2])
    ax4_3.set_box_aspect(0.8) 
    fig_func.tf_activity_lineplot(endo_scores_long, regulon_dict, 'Dbp', significance=[('m3', 'm12', 0.001), ('m3', 'm24', 0.001)],
                     legend=False, hue='cell_type_fine', palette=cell_type_fine_colors, ax=ax4_3)
    
    ax4_4 = fig.add_subplot(gs03[3])
    ax4_4.set_box_aspect(0.8) 
    fig_func.tf_activity_lineplot(endo_scores_long, regulon_dict, 'Elk4', significance=[('m3', 'm12', 0.01), ('m3', 'm24', 0.1)],
                     legend=False, hue='cell_type_fine', palette=cell_type_fine_colors, ax=ax4_4)
    
    # Fifth row
    gs04 = gs0[4].subgridspec(1, 1, width_ratios=[1])
    ax4 = fig.add_subplot(gs04[0])
    ax4.imshow(fig5d_l, aspect="auto")
    ax4.axis('off')

    # Subfigure labels
    ax1_1.text(-0.1, 1.1, "A", transform=ax1_1.transAxes, fontsize=20, fontweight='bold', va='top', ha='left')
    ax1_2.text(-0.1, 1.1, "B", transform=ax1_2.transAxes, fontsize=20, fontweight='bold', va='top', ha='left')
    ax2_1.text(-0.1, 1.3, "C", transform=ax2_1.transAxes, fontsize=20, fontweight='bold', va='top', ha='left')
    ax4.text(-0.1, 1.3, "D", transform=ax4_1.transAxes, fontsize=20, fontweight='bold', va='top', ha='left')

    # Title and final layout adjustments
    plt.suptitle("Figure 5. Gene regulatory networks during aging by single-cell \n regulatory network inference and clustering (SCENIC)")
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.25, wspace=0.35, top=0.92, bottom=0.01)
    #fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0.25, hspace=0.1)  # Adjust global spacing
    plt.show()

    

# %%
fig.savefig(os.path.join(FIGURE_PATH, "figure5.png"), bbox_inches='tight', dpi=300)

# %%
