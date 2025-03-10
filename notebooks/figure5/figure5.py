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
import matplotlib.image as mpimg  # Load images
import matplotlib.colors as mcolors
import matplotlib.cm as cm
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

# %%
fig5a = plt.imread(f"{FIGURE_PATH}/figure5/figure5a.png")
fig5b = plt.imread(f"{FIGURE_PATH}/figure5/figure5b.png")

# %% [markdown]
# # Main Figure

# %%
with plt.rc_context({"figure.figsize": (12, 18), "figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10, "xtick.labelsize": 10, 
                     "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 10, "figure.titlesize": 16}):
    
    fig = plt.figure()
    
    # Define the main grid with reduced spacing
    gs0 = gridspec.GridSpec(5, 1, figure=fig, height_ratios=[2, 1, 1, 1, 1], wspace=0.05)

    # First row
    gs00 = gs0[0].subgridspec(1, 2, width_ratios=[1, 1])
    ax1_1 = fig.add_subplot(gs00[0])
    ax1_1.imshow(fig5a)
    ax1_1.axis('off')
    
    ax1_2 = fig.add_subplot(gs00[1])
    ax1_2.imshow(fig5b)
    ax1_2.axis('off')

    # Second row 
    gs01 = gs0[1].subgridspec(1, 4, width_ratios=[1, 1, 1, 1])
    ax2_1 = fig.add_subplot(gs01[0])
    ax2_2 = fig.add_subplot(gs01[1])
    ax2_3 = fig.add_subplot(gs01[2])
    ax2_4 = fig.add_subplot(gs01[3])
   
    # Third row
    gs02 = gs0[2].subgridspec(1, 4, width_ratios=[1, 1, 1, 1])
    ax3_1 = fig.add_subplot(gs02[0])
    ax3_2 = fig.add_subplot(gs02[1])
    ax3_3 = fig.add_subplot(gs02[2])
    ax3_4 = fig.add_subplot(gs02[3])

    # Fourth row
    gs03 = gs0[3].subgridspec(1, 4, width_ratios=[1, 1, 1, 1])
    ax4_1 = fig.add_subplot(gs03[0])
    ax4_2 = fig.add_subplot(gs03[1])
    ax4_3 = fig.add_subplot(gs03[2])
    ax4_4 = fig.add_subplot(gs03[3])
    
    # Fifth row
    gs04 = gs0[4].subgridspec(1, 4, width_ratios=[1, 1, 1, 1])
    ax5_1 = fig.add_subplot(gs04[0])
    ax5_2 = fig.add_subplot(gs04[1])
    ax5_3 = fig.add_subplot(gs04[2])
    ax5_4 = fig.add_subplot(gs04[3])
    

    # Subfigure labels
    axes_heading = [ax1_1, ax1_2, ax2_1, ax3_1, ax4_1, ax5_1]
    for i, ax in enumerate(axes_heading, 1):
        ax.text(-0.1, 1.1, f"{chr(64+i)}", transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='left')
    
    # Title and final layout adjustments
    plt.suptitle("Figure 5. Gene regulatory networks during aging by single-cell \n regulatory network inference and clustering (SCENIC)")
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.25)
    plt.show()

    

# %%
