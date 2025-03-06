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

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# ## Main Figure

# %%
with plt.rc_context({"figure.figsize": (12, 18), "figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10,"xtick.labelsize": 10, 
                    "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 10, "figure.titlesize": 16}):
    fig = plt.figure()
    gs0 = gridspec.GridSpec(5, 1, figure=fig, height_ratios=[2, 2, 1, 1, 1])
    
    # first row
    gs00 = gs0[0].subgridspec(1, 4, width_ratios=[1, 1, 1, 1])
    ax1_1 = fig.add_subplot(gs00[0, 0])
    ax1_2 = fig.add_subplot(gs00[0, 1])
    ax1_3 = fig.add_subplot(gs00[0, 2])
    ax1_4 = fig.add_subplot(gs00[0, 3])
    
    # second row
    gs01 = gs0[1].subgridspec(1, 2, width_ratios=[1, 1.2])
    
    ax2_1 = fig.add_subplot(gs01[0, 0])
    ax2_2 = fig.add_subplot(gs01[0, 1])

    # third row
    gs02 = gs0[2].subgridspec(1, 2, width_ratios=[1, 1]) 
    ax3_1 = fig.add_subplot(gs02[0, 0])
    ax3_2 = fig.add_subplot(gs02[0, 1])
    
    # fourth row
    gs03 = gs0[3].subgridspec(1, 2, width_ratios=[1, 1])
    ax4_1 = fig.add_subplot(gs03[0, 0])
    ax4_2 = fig.add_subplot(gs03[0, 1])
    
    # fifth row
    gs04 = gs0[4].subgridspec(1, 2, width_ratios=[1, 1])
    ax5_1 = fig.add_subplot(gs04[0, 0])
    ax5_2 = fig.add_subplot(gs04[0, 1])
    
    axes_heading = [ax1_1, ax1_3, ax2_1, ax3_1, ax4_1, ax5_1]
    for i, ax in enumerate(axes_heading, 1):
        ax.text(-0.1, 1.1, f"{chr(64+i)}", transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='left') # subfigure labels
    

    plt.suptitle("Figure 4 Cell-cell communication during aging")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.tight_layout()
    plt.show()
    

# %%

# %%
