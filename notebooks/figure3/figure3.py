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
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.image as mpimg  # Load images
import matplotlib.gridspec as gridspec

import starbars

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
# ## Figure 3A

# %%
# Load image
image_path = os.path.join(FIGURE_PATH, "figure3a.png")  # Adjust path as needed
fig3a = mpimg.imread(image_path)

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
    ax1.imshow(fig3a)
    ax1.axis('off')
    
    ax2 = fig.add_subplot(gs00[0, 1])
    
    # second row
    gs01 = gs0[1].subgridspec(1, 5, width_ratios=[0, 1, 0, 1, 0])
    
    ax3 = fig.add_subplot(gs01[0, 1])
    
    ax4 = fig.add_subplot(gs01[0, 3])
    
    # third row
    gs02 = gs0[2].subgridspec(1, 4, width_ratios=[1, 1, 1, 3]) 
    
    ax5 = fig.add_subplot(gs02[0, 0])
    ax6 = fig.add_subplot(gs02[0, 1])
    ax7 = fig.add_subplot(gs02[0, 2])
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
