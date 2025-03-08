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

# %% [markdown]
# ### Figure 4A

# %%
fig4a = plt.imread(f"{FIGURE_PATH}/figure4/figure4a.png")

# %% [markdown]
# ### Figure 4B

# %%
fig4b_1 = plt.imread(f"{FIGURE_PATH}/figure4/figure4b_1.png")
fig4b_2 = plt.imread(f"{FIGURE_PATH}/figure4/figure4b_2.png")

# %% [markdown]
# ### Figure 4C

# %%
fig4c_1 = plt.imread(f"{FIGURE_PATH}/figure4/figure4c_1.png")
fig4c_2 = plt.imread(f"{FIGURE_PATH}/figure4/figure4c_2.png")

# %% [markdown]
# ### Figures 4D-F

# %%
fig4d_1 = plt.imread(f"{FIGURE_PATH}/figure4/figure4d_1.png")
fig4d_2 = plt.imread(f"{FIGURE_PATH}/figure4/figure4d_2.png")

fig4e_1 = plt.imread(f"{FIGURE_PATH}/figure4/figure4e_1.png")
fig4e_2 = plt.imread(f"{FIGURE_PATH}/figure4/figure4e_2.png")

fig4f_1 = plt.imread(f"{FIGURE_PATH}/figure4/figure4f_1.png")
fig4f_2 = plt.imread(f"{FIGURE_PATH}/figure4/figure4f_2.png")

# %% [markdown]
# # Main Figure

# %%
with plt.rc_context({"figure.figsize": (12, 18), "figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10, "xtick.labelsize": 10, 
                     "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 10, "figure.titlesize": 16}):
    
    fig = plt.figure()
    
    # Define the main grid with reduced spacing
    gs0 = gridspec.GridSpec(5, 1, figure=fig, height_ratios=[1.4, 2, 1, 1, 1], wspace=0.05)

    # First row with reduced subplot spacing
    gs00 = gs0[0].subgridspec(1, 3, width_ratios=[2, 1, 1], wspace=0.05)
    ax1_1 = fig.add_subplot(gs00[0, 0])
    ax1_1.imshow(fig4a, aspect='auto')
    ax1_1.axis('off')

    ax1_2 = fig.add_subplot(gs00[0, 1])
    ax1_2.imshow(fig4b_1, aspect='auto')
    ax1_2.axis('off')

    ax1_3 = fig.add_subplot(gs00[0, 2])
    ax1_3.imshow(fig4b_2, aspect='auto')
    ax1_3.axis('off')

    # Second row with reduced spacing
    gs01 = gs0[1].subgridspec(1, 2, width_ratios=[1, 1], wspace=0.05)
    ax2_1 = fig.add_subplot(gs01[0, 0])
    ax2_1.imshow(fig4c_1, aspect='auto')
    ax2_1.axis('off')

    ax2_2 = fig.add_subplot(gs01[0, 1])
    ax2_2.imshow(fig4c_2, aspect='auto')
    ax2_2.axis('off')

    # Third row with reduced spacing
    gs02 = gs0[2].subgridspec(1, 2, width_ratios=[1, 1], wspace=0.05)
    ax3_1 = fig.add_subplot(gs02[0, 0])
    ax3_1.imshow(fig4d_1, aspect='auto')
    ax3_1.axis('off')

    ax3_2 = fig.add_subplot(gs02[0, 1])
    ax3_2.imshow(fig4d_2, aspect='auto')
    ax3_2.axis('off')

    # Fourth row
    gs03 = gs0[3].subgridspec(1, 2, width_ratios=[1, 1], wspace=0.05)
    ax4_1 = fig.add_subplot(gs03[0, 0])
    ax4_1.imshow(fig4e_1, aspect='auto')
    ax4_1.axis('off')
    
    ax4_2 = fig.add_subplot(gs03[0, 1])
    ax4_2.imshow(fig4e_2, aspect='auto')
    ax4_2.axis('off')
    
    # Fifth row
    gs04 = gs0[4].subgridspec(1, 2, width_ratios=[1, 1], wspace=0.05)
    ax5_1 = fig.add_subplot(gs04[0, 0])
    ax5_1.imshow(fig4f_1, aspect='auto')
    ax5_1.axis('off')
    
    ax5_2 = fig.add_subplot(gs04[0, 1])
    ax5_2.imshow(fig4f_2, aspect='auto')
    ax5_2.axis('off')

    # Subfigure labels
    axes_heading = [ax1_1, ax1_2, ax2_1, ax3_1, ax4_1, ax5_1]
    for i, ax in enumerate(axes_heading, 1):
        ax.text(-0.1, 1.1, f"{chr(64+i)}", transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='left')

    ax3_1.text(0, 1.05, "Senescence", transform=ax3_1.transAxes, fontsize=12, va='top', ha='left')
    ax4_1.text(0, 1.05, "Inflammation", transform=ax4_1.transAxes, fontsize=12, va='top', ha='left')
    ax5_1.text(0, 1.05, "Fib.6", transform=ax5_1.transAxes, fontsize=12, va='top', ha='left')
    
    # Title and final layout adjustments
    plt.suptitle("Figure 4 Cell-cell communication during aging")
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0.05, hspace=0.1)  # Adjust global spacing
    plt.show()

    

# %%
fig.savefig(os.path.join(FIGURE_PATH, "figure4.png"), bbox_inches='tight', dpi=300)

# %%
