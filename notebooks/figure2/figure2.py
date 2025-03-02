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
# # Figure 2

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
import matplotlib.patheffects as path_effects
import starbars

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

# %%
adata = sc.read_h5ad(os.path.join(DATA_PATH, 'processed/scRNA_all.h5ad'))
geneset_scores = pd.read_csv(os.path.join(DATA_PATH, 'geneset_scores.csv'), index_col=0)

# %%
month_colors = {
    "m3": "#98df8a",
    'm12': "#fdd0a2",
    'm24': "#ff9896"
}

# %%
sasp_fib1 = geneset_scores[geneset_scores['cell_type_fine'] == 'Fib.1']
sasp_fib6 = geneset_scores[geneset_scores['cell_type_fine'] == 'Fib.6']
sasp_mc = geneset_scores[geneset_scores['cell_type_fine'] == 'MC.2']
sasp_neutrophil = geneset_scores[geneset_scores['cell_type_fine'] == 'Granulocyte/Neutrophil']


# %%
def add_median_labels(ax: plt.Axes, fmt: str = ".4f") -> None:
    """Add text labels to the median lines of a seaborn boxplot.

    Args:
        ax: plt.Axes, e.g., the return value of sns.boxplot()
        fmt: format string for the median value
    """
    lines = ax.get_lines()
    
    # Automatically find median lines (they are the only lines with a single X or Y value)
    median_lines = [line for line in lines if line.get_linestyle() == '-' and len(set(line.get_ydata())) == 1]

    for median in median_lines:
        x, y = (data.mean() for data in median.get_data())

        # Annotate median
        text = ax.text(x, y, f'{y:{fmt}}', ha='center', va='center',
                       fontweight='bold', fontsize=8, color='white')

        # Add contrast stroke
        text.set_path_effects([
            path_effects.Stroke(linewidth=2, foreground=median.get_color()),
            path_effects.Normal(),
        ])


# %% [markdown]
# ## Main Figure

# %%
with plt.rc_context({"figure.figsize": (13, 16), "figure.dpi": 150, "figure.frameon": True, 
                     "axes.labelsize": 10, "axes.titlesize": 10,"xtick.labelsize": 10, 
                    "ytick.labelsize": 10, "legend.fontsize": 8, "legend.title_fontsize": 10, "figure.titlesize": 16}):
    fig = plt.figure()
    gs0 = gridspec.GridSpec(4, 1, figure=fig, height_ratios=[1.4, 1, 1, 1])
    
    # first row
    gs00 = gs0[0].subgridspec(1, 2)
    
    ax1 = fig.add_subplot(gs00[0, 0])
    sns.histplot(geneset_scores['SASP'], bins=120, ax=ax1)
    for bar in ax1.patches:
        if bar.get_x() + bar.get_width() / 2 > 0.09:  # Bar center is above threshold
            bar.set_facecolor("#10518B")
        else:
            bar.set_facecolor("#BDD9F1")
    ax1.set_title("SASP Gene Set")
    ax1.set_xlabel("AUC histogram")
    ax1.set_ylabel("Frequency")
    
    ax2 = fig.add_subplot(gs00[0, 1])
    sc.pl.umap(adata, color='SASP_score', ax=ax2, show=False, title="", cmap="plasma", size=12, vmin=0.03, colorbar_loc=None)
    fig1.repel_umap_labels(
        adata,
        "cell_type",
        include=['Fibroblast', "Macrophage", "Endothelial"],
        ax=ax2,
        adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
        text_kwargs=dict(fontsize=12, weight='bold')
    )
    fig1.augur_colorbar(ax2, "AUC Score", label_fontsize=10, tick_fontsize=8, pad_size=0.05, size="1.5%")
    
    # second row
    gs01 = gs0[1].subgridspec(1, 3)
    
    ax3 = fig.add_subplot(gs01[0, 0])
    sns.kdeplot(geneset_scores, x="SASP", hue="month", fill=True, palette=month_colors, hue_order=["m24", "m12", "m3"], alpha=0.4, ax=ax3)
    sns.kdeplot(geneset_scores, x="SASP", hue='month', linewidth=1, fill=False, palette={'m3': 'black', 'm12': 'black', 'm24': 'black'}, ax=ax3)
    handles = [mpatches.Patch(color=month_colors[m], label=m) for m in ["m3", "m12", "m24"]]
    labels = ["3 months", "12 months", "24 months"]
    ax3.legend(handles, labels, title="Time point")
    ax3.set_xlabel("AUC Score")
    
    ax4 = fig.add_subplot(gs01[0, 1])
    sns.boxplot(sasp_fib1, x="month", y="SASP", hue='month', palette=month_colors, linecolor="black", linewidth=1.2, showfliers=False, showcaps=False, ax=ax4)
    sig = [('m3', 'm12', 0.001), ('m3', 'm24', 0.001), ('m12', 'm24', 0.01)]
    starbars.draw_annotation(sig, ax=ax4)
    ax4.set_title("SASP for Fib1.Cxcl1")
    add_median_labels(ax4)
    
    ax5 = fig.add_subplot(gs01[0, 2])
    sns.boxplot(sasp_fib6, x="month", y="SASP", hue='month', palette=month_colors, linecolor="black", linewidth=1.2, showfliers=False, showcaps=False, ax=ax5)
    sig = [('m3', 'm12', 1), ('m3', 'm24', 0.001), ('m12', 'm24', 0.001)]
    starbars.draw_annotation(sig, ax=ax5)
    ax5.set_title("SASP for Fib6.Erbb4")
    add_median_labels(ax5)
    
    # third row
    gs02 = gs0[2].subgridspec(1, 3) 
    
    ax6 = fig.add_subplot(gs02[0, 0])
    sns.boxplot(geneset_scores, x="month", y="SASP", hue='month', palette=month_colors, linecolor="black", linewidth=1.2, showfliers=False, showcaps=False, ax=ax6)
    sig = [('m3', 'm12', 0.001), ('m3', 'm24', 0.001), ('m12', 'm24', 1)]
    starbars.draw_annotation(sig, ax=ax6)
    add_median_labels(ax6)
   
    ax7 = fig.add_subplot(gs02[0, 1])
    sns.boxplot(sasp_mc, x="month", y="SASP", hue='month', palette=month_colors, linecolor="black", linewidth=1.2, showfliers=False, showcaps=False, ax=ax7)
    sig = [('m3', 'm12', 1), ('m3', 'm24', 1), ('m12', 'm24', 0.01)]
    starbars.draw_annotation(sig, ax=ax7)
    ax7.set_title("SASP for MC.Cd209a")
    add_median_labels(ax7)
    
    ax8 = fig.add_subplot(gs02[0, 2])
    sns.boxplot(sasp_neutrophil, x="month", y="SASP", hue='month', palette=month_colors, linecolor="black", linewidth=1.2, showfliers=False, showcaps=False, ax=ax8)
    sig = [('m3', 'm12', 0.001), ('m3', 'm24', 0.001), ('m12', 'm24', 1)]
    starbars.draw_annotation(sig, ax=ax8)
    ax8.set_title("SASP for Neutrophil")
    add_median_labels(ax8)
    
    # fourth row
    gs03 = gs0[3].subgridspec(1, 3) 
    
    ax9 = fig.add_subplot(gs03[0, 0])
    sns.boxplot(geneset_scores, x="month", y="Cardiac Fibrosis", hue='month', palette=month_colors, linecolor="black", linewidth=1.2, showfliers=False, showcaps=False, ax=ax9)
    sig = [('m3', 'm12', 0.001), ('m3', 'm24', 0.001), ('m12', 'm24', 0.001)]
    starbars.draw_annotation(sig, ax=ax9)
    ax9.set_title("Cardiac Fibrosis for All Clusters")
    add_median_labels(ax9)
    
    ax10 = fig.add_subplot(gs03[0, 1])
    sns.boxplot(geneset_scores, x="month", y="Heart Failure", hue='month', palette=month_colors, linecolor="black", linewidth=1.2, showfliers=False, showcaps=False, ax=ax10)
    sig = [('m3', 'm12', 0.001), ('m3', 'm24', 0.001), ('m12', 'm24', 0.01)]
    starbars.draw_annotation(sig, ax=ax10)
    ax10.set_title("Heart Failure for All Clusters")
    add_median_labels(ax10)
    
    ax11 = fig.add_subplot(gs03[0, 2])
    sns.boxplot(geneset_scores, x="month", y="Inflammation", hue='month', palette=month_colors, linecolor="black", linewidth=1.2, showfliers=False, showcaps=False, ax=ax11)
    sig = [('m3', 'm12', 0.001), ('m3', 'm24', 0.001), ('m12', 'm24', 0.001)]
    starbars.draw_annotation(sig, ax=ax11)
    ax11.set_title("Inflammation for All Clusters")
    add_median_labels(ax11)
    
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11]
    for i, ax in enumerate(axes, 1):
        ax.grid(False)
    
    axes_heading = [ax for ax in axes if ax not in [ax5, ax6, ax7, ax8]]
    for i, ax in enumerate(axes_heading, 1):
        ax.text(-0.1, 1.1, f"{chr(64+i)}", transform=ax.transAxes, fontsize=20, fontweight='bold', va='top', ha='left') # subfigure labels
    
    axes_boxplot = [ax for ax in axes if ax not in [ax1, ax2, ax3]]
    for ax in axes_boxplot:
        ax.set_xticklabels(["3 month", "12 month", "24 month"])
        ax.set_xlabel("")
        ax.set_ylabel("AUC Score")
        for spine in ax.spines.values():
            spine.set_visible(False)
    

    plt.suptitle("Figure 2. Gene signatures profiling of aging-susceptible of cardiac non-myocyte changes")
    fig.subplots_adjust(wspace=0.3)
    fig.tight_layout()
    plt.show()

# %%
fig.savefig(os.path.join(FIGURE_PATH, "figure2.png"), bbox_inches='tight', dpi=300)

# %%
