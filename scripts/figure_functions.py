import os
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as path_effects
import matplotlib.cm as cm
import plotly.graph_objects as go
import plotly.io as pio

from PIL import Image
from io import BytesIO
from adjustText import adjust_text
from typing import Dict, List, Optional, Set, Tuple
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib_venn import venn3_unweighted

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


def get_augur_colors():
    augur_colors = ["#27408B",  # royalblue4
                    "#FFF8DC",  # cornsilk1
                    "#FF8C00",  # DarkOrange
                    "#CD4F39"]  # tomato3

    # Create a custom colormap
    return mcolors.LinearSegmentedColormap.from_list("custom_cmap", augur_colors)


def get_celltype_colors():
    """
    Returns a dictionary of colors for each cell type in the dataset.
    """

    # colors of cell type on UMAP
    cell_type_fine_colors = {
        
        # B-Cell
        "B-Cell": "#1f77b4",

        # Endo types in distinct blue shades
        "Endo.1": "#2171b5", 
        "Endo.2": "#4292c6", 
        "Endo.3": "#6baed6", 
        "Endo.4": "#9ecae1", 
        "Endo.5": "#c6dbef", 
        "Endo.6": "#deebf7", 
        "Endo.7": "#08306b", 

        # Fib types in distinct orange shades
        "Fib.1": "#e6550d", 
        "Fib.2": "#fd8d3c", 
        "Fib.3": "#fdae6b", 
        "Fib.4": "#fdd0a2", 
        "Fib.5": "#feedde", 
        "Fib.6": "#a63603", 

        # MC types in similar colors
        "MC.1": "#31a354", 
        "MC.2": "#74c476", 
        "MC.3": "#a1d99b", 
        "MC/B-Cell": "#c7e9c0", 

        # Others
        "Neutrophil": "#9467bd",
        "Smooth Muscle": "#17becf",
        "T-Cell": "#aec7e8"
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
    
    return cell_type_colors, cell_type_fine_colors

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
        
def venn3_custom(set1: set, set2: set, set3: set, labels: Tuple[str], title: str = "", 
                 cmap: str = "Reds", ax: Optional[plt.Axes] = None, normalize_range: Tuple[int] = (0, 100),
                 set_label_size: int = 12, subset_label_size: int = 10, title_size: int = 12) -> plt.Axes:
    """Create a custom Venn diagram for three sets."""
    subset_sizes = {
        "100": len(set1 - set2 - set3),
        "010": len(set2 - set1 - set3),
        "001": len(set3 - set1 - set2), 
        "110": len(set1 & set2 - set3),  
        "101": len(set1 & set3 - set2),  
        "011": len(set2 & set3 - set1), 
        "111": len(set1 & set2 & set3),  
    }
    cmap = cm.get_cmap(cmap)
    norm = plt.Normalize(vmin=normalize_range[0], vmax=normalize_range[1])
    venn_color_dict = {key: cmap(norm(value)) for key, value in subset_sizes.items()}
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    
    # Draw Venn diagram
    venn = venn3_unweighted([set1, set2, set3], set_labels=labels)
    # 🔹 Apply gradient fill to each region
    for subset in ["100", "010", "001", "110", "101", "011", "111"]:
        patch = venn.get_patch_by_id(subset)
        if patch:
            patch.set_color(venn_color_dict[subset])
            patch.set_alpha(0.6)
            patch.set_edgecolor("black")
            
    for label in venn.set_labels:
        if label:
            label.set_fontsize(set_label_size)

    for label in venn.subset_labels:
        if label:
            label.set_fontsize(subset_label_size)
            
    ax.set_title(title, fontsize=title_size)

    return ax