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


