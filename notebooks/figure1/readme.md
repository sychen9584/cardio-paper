## Figure 1: Implementation Notes

> **Note**: `.py` files are converted from Jupyter notebooks using **jupytext** prior to commit.

### **Panel Overview**
| Panel | Description | Visualization | Tool/Function |
|-------|-------------|---------------|---------------|
| **1A** | Schematic workflow | *Not reproduced* | – |
| **1B** | UMAP from scRNA-seq | Dimensionality reduction | `scanpy.pl.umap()` |
| **1C** | Sankey diagram: coarse → fine annotations | Cell type transitions | `plotly.Sankey()` |
| **1D** | UMAP from scATAC-seq | Dimensionality reduction | `scanpy.pl.umap()` |
| **1E** | Stacked barplot of coarse annotations | Cell type proportions | `matplotlib.pyplot.bar()` |
| **1F** | Stacked barplot of fine annotations | Cell type proportions | `matplotlib.pyplot.bar()` |
| **1G** | Augur: perturbation AUC scores on UMAPs | Perturbation analysis | `pertpy` → *augur_analysis.ipynb* |
| **1H** | Scatterplot: AUGUR scores (3–12 mo vs. 12–24 mo) | Temporal comparison | `seaborn.scatterplot()` |
