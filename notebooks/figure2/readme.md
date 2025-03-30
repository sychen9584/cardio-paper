## Figure 2: Implementation Notes

> **Note**: `.py` files are converted from Jupyter notebooks using **jupytext** prior to commit.

### **Panel Overview**
| Panel | Description | Visualization | Tool/Function |
|-------|-------------|---------------|---------------|
| **2A** | Histogram of SASP gene set AUC scores across all cells | Gene signature scoring | `decoupler` → *aucell.ipynb*, `seaborn.histplot()` |
| **2B** | UMAP of SASP scores across cells | Dimensionality reduction | `scanpy.pl.umap()` |
| **2C** | KDE + boxplots of SASP scores by age group | Distribution & comparison | `seaborn.kdeplot()`, `seaborn.boxplot()` |
| **2D** | SASP scores in selected cell types | Grouped comparison | `seaborn.boxplot()` |
| **2E–2G** | **Pathology-associated scores** (Cardiac fibrosis, Heart failure, Inflammation) across cell types | Gene signature scoring | `seaborn.boxplot()` |
