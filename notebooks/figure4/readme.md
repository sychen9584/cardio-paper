## Figure 4: Implementation Notes

> **Note1**: `.py` files are converted from Jupyter notebooks using **jupytext** prior to commit.
> **Note 2**: The original **R implementation of CellChat** was used instead of **liana** in Python because:  
> • The Python version produced very different interaction results, relying more heavily on heuristics.  
> • liana's CellChat does not offer the visualization functions needed to replicate Figure 4 panels from the original publication.

### **Panel Overview**
| Panel | Description | Visualization | Tool/Function |
|-------|-------------|---------------|---------------|
| **4A** | Barplots showing total number of inferred interactions and overall interaction strength across age groups (3, 12, 24 months) | Global communication metrics | `CellChat::computeNetStats()`, `ggplot2::geom_bar()` |
| **4B** | Circle plots of global intercellular signaling networks for 3 and 24 months | Network topology | `CellChat::netVisual_circle()` |
| **4C** | Heatmaps showing differential interaction strength between sender and receiver cell types across timepoints (3 vs 12 mo, 12 vs 24 mo) | Pairwise signaling comparison | `CellChat::netVisual_heatmap()` |
| **4D** | Heatmaps of senescence-related pathways (e.g., IGF, TGFβ) over time | Pathway-specific dynamics | `CellChat::netVisual_heatmap()` |
| **4E** | Heatmaps of inflammation-related signaling pathways (e.g., CCL, IL1) | Immune signaling variation | `CellChat::netVisual_heatmap()` |
| **4F** | Heatmaps of fibroblast-related signaling pathways (e.g., PTN, MK) | Fibroblast-centric communication | `CellChat::netVisual_heatmap()` |