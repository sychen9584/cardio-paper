## Figure 3: Implementation Notes

> **Note**: `.py` files are converted from Jupyter notebooks using **jupytext** prior to commit.

### **Panel Overview**
| Panel | Description | Visualization | Tool/Function |
|-------|-------------|---------------|---------------|
| **3A** | Heatmap of marker genes across cell types | Cell type markers | `pycomplexheatmap.ClusterMapPlotter()`, `HeatmapAnnotation()` |
| **3B** | Barplot of enriched GO terms per cell type | Gene ontology enrichment | `gprofiler.profile()`, `seaborn.barplot()` |
| **3C** | Venn diagram of DEGs across fibroblasts, macrophages, endothelial cells (3 vs 12 mo, 12 vs 24 mo) | DEG overlap | `matplotlib_venn.venn3()` |
| **3D** | K-means clustering heatmap of DEG fold-changes in key cell types | Unsupervised clustering | `sklearn.cluster.KMeans()`, `pycomplexheatmap.ClusterMapPlotter()` |
| **3E** | GO term enrichment for select K-means clusters | Cluster-level functional annotation | `gprofiler.profile()`, `seaborn.barplot()` |
