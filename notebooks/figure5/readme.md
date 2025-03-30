## Figure 5: Implementation Notes

> **Note**: `.py` files are converted from Jupyter notebooks using **jupytext** prior to commit.

### **Panel Overview**
| Panel | Description | Visualization | Tool/Function |
|-------|-------------|---------------|---------------|
| **5A** | Heatmap of regulon activity across cell types for selected transcription factors | Global TF activity landscape | `pycomplexheatmap.ClusterMapPlotter()` |
| **5B** | Heatmap of regulon activity over time for selected TFs in key cell types | TF activity dynamics | `pycomplexheatmap.ClusterMapPlotter()` |
| **5C** | Line plots showing increased regulon activity over aging (m3 → m24) for upregulated TFs + GO enrichment of their target genes | Temporal activity trends + functional relevance | `seaborn.lineplot()`, `gprofiler.profile()` |
| **5D** | Line plots showing decreased regulon activity over aging for downregulated TFs + GO enrichment of their target genes | Temporal activity decline + functional relevance | `seaborn.lineplot()`, `gprofiler.profile()` |