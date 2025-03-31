## Figure 6: Implementation Notes

> **Note**: `.py` files are converted from Jupyter notebooks using **jupytext** prior to commit. 

### **Panel Overview**
| Panel | Description | Visualization | Tool/Function |
|-------|-------------|---------------|---------------|
| **6A** | Pie charts showing genomic distribution of accessible peaks by cell type | Peak annotation breakdown | `annotatePeaks.pl`, `matplotlib.pyplot.pie()` |
| **6B** | Venn diagram showing overlap of differentially accessible regions (DARs) among fibroblasts, endothelial cells, and macrophages | Shared vs. cell type-specific DARs | `matplotlib_venn.venn3()` |
| **6C** | Sequence logos of top enriched TF motifs in DARs | TF motif discovery | `findMotifsGenome.pl`, `cairosvg.svg2png()` |
| **6D** | Barplots showing average chromVAR deviation scores (TF activity) across aging timepoints (3, 12, 24 months) in each cell type | Motif accessibility dynamics | `pychromVAR.compute_deviations()`, `seaborn.catplot()` |