# Replicating "Decoding aging in the heart via single cell dual omics of non-cardiomyocytes" in Python 🐍

This project aims to replicate the analysis from "**[Replicating "Decoding aging in the heart via single cell dual omics of non-cardiomyocytes](https://www.cell.com/iscience/fulltext/S2589-0042(24)02696-8)**" using Python instead of R.
The original study applied single-cell RNA-seq (scRNA-seq) and single-cell ATAC-seq (scATAC-seq) analysis in R, utilizing packages such as [Seurat](https://satijalab.org/seurat/) and [Signac](https://stuartlab.org/signac/).
This repository provides a step-by-step Python implementation using the [scverse](https://scverse.org/) ecosystem, which consists of Scanpy, muon, and other related libraries.

## Highlights from the Song et al.
- Single cell dual-omics profiling reveals nonmyocyte heterogeneity in heart aging
- Aging non-cardiomyocytes show cell-type-specific transcriptomics and epigenomic shifts
- Senescence-associated secretory phenotype (SASP) elevates in macrophages and fibroblasts
- Fibroblast subcluster with ERBB4 expression shows unique aging impact and cell interactions
  
<p align="center">
<img src="https://github.com/user-attachments/assets/7da3ed2e-3b76-4ebd-9349-85246e3f3ce0" width=50% height=50%>
</p>

## Goals of this project
- Reproduce key figures and results from the original paper
- Compare results between R-based and Python-based workflows
- Explore Python alternatives for single-cell genomics analysis

## Progress
- [x] Data download
- [x] scRNA-seq processing
- [x] scATAC-seq processing
- [x] [Figure 1: Integrated single-cell transcriptomic and epigenomic anlaysis of no-cardiac (non-CM) cells during aging](figures/figure1.png)
- [X] [Figure 2: Gene signatures profiling of aging-susceptible of cardiac non-myocyte changes](figures/figure2.png)
- [x] [Figure 3: Differentially expressed genes and cluster analysis during aging](figures/figure3.png)
- [X] [Figure 4: Cell-cell communication during aging](figures/figure4.png)
- [X] [Figure 5: Gene regulatory networks during aging by single-cell regulatory network inference and clustering (SCENIC)](figures/figure5.png)
- [X] [Figure 6: Integrated single-cell ATAC seq epigenomic analysis of non-CM cells during aging](figures/figure6.png)
- [X] Implementation notes for each figure (README.md files in /notebooks/)
- [ ] Final summary write up
      
## Comparison of tech stack

| Task  | R | Python  |
| ------------- | ------------- | ------------- |
| scRNA-seq data structure + processing  | Seurat | Scanpy |
| scATAC-seq data structure + processing  | Signac | Muon + Custom codes |
| Label transfer | Signac | scANVI | 
| General Vis. | ggplot2 | Matplotlib + Seaborn |
| Heatmap Vis. | ComplexHeatmap | pyComplexHeatmap |
| Venn diagram | VennDiagram | matplotlib-venn |
| Perturbation analysis | Augur | pertpy |
| Gene set enrichment | AUCell | decoupler |
| Gene Ontology | clusterProfiler | gprofiler-official |
| Cell-cell communication | CellChat | liana |
| GRN inference | SCENIC | pySCENIC |
| TF motif analysis | Signac | HOMER + pychromVAR |

## References
Song, Yiran, et al. "Decoding aging in the heart via single cell dual omics of non-cardiomyocytes." iScience 27.12 (2024).
