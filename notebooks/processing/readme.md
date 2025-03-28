**Note**: The .py files are coverted from .ipynb Jupyter notebooks using *jupytext* prior to commit.
# Data download
- **Data sources**: [GSE277702](https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE277702) for cell-ranger processed data files and [zenodo](https://zenodo.org/records/14888193) for scATAC-seq fragment files and R source scripts
- **data_download.py**: Downloading raw data files from GEO repository
-  **generate_feature_tsv.py**: Recreating features.tsv files from reference genome annotations due to them not being included in GEO repository
# scRNA-seq
## Single sample
- **scRNA_pilot_workflow.py**: Example notebook for processing a single scRNA-seq sample
- **/scripts/preprocessing.py**: Wrapper functions for preprocessing a single sample.

> sc -> Scanpy
> ad -> Anndata
> dc -> decoupler

```mermaid
graph TD
A(Load in adata object) -- sc.read_10x_mtx --> B(QC filters)
B -- sc.pp.calculate_qc_metrics <br> sc.pp.filter_cells --> C(Doublet Removal)
C -- sc.pp.scrublet <br> sc.pl.scrublet_score_distribution --> D(Normalization)
D -- sc.pp.normalize_total <br> sc.pp.lop1p --> E(Preprocessed sample)
```
## Multi samples
- **scRNA_processing.py**: Notebook for merging all scRNA-seq samples and performing downstream tasks like dimensional reduction, clustering, cell type annotations, and identification of DEGs.  
```mermaid
graph TD
A1(3month sample1) --> B1(Processed 3MS1)
A2(12month sample2) -- preprocess using workflow detailed above for single sample --> B2(Processed 12MS2)
A3(24month sample1) --> B3(Processed 24MS1)
B1 --> C(Merged adata object)
B2 -- ad.experimental.concat_on_disk --> C
B3 --> C
C --> D(Identify highly variable genes)
D -- sc.pp.highly_variable_genes --> E(Scale data)
E -- sc.pp.scale --> F("Principal component analysis (PCA)")
F -- sc.pp.pca --> G(Construct neighborhood graph)
G -- sc.pp.neighbors --> H(UMAP visualization)
G -- sc.pp.neighbors --> I(Leiden clustering)
H -- sc.tl.umap <br> sc.pl.umap --> J(Cell type annotation)
I -- sc.tl.leiden --> J
J -- dc.run_mlm --> K(Differentially expressed genes)
K -- sc.tl.rank_genes_group <br> sc.tl.filter_rank_genes_group --> L(Save as h5ad file)
```
# scATAC-seq
## Single sample
- **scATAC_pilot_workflow.py**: Example notebook for processing a single scATAC-seq sample
- **/scripts/scatac_preprocessing.py**: Wrapper functions for preprocessing a single scATAC-seq sample.
- **/scripts/scDblFinder_script.R**: R script for running scDblFinder

> sc -> Scanpy <br>
> ad -> Anndata <br>
> mu -> muon <br>
> ac -> mu.atac

```mermaid
graph TD
A1("matrix (.mtx)") --> B(adata object)
A2("peaks (.bed)") -- ad.AnnData --> B
A3("barcodes (.tsv)") --> B
B --> C1(Add fragment files)
B --> C2(Add peak annotations)
C1 -- ac.tools.locate_fragments --> D(Doublet Removal)
C2 -- "annotatePeaks.pl (HOMER) <br> ac.tools.add_peak_annotation" --> D
D -- "scDblFinder (R)" --> E(Quality Control Metrics)
E --> F1(General QC)
E --> F2(Nucleosome Signal)
E --> F3(TSS Enrichment)
F1 -- sc.pp.calculate_qc_metrics --> G(Filtering)
F2 -- ac.tl.nucleosome_signal <br> mu.pl.histogram --> G
F3 -- ac.tl.tss_enrichment <br> ac.pl.tss_enrichment --> G
G -- mu.pp.filter_obs --> H(Preprocessed sample)
```

## Multi samples