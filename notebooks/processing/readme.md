**Note**: The .py files are coverted from .ipynb Jupyter notebooks using *jupytext* prior to commit.
# Data download
- **Data sources**: [GSE277702](https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE277702) for cell-ranger processed data files and [zenodo](https://zenodo.org/records/14888193) for scATAC-seq fragment files and R source scripts
- **data_download.py**: Downloading raw data files from GEO repository
-  **generate_feature_tsv.py**: Recreating features.tsv files from reference genome annotations due to them not being included in GEO repository
# scRNA-seq
## Single sample
- **scRNA_pilot_workflow.py**: Example notebook for processing a single scRNA-seq sample
- **/scripts/preprocessing.py**: Wrapper functions for preprocessing a single sample.

> sc -> Scanpy <br>
> ad -> Anndata <br>
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
> ac -> mu.atac <br>
> scvi -> scvi-tools

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
- **scATAC_processing.py**: Notebook for merging all scATAC-seq samples and performing downstream tasks like dimensional reduction, clustering, and cell type annotation transfers.
- **scanvi_scatac_annot.py**: Notebook for annotating cell types in scATAC-seq data by learning and transfering labels from scRNA-seq data using scVI and scANVI packages. This notebook was ran on Google Colab to take advantage of GPU computing.
- **/scripts/merge_fragment_files.sh**: Bash script for merging all fragment files from individual samples into one prior to merging their adata objects.

### Merging and visualization of scATAC-seq data
```mermaid
graph TD
A1(3month sample1) --> B1(Processed 3MS1)
A2(12month sample2) -- preprocess using workflow detailed above for single sample --> B2(Processed 12MS2)
A3(24month sample1) --> B3(Processed 24MS1)
B1 --> C1(Unified 3MS1)
B2 -- Create a unified peak set for all samples using pyranges  --> C2(Unified 12MS2)
B3 --> C3(Unified 24MS1)
C1 --> D
C2 -- ad.experimental.concat_on_disk --> D(Merged adata object)
C3 --> D
D --> E(TF-IDF normalization)
E -- ac.pp.tfidf --> F("Latent semantic indexing embedding (LSI)")
F -- sc.tl.lsi --> G(Harmony batch correction)
G -- sc.external.pp.harmony_integrate --> H(Construct neighborhood graph)
H -- sc.pp.neighbors --> I1(UMAP visualization)
H -- sc.pp.neighbors --> I2(Leiden clustering)
I1 -- sc.tl.umap <br> sc.pl.umap --> J
I2 -- sc.tl.leiden --> J(Save as h5ad file) 
```
### Computation of gene activity matrix and cell type label transfer from scRNA-seq data
Basically transforming scATAC-seq data from *cell X peak matrix* to *cell X gene matrix* so it can be compared directly to scRNA-seq data for label transfer. 
```mermaid
graph TD
A(Merged adata object) --> B(Compute gene activity matrix)
B -- Custom function in scatac_preprocessing.py <br> that group counts in unified peaks to their nearest genes <br> and return a new adata object --> C(Gene activity adata object)
C --> D(Normalization)
D -- sc.pp.normalize_total <br> sc.pp.log1p --> E1(scATAC-seq object)
E1 --> F(Subset to only HVGs and concat as one adata object)
E2(scRNA-seq data object) --> F
F -- ad.concat --> G(Train scVI model jointly on scRNA and scATAC data)
G -- "vae=scvi.model.SCVI() <br> vae.train()" --> H(Train scANVI model based on scVI weights)
H -- "scanvi=scvi.model.SCANVI.from_scvi_model <br> scanvi.train()" --> I(Predict cell type labels for scATAC-seq cells)
I -- "scanvi.predict()" --> J(Save labels to h5ad)
```
