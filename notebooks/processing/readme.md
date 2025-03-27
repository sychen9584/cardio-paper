**Note**: The .py files are coverted from .ipynb Jupyter notebooks using *jupytext* prior to commit.
# Data download
- **Data sources**: [GSE277702](https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE277702) for cell-ranger processed data files and [zenodo](https://zenodo.org/records/14888193) for scATAC-seq fragment files and R source scripts
- **data_download.py**: Downloading raw data files from GEO repository
-  **generate_feature_tsv.py**: Recreating features.tsv files from reference genome annotations due to them not being included in GEO repository
# scRNA-seq
## Single sample
- **scRNA_pilot_workflow.py**: Example notebook for processing a single scRNA-seq sampel
- **/scripts/preprocessing.py**: Wrapper functions for processing a single sample.

```mermaid
graph TD
A(Load in adata object) -- sc.read_10x_mtx --> B(QC filters)
B -- sc.pp.calculate_qc_metrics <br> sc.pp.filter_cells --> C(Doublet Removal)
C -- sc.pp.scrublet <br> sc.pl.scrublet_score_distribution --> D(Normalization)
```

## Multi sample

# scATAC-seq
## Single sample
## Multi sample
