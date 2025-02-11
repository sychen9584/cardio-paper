import os
import logging
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import gc

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=True,
)


def preprocess_adata(data_path: str,
                     sample_name: str,
                     figure_path: str,
                     min_genes: int = 1500,
                     max_genes: int = 7500,
                     mt_threshold: int = 10) -> None:
    
    logging.info(f"Preprocessing {sample_name}...")
    
    # set up figure directory saving path
    sample_figure_path = os.path.join(figure_path, sample_name, 'preprocess')
    os.makedirs(sample_figure_path, exist_ok=True)
    
    # load data
    adata = load_data(data_path, sample_name)
    
    # QC metrics
    calculate_qc_metrics(adata, sample_figure_path)
    
    # Filter out cells with <1500 and >7500 detected genes
    adata = filter_cells(adata, min_genes, max_genes, mt_threshold)
    
    # doublet detection and removal
    logging.info(f"Running Scrublet for doublet detection...")
    sc.pp.scrublet(adata)
    sc.pl.scrublet_score_distribution(adata, show=False)
    save_figure("scrublet_score_distribution.png", sample_figure_path)
    
    adata = adata[adata.obs['predicted_doublet'] == False]
    
    # Normalizing to 10K counts
    normalize_and_log_transform(adata)
    
    # save preprocessed data
    save_preprocessed_data(adata, data_path, sample_name)
    
    # remove object in environment to save RAM space
    del adata
    gc.collect()
  
  
def load_data(data_path: str, sample_name: str) -> sc.AnnData:
    """Load the 10X Genomics data into an AnnData object."""
    logging.info(f"Loading data for {sample_name}...")
    adata = sc.read_10x_mtx(path=os.path.join(data_path, "raw", sample_name))
    return adata


def calculate_qc_metrics(adata: sc.AnnData, figure_path: str) -> None:
    """Calculate QC metrics and create QC plots."""
    logging.info("Calculating QC metrics...")
    
    adata.var["mt"] = adata.var_names.str.startswith("mt-")  # mitochondrial genes
    
    if not any(adata.var_names.str.startswith("mt-")):
        logging.warning("No mitochondrial genes detected in var_names.")
        
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)

    # QC Violin Plot
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show=False)
    save_figure("qc_violin_plot.png", figure_path)

    # Scatter Plot
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
    save_figure("qc_scatter_plot.png", figure_path)
    
    
def filter_cells(adata: sc.AnnData, min_genes: int, max_genes: int, mt_threshold: int) -> sc.AnnData:
    """Filter out low-quality cells based on gene counts and mitochondrial content."""
    logging.info("Filtering low-quality cells...")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    adata = adata[adata.obs['pct_counts_mt'] < mt_threshold]
    return adata


def normalize_and_log_transform(adata: sc.AnnData) -> None:
    """Normalize data to 10,000 counts and log-transform."""
    logging.info("Normalizing and log-transforming data...")
    
    adata.layers["counts"] = adata.X.copy()
    adata.layers["counts"] = adata.layers["counts"].tocsr()
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    adata.layers['lognormalized'] = adata.X.copy()
    adata.layers["lognormalized"] = adata.layers["lognormalized"].tocsr()
    
    adata.X = adata.X.tocsr()
    adata.raw = adata.copy()


def identify_highly_variable_genes(adata: sc.AnnData, figure_path: str, hvg_params: dict = None) -> None:
    """Identify highly variable genes and save a plot."""
    logging.info("Identifying highly variable genes...")
    hvg_params = hvg_params or {"n_top_genes": 2000, "flavor": "seurat_v3", 'layer': 'counts'}
    sc.pp.highly_variable_genes(adata, **hvg_params)
    
    sc.pl.highly_variable_genes(adata, show=False)
    save_figure("highly_variable_genes.png", figure_path)


def save_preprocessed_data(adata: sc.AnnData, data_path: str, sample_name: str) -> None:
    """Save the preprocessed AnnData object to disk."""
    logging.info(f"Saving preprocessed data for {sample_name}...")
    processed_path = os.path.join(data_path, "preprocessed", sample_name)
    os.makedirs(processed_path, exist_ok=True)
        
    #if "counts" in adata.layers:
        #del adata.layers["counts"]
        
    adata.write_h5ad(os.path.join(processed_path, f"{sample_name}_preprocessed.h5ad"))
    logging.info(f"Preprocessed data saved to {processed_path}.")


def save_figure(fig_name, figure_path):
    """Save the current Matplotlib figure."""
    plt.savefig(os.path.join(figure_path, fig_name), dpi=300, bbox_inches="tight")
    plt.close()