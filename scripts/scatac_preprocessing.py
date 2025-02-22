from typing import List, Union, Optional, Callable, Iterable
import os
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import anndata
import subprocess
import muon as mu
from muon import atac as ac

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=True,
)

def preprocess_atac_data(data_path: str,
                         sample_name: str,
                         figure_path: str,
                         gene_annot_file: str,
                         script_path: str,
                         nuc_signal_threshold: int = 4,
                         tss_threshold: int = 3,
                         nuc_params: dict = None,
                         tss_params: dict = None) -> None:
    """
    Preprocess ATAC-seq data for downstream analysis.
    """
    logging.info(f"Preprocessing {sample_name}...")
    
    # set up figure directory saving path
    sample_figure_path = os.path.join(figure_path, sample_name, 'preprocess')
    os.makedirs(sample_figure_path, exist_ok=True)
    
    # load data
    logging.info("Loading ATAC-seq data from directory...")
    adata = _load_atac_adata(data_path, sample_name)
    
    # locate fragments file
    logging.info("Locating fragments...")
    ac.tools.locate_fragments(adata, fragments=os.path.join(data_path, sample_name, "fragments.tsv.gz"))
        
    # HOMER annotation
    logging.info("Annotating peaks with HOMER...")
    homer_df = _homer_annotation(data_path, sample_name)
    ac.tools.add_peak_annotation(adata, annotation=homer_df)
    
    # Save raw data as h5ad for scDblFinder
    adata.write_h5ad(os.path.join(data_path, sample_name, "raw.h5ad"))
    
    # Run scDblFinder and append doublet scores to adata
    logging.info("Running scDblFinder for doublet detection...")
    _run_scdblfinder(script_path, data_path, sample_name)
    dbl_scores = pd.read_csv(os.path.join(data_path, sample_name, "doublet_scores.csv")).set_index('barcode')
    adata.obs['dbl_score'] = dbl_scores['doublet_score']
    adata.obs['dbl_class'] = dbl_scores['doublet_class']
    
    # Calculate QC metrics
    logging.info("Calculating scATAC-seq QC metrics...")
    _calculate_qc_metrics(adata, gene_annot_file, sample_figure_path, nuc_signal_threshold, tss_threshold, nuc_params, tss_params)
    
    # Save preprocessed data
    logging.info("Saving adata with QC metrics...")
    adata.write_h5ad(os.path.join(data_path, sample_name, "qc_metrics.h5ad"))
    

def _load_atac_adata(data_path: str, sample_name: str) -> anndata.AnnData:
    """
    Load ATAC-seq data from the specified directory and return an AnnData object.
    """
    # Define the directory containing the ATAC-seq data
    matrix_path = os.path.join(data_path, sample_name, "matrix.mtx.gz")

    # Load the sparse peak-barcode matrix
    mat = scipy.io.mmread(matrix_path).T.tocsc()  # Convert to CSC format for efficiency

    # Load peak information (BED file)
    peaks_path = os.path.join(data_path, sample_name, "peaks.bed.gz")
    peaks = pd.read_csv(peaks_path, sep="\t", header=None, names=["chrom", "start", "end"])

    # Load barcode metadata
    barcodes_path = os.path.join(data_path, sample_name, "barcodes.tsv.gz")
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, names=["barcode"])

    # Convert barcodes and peaks into the required format
    barcodes.index = barcodes["barcode"]  # Set barcodes as index
    peaks["peak_name"] = peaks['chrom'] + ":" + peaks['start'].astype(str) + "-" + peaks['end'].astype(str)
    peaks.index = peaks["peak_name"]  # Set peaks as index

    # Create AnnData object for ATAC-seq counts
    adata_atac = anndata.AnnData(
        X=mat,  # Sparse matrix
        obs=barcodes,  # Barcodes (cells)
        var=peaks  # Peaks (features)
    )
    
    return adata_atac

def _homer_annotation(data_path: str, sample_name: str) -> pd.DataFrame:
    """
    Annotate peaks with HOMER and wrangle the output into a 10X tsv format.
    """
    dir_path = os.path.join(data_path, sample_name)
    bed_file = os.path.join(dir_path, "peaks.bed.gz")
    bed_unzipped_file = os.path.join(dir_path, "peaks.bed")
    annot_file = os.path.join(dir_path, "annotated_peaks.txt")
    output_file = os.path.join(dir_path, "peak_annotations.tsv")
    
    # Check if output file already exists
    if os.path.exists(output_file):
        logging.info(f"Output file {output_file} already exists. Loading the file.")
        return pd.read_csv(output_file, sep="\t")
    
    # Ensure input file exists
    if not os.path.exists(bed_file):
        raise FileNotFoundError(f"Input file {bed_file} does not exist.")
    

    # HOMER Bash script
    bash_script = f"""
    echo "Annotating peaks with HOMER..." && \
    gunzip -c {bed_file} > {bed_unzipped_file} && \
    annotatePeaks.pl {bed_unzipped_file} mm10 > {annot_file} && \
    """
    # Run the script
    result = subprocess.run(bash_script, shell=True, text=True, capture_output=True)

    # Check for errors
    if result.returncode != 0:
        raise RuntimeError(f"Error in HOMER annotation: {result.stderr}")

    # Load HOMER annotation output and wrangle into 10X tsv format
    homer_df = pd.read_csv(annot_file, sep="\t")
    # Ensure required columns exist dynamically
    required_columns = ["Chr", "Start", "End", "Nearest Ensembl", "Gene Name", "Distance to TSS", "Annotation"]
    missing_columns = [col for col in required_columns if col not in homer_df.columns]
    if missing_columns:
        raise ValueError(f"Missing expected columns in HOMER output: {missing_columns}")
    
    homer_df = homer_df.iloc[:, required_columns]
    homer_df.columns = ['chrom', 'start', 'end', 'gene_id', 'gene', 'distance', 'peak_type']
    homer_df['peak_type'] = homer_df['peak_type'].str.extract(r'([^\(]+)').squeeze().str.strip()
    
    # Standardize peak type classifications
    peak_type_map = {
        'intron': 'distal',
        'exon': 'distal',
        'promoter-TSS': 'promoter',
        "5' UTR": 'distal',
        'non-coding': 'distal',
        'TTS': 'distal',
        "3' UTR": 'distal'
    }
    homer_df['peak_type'] = homer_df['peak_type'].replace(peak_type_map).fillna('intergenic')
    homer_df['start'] -= 1 # Convert to 0-based indexing
    
    homer_df.to_csv(output_file, sep="\t", index=False)
    logging.info(f"Annotated peaks saved to {output_file}")

    return homer_df
        
        
def _run_scdblfinder(script_path: str, data_path: str, sample_name: str) -> pd.DataFrame:
    '''
    Run scDBLFinder on raw ATAC-seq data and save the output to a file.
    '''
    if not os.path.exists(script_path):
        raise FileNotFoundError(f"R script not found: {script_path}")
    
    input_h5ad = os.path.join(data_path, sample_name, "raw.h5ad")
    
    command = ['Rscript', script_path, input_h5ad]
    logging.info(command)
    result = subprocess.run(command, text=True, capture_output=True)
    
    if result.returncode == 0:
        logging.info("scDblFinder executed successfully!")
    else:
        raise RuntimeError(f"Error running scDblFinder R script: {result.stderr}")
    
def _calculate_qc_metrics(adata: anndata.AnnData, 
                          gene_annot_file: str, 
                          figure_path: str,
                          nuc_signal_threshold: int = 4,
                          tss_threshold: int = 3,
                          nuc_params: dict = None,
                          tss_params: dict = None) -> None:
    
    """
    Calculate scATAC-seq QC metrics and generate diagnostic plots.
    """
    default_n = 10e3*adata.n_obs
    nuc_params = nuc_params or {"n": default_n}
    tss_params = tss_params or {"n_tss": 3000, "random_state": 42, 'layer': 'counts', 'extend_upstream': 2000, 'extend_downstream': 2000}
    
    # Calculate general qc metrics using scanpy
    logging.info('--- General QC metrics')
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    # log-transform total counts and add as column
    adata.obs["log_total_fragment_counts"] = np.log10(adata.obs["total_fragment_counts"])

    # Rename columns
    adata.obs.rename(
        columns={
            "n_genes_by_counts": "n_features_per_cell",
            "total_counts": "total_fragment_counts",
        },
        inplace=True,
    )
    
    # QC Violin Plot
    sc.pl.violin(adata, ['total_fragment_counts', 'n_features_per_cell'], jitter=0.4, multi_panel=True, show=False)
    _save_figure("qc_violin_plot.png", figure_path)
    
    # Calculate nucleosome signal distribution across cells
    logging.info('--- Nucleosome signal')
    ac.tl.nucleosome_signal(adata, **nuc_params)
    mu.pl.histogram(adata, "nucleosome_signal", show=False)
    _save_figure("nucleosome_signal", figure_path)
    
    # Add group labels for above and below the nucleosome signal threshold
    adata.obs["nuc_signal_filter"] = [
        "NS_FAIL" if ns > nuc_signal_threshold else "NS_PASS"
        for ns in adata.obs["nucleosome_signal"]
    ]
    
    # TSS enrichment
    logging.info('--- TSS enrichment')
    gene_intervals = pd.read_csv(gene_annot_file, sep="\t")
    tss = tss_enrichment(adata, features=gene_intervals, **tss_params)
    
    fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))
    # histogram of the TSS scores
    p1 = sns.histplot(adata.obs, x="tss_score", ax=axs[0])
    p1.set_title("Full range")

    p2 = sns.histplot(
        adata.obs,
        x="tss_score",
        binrange=(0, adata.obs["tss_score"].quantile(0.95)),
        ax=axs[1],
    )
    p2.axvline(x=3)
    p2.set_title("Up to 95% percentile")
    plt.suptitle("Distribution of the TSS score")
    plt.tight_layout()
    _save_figure("tss_score_distribution.png", figure_path)
    
    tss.obs["tss_filter"] = [
        "TSS_FAIL" if score < tss_threshold else "TSS_PASS"
        for score in adata.obs["tss_score"]
    ]

    # Temporarily set different color palette
    sns.set_palette(palette="Set1")
    ac.pl.tss_enrichment(tss, color="tss_filter", show=False)
    _save_figure("tss_enrichment.png", figure_path)
    # reset color palette
    sns.set_palette(palette="tab10")

    
def _save_figure(fig_name, figure_path):
    """Save the current Matplotlib figure."""
    plt.savefig(os.path.join(figure_path, fig_name), dpi=300, bbox_inches="tight")
    plt.close()


def tss_enrichment(
    data: Union[anndata.AnnData, mu.MuData],
    features: Optional[pd.DataFrame] = None,
    extend_upstream: int = 1000,
    extend_downstream: int = 1000,
    n_tss: int = 2000,
    return_tss: bool = True,
    random_state=None,
    barcodes: Optional[str] = None,
):
    """
    Calculate TSS enrichment according to ENCODE guidelines. Adds a column tss_score to the .obs DataFrame and
    optionally returns a tss score object.

    Parameters
    ----------
    data
        AnnData object with peak counts or multimodal MuData object with 'atac' modality.
    features
        A DataFrame with feature annotation, e.g. genes.
        Annotation has to contain columns: Chromosome, Start, End.
    extend_upsteam
        Number of nucleotides to extend every gene upstream (2000 by default to extend gene coordinates to promoter regions)
    extend_downstream
        Number of nucleotides to extend every gene downstream (0 by default)
    n_tss
        How many randomly chosen TSS sites to pile up. The fewer the faster. Default: 2000.
    return_tss
        Whether to return the TSS pileup matrix. Needed for enrichment plots.
    random_state : int, array-like, BitGenerator, np.random.RandomState, optional
        Argument passed to pandas.DataFrame.sample() for sampling features.
    barcodes
        Column name in the .obs of the AnnData
        with barcodes corresponding to the ones in the fragments file.

    Returns
    ----------
    AnnData
        AnnData object with a 'tss_score' column in the .obs slot.


    """
    if isinstance(data, AnnData):
        adata = data
    elif isinstance(data, MuData) and "atac" in data.mod:
        adata = data.mod["atac"]
    else:
        raise TypeError("Expected AnnData or MuData object with 'atac' modality")

    if features is None:
        # Try to gene gene annotation in the data.mod['rna']
        if (
            isinstance(data, MuData)
            and "rna" in data.mod
            and "interval" in data.mod["rna"].var.columns
        ):
            features = get_gene_annotation_from_rna(data)
        else:
            raise ValueError(
                "Argument features is required. It should be a BED-like DataFrame with gene coordinates and names."
            )

    if features.shape[0] > n_tss:
        # Only use n_tss randomly chosen sites to make function faster
        features = features.sample(n=n_tss, random_state=random_state)

    # Pile up tss regions
    tss_pileup = _tss_pileup(
        adata,
        features,
        extend_upstream=extend_upstream,
        extend_downstream=extend_downstream,
        barcodes=barcodes,
    )

    flank_means, center_means = _calculate_tss_score(data=tss_pileup)

    tss_pileup.X = tss_pileup.X / flank_means[:, None]

    tss_scores = center_means / flank_means

    adata.obs["tss_score"] = tss_scores
    tss_pileup.obs["tss_score"] = tss_scores

    if isinstance(data, AnnData):
        logging.info('Added a "tss_score" column to the .obs slot of the AnnData object')
    else:
        logging.info("Added a \"tss_score\" column to the .obs slot of tof the 'atac' modality")

    if return_tss:
        return tss_pileup

    
def _tss_pileup(
    adata: anndata.AnnData,
    features: pd.DataFrame,
    extend_upstream: int = 1000,
    extend_downstream: int = 1000,
    barcodes: Optional[str] = None,
) -> anndata.AnnData:
    """
    Pile up reads in TSS regions while accounting for gene strand orientation.
    """

    if "files" not in adata.uns or "fragments" not in adata.uns["files"]:
        raise KeyError("No fragments file found. Run `muon.atac.tl.locate_fragments` first.")

    try:
        import pysam
    except ImportError:
        raise ImportError(
            "pysam is required to work with the fragments file. Install it via `pip install pysam`."
        )

    n = adata.n_obs
    n_features = extend_downstream + extend_upstream + 1

    # Dictionary to map barcodes to indices
    if barcodes and barcodes in adata.obs.columns:
        d = {k: v for k, v in zip(adata.obs.loc[:, barcodes], range(n))}
    else:
        d = {k: v for k, v in zip(adata.obs.index, range(n))}

    # Matrix to store TSS pileup
    mx = np.zeros((n, n_features), dtype=int)

    fragments = pysam.TabixFile(adata.uns["files"]["fragments"], parser=pysam.asBed())

    # Ensure only valid chromosomes are used
    chromosomes = fragments.contigs
    features = features[features["Chromosome"].isin(chromosomes)]

    for i in range(features.shape[0]):  # Iterate over TSS sites
        f = features.iloc[i]
        
        # **Fix: Adjust for strand**
        if f["Strand"] == "+":  # Forward strand
            tss_start = f.Start - extend_upstream
            tss_end = f.Start + extend_downstream
        else:  # Reverse strand (`-`)
            tss_start = f.End - extend_downstream
            tss_end = f.End + extend_upstream

        # Fetch overlapping fragments
        for fr in fragments.fetch(f.Chromosome, tss_start, tss_end):
            try:
                rowind = d[fr.name]  # Cell barcode
                score = int(fr.score)  # Number of cuts per fragment
                colind_start = max(fr.start - tss_start, 0)
                colind_end = min(fr.end - tss_start, n_features)
                mx[rowind, colind_start:colind_end] += score
            except:
                pass  # Ignore missing barcodes

    fragments.close()

    # Annotate positions
    anno = pd.DataFrame(
        {"TSS_position": range(-extend_upstream, extend_downstream + 1)},
    )
    anno.index = anno.index.astype(str)

    return anndata.AnnData(X=mx, obs=adata.obs, var=anno, dtype=int)


def _calculate_tss_score(data: anndata.AnnData, flank_size: int = 100, center_size: int = 1001):
    """
    Calculate TSS enrichment scores (defined by ENCODE) for each cell.

    Parameters
    ----------
    data
        AnnData object with TSS positons as generated by tss_pileup.
    flank_size
        Number of nucleotides in the flank on either side of the region (ENCODE standard: 100bp).
    center_size
        Number of nucleotides in the center on either side of the region (ENCODE standard: 1001bp).
    """
    region_size = data.X.shape[1]

    if center_size > region_size:
        raise ValueError(
            f"center_size ({center_size}) must smaller than the piled up region ({region_size})."
        )

    if center_size % 2 == 0:
        raise ValueError(f"center_size must be an uneven number, but is {center_size}.")

    # Calculate flank means
    flanks = np.hstack((data.X[:, :flank_size], data.X[:, -flank_size:]))
    flank_means = flanks.mean(axis=1)

    # Replace 0 means with population average (to not have 0 division after)
    flank_means[flank_means == 0] = flank_means.mean()

    # Calculate center means
    center_dist = (region_size - center_size) // 2  # distance from the edge of data region
    centers = data.X[:, center_dist:-center_dist]
    center_means = centers.mean(axis=1)

    return flank_means, center_means


def get_gene_annotation_from_rna(data: Union[anndata.AnnData, mu.MuData]) -> pd.DataFrame:
    """
    Get data frame with start and end positions from interval
    column of the 'rna' layers .var.

    Parameters
    ----------
    mdata: MuData
            MuData object
    """

    if isinstance(data, AnnData):
        adata = data
    elif isinstance(data, MuData) and "rna" in data.mod:
        adata = data.mod["rna"]
    else:
        raise TypeError("Expected AnnData or MuData object with 'rna' modality")

    if "interval" in adata.var.columns:
        features = pd.DataFrame([s.replace(":", "-", 1).split("-") for s in adata.var.interval])
        features.columns = ["Chromosome", "Start", "End"]
        features["gene_id"] = adata.var.gene_ids.values
        features["gene_name"] = adata.var.index.values
        features.index = adata.var.index
        # Remove genes with no coordinates indicated
        features = features.loc[~features.Start.isnull()]
        features.Start = features.Start.astype(int)
        features.End = features.End.astype(int)
    else:
        raise ValueError(".var object does not have a column named interval")
    return features
