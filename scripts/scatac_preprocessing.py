from typing import List, Union, Optional, Callable, Iterable
import numpy as np
import pandas as pd
from scanpy import logging
from anndata import AnnData
from mudata import MuData


def tss_enrichment(
    data: Union[AnnData, MuData],
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
    adata: AnnData,
    features: pd.DataFrame,
    extend_upstream: int = 1000,
    extend_downstream: int = 1000,
    barcodes: Optional[str] = None,
) -> AnnData:
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

    return AnnData(X=mx, obs=adata.obs, var=anno, dtype=int)


def _calculate_tss_score(data: AnnData, flank_size: int = 100, center_size: int = 1001):
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


def get_gene_annotation_from_rna(data: Union[AnnData, MuData]) -> pd.DataFrame:
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
