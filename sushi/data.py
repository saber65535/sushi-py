"""
sushi.data -- load example datasets (ported from R Sushi 1.32.0).

The original R Sushi package ships 14 example datasets in .rda files.
This module loads them from CSV/npz files in sushi/data/ (created from
the original .rda files via rdata parser at install time).

Usage:
    >>> import sushi
    >>> dnase = sushi.data.Sushi_DNaseI_bedgraph()   # returns pandas DataFrame
    >>> hic = sushi.data.Sushi_HiC_matrix()          # returns dict w/ matrix + positions

Each loader function reads from sushi/data/<name>.csv (or .npz for Hi-C
matrix) the first time it's called and caches the result.
"""

import os
import functools
import numpy as np
import pandas as pd

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@functools.lru_cache(maxsize=None)
def _load_csv(name: str) -> pd.DataFrame:
    """Load a CSV dataset, cached after first read."""
    path = os.path.join(_DATA_DIR, name + ".csv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Dataset {name!r} not found at {path}. "
            f"Did the package data/ folder get copied correctly?"
        )
    return pd.read_csv(path)


@functools.lru_cache(maxsize=None)
def _load_npz(name: str) -> dict:
    """Load a .npz dataset, cached after first read."""
    path = os.path.join(_DATA_DIR, name + ".npz")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Dataset {name!r} not found at {path}")
    z = np.load(path, allow_pickle=False)
    return {k: z[k] for k in z.files}


# === Dataset loaders ===
def Sushi_DNaseI_bedgraph() -> pd.DataFrame:
    """DNaseI hypersensitivity signal track on chr11 (hg18).

    Format: bedGraph = [chrom, start, end, value]
    Rows: 5,714
    """
    return _load_csv("Sushi_DNaseI.bedgraph")


def Sushi_ChIPSeq_CTCF_bedgraph() -> pd.DataFrame:
    """CTCF ChIP-seq signal track on chr11 (hg18).

    Rows: 27,576
    """
    return _load_csv("Sushi_ChIPSeq_CTCF.bedgraph")


def Sushi_ChIPSeq_pol2_bed() -> pd.DataFrame:
    """Pol2 ChIP-seq aligned reads (hg18).

    6-column BED (chrom, start, end, name, score, strand).
    Rows: 208
    """
    return _load_csv("Sushi_ChIPSeq_pol2.bed")


def Sushi_ChIPSeq_pol2_bedgraph() -> pd.DataFrame:
    """Pol2 ChIP-seq signal track on chr11 (hg18).

    Rows: 7,528
    """
    return _load_csv("Sushi_ChIPSeq_pol2.bedgraph")


def Sushi_ChIPExo_CTCF_bedgraph() -> pd.DataFrame:
    """CTCF ChIP-exo signal track (hg18).

    Rows: 6,292
    """
    return _load_csv("Sushi_ChIPExo_CTCF.bedgraph")


def Sushi_ChIPSeq_severalfactors_bed() -> pd.DataFrame:
    """Multiple ChIP-seq factor binding sites on chr15 (hg18).

    7-column BED with extra `row` and `name` columns.
    Rows: 130
    """
    return _load_csv("Sushi_ChIPSeq_severalfactors.bed")


def Sushi_RNASeq_K562_bedgraph() -> pd.DataFrame:
    """RNA-Seq K562 signal track (hg18).

    Rows: 1,073
    """
    return _load_csv("Sushi_RNASeq_K562.bedgraph")


def Sushi_5C_bedpe() -> pd.DataFrame:
    """5C interactions in K562/HeLa/GM12878 (hg18).

    11-column BEDPE with samplenumber (1/2/3) for cell line.
    Rows: 4,787
    """
    return _load_csv("Sushi_5C.bedpe")


def Sushi_ChIAPET_pol2_bedpe() -> pd.DataFrame:
    """Pol2 ChIA-PET interactions in K562 (hg18).

    Rows: 48,634
    """
    return _load_csv("Sushi_ChIAPET_pol2.bedpe")


def Sushi_HiC_matrix() -> dict:
    """Hi-C interaction matrix on chr11 (hg18, Dixon et al. 2012).

    Returns a dict with:
        matrix: 114 x 114 numpy array of interaction scores
        positions: 114-element array of genomic positions (in bp)
    """
    npz_path = os.path.join(_DATA_DIR, "Sushi_HiC.matrix.npz")
    csv_path = os.path.join(_DATA_DIR, "Sushi_HiC.matrix.csv")
    if os.path.exists(npz_path):
        return _load_npz("Sushi_HiC.matrix")
    elif os.path.exists(csv_path):
        # The CSV was saved from the rdata-parsed DataFrame with index=True
        # and header=True, so:
        #   - First column = genomic positions (rownames)
        #   - First row    = genomic positions (colnames)
        #   - Body         = 114 x 114 square interaction matrix
        df = pd.read_csv(csv_path, index_col=0)
        positions = np.asarray(df.index, dtype=float)
        return {"matrix": df.values, "positions": positions,
                "col_positions": np.asarray(df.columns, dtype=float)}
    else:
        raise FileNotFoundError(
            f"Neither Sushi_HiC.matrix.npz nor .csv found in {_DATA_DIR}")


def Sushi_GWAS_bed() -> pd.DataFrame:
    """GWAS blood-pressure variants (hg18, Ehret et al. 2011).

    6-column BED with p-value in col 5.
    Rows: 32,760
    """
    return _load_csv("Sushi_GWAS.bed")


def Sushi_genes_bed() -> pd.DataFrame:
    """Human gene annotations on chr15 (hg18).

    Rows: 5
    """
    return _load_csv("Sushi_genes.bed")


def Sushi_transcripts_bed() -> pd.DataFrame:
    """Human transcript annotations on chr15 (hg18) with FPKM scores.

    Rows: 143
    """
    return _load_csv("Sushi_transcripts.bed")


def Sushi_hg18_genome() -> pd.DataFrame:
    """hg18 (NCBI36) chromosome lengths.

    Returns a DataFrame with columns: chrom, size.
    Rows: 22
    """
    df = _load_csv("Sushi_hg18_genome")
    return df.rename(columns={"V1": "chrom", "V2": "size"})


# Expose all loaders in __all__
__all__ = [
    "Sushi_DNaseI_bedgraph",
    "Sushi_ChIPSeq_CTCF_bedgraph",
    "Sushi_ChIPSeq_pol2_bed",
    "Sushi_ChIPSeq_pol2_bedgraph",
    "Sushi_ChIPExo_CTCF_bedgraph",
    "Sushi_ChIPSeq_severalfactors_bed",
    "Sushi_RNASeq_K562_bedgraph",
    "Sushi_5C_bedpe",
    "Sushi_ChIAPET_pol2_bedpe",
    "Sushi_HiC_matrix",
    "Sushi_GWAS_bed",
    "Sushi_genes_bed",
    "Sushi_transcripts_bed",
    "Sushi_hg18_genome",
]
