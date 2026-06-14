"""
Sushi Python Version -- A Python port of Sushi (Phanstiel et al., 2015).

Sushi is a flexible R/Bioconductor package for publication-quality
visualization of genomic data (BED, BEDGraph, BEDPE, Hi-C, gene models,
Manhattan plots). This is a Python implementation that uses matplotlib
instead of base-R graphics.

Modules:
    sushi._helpers      -- pure helpers (SushiColors, maptocolors, opaque, ...)
    sushi._axes         -- labelgenome, addlegend, zoombox, zoomsregion
    sushi.plotters      -- 6 plot functions (plotBedgraph, plotBed,
                           plotBedpe, plotHic, plotGenes, plotManhattan)
    sushi.data          -- load example datasets

Public API (also available via top-level imports):
    plotBedgraph, plotBed, plotBedpe, plotHic, plotGenes, plotManhattan
    labelgenome, addlegend, zoombox, zoomsregion, labelplot
    SushiColors, maptocolors, maptolwd, opaque
    sortChrom, chromOffsets, convertstrandinfo

Conventions:
    * Chromosome names are Python strings.
    * All genomic coordinates are 1-based inclusive (BED convention).
    * Each plot function takes a matplotlib Axes (caller's responsibility
      to call fig.subplots_adjust(...) and savefig).
    * Output paths are not handled by the plotters; the caller decides.

Version: 0.1.16 (2026-06-14)
"""

__version__ = "0.1.16"

from ._helpers import (
    SushiColors,
    maptocolors,
    maptolwd,
    opaque,
    convertstrandinfo,
    sortChrom,
    chromOffsets,
    labelplot,
)
from ._axes import (
    labelgenome,
    addlegend,
    zoombox,
    zoomsregion,
)
from .plotters import (
    plotBedgraph,
    plotBed,
    plotBedpe,
    plotHic,
    plotGenes,
    plotManhattan,
)

# data loading
from . import data

__all__ = [
    # version
    "__version__",
    # helpers
    "SushiColors",
    "maptocolors",
    "maptolwd",
    "opaque",
    "convertstrandinfo",
    "sortChrom",
    "chromOffsets",
    "labelplot",
    # axes
    "labelgenome",
    "addlegend",
    "zoombox",
    "zoomsregion",
    # plotters
    "plotBedgraph",
    "plotBed",
    "plotBedpe",
    "plotHic",
    "plotGenes",
    "plotManhattan",
    # data
    "data",
]
