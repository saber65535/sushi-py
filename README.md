# Sushi_Python_Version

**A Python port of [Sushi](https://github.com/PhanstielLab/Sushi) (Phanstiel et al., 2015, Bioinformatics) — flexible, quantitative, publication-quality genomic visualizations.**

Sushi R package is a Bioconductor package for plotting BED/BEDGraph/BEDPE/Hi-C/gene
data in the style of multi-panel paper figures. This Python version reproduces
the full public API (16 functions + 14 example datasets) using **matplotlib**.

## Status

| Component | Status |
| --- | --- |
| `plotBedgraph` (with colorbycol gradient, overlay, flip) | ✅ |
| `plotBed` (region / circles / density / splitstrand) | ✅ |
| `plotBedpe` (loops / ribbons / lines) | ✅ |
| `plotHic` (lower-triangle heatmap) | ✅ |
| `plotGenes` (box / arrow, colorby FPKM) | ✅ |
| `plotManhattan` (single + multi-chrom) | ✅ |
| `labelgenome`, `addlegend`, `zoombox`, `zoomsregion`, `labelplot` | ✅ |
| `SushiColors`, `maptocolors`, `maptolwd`, `opaque` | ✅ |
| `sortChrom`, `chromOffsets`, `convertstrandinfo` | ✅ |
| 14 example datasets ported from .rda → CSV | ✅ |
| Full vignette reproducing R `inst/doc/Sushi.R` (15 examples) | ✅ 15/15 pass |
| PyPI / pip install | ⏳ see `docs/` for install |

Tested with **Python 3.10**, **matplotlib 3.10**, **numpy 1.26**, **pandas 2.2**, **scipy 1.15**.

## Quick start

```python
import matplotlib.pyplot as plt
import sushi
from sushi import plotBedgraph, labelgenome, SushiColors

# Load an example bedGraph track (DNaseI hypersensitivity)
bg = sushi.data.Sushi_DNaseI_bedgraph()      # pandas DataFrame (chrom, start, end, value)

# Plot it on a fresh figure
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)             # leave space for axis labels
plotBedgraph(ax, bg, "chr11", 1650000, 2350000,
             colorbycol=SushiColors(5))      # gradient fill
labelgenome(ax, "chr11", 1650000, 2350000,
             n=4, scale="Mb", side=1)
fig.savefig("dnasei.png", dpi=120, bbox_inches="tight")
```

Output:
```
[ a DNaseI signal track on chr11:1.65-2.35Mb, colored with SushiColors(5) ]
```

## API at a glance

### Top-level imports
```python
from sushi import (
    # plotters
    plotBedgraph, plotBed, plotBedpe, plotHic, plotGenes, plotManhattan,
    # axes / overlays
    labelgenome, addlegend, zoombox, zoomsregion, labelplot,
    # helpers
    SushiColors, maptocolors, maptolwd, opaque,
    sortChrom, chromOffsets, convertstrandinfo,
    # data
    data,                            # submodule
)
```

### R vs Python: signature mapping

| R (`Sushi::`) | Python (`sushi.`) | Notes |
| --- | --- | --- |
| `plotBedgraph(signal, chrom, cs, ce, ...)` | `plotBedgraph(ax, signal, chrom, cs, ce, ...)` | Python takes **ax** as first arg; R creates a plot in place |
| `plotBed(beddata, chrom, cs, ce, type=...)` | `plotBed(ax, beddata, chrom, cs, ce, type=...)` | same |
| `plotBedpe(bedpe, chrom, cs, ce, plottype=...)` | `plotBedpe(ax, bedpe, ...)` | same |
| `plotHic(hic, chrom, cs, ce, max_y, zrange, palette, flip)` | `plotHic(ax, hic, ...)` | returns `[zrange, palette]` for use with addlegend |
| `plotGenes(geneinfo, chrom, cs, ce, types, ...)` | `plotGenes(ax, geneinfo, ...)` | same |
| `plotManhattan(bedfile, pvalues, genome, col)` | `plotManhattan(ax, bedfile, pvalues, genome, col)` | same |
| `labelgenome(...)` | `labelgenome(ax, ...)` | same |
| `addlegend(range, palette, title, ...)` | `addlegend(ax, range_val=..., palette=..., ...)` | range renamed to range_val (R `range` is reserved in Python) |
| `SushiColors(N)` | `SushiColors(N)` | same — returns a callable palette |
| `maptocolors(vec, col, num, range)` | `maptocolors(vec, col, num, range=)` | same — `range=` keyword preserved |

R function calls `par(mfrow=...)` and `layout(...)` are replaced by Python
`fig, axes = plt.subplots(...)` + `fig.add_gridspec(...)`. The caller
controls figure layout.

## Documentation

- **[docs/QUICKSTART.md](docs/QUICKSTART.md)** — 5-minute tutorial with runnable examples
- **[docs/API.md](docs/API.md)** — full parameter reference for every function
- **[docs/CHANGELOG.md](docs/CHANGELOG.md)** — v0.1.0 release notes + R vs Python diff
- **[docs/examples/](docs/examples/)** — one .py file per vignette example
- **[docs/VIGNETTE.md](docs/VIGNETTE.md)** — side-by-side R vignette vs Python output

## Repo layout

```
Sushi_Python_Version/
├── sushi/                       # the Python package
│   ├── __init__.py             # public API re-exports
│   ├── _helpers.py             # pure helpers (no matplotlib)
│   ├── _axes.py                # labelgenome, addlegend, zoombox, zoomsregion, labelplot
│   ├── plotters.py             # the 6 plot functions
│   ├── data.py                 # dataset loaders (lru_cache on CSV reads)
│   ├── data/                   # 14 CSV files + Sushi_HiC.matrix.csv (port from .rda)
│   └── vignette.py             # reproduces R inst/doc/Sushi.R 15 examples
│
├── docs/                        # human-facing docs (for the user and for future agents)
│   ├── README.md               # this file (also at top-level for visibility)
│   ├── QUICKSTART.md
│   ├── API.md
│   ├── CHANGELOG.md
│   ├── VIGNETTE.md
│   └── examples/
│       └── *.py                # one per vignette example
│
├── output/                      # scratch directory for ad-hoc outputs
└── README.md                    # you are here
```

## What is ported vs not

| Sushi R feature | Python port | Notes |
| --- | --- | --- |
| BED / BEDGraph / BEDPE / Hi-C / genes / GWAS plotting | ✅ | same algorithms |
| `par(mfrow=...)` multi-panel layout | ✅ via matplotlib gridspec | caller responsibility |
| `sushiColors` palettes (2..7 stops) | ✅ | discrete stops, no interpolation |
| `maptocolors`, `maptolwd`, `opaque`, `labelplot`, `sortChrom`, `chromOffsets`, `convertstrandinfo` | ✅ | faithful |
| 14 example datasets | ✅ (CSV format) | data was in `.rda`; converted via the `rdata` Python library |
| `labelgenome`, `addlegend`, `zoombox`, `zoomsregion` | ✅ | — |
| `biomaRt` gene fetch (used in vignette example 33-35) | ❌ | out of scope; users fetch their own gene lists |
| `pdf()` / `png()` device drivers | ❌ | use `matplotlib.pyplot.savefig` directly |
| `sushiR::plotBedpe(..., method="interaction3D")` | ❌ | not in upstream R package either |
| `zoombox` 3D perspective | ❌ | 2D only, matches R `zoombox` |

## Credits

Ported from [Sushi 1.32.0](https://github.com/PhanstielLab/Sushi) by Douglas H Phanstiel et al., GPL (>= 2).
Original paper: *Sushi.r: flexible, quantitative, and integrative genomic visualizations for publication-quality multi-panel figures*, Bioinformatics (2015).

This port was developed in June 2026 by Qianli on the 师大服务器 (`10.68.162.201`).
