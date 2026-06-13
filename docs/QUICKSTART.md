# Sushi Python — Quick Start

A 5-minute tour of every plot function. Each block is a complete, runnable
script. Save any block as `foo.py` and run it; output goes to `/tmp/foo.png`.

## Setup

```bash
# Add the package to PYTHONPATH
export PYTHONPATH=/datapool/life-zhanghk/qianli/39-Sushi_Python:$PYTHONPATH
python -c "import sushi; print(sushi.__version__)"   # → 0.1.0
```

## 1. Plot a bedGraph with a gradient fill

```python
import matplotlib.pyplot as plt
import sushi
from sushi import plotBedgraph, labelgenome, SushiColors

bg = sushi.data.Sushi_DNaseI_bedgraph()      # DNaseI hypersensitivity on chr11

fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedgraph(ax, bg, "chr11", 1650000, 2350000, colorbycol=SushiColors(5))
labelgenome(ax, "chr11", 1650000, 2350000, n=4, scale="Mb", side=1)
fig.savefig("/tmp/dnasei.png", dpi=120, bbox_inches="tight")
```

Result: a DNaseI signal track filled with the 5-color Sushi palette (blue → purple → red → orange → yellow).

## 2. Overlay two tracks

```python
ctcf = sushi.data.Sushi_ChIPSeq_CTCF_bedgraph()
dnase = sushi.data.Sushi_DNaseI_bedgraph()

fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedgraph(ax, ctcf,  "chr11", 1955000, 1960000,
             transparency=0.5, color=SushiColors(2)(2)[0])
plotBedgraph(ax, dnase, "chr11", 1955000, 1960000,
             transparency=0.5, color=SushiColors(2)(2)[1],
             overlay=True, rescaleoverlay=True)
labelgenome(ax, "chr11", 1955000, 1960000, n=3, scale="Kb", side=1)
fig.savefig("/tmp/overlay.png", dpi=120, bbox_inches="tight")
```

`overlay=True` says "don''t create a new empty plot — draw on top of the existing one". `rescaleoverlay=True` rescales the second track to share the y-range.

## 3. Plot ChIP-seq peaks as colored regions

```python
import sushi
from sushi import plotBed, labelgenome, SushiColors, maptocolors

bed = sushi.data.Sushi_ChIPSeq_severalfactors_bed()    # 130 rows, 7 cols
bed["color"] = maptocolors(bed["row"].tolist(), col=SushiColors(6))

fig, ax = plt.subplots(figsize=(10, 5))
fig.subplots_adjust(bottom=0.22)
plotBed(ax, bed, "chr15", 72800000, 73100000,
        type="region",
        color=bed["color"].tolist(),
        rownumber=bed["row"].tolist(), row="given",
        plotbg="#F2F2F2",
        rowlabels=bed["name"].drop_duplicates().tolist()[:11],
        rowlabelcex=0.7)
labelgenome(ax, "chr15", 72800000, 73100000, n=3, scale="Mb", side=1)
fig.savefig("/tmp/several_factors.png", dpi=120, bbox_inches="tight")
```

Three `type` values are supported: `"region"` (filled rectangles), `"circles"` (one circle per row, like WashU browser), `"density"` (vertical gradient heatmap per row).

## 4. Hi-C lower-triangle heatmap

```python
import numpy as np, pandas as pd
import sushi

hic = sushi.data.Sushi_HiC_matrix()       # {matrix: (114,113), positions: (114,)}
pos = hic["positions"].astype(int)
m = hic["matrix"]
# R stores lower-triangular; pad to a square so plotHic is happy
full = np.full((m.shape[0], m.shape[0]), np.nan)
for i in range(m.shape[0]):
    for j in range(i, m.shape[1]):
        full[i, j] = m[i, j]
hic_df = pd.DataFrame(full, index=pos, columns=pos)

fig, ax = plt.subplots(figsize=(10, 6))
fig.subplots_adjust(bottom=0.22, right=0.82)
phic = sushi.plotHic(ax, hic_df, "chr11", 500000, 5050000,
                     max_y=20, zrange=(0, 28), palette=SushiColors(7))
sushi.addlegend(ax, range_val=phic[0], palette=phic[1], title="score",
                side="right", bottominset=0.4, topinset=0, xoffset=-0.035,
                labelside="left", width=0.025, title_offset=0.035)
sushi.labelgenome(ax, "chr11", 500000, 5050000, n=4, scale="Mb",
                  edgeblankfraction=0.20, side=1)
fig.savefig("/tmp/hic.png", dpi=120, bbox_inches="tight")
```

`plotHic` returns `[zrange, palette]` — pass both to `addlegend` to get a matching color bar.

## 5. Plot gene / transcript models

```python
genes = sushi.data.Sushi_genes_bed()         # 5 rows on chr15
trans = sushi.data.Sushi_transcripts_bed()   # 143 rows on chr15

# (a) Arrow style (one big gene)
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
sushi.plotGenes(ax, genes, "chr15", 72998000, 73020000,
                types=["exon"] * len(genes),
                maxrows=1, bheight=0.2, plotgenetype="arrow",
                bentline=False, labeloffset=0.4, fontsize=1.2,
                arrowlength=0.025, labeltext=True)
sushi.labelgenome(ax, "chr15", 72998000, 73020000, n=3, scale="Mb", side=1)
fig.savefig("/tmp/genes.png", dpi=120, bbox_inches="tight")

# (b) Box style + colorby FPKM
fig, ax = plt.subplots(figsize=(10, 5))
fig.subplots_adjust(bottom=0.22)
pg = sushi.plotGenes(ax, trans, "chr15", 72965000, 72990000,
                     types=trans["type"].tolist(),
                     colorby=np.log10(trans["score"].values + 0.001).tolist(),
                     colorbycol=SushiColors(5), colorbyrange=(0, 1.0),
                     labeltext=True, maxrows=50, bheight=0.4,
                     plotgenetype="box")
sushi.labelgenome(ax, "chr15", 72965000, 72990000, n=3, scale="Mb", side=1)
fig.savefig("/tmp/transcripts.png", dpi=120, bbox_inches="tight")
```

## 6. Manhattan plot

```python
gwas = sushi.data.Sushi_GWAS_bed()      # 32760 SNPs across the genome
genome = sushi.data.Sushi_hg18_genome() # 22 chrom sizes

fig, ax = plt.subplots(figsize=(12, 4))
fig.subplots_adjust(bottom=0.22)
sushi.plotManhattan(ax, gwas, pvalues=gwas.iloc[:, 4].values,
                    genome=genome, col=SushiColors(6), space=0.01)
sushi.labelgenome(ax, "chr1", 0, 1, genome=genome, n=4, scale="Mb",
                  edgeblankfraction=0.20, space=0.01)
fig.savefig("/tmp/manhattan.png", dpi=120, bbox_inches="tight")
```

## 7. BEDPE interactions as loops / ribbons

```python
bedpe = sushi.data.Sushi_5C_bedpe()      # 4787 5C interactions

# (a) parabolic loops
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
sushi.plotBedpe(ax, bedpe, "chr11", 1650000, 2350000,
                heights=bedpe["score"].tolist(),
                plottype="loops",
                colorby=bedpe["samplenumber"].tolist(),
                colorbycol=SushiColors(3))
sushi.labelgenome(ax, "chr11", 1650000, 2350000, n=3, scale="Mb", side=1)
fig.savefig("/tmp/loops.png", dpi=120, bbox_inches="tight")

# (b) filled ribbons
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
sushi.plotBedpe(ax, bedpe, "chr11", 1650000, 2350000,
                heights=bedpe["score"].tolist(),
                plottype="ribbons",
                colorby=bedpe["samplenumber"].tolist(),
                colorbycol=SushiColors(3), border="black")
sushi.labelgenome(ax, "chr11", 1650000, 2350000, n=3, scale="Mb", side=1)
fig.savefig("/tmp/ribbons.png", dpi=120, bbox_inches="tight")
```

`plottype="lines"` is a third option that draws each end as a colored rectangle connected by a horizontal segment.

## 8. Zoom-in figure (1 wide + 2 narrow panels)

```python
bg = sushi.data.Sushi_DNaseI_bedgraph()
fig = plt.figure(figsize=(10, 8))
gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])
ax_top = fig.add_subplot(gs[0, :])
ax_b1  = fig.add_subplot(gs[1, 0])
ax_b2  = fig.add_subplot(gs[1, 1])
fig.subplots_adjust(bottom=0.10, hspace=0.5, left=0.10, right=0.95)

chrom, cs, ce = "chr11", 1900000, 2350000
sushi.plotBedgraph(ax_top, bg, chrom, cs, ce, colorbycol=SushiColors(5))
sushi.labelgenome(ax_top, chrom, cs, ce, n=4, scale="Mb", side=1)
sushi.zoomsregion(ax_top, [1955000, 1960000],
                  extend=(0.01, 0.13), wideextend=0.05, offsets=(0, 0.580))
sushi.zoomsregion(ax_top, [2279000, 2284000],
                  extend=(0.01, 0.13), wideextend=0.05, offsets=(0.580, 0))

zr1 = (1955000, 1960000)
sushi.plotBedgraph(ax_b1, bg, chrom, zr1[0], zr1[1], colorbycol=SushiColors(5))
sushi.labelgenome(ax_b1, chrom, zr1[0], zr1[1], n=4, scale="Kb",
                  edgeblankfraction=0.20)
sushi.zoombox(ax_b1, zoomregion=None, lwd=1, col="black")

zr2 = (2279000, 2284000)
sushi.plotBedgraph(ax_b2, bg, chrom, zr2[0], zr2[1], colorbycol=SushiColors(5))
sushi.labelgenome(ax_b2, chrom, zr2[0], zr2[1], n=4, scale="Kb",
                  edgeblankfraction=0.20)
sushi.zoombox(ax_b2, zoomregion=None, lwd=1, col="black")
fig.savefig("/tmp/zoom.png", dpi=120, bbox_inches="tight")
```

`zoomsregion` draws the trapezoid connector on the wide panel; `zoombox` draws the indicator on each zoomed-in panel.

## What''s next?

- See **[API.md](API.md)** for the full parameter reference
- See **[VIGNETTE.md](VIGNETTE.md)** for the side-by-side R-vs-Python port of all 15 vignette examples
- See **[examples/](examples/)** for one self-contained .py per example
- Run `python sushi/vignette.py` to regenerate all 15 figures into `sushi/vignette_output/`
