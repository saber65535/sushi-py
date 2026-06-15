# R ↔ Python Vignette Comparison

This document pairs each example from the R package''s `inst/doc/Sushi.Rnw`
vignette with its Python equivalent. The R code is the original; the
Python code is what `sushi/vignette.py` runs. All 15 outputs are saved to
`sushi/vignette_output/`.

The Python side differs from R in three ways:

1. The caller manages the figure (`plt.subplots`, `fig.subplots_adjust`),
   not `par(mfrow=...)`.
2. Every plot function takes an `Axes` as its first argument.
3. `addlegend(range, ...)` is renamed to `addlegend(ax, range_val=..., ...)`
   to avoid shadowing Python''s built-in `range`.

---

## Example 1 — `plotBedgraph` with gradient

R:
```r
data(Sushi_DNaseI.bedgraph)
chrom       <- "chr11"
chromstart  <- 1650000
chromend    <- 2350000
plotBedgraph(Sushi_DNaseI.bedgraph, chrom, chromstart, chromend,
             colorbycol = SushiColors(5))
labelgenome(chrom, chromstart, chromend, n=4, scale="Mb", side=1)
```

Python:
```python
bg = sushi.data.Sushi_DNaseI_bedgraph()
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedgraph(ax, bg, "chr11", 1650000, 2350000, colorbycol=SushiColors(5))
labelgenome(ax, "chr11", 1650000, 2350000, n=4, scale="Mb", side=1)
fig.savefig("01_dnasei_bedgraph.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/01_dnasei_bedgraph.png`

---

## Example 2 — Overlay two bedGraph tracks

R:
```r
plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph, chrom, chromstart, chromend,
             transparency=.50, flip=FALSE, color="blue")
plotBedgraph(Sushi_DNaseI.bedgraph, chrom, chromstart, chromend,
             transparency=.50, flip=FALSE, color="#E5001B",
             overlay=TRUE, rescaleoverlay=TRUE)
labelgenome(chrom, chromstart, chromend, n=3, scale="Kb", side=1)
legend(..., fill=opaque(c("blue","#E5001B")), ...)
```

Python:
```python
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedgraph(ax, ctcf,  "chr11", 1955000, 1960000, transparency=0.5,
             color=SushiColors(2)(2)[0])
plotBedgraph(ax, dnase, "chr11", 1955000, 1960000, transparency=0.5,
             color=SushiColors(2)(2)[1], overlay=True, rescaleoverlay=True)
labelgenome(ax, "chr11", 1955000, 1960000, n=3, scale="Kb", side=1)
# Add legend manually if needed (matplotlib Axes.legend)
```

Output: `sushi/vignette_output/02_overlay_dnasei_ctcf.png`

---

## Example 3 — Flip one of two tracks

R:
```r
par(mfrow=c(2,1), mar=c(1,4,1,1))
plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph, chrom, chromstart, chromend,
             transparency=.50, color="blue")
# y-axis ticks added manually in R
plotBedgraph(Sushi_DNaseI.bedgraph, chrom, chromstart, chromend,
             transparency=.50, flip=TRUE, color="#E5001B")
labelgenome(chrom, chromstart, chromend, side=3, n=3, scale="Kb")
```

Python:
```python
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedgraph(ax, ctcf,  "chr11", 1955000, 1960000, transparency=0.5,
             color=SushiColors(2)(2)[0])
plotBedgraph(ax, dnase, "chr11", 1955000, 1960000, transparency=0.5,
             flip=True, color=SushiColors(2)(2)[1])
labelgenome(ax, "chr11", 1955000, 1960000, side=3, n=3, scale="Kb")
fig.savefig("03_flip_dnasei_ctcf.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/03_flip_dnasei_ctcf.png`

> Multi-panel `par(mfrow=c(2,1))` is replaced by `plt.subplots(2, 1)` in the
> Python vignette. For more elaborate layouts see Example 15 (zoom).

---

## Example 4 — `plotHic` with `addlegend`

R:
```r
data(Sushi_HiC.matrix)
phic = plotHic(Sushi_HiC.matrix, chrom, chromstart, chromend,
               max_y = 20, zrange=c(0, 28), palette = SushiColors(7))
addlegend(phic[[1]], palette=phic[[2]], title="score", side="right",
          bottominset=0.4, topinset=0, xoffset=-.035, labelside="left",
          width=0.025, title.offset=0.035)
labelgenome(chrom, chromstart, chromend, n=4, scale="Mb",
            edgeblankfraction=0.20)
```

Python:
```python
import numpy as np, pandas as pd
hic = sushi.data.Sushi_HiC_matrix()
pos = hic["positions"].astype(int)
m = hic["matrix"]
# R stores lower-triangular; pad to a square for pandas + plotHic
full = np.full((m.shape[0], m.shape[0]), np.nan)
for i in range(m.shape[0]):
    for j in range(i, m.shape[1]):
        full[i, j] = m[i, j]
hic_df = pd.DataFrame(full, index=pos, columns=pos)

fig, ax = plt.subplots(figsize=(10, 6))
fig.subplots_adjust(bottom=0.22, right=0.82)
phic = plotHic(ax, hic_df, "chr11", 500000, 5050000, max_y=20,
               zrange=(0, 28), palette=SushiColors(7))
addlegend(ax, range_val=phic[0], palette=phic[1], title="score",
          side="right", bottominset=0.4, topinset=0, xoffset=-0.035,
          labelside="left", width=0.025, title_offset=0.035)
labelgenome(ax, "chr11", 500000, 5050000, n=4, scale="Mb",
            edgeblankfraction=0.20, side=1)
fig.savefig("04_hic_chr11.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/04_hic_chr11.png`

---

## Example 5 — `plotBedpe` loops

R:
```r
data(Sushi_5C.bedpe)
chrom            <- "chr11"
chromstart       <- 1650000
chromend         <- 2350000
pbpe = plotBedpe(Sushi_5C.bedpe, chrom, chromstart, chromend,
                 heights = Sushi_5C.bedpe$score, plottype = "loops",
                 colorby = Sushi_5C.bedpe$samplenumber,
                 colorbycol = SushiColors(3))
labelgenome(chrom, chromstart, chromend, n=3, scale="Mb", side=1)
legend("topright", legend=c("K562","HeLa","GM12878"),
       col=SushiColors(3)(3), pch=19, bty="n", text.font=2)
```

Python:
```python
bedpe = sushi.data.Sushi_5C_bedpe()
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedpe(ax, bedpe, "chr11", 1650000, 2350000,
          heights=bedpe["score"].tolist(), plottype="loops",
          colorby=bedpe["samplenumber"].tolist(),
          colorbycol=SushiColors(3))
labelgenome(ax, "chr11", 1650000, 2350000, n=3, scale="Mb", side=1)
# Manual legend (matplotlib)
from matplotlib.lines import Line2D
handles = [Line2D([0],[0], marker="o", linestyle="",
                   markerfacecolor=c, markeredgecolor=c, label=l)
           for c, l in zip(SushiColors(3)(3), ["K562","HeLa","GM12878"])]
ax.legend(handles=handles, loc="upper right", frameon=False)
fig.savefig("05_bedpe_5C_loops.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/05_bedpe_5C_loops.png`

---

## Example 6 — `plotBedpe` ribbons

R:
```r
plotBedpe(Sushi_5C.bedpe, chrom, chromstart, chromend,
          heights=Sushi_5C.bedpe$score, plottype="ribbons",
          colorby=Sushi_5C.bedpe$samplenumber,
          colorbycol=SushiColors(3), border="black")
```

Python:
```python
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedpe(ax, bedpe, "chr11", 1650000, 2350000,
          heights=bedpe["score"].tolist(), plottype="ribbons",
          colorby=bedpe["samplenumber"].tolist(),
          colorbycol=SushiColors(3), border="black")
labelgenome(ax, "chr11", 1650000, 2350000, n=3, scale="Mb", side=1)
fig.savefig("06_bedpe_5C_ribbons.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/06_bedpe_5C_ribbons.png`

---

## Example 7 — `plotBed` region mode

R:
```r
chrom       <- "chr11"
chromstart  <- 2281200
chromend    <- 2282200
plotBed(beddata=Sushi_ChIPSeq_pol2.bed, chrom, chromstart, chromend,
        colorby=Sushi_ChIPSeq_pol2.bed$strand,
        colorbycol=SushiColors(2), row="auto", wiggle=0.001)
```

Python:
```python
bed = sushi.data.Sushi_ChIPSeq_pol2_bed()
bed = bed[bed["chrom"] == "chr11"].reset_index(drop=True)   # pre-filter
cb = sushi.convertstrandinfo(bed["strand"].tolist())  # strand is numeric (-1/1) in the data
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22, left=0.12)
plotBed(ax, bed, "chr11", 2281200, 2282200, colorby=cb,
        colorbycol=SushiColors(2), row="auto", wiggle=0.001)
labelgenome(ax, "chr11", 2281200, 2282200, n=2, scale="Kb", side=1)
fig.savefig("07_bed_pol2_region.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/07_bed_pol2_region.png`

> The Python vignette pre-filters `bed` to the chromosome region before
> passing it to `plotBed` because the auto-filter inside `plotBed` does not
> also trim the `colorby` / `rownumber` arrays to the filtered length.
> Pre-filtering keeps the two array lengths in sync. This is a known
> ergonomic quirk of the Python port; future versions may add automatic
> pre-filtering.

---

## Example 8 — `plotBed` circles

Same pattern as Example 7, with `type="circles"`. See `sushi/vignette.py`.

Output: `sushi/vignette_output/08_bed_pol2_circles.png`

---

## Example 9 — `plotBed` splitstrand

Same as Example 7 with `splitstrand=True`. See `sushi/vignette.py`.

Output: `sushi/vignette_output/09_bed_pol2_splitstrand.png`

---

## Example 10 — `plotBed` several factors, circles

R:
```r
chrom            = "chr15"
chromstart       = 72800000
chromend         = 73100000
Sushi_ChIPSeq_severalfactors.bed$color =
   maptocolors(Sushi_ChIPSeq_severalfactors.bed$row, col=SushiColors(6))
plotBed(beddata=Sushi_ChIPSeq_severalfactors.bed, chrom, chromstart, chromend,
        rownumber=Sushi_ChIPSeq_severalfactors.bed$row, type="circles",
        color=Sushi_ChIPSeq_severalfactors.bed$color, row="given",
        plotbg="grey95",
        rowlabels=unique(Sushi_ChIPSeq_severalfactors.bed$name),
        rowlabelcol=unique(Sushi_ChIPSeq_severalfactors.bed$color),
        rowlabelcex=0.75)
```

Python:
```python
bed_full = sushi.data.Sushi_ChIPSeq_severalfactors_bed()
bed = bed_full[(bed_full["chrom"] == "chr15") &
               (((bed_full["start"] > 72800000) & (bed_full["start"] < 73100000)) |
                ((bed_full["end"] > 72800000) & (bed_full["end"] < 73100000)))].reset_index(drop=True)
bed["color"] = sushi.maptocolors(bed["row"].tolist(), col=SushiColors(6))
fig, ax = plt.subplots(figsize=(10, 5))
fig.subplots_adjust(bottom=0.22)
plotBed(ax, bed, "chr15", 72800000, 73100000, type="circles",
        color=bed["color"].tolist(),
        rownumber=bed["row"].tolist(), row="given",
        plotbg="#F2F2F2",
        rowlabels=bed_full["name"].drop_duplicates().tolist()[:11],
        rowlabelcex=0.7)
labelgenome(ax, "chr15", 72800000, 73100000, n=3, scale="Mb", side=1)
fig.savefig("10_severalfactors_circles.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/10_severalfactors_circles.png`

---

## Example 11 — `plotBed` several factors, region

Same as Example 10 with `type="region"`. See `sushi/vignette.py`.

Output: `sushi/vignette_output/11_severalfactors_region.png`

---

## Example 12 — `plotManhattan` (multi-chrom)

R:
```r
data(Sushi_GWAS.bed); data(Sushi_hg18_genome)
plotManhattan(bedfile=Sushi_GWAS.bed, pvalues=Sushi_GWAS.bed[,5],
              genome=Sushi_hg18_genome, col=SushiColors(6),
              cex=0.75, space=0.05)
labelgenome(genome=Sushi_hg18_genome, n=4, scale="Mb",
            edgeblankfraction=0.20, cex.axis=.5, space=0.05)
axis(side=2, las=2, tcl=.2)
mtext("log10(P)", side=2, line=1.75, cex=1, font=2)
```

Python:
```python
gwas = sushi.data.Sushi_GWAS_bed()
genome = sushi.data.Sushi_hg18_genome()
fig, ax = plt.subplots(figsize=(12, 4))
fig.subplots_adjust(bottom=0.22)
plotManhattan(ax, gwas, pvalues=gwas.iloc[:, 4].values,
              genome=genome, col=SushiColors(6), space=0.01)
# In multi-chrom mode, labelgenome needs a single chrom/region
# to draw the scale label; we use chr1 as a placeholder.
labelgenome(ax, "chr1", 0, 1, genome=genome, n=4, scale="Mb",
            edgeblankfraction=0.20, space=0.01)
fig.savefig("12_manhattan_gwas.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/12_manhattan_gwas.png`

---

## Example 13 — `plotGenes` (single gene, arrow)

R:
```r
pg = plotGenes(Sushi_genes.bed, chrom, chromstart, chromend,
               types=Sushi_genes.bed$type, maxrows=1, bheight=0.2,
               plotgenetype="arrow", bentline=FALSE,
               labeloffset=.4, fontsize=1.2, arrowlength=0.025,
               labeltext=TRUE)
```

Python:
```python
genes = sushi.data.Sushi_genes_bed()
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotGenes(ax, genes, "chr15", 72998000, 73020000,
          types=["exon"] * len(genes),
          maxrows=1, bheight=0.2, plotgenetype="arrow",
          bentline=False, labeloffset=0.4, fontsize=1.2,
          arrowlength=0.025, labeltext=True)
labelgenome(ax, "chr15", 72998000, 73020000, n=3, scale="Mb", side=1)
fig.savefig("13_genes_chr15.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/13_genes_chr15.png`

---

## Example 14 — `plotGenes` (transcripts, box, colorby FPKM)

R:
```r
pg = plotGenes(Sushi_transcripts.bed, chrom, chromstart, chromend,
               types=Sushi_transcripts.bed$type,
               colorby=log10(Sushi_transcripts.bed$score + 0.001),
               colorbycol=SushiColors(5), colorbyrange=c(0, 1.0),
               labeltext=TRUE, maxrows=50, bheight=0.4, plotgenetype="box")
addlegend(pg[[1]], palette=pg[[2]], title="log10(FPKM)", side="right",
          bottominset=0.4, topinset=0, xoffset=-.035, labelside="left",
          width=0.025, title.offset=0.055)
```

Python:
```python
import numpy as np
trans = sushi.data.Sushi_transcripts_bed()
fig, ax = plt.subplots(figsize=(10, 5))
fig.subplots_adjust(bottom=0.22, right=0.82)
pg = plotGenes(ax, trans, "chr15", 72965000, 72990000,
               types=trans["type"].tolist(),
               colorby=np.log10(trans["score"].values + 0.001).tolist(),
               colorbycol=SushiColors(5), colorbyrange=(0, 1.0),
               labeltext=True, maxrows=50, bheight=0.4,
               plotgenetype="box")
labelgenome(ax, "chr15", 72965000, 72990000, n=3, scale="Mb", side=1)
addlegend(ax, range_val=pg[0], palette=pg[1], title="log10(FPKM)",
          side="right", bottominset=0.4, topinset=0, xoffset=-0.035,
          labelside="left", width=0.025, title_offset=0.055)
fig.savefig("14_transcripts_colorby.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/14_transcripts_colorby.png`

---

## Example 15 — Zoom-in (1 wide + 2 narrow)

R:
```r
layout(matrix(c(1,1, 2,3), 2, 2, byrow=TRUE))
par(mar=c(3,4,2,1))
chrom            = "chr11"; chromstart=1900000; chromend=2350000
plotBedgraph(Sushi_DNaseI.bedgraph, chrom, chromstart, chromend,
             colorbycol=SushiColors(5))
labelgenome(...)
zoomsregion(zoomregion1, extend=c(0.01,0.13), wideextend=0.05,
            offsets=c(0, 0.580))
zoomsregion(zoomregion2, extend=c(0.01,0.13), wideextend=0.05,
            offsets=c(0.580, 0))
plotBedgraph(Sushi_DNaseI.bedgraph, chrom, zoomregion1[1], zoomregion1[2],
             colorbycol=SushiColors(5))
zoombox(...)
```

Python:
```python
fig = plt.figure(figsize=(10, 8))
gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])
ax_top = fig.add_subplot(gs[0, :])
ax_b1  = fig.add_subplot(gs[1, 0])
ax_b2  = fig.add_subplot(gs[1, 1])
fig.subplots_adjust(bottom=0.10, hspace=0.5, left=0.10, right=0.95)

chrom, cs, ce = "chr11", 1900000, 2350000
plotBedgraph(ax_top, bg, chrom, cs, ce, colorbycol=SushiColors(5))
labelgenome(ax_top, chrom, cs, ce, n=4, scale="Mb", side=1)
zoomsregion(ax_top, [1955000, 1960000],
            extend=(0.01, 0.13), wideextend=0.05, offsets=(0, 0.580))
zoomsregion(ax_top, [2279000, 2284000],
            extend=(0.01, 0.13), wideextend=0.05, offsets=(0.580, 0))

zr1 = (1955000, 1960000)
plotBedgraph(ax_b1, bg, chrom, zr1[0], zr1[1], colorbycol=SushiColors(5))
labelgenome(ax_b1, chrom, zr1[0], zr1[1], n=4, scale="Kb",
            edgeblankfraction=0.20)
zoombox(ax_b1, zoomregion=None, lwd=1, col="black")

zr2 = (2279000, 2284000)
plotBedgraph(ax_b2, bg, chrom, zr2[0], zr2[1], colorbycol=SushiColors(5))
labelgenome(ax_b2, chrom, zr2[0], zr2[1], n=4, scale="Kb",
            edgeblankfraction=0.20)
zoombox(ax_b2, zoomregion=None, lwd=1, col="black")
fig.savefig("15_zoom_in_out.png", dpi=100, bbox_inches="tight")
```

Output: `sushi/vignette_output/15_zoom_in_out.png`

---

## Run all 15 at once

```bash
cd /datapool/life-zhanghk/qianli/39-Sushi_Python
python sushi/vignette.py
ls sushi/vignette_output/
# → 15 PNG files, named 01_*.png .. 15_*.png
```
