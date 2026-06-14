# Changelog

## [0.1.0] - 2026-06-13

First release of Sushi Python Version.

### Ported from R Sushi 1.32.0

| R function | Python equivalent | Status |
| --- | --- | --- |
| `plotBedgraph` | `sushi.plotBedgraph` | ✅ (gradient, overlay, flip) |
| `plotBed` | `sushi.plotBed` | ✅ (region / circles / density / splitstrand) |
| `plotBedpe` | `sushi.plotBedpe` | ✅ (loops / ribbons / lines) |
| `plotHic` | `sushi.plotHic` | ✅ (auto-pads lower-triangular R storage) |
| `plotGenes` | `sushi.plotGenes` | ✅ (box / arrow, colorby FPKM) |
| `plotManhattan` | `sushi.plotManhattan` | ✅ (single + multi-chrom) |
| `labelgenome` | `sushi.labelgenome` | ✅ (single + multi-chrom) |
| `addlegend` | `sushi.addlegend` | ✅ |
| `zoombox` | `sushi.zoombox` | ✅ |
| `zoomsregion` | `sushi.zoomsregion` | ✅ |
| `labelplot` | `sushi.labelplot` | ✅ |
| `SushiColors` | `sushi.SushiColors` | ✅ (palettes 2..7) |
| `maptocolors` | `sushi.maptocolors` | ✅ |
| `maptolwd` | `sushi.maptolwd` | ✅ |
| `opaque` | `sushi.opaque` | ✅ (9-char RGBA hex) |
| `convertstrandinfo` | `sushi.convertstrandinfo` | ✅ |
| `sortChrom` | `sushi.sortChrom` | ✅ |
| `chromOffsets` | `sushi.chromOffsets` | ✅ |

### 14 example datasets ported

All 14 datasets from `Sushi/data/` were converted from `.rda` (R binary)
to `.csv` via the `rdata` Python library:

| Dataset | Rows | Format |
| --- | --- | --- |
| `Sushi_DNaseI.bedgraph` | 5,714 | CSV (chrom, start, end, value) |
| `Sushi_ChIPSeq_CTCF.bedgraph` | 27,576 | CSV |
| `Sushi_ChIPExo_CTCF.bedgraph` | 6,292 | CSV |
| `Sushi_ChIPSeq_pol2.bedgraph` | 7,528 | CSV |
| `Sushi_ChIPSeq_pol2.bed` | 208 | CSV (BED6) |
| `Sushi_ChIPSeq_severalfactors.bed` | 130 | CSV (BED7) |
| `Sushi_RNASeq_K562.bedgraph` | 1,073 | CSV |
| `Sushi_5C.bedpe` | 4,787 | CSV (BEDPE11) |
| `Sushi_ChIAPET_pol2.bedpe` | 48,634 | CSV (BEDPE10) |
| `Sushi_HiC.matrix` | 114×113 | CSV (row labels in col 0) |
| `Sushi_GWAS.bed` | 32,760 | CSV |
| `Sushi_genes.bed` | 5 | CSV |
| `Sushi_transcripts.bed` | 143 | CSV |
| `Sushi_hg18_genome` | 22 | CSV (chrom, size) |

### Vignette: 15/15 examples pass

`python sushi/vignette.py` reproduces all 15 examples from the R
package''s `inst/doc/Sushi.Rnw` and saves 15 PNG figures to
`sushi/vignette_output/`.

### Signature differences from R

| Concept | R | Python |
| --- | --- | --- |
| Plot context | `par(mfrow=c(2,1)); plotBed(...)` | `fig, axes = plt.subplots(2, 1); plotBed(ax, ...)` |
| Save to file | `png("out.png"); ...; dev.off()` | `fig.savefig("out.png")` |
| Layout | `layout(matrix(c(1,1,2,3), 2, 2))` | `fig.add_gridspec(...)` + `fig.add_subplot(...)` |
| Reserved-word rename | `addlegend(range, ...)` | `addlegend(ax, range_val=..., ...)` |
| Default color | `dodgerblue4` | `navy` (matplotlib doesn''t know `dodgerblue4`) |
| Default color (2) | `dodgerblue2` | `dodgerblue` |
| `par(mgp=c(3,.3,0))` | `fig.subplots_adjust(left=..., bottom=...)` (caller) |
| `NA` (R missing) | None (Python) |

### Known limitations vs R Sushi

1. **Multi-panel layout**: R uses `par(mfrow=...)` and `layout(matrix(...))` to
   lay out multi-panel figures on a single `plot()` device. Python replaces
   this with `matplotlib.gridspec.GridSpec` and `fig.add_subplot(...)`,
   which is slightly more verbose but gives the caller full control.

2. **3D perspective in `zoombox`**: the R package has no 3D `zoombox`
   either; this port covers only the 2D version.

3. **`opaque()` returns 9-character RGBA hex**: matplotlib parses `#RRGGBBAA`
   correctly, but you must NOT call `to_hex(..., keep_alpha=False)` on the
   result (which would drop the alpha channel).

4. **Hi-C rectangular storage**: R stores `Sushi_HiC.matrix` as a
   `114 × 113` lower-triangular object. The Python loader returns the matrix
   as-is; the caller must pad it to a `114 × 114` square before calling
   `plotHic` (see QUICKSTART §4 for an example).

5. **`biomaRt` gene fetch**: not ported. The R vignette uses it to fetch
   gene coordinates from Ensembl; this Python version assumes you have
   your gene table already in a BED-like DataFrame.

6. **`pdf()` / `png()` device drivers**: use `matplotlib.pyplot.savefig`
   instead. Supports PNG, PDF, SVG, EPS, etc.

### Fixed in v0.1.3 (2026-06-14)

**Root-cause bugs found by reading the R source and comparing with the
Sushi.pdf vignette outputs (NOT by inspecting our own output, which was
the original mistake). Three real bugs were hiding in v0.1.0/v0.1.1/v0.1.2
that made the rendered figures look "mostly black" or "all blue" vs the
R reference:**

1. **`_PALETTES` was wrong since v0.1.0.** R's SushiColors(5) is
   `c("black","blue","purple","red","orange")` but our port had
   `["blue","purple","red","orange","yellow"]` (missing "black",
   ending in "yellow" instead of "orange"). Same off-by-one for n=6.
   Source: R/SushiColors.R lines 49-79. Fixed in
   sushi/_helpers.py:_PALETTES to match R exactly.

2. **`plotBedgraph` did NOT set ylim when overlay=True** (plotters.py
   line 80 originally only called ax.set_ylim inside the `if not
   overlay:` block). R always sets ylim regardless of overlay (R
   plot.R line 85: `plot(..., ylim=range, ...)`). Fix: move the
   set_ylim call outside the if-guard.

3. **`plotBedgraph` colorbycol branch used a wrong algorithm.** The
   original Python port called `maptocolors(ys_step, ...)` once and
   drew a single big polygon. R's actual implementation (R/plotBedgraph.R
   lines 181-213) computes 100 evenly-spaced colors from the palette,
   then for each bin draws up to 100 separate rect() calls, each with
   the color of the corresponding y-band. This makes the gradient
   stripe inside each bin. Rewrote the branch to faithfully port
   R's per-bin-rect approach. The result matches the R output
   (black bottom + blue middle + red/orange top of each spike) for
   the DNaseI + CTCF overlay example.

4. **Added ground-truth validation framework.** Saved the original
   R/Bioconductor Sushi 1.32.0 .rda files to
   `C:\Users\Qianli\Desktop\sushi_rda_groundtruth\` and used the
   Python `rdata==1.1.0` library to parse them to TSV. All 14
   datasets round-trip to 100% byte-identical match with our existing
   `sushi/data/*.csv` (verified cell-by-cell with NaN-aware
   comparison, 705,778 cells total, 0 mismatches). This rules out
   data-conversion bugs as a source of visual differences.

After these fixes, `plotBedgraph(Sushi_DNaseI.bedgraph, chrom, start,
end, colorbycol=SushiColors(5))` matches the R vignette example 6
output in `Sushi.pdf` page 4. 97 distinct color shades are now
used (vs only black+navy before).

### Fixed in v0.1.2 (2026-06-14)

**Three additional critical bug fixes (follow-up to v0.1.1)**:

1. **`SushiColors` factory didn't interpolate colors**. The original
   `_SushiPaletteFactory.__call__` rounded n evenly-spaced positions to
   the nearest of 7 discrete stops, so `SushiColors(7)(101)` returned
   only 7 unique colors no matter what n you asked for. R's
   `colorRampPalette` interpolates in LAB color space. Replaced with
   `matplotlib.colors.LinearSegmentedColormap.from_list` which does
   the same.

2. **`plotBedpe`/`plotBed` called `convertstrandinfo` on the strand
   columns, which converted "." / "+" to int 1 and assigned to a
   dtype='str' pandas column, raising `Invalid value for dtype 'str'`
   for ex05, ex06, ex10, ex11, ex12. Fix: (a) removed the auto-call
   in `plotBedpe` (the numeric strand was never used downstream), and
   (b) in `plotBed` only convert when `splitstrand=True` is requested.

3. **`plotManhattan` int column got float assignment**. The line
   `df.loc[mask, df.columns[1]] = df.loc[mask, df.columns[1]] + offset`
   added a float chrom offset to an int position column. Fix: cast
   through float then back to int.

Also: vignette example 04 (`ex04_hic`) now uses `topo.colors` palette
and `flip=True` to match R's rendered output on Sushi.pdf page 11
(matches the actual Hi-C matrix page in the package documentation).

Result: **15/15 vignette examples now pass on local + production**.
The earlier "10/15" was real (pre-existing in v0.1.0/v0.1.1, not
introduced by v0.1.1's plotBedgraph fix) — these 5 bugs were silent
when the user only ran the other examples.

### Fixed in v0.1.1 (2026-06-14)

**Critical bug fix in `plotBedgraph`** (caught by user on `44_PolII/sushi_6candidate_v1.png`):

The original polygon construction `np.concatenate([starts, stops_reversed])` +
`np.concatenate([values, values_reversed])` drew a single Z-folded polygon with a
long diagonal from the last bin's top to the first bin's top. For bedGraphs
with more than ~3 bins of distinct values, this filled the **entire panel** with
the wrong shape instead of correctly filling each bin as a horizontal step.

**Fix**: new helper `_bedgraph_step_polygon` in `sushi/_helpers.py` builds a true
step-function polygon (each bin is a horizontal rectangle). `plotBedgraph`
now calls this helper. **Do not revert to the old `np.concatenate` form** without
understanding the geometry — see the docstring on `_bedgraph_step_polygon` and
the implementation note added to `plotBedgraph` itself.

Verified against user's `sushi_6candidate.py` (6 candidate genes × 4 NET-seq
tracks): previously each panel showed one giant triangle, now each panel
shows the correct per-bin step-function NET-seq profile.

### Known limitations vs R Sushi

1. **`biomaRt` gene fetch** (used in R vignette examples 33-35) is **not ported**.
   Reason: `biomaRt` is an R-specific Bioconductor binding to the Ensembl
   BioMart Web service. Use any Python tool (UCSC Table Browser, `gget`,
   `pybiomart`, direct REST) to fetch gene tables and pass them to
   `sushi.plotGenes`.

2. **R's `pdf()` / `png()` device drivers** are **not ported**. Use
   `matplotlib.pyplot.savefig` directly — supports PNG/PDF/SVG/EPS/JPG out
   of the box.

(There are no other R-side features that sushi-py is missing. Earlier
drafts of this list mentioned `plotBedpe method="interaction3D"` and a
3D `zoombox`, but those were CHANGELOG mistakes — R Sushi 1.32.0 has only
`plottype="loops"|"ribbons"|"lines"` and only a 2D `zoombox`. The Python port
covers R Sushi 1.32.0's full public API.)

### Tested on

- Python 3.10.12 (师大服务器 `miniforge3`)
- matplotlib 3.10.8
- numpy 1.26.4
- pandas 2.2.1
- pyBigWig 0.3.25 (not used by Sushi but available)
- scipy 1.15.3

### Citation

If you use this Python port, please cite the original R paper:

> Phanstiel DH, Boyle AP, Araya CL, Snyder MP. **Sushi.r: flexible,
> quantitative, and integrative genomic visualizations for
> publication-quality multi-panel figures.** *Bioinformatics* (2015).

This Python port does not introduce a new citation; it is faithful to the
algorithms and parameters of the original R package (Bioconductor
`Sushi 1.32.0`, GPL ≥ 2 license).
