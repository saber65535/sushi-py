# Sushi Python — API Reference

This document lists every public function with its full parameter reference.
For runnable examples see [QUICKSTART.md](QUICKSTART.md) or [examples/](examples/).

## Plotting helpers

All plotters take a `matplotlib.axes.Axes` as the first argument and modify
it in place. The caller is responsible for `plt.subplots(...)`,
`fig.subplots_adjust(...)` and `fig.savefig(...)`.

---

### `sushi.plotBedgraph(ax, signal, chrom, chromstart, chromend, ...)`

Plot a bedGraph signal track.

**R signature**: `plotBedgraph(signal, chrom, chromstart, chromend, ...)`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `ax` | `Axes` | (required) | Target matplotlib axes |
| `signal` | `pd.DataFrame` or `np.ndarray` | (required) | bedGraph = (chrom, start, stop, value) |
| `chrom` | `str` | (required) | Chromosome name to plot (e.g. "chr11") |
| `chromstart`, `chromend` | `int` | (required) | 1-based inclusive coords |
| `range_` | `(min, max)` | `None` | Y-axis range. If None, auto-scales to `ymax * max(value)` |
| `color` | `str` | `SushiColors(2)(2)[0]` | Fill color (hex or named) |
| `lwd` | `float` | `1.0` | Line width for outline |
| `linecolor` | `str` | same as `color` | Outline color; `NA` in R = no outline |
| `addscale` | `bool` | `False` | Add a small scale label upper-right |
| `overlay` | `bool` | `False` | If True, draw on top of an existing axes (no empty plot creation) |
| `rescaleoverlay` | `bool` | `False` | Rescale overlay values to share y-range |
| `transparency` | `float` | `1.0` | Alpha for the fill (0=transparent, 1=opaque) |
| `flip` | `bool` | `False` | Negate values and draw inverted |
| `ymax` | `float` | `1.04` | Fraction of max value to set as upper bound |
| `colorbycol` | `callable` | `None` | If set, draw a vertical gradient under the curve (e.g. `SushiColors(5)`) |

**Returns**: the `ax` (for chaining).

---

### `sushi.plotBed(ax, beddata, chrom, chromstart, chromend, type='region', ...)`

Plot BED-format genomic elements as **region**, **circles**, or **density** heatmap.

**R signature**: `plotBed(beddata, chrom, chromstart, chromend, type, ...)`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `beddata` | `pd.DataFrame` or `np.ndarray` | (required) | BED columns (chrom, start, stop, [name, score, strand], ...) |
| `type` | `str` | `"region"` | One of `"region"`, `"circles"`, `"density"` |
| `colorby` | list of float | `None` | Numeric vector to scale colors by |
| `colorbycol` | `callable` | `None` | Palette function (required when `colorby` is set) |
| `colorbyrange` | `(min, max)` | `None` | Value range for color mapping |
| `rownumber` | list of int | `None` | Per-row assignment. If `None`, defaults to `1` for all |
| `row` | `str` | `"auto"` | `"auto"` packs to fit; `"given"`/`"supplied"` uses `rownumber` as-is |
| `height` | `float` | `0.4` | Row height (for region/circles) |
| `plotbg` | `str` | `"white"` | Background fill of the panel |
| `wiggle` | `float` | `0.02` | Fraction of plot width left blank on either side of each element |
| `splitstrand` | `bool` | `False` | Reverse-strand elements drawn below the x-axis (needs col 6 = strand) |
| `numbins` | `int` | `200` | Bins for density mode |
| `binsmoothing` | `int` | `10` | Bin-smoothing window |
| `palettes` | list of `callable` | `SushiColors(7)` | One palette per row (density mode) |
| `rowlabels` | list of str | `None` | Y-axis labels per row (drawn outside the panel) |
| `rowlabelcol` | `str` or list | `"dodgerblue"` | Color for row labels |
| `rowlabelfont` | `int` | `2` (bold) | matplotlib font weight |
| `rowlabelcex` | `float` | `1.0` | Row label size |
| `maxrows` | `int` | `1000000` | Maximum rows to assign in auto-pack |
| `color` | `str` | `"navy"` | Default color when `colorby` is `None` |

**Returns**: the `ax`.

---

### `sushi.plotBedpe(ax, bedpedata, chrom, chromstart, chromend, heights, ...)`

Plot BEDPE paired-end data as **loops**, **ribbons**, or **lines**.

**R signature**: `plotBedpe(bedpedata, chrom, chromstart, chromend, heights, ...)`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `bedpedata` | `pd.DataFrame` or `np.ndarray` | (required) | BEDPE = (chrom1, start1, stop1, chrom2, start2, stop2, ...) |
| `heights` | list of float | (required) | Per-row apex heights for loops/ribbons |
| `color` | `str` or list | `"black"` | Single color or per-row colors |
| `colorby` | list of float | `None` | Numeric vector to scale colors by |
| `colorbycol` | `callable` | `None` | Palette for colorby |
| `colorbyrange` | `(min, max)` | `None` | Range for colorby |
| `border` | `str` | `None` | Border color for ribbons (default: same as fill) |
| `lwdby` | list of float | `None` | Per-row line widths |
| `lwdrange` | `(min, max)` | `(1, 5)` | Range mapped from `lwdby` |
| `offset` | `float` | `0` | Vertical offset from x-axis |
| `flip` | `bool` | `False` | Negate heights (draw below axis) |
| `lwd` | `float` | `1.0` | Default line width |
| `plottype` | `str` | `"loops"` | One of `"loops"`, `"ribbons"`, `"lines"` |
| `maxrows` | `int` | `10000` | Max rows (only for `lines` plottype) |
| `height` | `float` | `0.3` | Box height for `lines` plottype |
| `ymax` | `float` | `1.04` | Fraction of max height for y-axis ceiling |

**Returns**: the `ax`.

---

### `sushi.plotHic(ax, hicdata, chrom, chromstart, chromend, max_y=30, zrange=None, palette=None, flip=False)`

Plot a Hi-C interaction matrix as a lower-triangle heatmap.

**R signature**: `plotHic(hicdata, chrom, chromstart, chromend, max_y, zrange, palette, flip)`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `hicdata` | `pd.DataFrame` | (required) | Square matrix. Row and column labels should be genomic positions. If you have R''s lower-triangular storage (shape `N, N-1`), build a square matrix yourself (see QUICKSTART §4). |
| `max_y` | `float` | `30.0` | Maximum bin distance on the y-axis |
| `zrange` | `(min, max)` | `None` | Score range. `None` = auto. Values outside are clipped. |
| `palette` | `callable` | `SushiColors(7)` | Color palette function |
| `flip` | `bool` | `False` | Draw inverted (max_y down to 0) |

**Returns**: `[zrange, palette]` — pass to `addlegend` for a matching color bar.

---

### `sushi.plotGenes(ax, geneinfo, chrom, chromstart, chromend, types, ...)`

Plot a gene / transcript model.

**R signature**: `plotGenes(geneinfo, chrom, chromstart, chromend, types, ...)`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `geneinfo` | `pd.DataFrame` | (required) | (chrom, start, stop, name, score, strand, ...) |
| `types` | list of str | (required) | Per-row `exon` / `utr` / `cds` classification. Determines box height: UTR = `bheight/2`, exon/cds = `bheight`. |
| `colorby` | list of float | `None` | Numeric vector to scale colors by |
| `colorbycol` | `callable` | `None` | Palette |
| `colorbyrange` | `(min, max)` | `None` | Range |
| `labeltext` | `bool` | `True` | Show gene labels above the box |
| `labeloffset` | `float` | `0.4` | Y-offset of labels |
| `fontsize` | `float` | `0.7` | Label font size |
| `fonttype` | `int` | `2` (bold) | matplotlib font weight |
| `labelat` | `str` | `"middle"` | `"start"`, `"end"`, or `"middle"` |
| `arrowlength` | `float` | `0.025` | Arrow length (relative to bp range) |
| `bheight` | `float` | `0.2` | Box height |
| `lheight` | `float` | `0.05` | Bent-line peak height |
| `bentline` | `bool` | `True` | Draw bent lines between exons |
| `plotgenetype` | `str` | `"box"` | `"box"` or `"arrow"` |
| `packrow` | `bool` | `True` | Auto-pack to fit non-overlapping rows |
| `maxrows` | `int` | `10000` | Max rows |
| `color` | `str` | `"navy"` | Default color |
| `wigglefactor` | `float` | `0.02` | Fraction of plot width for wiggle |

**Returns**: `[colorbyrange, colorbycol]` (for `addlegend`).

---

### `sushi.plotManhattan(ax, bedfile, pvalues=None, chrom=None, chromstart=None, chromend=None, genome=None, col=SushiColors(5), space=0.01, ymax=1.04)`

Plot a Manhattan-style `-log10(pvalue)` scatter.

**R signature**: `plotManhattan(bedfile, chrom, chromstart, chromend, pvalues, genome, col, space, ymax)`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `bedfile` | `pd.DataFrame` | (required) | (chrom, start, name, score, pvalue, strand, ...) or `pvalues` supplied separately |
| `pvalues` | list of float | `None` | p-values to plot as `-log10(p)`. If `None`, looks for `pvalue` column, then col 4 |
| `chrom`, `chromstart`, `chromend` | `str, int, int` | `None` | Required for single-chrom mode (when `genome=None`) |
| `genome` | `pd.DataFrame` | `None` | (chrom, size) DataFrame. If provided, plot multi-chromosome with cumulative offsets |
| `col` | `callable` or list | `SushiColors(5)` | Palette or list of colors |
| `space` | `float` | `0.01` | Fraction of total length as gap between chromosomes |
| `ymax` | `float` | `1.04` | Fraction of max y to set as y-axis ceiling |

**Returns**: `None` (draws in place).

---

## Axes / overlay helpers

### `sushi.labelgenome(ax, chrom, chromstart, chromend, ...)`

Draw genomic coordinate axis ticks on an axes.

**R signature**: `labelgenome(chrom, chromstart, chromend, ...)`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `chrom` | `str` | (required for single-chrom) | Chromosome name |
| `chromstart`, `chromend` | `int` | (required for single-chrom) | Region bounds |
| `genome` | `pd.DataFrame` | `None` | If provided, multi-chrom mode (chrom names drawn as tick labels at chrom centers) |
| `space` | `float` | `0.01` | Gap between chroms (multi-chrom mode) |
| `scale` | `str` | `"bp"` | `"bp"`, `"Kb"`, or `"Mb"` |
| `side` | `int` | `1` | matplotlib axes side (1=bottom, 3=top) |
| `scipen` | `int` | `20` | Higher values decrease the chance of scientific notation (kept for R compat; matplotlib uses a different formatter) |
| `n` | `int` | `5` | Desired number of tick labels |
| `chromfont`, `chromcex`, `chromline`, `chromadjust` | — | `2, 1.0, 0.5, 0.015` | Chromosome-name label styling |
| `scalefont`, `scalecex`, `scaleline`, `scaleadjust` | — | `2, 1.0, 0.5, 0.985` | Scale label styling |
| `line` | `float` | `0.18` | Vertical offset of tick labels |
| `edgeblankfraction` | `float` | `0.10` | Fraction of x-range left blank on each side for chrom/scale labels |

**Caller must do `fig.subplots_adjust(bottom=0.18)` (or similar) so the
labels are not clipped by the savefig boundary.**

---

### `sushi.addlegend(ax, range_val=None, title='', palette=None, ...)`

Draw a color-scaled legend on the right/left of an axes.

**R signature**: `addlegend(range, title, palette, ...)`

> Note: Python keyword `range` shadows the built-in; renamed to `range_val`.

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `range_val` | `(min, max)` | `None` (required via kwarg) | Color mapping range |
| `title` | `str` | `""` | Legend title (vertical text) |
| `palette` | `callable` | `SushiColors(7)` | Palette function |
| `labels_digits` | `int` | `1` | Decimal digits on tick labels |
| `side` | `str` | `"right"` | `"right"` or `"left"` |
| `labelside` | `str` | `"left"` | `"left"` or `"right"` |
| `xoffset` | `float` | `0.1` | Fraction of plot width to offset legend from axes |
| `width` | `float` | `0.05` | Fraction of plot width |
| `bottominset`, `topinset` | `float` | `0.025, 0.025` | Inset from bottom/top of axes |
| `tick_num` | `int` | `5` | Number of tick marks |
| `tick_length` | `float` | `0.01` | Tick mark length (fraction of plot width) |
| `txt_font`, `txt_cex` | `int, float` | `1, 0.75` | Tick label font / size |
| `title_offset`, `title_font`, `title_cex` | `float, int, float` | `0.05, 2, 1.0` | Title styling |

---

### `sushi.zoombox(ax, zoomregion=None, ...)`

Draw a zoom indicator above an axes. Use on the **second** plot of a zoom-in.

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `zoomregion` | `(start, stop)` | `None` | If provided and `passthrough=False`, the top of the box has a gap |
| `lty`, `lwd`, `col` | matplotlib line kwargs | `"--"`, `1.0`, `"black"` | Style |
| `topextend` | `float` | `2.0` | Fraction of y-range to extend lines above the panel top |
| `passthrough` | `bool` | `False` | If True, extend lines below the panel bottom too |

---

### `sushi.zoomsregion(ax, region, chrom=None, genome=None, ...)`

Draw a trapezoid connector from a parent panel to a zoomed panel. Use on the **first** plot of a zoom-in.

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `region` | `(start, stop)` | (required) | Region the zoom covers |
| `chrom` | `str` | `None` | Required when `genome` is provided |
| `genome` | `pd.DataFrame` | `None` | Multi-chrom mode requires this |
| `padding` | `float` | `0.005` | Min size of the zoom region (fraction of panel width); expanded symmetrically if too narrow |
| `col` | `str` | `"white"` | Fill color (use `"white"` for invisible interior) |
| `zoomborder` | `str` | `"black"` | Border color |
| `lty`, `lwd` | matplotlib line kwargs | `"--"`, `1.0` | Border style |
| `extend` | `float` or `(float, float)` | `0.0` | Vertical extension above/below panel (fraction of y-range) |
| `wideextend` | `float` | `0.1` | How far below the panel the wide part of the trapezoid starts |
| `offsets` | `(float, float)` | `(0, 0)` | Offsets to left/right sides of the wide portion |
| `highlight` | `bool` | `False` | If True, draw a highlight box instead of a zoom connector |

---

### `sushi.labelplot(ax, letter=None, title=None, ...)`

Add a paper letter + title overlay to the top of an axes.

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `letter` | `str` | `None` | e.g. `"A"`, `"A)"`, `"1"` |
| `title` | `str` | `None` | Title text |
| `letteradj`, `titleadj` | `float` | `-0.05, 0.0` | matplotlib `axes fraction` coords |
| `letterfont`, `titlefont` | `int` | `2, 2` (bold) | matplotlib font weight |
| `lettercex`, `titlecex` | `float` | `1.2, 1.0` | Font size |
| `letterline`, `titleline` | `float` | `0.5, 0.5` | Y offset (axes fraction) |
| `lettercol`, `titlecol` | `str` | `"black", "black"` | Color |

---

## Helper functions

### `sushi.SushiColors(palette=2)`

| `palette` value | Returns |
| --- | --- |
| `2` | 2-color palette `["blue", "#E5001B"]` (K562 blue, ENCODE red) |
| `3` | `["blue", "#E5001B", "orange"]` |
| `4` | `["blue", "#5900E5", "#E5001B", "orange"]` |
| `5` | `["blue", "#5900E5", "#E5001B", "orange", "yellow"]` |
| `6` | `["blue", "#5900E5", "#E5001B", "orange", "yellow"]` |
| `7` | `["black", "blue", "#5900E5", "#E5001B", "orange", "yellow", "white"]` |
| `"list"` | `None` (use `_PALETTES` to inspect) |

The returned object is callable: `palette(N)` returns a list of N hex color strings.

### `sushi.maptocolors(vec, col, num=100, range=None)`

Map a numeric vector to a discrete color palette (R `maptocolors` semantics).

- `vec`: numeric values
- `col`: palette function (e.g. `SushiColors(7)`)
- `num`: number of bins
- `range`: optional `(min, max)` for fixed-range mapping; values outside are clipped

Returns: list of hex color strings, one per element.

### `sushi.maptolwd(lwdby, range=(1, 5))`

Map a numeric vector to line widths in the inclusive `range`.

### `sushi.opaque(color, transparency=0.5)`

Apply alpha channel to one or many colors. Returns RGBA hex (`#RRGGBBAA`, 9 chars). **Important**: pass through matplotlib before any other color processing.

### `sushi.convertstrandinfo(strandvector)`

Convert `["+", "-", ...]` to `[1, -1, ...]`. Numeric vectors are returned unchanged.

### `sushi.sortChrom(genome)`

Sort `(chrom, size)` table as chr1, chr2, ..., chrX, chrY, chrM (numeric first).

### `sushi.chromOffsets(genome, space=0.01)`

Compute per-chromosome offsets (start, stop) for multi-chrom Manhattan plots. Returns DataFrame with columns `chrom, size, start, stop`.

---

## Data loaders

All loaders return cached results (`functools.lru_cache`). First call
reads from CSV; subsequent calls are O(1).

| Function | Returns | Rows |
| --- | --- | --- |
| `sushi.data.Sushi_DNaseI_bedgraph()` | DataFrame (chrom, start, end, value) | 5,714 |
| `sushi.data.Sushi_ChIPSeq_CTCF_bedgraph()` | DataFrame | 27,576 |
| `sushi.data.Sushi_ChIPSeq_pol2_bed()` | DataFrame (BED6) | 208 |
| `sushi.data.Sushi_ChIPSeq_pol2_bedgraph()` | DataFrame | 7,528 |
| `sushi.data.Sushi_ChIPExo_CTCF_bedgraph()` | DataFrame | 6,292 |
| `sushi.data.Sushi_ChIPSeq_severalfactors_bed()` | DataFrame (BED7+) | 130 |
| `sushi.data.Sushi_RNASeq_K562_bedgraph()` | DataFrame | 1,073 |
| `sushi.data.Sushi_5C_bedpe()` | DataFrame (BEDPE11) | 4,787 |
| `sushi.data.Sushi_ChIAPET_pol2_bedpe()` | DataFrame (BEDPE10) | 48,634 |
| `sushi.data.Sushi_HiC_matrix()` | dict `{matrix, positions}` | 114×113 |
| `sushi.data.Sushi_GWAS_bed()` | DataFrame | 32,760 |
| `sushi.data.Sushi_genes_bed()` | DataFrame | 5 |
| `sushi.data.Sushi_transcripts_bed()` | DataFrame (with FPKM) | 143 |
| `sushi.data.Sushi_hg18_genome()` | DataFrame (chrom, size) | 22 |
