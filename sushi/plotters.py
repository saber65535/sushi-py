"""
sushi.plotters -- all 6 plot functions in one module.

Public API:
    plotBedgraph, plotBed, plotBedpe, plotHic, plotGenes, plotManhattan

Each function takes a matplotlib Axes as the first argument and modifies it
in place. See each function's docstring for parameters.
"""

from __future__ import annotations

from typing import Optional, Sequence, Callable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgba, to_hex
from matplotlib.patches import Polygon

from ._helpers import (
    SushiColors, maptocolors, maptolwd, opaque,
    convertstrandinfo, sortChrom, chromOffsets,
    _checkrow, _downsample_bedgraph,
    _bedgraph_step_polygon,
)


# ============================================================================
# plotBedgraph
# ============================================================================
def plotBedgraph(ax, signal, chrom: str, chromstart: int, chromend: int,
                 range_=None, color=None,
                 lwd: float = 1.0, linecolor=None,
                 addscale: bool = False,
                 overlay: bool = False, rescaleoverlay: bool = False,
                 transparency: float = 1.0,
                 flip: bool = False,
                 ymax: float = 1.04,
                 colorbycol=None,
                 **plot_kwargs):
    """Plot a bedGraph signal track.

    R signature:
        plotBedgraph(signal, chrom, chromstart, chromend, ...)

    Implementation note (2026-06-14 fix):
        The bedgraph polygon is a true step function (each bin spans its
        full [start, stop] at constant value), built via the helper
        ``_bedgraph_step_polygon``. The original (buggy) implementation
        concatenated ``[starts, stops_reversed]`` + ``[values, values_reversed]``
        which drew a long diagonal from the last bin's top to the first
        bin's top, producing the visual artifact of one giant triangle
        per panel. Do not revert this to the old ``np.concatenate`` form
        without understanding why it was wrong.


    See docs/API.md for full parameter documentation.
    """
    if color is None:
        color = SushiColors(2)(2)[0]
    if linecolor is None:
        linecolor = color
    if overlay:
        colorbycol = None

    # filter region
    if isinstance(signal, pd.DataFrame):
        chroms = signal.iloc[:, 0].astype(str).values
        starts = signal.iloc[:, 1].values
        stops = signal.iloc[:, 2].values
        vals = signal.iloc[:, 3].values
    else:
        chroms = np.asarray([str(r[0]) for r in signal])
        starts = np.asarray([r[1] for r in signal])
        stops = np.asarray([r[2] for r in signal])
        vals = np.asarray([r[3] for r in signal])

    mask = ((chroms == chrom) &
            ((starts > chromstart) & (starts < chromend) |
             (stops > chromstart) & (stops < chromend)))
    if mask.any():
        signaltrack = np.column_stack([starts[mask], stops[mask], vals[mask]]).astype(float)
    else:
        signaltrack = np.zeros((0, 3))

    if overlay and len(signaltrack) < 2:
        return ax
    if len(signaltrack) < 2:
        if range_ is None:
            range_ = (0, 1)
        ax.set_xlim(chromstart, chromend)
        ax.set_ylim(range_)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        return ax

    signaltrack = _downsample_bedgraph(signaltrack, max_points=8000)

    if flip:
        signaltrack = signaltrack.copy()
        signaltrack[:, 2] = -signaltrack[:, 2]

    if range_ is None:
        if not flip:
            range_ = (0, ymax * float(np.max(signaltrack[:, 2])))
        else:
            range_ = (ymax * float(np.min(signaltrack[:, 2])), 0)

    if not overlay:
        ax.set_xlim(chromstart, chromend)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
    # Set ylim always (matches R which always calls plot(..., ylim=range, ...))
    # regardless of overlay. The rescaleoverlay block below may override this.
    if not overlay:
        ax.set_ylim(range_)
    else:
        # In overlay mode, R still sets the ylim to the global range so
        # multiple overlays use the same scale.
        ax.set_ylim(range_)
    if not overlay:
        ax.autoscale(enable=False)

    if rescaleoverlay:
        _, xmax = ax.get_xlim()
        ymin, ymax_ = ax.get_ylim()
        if not flip:
            signaltrack[:, 2] = ymax_ * signaltrack[:, 2] / max(1e-9, np.max(np.abs(signaltrack[:, 2])))
        else:
            signaltrack[:, 2] = abs(ymin) * signaltrack[:, 2] / max(1e-9, np.max(np.abs(signaltrack[:, 2])))

    rgba = to_rgba(color, alpha=transparency)
    fill_color = to_hex(rgba)

    if colorbycol is None:
        # No gradient color -- fill the step-function polygon in a single
        # solid color (R does this with `polygon(signaltrack, col=finalcolor)`
        # where signaltrack is the (n, 2) x/value matrix).
        poly_x, poly_y = _bedgraph_step_polygon(signaltrack, baseline=range_[0])
        ax.fill(poly_x, poly_y, color=fill_color, linewidth=0)
        # Outline = top edge of step function only (no baseline retrace)
        starts = signaltrack[:, 0]
        stops = signaltrack[:, 1]
        values = signaltrack[:, 2]
        outline_x = np.concatenate([starts, [stops[-1]]])
        outline_y = np.concatenate([values, [values[-1]]])
        ax.plot(outline_x, outline_y, color=linecolor, linewidth=lwd, clip_on=False)
    else:
        # Gradient color background (R plotBedgraph with colorbycol).
        #
        # R's actual implementation (R/plotBedgraph.R lines 171-213):
        #   bgcol = maptocolors((1:100), colorbycol)  # 100 evenly-spaced colors
        #   tops = seq(plotbot, plottop, length.out=101)[2:101]  # 100 top y values
        #   bots = seq(plotbot, plottop, length.out=101)[1:100]  # 100 bottom y values
        #   for i in (3:(nrow(signaltrack)-1)) {  # each bin (skipping
        #       # the first row that R prepends as a baseline, and the last)
        #     ybots = bots[which(bots <= signaltrack[i, 2])]  # rect bottoms
        #     ytops = tops[which(tops <  signaltrack[i, 2])]  # rect tops
        #     ytops = c(ytops, signaltrack[i, 2])  # last top = actual value
        #     rect(xleft=signaltrack[i-1,1], ybottom=ybots,
        #          xright=signaltrack[i,1], ytop=ytops, col=xybgcol)
        #   }
        #
        # R's ylim is set to (plotbot, plottop) = (0, 1.04*max(values)) by
        # default. For each bin, it draws a column of 100 rects from
        # plotbot up to signaltrack[i, 2], with each rect's color
        # determined by its y position. R does NOT clip bins to plotbot;
        # the last rect's ytop is the bin value itself (which can exceed
        # plottop, making it use the last color = SushiColors(n)(100) = end
        # of the gradient).
        if not overlay:
            plotbot, plottop = range_
        else:
            plotbot, plottop = ax.get_ylim()
        if not flip:
            tops = np.linspace(plotbot, plottop, 101)[1:101]   # 100 tops
            bots = np.linspace(plotbot, plottop, 101)[:100]     # 100 bots
        else:
            tops_rev = np.linspace(plotbot, plottop, 101)
            tops = tops_rev[:100][::-1]  # 100 tops, descending
            bots = tops_rev[1:101][::-1]  # 100 bots, descending
        bgcol = colorbycol(100)  # 100 evenly-spaced colors
        # R prepends a baseline row: c(min(start), -0.00001) so the loop
        # starts at row 2 (i=2 in 1-based) and the bin value is at
        # signaltrack[i, 1]. Our 3-col matrix has the same row order.
        # We start at i=1 (Python) = second row in our 0-indexed matrix.
        for i in range(1, len(signaltrack)):
            bin_value = signaltrack[i, 2]
            # R clips to plotbot..plottop for the band lookup, then uses
            # the actual value for the last top.
            if not flip:
                bot_mask = bots <= bin_value
                top_mask = tops < bin_value
            else:
                bot_mask = bots >= bin_value
                top_mask = tops > bin_value
            ybots = bots[bot_mask]
            ytops = tops[top_mask]
            ytops = np.concatenate([ytops, [bin_value]]) if ytops.size else np.array([bin_value])
            if ybots.size == 0:
                continue
            xleft = signaltrack[i - 1, 0]
            xright = signaltrack[i, 0]   # NOTE: in R's signaltrack (n,2),
                                        # col 1 is value; col 0 is x. In our
                                        # (n, 3) matrix, cols are (start, stop,
                                        # value). R's i-th bin xrange is
                                        # signaltrack[i-1, 1] to signaltrack[i, 1].
                                        # We approximate with i-1's start and i's
                                        # start (which equal the bin boundary).
            xs = [xleft, xright, xright, xleft]
            for j in range(len(ytops)):
                ybottom = ybots[j] if j < len(ybots) else ybots[-1]
                ytop = ytops[j]
                color_idx = min(j, len(bgcol) - 1)
                ax.fill(xs + [xs[0]], [ybottom, ybottom, ytop, ytop, ybottom],
                        color=bgcol[color_idx], linewidth=0, edgecolor="none")
        # Outline on top
        starts = signaltrack[:, 0]
        stops = signaltrack[:, 1]
        values = signaltrack[:, 2]
        outline_x = np.concatenate([starts, [stops[-1]]])
        outline_y = np.concatenate([values, [values[-1]]])
        ax.plot(outline_x, outline_y, color=linecolor, linewidth=lwd, clip_on=False)

    if addscale:
        ax.text(1.0, 1.02, f"{range_[0]}-{range_[1]}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8)

    return ax


# ============================================================================
# plotBed
# ============================================================================
def plotBed(ax, beddata, chrom: str, chromstart: int, chromend: int,
            type: str = "region",
            colorby=None, colorbycol=None, colorbyrange=None,
            rownumber=None,
            row: str = "auto",
            height: float = 0.4,
            plotbg: str = "white",
            wiggle: float = 0.02,
            splitstrand: bool = False,
            numbins: int = 200, binsmoothing: int = 10,
            palettes=None,
            rowlabels=None,
            rowlabelcol="dodgerblue",
            rowlabelfont: int = 2,
            rowlabelcex: float = 1.0,
            maxrows: int = 1000000,
            color="navy",
            **plot_kwargs):
    """Plot BED-format data as regions, circles, or density heatmap."""
    if palettes is None:
        palettes = SushiColors(7)
    # If palettes is a SushiColors factory (callable), expand it to a list of
    # a single 7-color palette. Otherwise, expect palettes to already be a list
    # of palettes (R's list(SushiColors(7)) pattern).
    if callable(palettes) and not isinstance(palettes, (list, tuple)):
        palettes = [palettes(7)]

    if isinstance(beddata, pd.DataFrame):
        df = beddata.copy()
    else:
        df = pd.DataFrame(beddata)

    if df.shape[1] >= 6:
        col5 = df.iloc[:, 5]
        try:
            int(col5.iloc[0])
        except (TypeError, ValueError):
            # Only convert strand to numeric if splitstrand=True needs
            # the +/- for sorting. Otherwise keep as string to avoid
            # "Invalid value for dtype 'str'" errors downstream.
            if splitstrand:
                df.iloc[:, 5] = convertstrandinfo(col5.tolist())

    df = df[(df.iloc[:, 0].astype(str) == chrom) &
            (((df.iloc[:, 1] > chromstart) & (df.iloc[:, 1] < chromend)) |
             ((df.iloc[:, 2] > chromstart) & (df.iloc[:, 2] < chromend)))].reset_index(drop=True)

    # Caller-supplied colorby/rownumber lists match the ORIGINAL df
    # length, but df has been filtered to the requested region. We
    # must truncate these lists to match the filtered df length, otherwise
    # pandas raises "Length of values does not match length of index".
    # R's plotBed doesn't have this issue because it works with R vectors
    # that are subset together with the bed data.
    if colorby is not None and len(colorby) > len(df):
        colorby = list(colorby)[:len(df)]
    if rownumber is not None and len(rownumber) > len(df):
        rownumber = list(rownumber)[:len(df)]

    if colorby is not None:
        df["plotcolor"] = maptocolors([float(c) for c in colorby], colorbycol,
                                       num=100, rng=colorbyrange)
    else:
        df["plotcolor"] = color
    if rownumber is not None:
        df["plotrow"] = rownumber
    else:
        df["plotrow"] = 1

    wiggle_abs = abs(chromend - chromstart) * wiggle

    def _local_checkrow(data, alldata, maxrows, wiggle):
        thestart = min(float(data[1]), float(data[2])) - wiggle
        thestop = max(float(data[1]), float(data[2])) + wiggle
        for row in range(1, maxrows + 1):
            occ = alldata[alldata["plotrow"] == row]
            if occ.empty:
                return row
            occ_s = occ.iloc[:, 1].values
            occ_e = occ.iloc[:, 2].values
            overlap = (
                ((thestart >= occ_s) & (thestart <= occ_e)) |
                ((thestop >= occ_s) & (thestop <= occ_e)) |
                ((thestart <= occ_s) & (thestop >= occ_e))
            )
            if not overlap.any():
                return row
        return None

    if row == "auto":
        df = df.sample(frac=1).reset_index(drop=True)
        df["plotrow"] = 0
        if splitstrand and df.shape[1] >= 6:
            plus = df[df.iloc[:, 5] > 0].copy()
            minus = df[df.iloc[:, 5] < 0].copy()
            for i in range(len(plus)):
                plus.at[plus.index[i], "plotrow"] = _local_checkrow(
                    plus.iloc[i].values, plus, maxrows, wiggle_abs)
            for i in range(len(minus)):
                minus.at[minus.index[i], "plotrow"] = _local_checkrow(
                    minus.iloc[i].values, minus, maxrows, wiggle_abs)
            minus["plotrow"] = -1 * minus["plotrow"]
            df = pd.concat([plus, minus], ignore_index=True)
        else:
            for i in range(len(df)):
                df.at[df.index[i], "plotrow"] = _local_checkrow(
                    df.iloc[i].values, df, maxrows, wiggle_abs)

    numberofrows = max(df["plotrow"]) if len(df) > 0 else 0
    if len(df) == 0:
        print("Warning: No bed elements in selected range")

    weight = {1: "normal", 2: "bold", 3: "italic", 4: "bold italic"}

    if type == "region":
        y_low = (min(df["plotrow"]) - height * 1.2) if len(df) > 0 else 0
        y_high = (min(maxrows, max(df["plotrow"]) if len(df) > 0 else 1) + height * 1.2)
        ax.set_xlim(chromstart, chromend)
        ax.set_ylim(y_low, y_high)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.add_patch(plt.Rectangle((chromstart, y_low),
                                    chromend - chromstart, y_high - y_low,
                                    facecolor=plotbg, edgecolor="none", zorder=0))
        if splitstrand:
            ax.axhline(0, color="black", linewidth=0.5)
        for i in range(len(df)):
            line = df["plotrow"].iloc[i]
            col = df["plotcolor"].iloc[i] if isinstance(df["plotcolor"].iloc[i], str) else "navy"
            ax.add_patch(plt.Rectangle(
                (df.iloc[i, 1], line - height),
                df.iloc[i, 2] - df.iloc[i, 1],
                2 * height,
                facecolor=col, edgecolor=col, linewidth=0.5, zorder=2))

    elif type == "circles":
        y_low = (min(df["plotrow"]) - height * 1.2) if len(df) > 0 else 0
        y_high = (min(maxrows, max(df["plotrow"]) if len(df) > 0 else 1) + height * 1.2)
        ax.set_xlim(chromstart, chromend)
        ax.set_ylim(y_low, y_high)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.add_patch(plt.Rectangle((chromstart, y_low),
                                    chromend - chromstart, y_high - y_low,
                                    facecolor=plotbg, edgecolor="none", zorder=0))
        if len(df) > 0:
            x_centers = (df.iloc[:, 1].values + df.iloc[:, 2].values) / 2.0
            colors = [df["plotcolor"].iloc[i] if isinstance(df["plotcolor"].iloc[i], str) else "navy"
                       for i in range(len(df))]
            ax.scatter(x_centers, df["plotrow"].values, s=height * 200,
                       c=colors, edgecolors="none", zorder=2)

    elif type == "density":
        # R uses a single row by default for density, but if row="given" and
        # rownumber is supplied, it iterates rows.
        if row == "given" and rownumber is not None:
            # group df by rownumber, render each row with its own palette
            unique_rows = sorted(set(rownumber))
            numberofrows = len(unique_rows)
            row_to_palette = {r: palettes[i] if i < len(palettes) else palettes[0]
                              for i, r in enumerate(unique_rows)}
        else:
            unique_rows = [1]
            numberofrows = 1
            row_to_palette = {1: palettes[0] if palettes else SushiColors(7)(7)}
        if len(palettes) < numberofrows:
            palettes = [palettes[0]] * numberofrows
        ax.set_xlim(chromstart, chromend)
        ax.set_ylim(0.5, numberofrows + 0.5)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        bins = np.linspace(chromstart, chromend, numbins + 1)
        for row_idx, row_num in enumerate(unique_rows, start=1):
            if row == "given" and rownumber is not None:
                mask = (np.asarray(rownumber) == row_num)
                sub = df[mask]
            else:
                sub = df
            if len(sub) == 0:
                continue
            x_centers = np.clip((sub.iloc[:, 1].values + sub.iloc[:, 2].values) / 2.0,
                                chromstart, chromend)
            counts, _ = np.histogram(x_centers, bins=bins)
            smooth = np.convolve(counts, np.ones(binsmoothing), mode="same")
            # The palette is the row's palette, expanded to num+1 colors
            row_palette = row_to_palette[row_num]
            if callable(row_palette):
                colors = maptocolors(smooth.tolist(), row_palette, num=100)
            else:
                # row_palette is already a list of colors
                colors = maptocolors(smooth.tolist(), row_palette, num=100)
            polys = []
            for i in range(len(bins) - 1):
                polys.append([[bins[i], row_idx - 0.5], [bins[i + 1], row_idx - 0.5],
                              [bins[i + 1], row_idx + 0.5], [bins[i], row_idx + 0.5]])
            pc = PolyCollection(polys, closed=True, facecolors=colors, edgecolors="none")
            ax.add_collection(pc)

    if rowlabels is not None and len(df) > 0:
        for i, lbl in enumerate(rowlabels[:numberofrows]):
            y_pos = i + 1
            c = rowlabelcol[i] if isinstance(rowlabelcol, list) and i < len(rowlabelcol) else rowlabelcol
            ax.annotate(str(lbl), xy=(chromstart - (chromend - chromstart) * 0.005, y_pos),
                        xycoords=("data", "data"),
                        ha="right", va="center",
                        fontsize=10 * rowlabelcex,
                        fontweight=weight.get(rowlabelfont, "normal"),
                        color=c, annotation_clip=False)
    return ax


# ============================================================================
# plotBedpe
# ============================================================================
def _plot_loop(ax, x1, x2, height, offset=0, flip=False, color="black", lwd=1.0):
    posneg = -1 if flip else 1
    if flip:
        offset = -offset
    x1, x2 = min(x1, x2), max(x1, x2)
    if x1 == x2:
        return
    mid = (x1 + x2) / 2.0
    a = height / ((mid - x1) * (mid - x2))
    xs = np.linspace(x1, x2, 200)
    ys = offset + posneg * a * (xs - x1) * (xs - x2)
    ax.plot(xs, ys, color=color, linewidth=lwd, solid_capstyle="butt")


def _plot_ribbon(ax, s1, e1, s2, e2, height, offset=0, flip=False,
                  color="black", border="black"):
    posneg = -1 if flip else 1
    if flip:
        offset = -offset
    pairs = [(s1, e2), (e1, s2)]
    if abs(s2 - s1) > 0 and abs(e2 - e1) > 0:
        xoffset = abs(abs(s2 - s1) - abs(e2 - e1)) / max(1e-9, (e2 - s1))
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        xrange = abs(xmax - xmin)
        yrange = abs(ymax - ymin)
        halfwidth = xoffset * yrange / 2.0
    else:
        halfwidth = 0
    verts_x, verts_y = [], []
    for idx, (xa, xb) in enumerate(pairs):
        xa, xb = min(xa, xb), max(xa, xb)
        if xa == xb:
            continue
        mid = (xa + xb) / 2.0
        h = height + halfwidth if idx == 0 else height - halfwidth
        a = h / ((mid - xa) * (mid - xb))
        if idx == 0:
            xs = np.linspace(xa, xb, 50)
        else:
            xs = np.linspace(xb, xa, 50)
        verts_x.extend(xs.tolist())
        for x in xs:
            verts_y.append(offset + posneg * a * (x - xa) * (x - xb))
    if verts_x:
        ax.fill(verts_x, verts_y, color=color, edgecolor=border, linewidth=0.5)


def plotBedpe(ax, bedpedata, chrom: str, chromstart: int, chromend: int,
              heights,
              color="black",
              colorby=None, colorbycol=None, colorbyrange=None,
              border=None,
              lwdby=None, lwdrange=(1, 5),
              offset: float = 0,
              flip: bool = False,
              lwd: float = 1.0,
              plottype: str = "loops",
              maxrows: int = 10000,
              height: float = 0.3,
              ymax: float = 1.04,
              **plot_kwargs):
    """Plot BEDPE interactions as loops, ribbons, or lines."""
    if isinstance(bedpedata, pd.DataFrame):
        df = bedpedata.copy()
    else:
        df = pd.DataFrame(bedpedata)

    # Note: in v0.1.0/v0.1.1, this function called convertstrandinfo on
    # cols 8/9, but the converted numeric strands were never used
    # (R's plotBedpe also only uses start1/start2/stop1/stop2 for
    # geometry, not strand1/strand2 for direction). Removing the call
    # lets the strand columns stay as strings (preserving dtype 'str'),
    # which is what plotBedpe and downstream callers expect.
    # Callers who need numeric strand (e.g. plotGenes) call
    # convertstrandinfo themselves on a column they have already
    # selected.

    if len(df) == 0:
        ax.set_xlim(chromstart, chromend)
        ax.set_ylim(0, 1)
        for spine in ax.spines.values():
            spine.set_visible(False)
        return ax

    df = df.iloc[:, :6].copy()
    df.columns = ["chrom1", "start1", "stop1", "chrom2", "start2", "stop2"]
    df["start1"] = df[["start1", "stop1"]].min(axis=1)
    df["stop1"] = df[["start1", "stop1"]].max(axis=1)
    df["start2"] = df[["start2", "stop2"]].min(axis=1)
    df["stop2"] = df[["start2", "stop2"]].max(axis=1)

    if plottype in ("loops", "ribbons"):
        df["heights"] = list(heights)
    df["color"] = color
    if colorby is not None:
        df["colorbyvalue"] = list(colorby)
    if isinstance(color, (list, tuple, np.ndarray)):
        df["color"] = list(color)
    df["lwd"] = lwd
    if lwdby is not None:
        df["lwdby"] = list(lwdby)
    df["pos1"] = (df["start1"] + df["stop1"]) / 2.0
    df["pos2"] = (df["start2"] + df["stop2"]) / 2.0

    df = df[(df["chrom1"].astype(str) == chrom) &
            (df["pos1"] > chromstart) & (df["pos1"] < chromend) &
            (df["chrom2"].astype(str) == chrom) &
            (df["pos2"] > chromstart) & (df["pos2"] < chromend)].reset_index(drop=True)

    if len(df) == 0:
        ax.set_xlim(chromstart, chromend)
        ax.set_ylim(0, 1)
        for spine in ax.spines.values():
            spine.set_visible(False)
        return ax

    if colorby is not None:
        df["color"] = maptocolors(df["colorbyvalue"].astype(float).tolist(),
                                    colorbycol, num=100, rng=colorbyrange)
        if colorbyrange is None:
            colorbyrange = (float(np.nanmin(df["colorbyvalue"])),
                             float(np.nanmax(df["colorbyvalue"])))

    if lwdby is not None:
        df["lwd"] = maptolwd(df["lwdby"].tolist(), range=lwdrange)

    if plottype == "lines":
        df["distance"] = (df["pos2"] - df["pos1"]).abs()
        df = df.sort_values("distance", ascending=False).reset_index(drop=True)
        wiggle = abs(chromend - chromstart) * 0.04
        df["plotrow"] = 0

        def _local_checkrow_bp(data, alldata, maxrows, wiggle):
            # data is [start1, stop1] (2 elements). Was incorrectly using
            # data[1]/data[2] which would IndexError on the original caller.
            # Fixed in v0.1.4: use data[0] (start1) and data[1] (stop1).
            thestart = min(float(data[0]), float(data[1])) - wiggle
            thestop = max(float(data[0]), float(data[1])) + wiggle
            for row in range(1, maxrows + 1):
                occ = alldata[alldata["plotrow"] == row]
                if occ.empty:
                    return row
                starts = np.concatenate([occ["start1"].values, occ["start2"].values])
                stops = np.concatenate([occ["stop1"].values, occ["stop2"].values])
                overlap = (
                    ((thestart >= starts) & (thestart <= stops)) |
                    ((thestop >= starts) & (thestop <= stops)) |
                    ((thestart <= starts) & (thestop >= stops))
                )
                if not overlap.any():
                    return row
            return None

        for i in range(len(df)):
            df.at[i, "plotrow"] = _local_checkrow_bp(
                [df.iloc[i]["start1"], df.iloc[i]["stop1"]],
                df, maxrows, wiggle)
        toprow = min(maxrows, max(df["plotrow"]) if len(df) else 1)
        if not flip:
            ax.set_ylim(0, toprow + 1)
        else:
            ax.set_ylim(-toprow - 1, 0)
        ax.set_xlim(chromstart, chromend)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        for i in range(len(df)):
            line = df["plotrow"].iloc[i] if not flip else -df["plotrow"].iloc[i]
            h_signed = -height if flip else height
            col = df["color"].iloc[i] if isinstance(df["color"].iloc[i], str) else "black"
            lwdi = df["lwd"].iloc[i] if isinstance(df["lwd"].iloc[i], (int, float)) else 1.0
            x0 = min(df["start1"].iloc[i], df["stop1"].iloc[i])
            x1 = max(df["start1"].iloc[i], df["stop1"].iloc[i])
            ax.add_patch(plt.Rectangle((x0, line - h_signed),
                                        x1 - x0, 2 * h_signed,
                                        facecolor=col, edgecolor=col, linewidth=lwdi))
            x0 = min(df["start2"].iloc[i], df["stop2"].iloc[i])
            x1 = max(df["start2"].iloc[i], df["stop2"].iloc[i])
            ax.add_patch(plt.Rectangle((x0, line - h_signed),
                                        x1 - x0, 2 * h_signed,
                                        facecolor=col, edgecolor=col, linewidth=lwdi))
            ax.plot([df["stop1"].iloc[i], df["start2"].iloc[i]],
                    [line, line], color=col, marker="", linewidth=lwdi, solid_capstyle="butt")
        return ax

    if plottype in ("loops", "ribbons"):
        if not flip:
            ax.set_ylim(0, (max(df["heights"]) + offset) * ymax)
        else:
            ax.set_ylim(-max(df["heights"] - offset) * ymax, 0)
        ax.set_xlim(chromstart, chromend)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        for i in range(len(df)):
            x1 = df["pos1"].iloc[i]
            x2 = df["pos2"].iloc[i]
            h = df["heights"].iloc[i]
            col = df["color"].iloc[i] if isinstance(df["color"].iloc[i], str) else "black"
            lwdi = df["lwd"].iloc[i] if isinstance(df["lwd"].iloc[i], (int, float)) else 1.0
            if plottype == "loops":
                _plot_loop(ax, x1, x2, h, offset=offset, flip=flip,
                           color=col, lwd=lwdi)
            else:
                bordercol = border if border is not None else col
                s1 = min(df["start1"].iloc[i], df["start2"].iloc[i])
                s2 = max(df["start1"].iloc[i], df["start2"].iloc[i])
                e1 = min(df["stop1"].iloc[i], df["stop2"].iloc[i])
                e2 = max(df["stop1"].iloc[i], df["stop2"].iloc[i])
                _plot_ribbon(ax, s1, e1, s2, e2, h, offset=offset, flip=flip,
                             color=col, border=bordercol)
        return ax

    return ax


# ============================================================================
# plotHic
# ============================================================================
def plotHic(ax, hicdata, chrom: str, chromstart: int, chromend: int,
            max_y: float = 30.0,
            zrange=None,
            palette=None,
            flip: bool = False):
    """Plot a Hi-C interaction matrix as a lower-triangle heatmap."""
    if palette is None:
        palette = SushiColors(7)
    if hasattr(hicdata, "index") and hasattr(hicdata, "columns"):
        rows = np.asarray(hicdata.index)
        cols = np.asarray(hicdata.columns)
        # If the DataFrame is rectangular (R's lower-triangular storage),
        # pad with NaN to make it square before any further processing.
        if len(rows) != len(cols):
            n_rows, n_cols = hicdata.shape
            full = _np.full((n_rows, n_rows), _np.nan)
            for i in range(n_rows):
                for j in range(i, n_cols):
                    full[i, j] = hicdata.iloc[i, j]
            cols = rows  # use rows for both axes
            hicdata = __import__("pandas").DataFrame(full, index=rows, columns=rows)
        rows_mask = (rows >= chromstart) & (rows <= chromend)
        cols_mask = (cols >= chromstart) & (cols <= chromend)
        hicregion = hicdata.iloc[rows_mask, cols_mask].values
    else:
        hicregion = np.asarray(hicdata)
        rows_mask = cols_mask = slice(None)
    # If hicregion is rectangular (e.g. R's lower-triangular storage where
    # the last row/col is shorter), use the smallest dimension for nbins so
    # that we never index out of bounds. R's plotHic handles this naturally.
    n_rows, n_cols = hicregion.shape
    nbins = min(n_rows, n_cols)
    stepsize = abs(chromstart - chromend) / (2 * nbins) if nbins else 0
    if zrange is None:
        min_z = float(np.nanmin(hicregion))
        max_z = float(np.nanmax(hicregion))
    else:
        min_z, max_z = zrange
    # Replace NaN with min_z so maptocolors can bucket them safely.
    # Then post-replace cells in the color matrix that are NaN in the original
    # hicregion with a transparent color so the upper-triangle doesn't render.
    nan_mask = np.isnan(hicregion)
    flat_for_color = np.where(nan_mask, min_z, hicregion).flatten()  # don't NaN-clobber for maptocolors
    hicmcol_flat = maptocolors(flat_for_color.tolist(), palette, num=100,
                                 rng=(min_z, max_z))
    hicmcol = np.array(hicmcol_flat).reshape(hicregion.shape)
    if not flip:
        ax.set_xlim(chromstart, chromend)
        ax.set_ylim(0, max_y)
    else:
        ax.set_xlim(chromstart, chromend)
        ax.set_ylim(-max_y, 0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    for rownum in range(nbins):
        y = -0.5 if not flip else 0.5
        x = chromstart + (rownum * 2 * stepsize) - (stepsize * 2)
        for colnum in range(rownum, nbins):
            x = x + stepsize
            y = (y + 0.5) if not flip else (y - 0.5)
            xs = [x - stepsize, x, x + stepsize, x, x - stepsize]
            ys = [y, y + 0.5, y, y - 0.5, y]
            # Skip upper-triangle cells (where hicdata was NaN)
            if nan_mask[rownum, colnum]:
                continue
            cell_color = hicmcol[rownum, colnum]
            ax.fill(xs, ys, color=cell_color, edgecolor="none")
    return [min_z, max_z], palette


# ============================================================================
# plotGenes
# ============================================================================
def plotGenes(ax, geneinfo, chrom: str, chromstart: int, chromend: int,
              types,
              colorby=None, colorbycol=None, colorbyrange=None,
              labeltext: bool = True, labeloffset: float = 0.4,
              fontsize: float = 0.7, fonttype: int = 2,
              labelat: str = "middle",
              arrowlength: float = 0.025,
              bheight: float = 0.2,
              lheight: float = 0.05,
              bentline: bool = True,
              plotgenetype: str = "box",
              packrow: bool = True,
              maxrows: int = 10000,
              color="navy",
              wigglefactor: float = 0.02,
              **plot_kwargs):
    """Plot a gene/transcript model."""
    if isinstance(geneinfo, pd.DataFrame):
        df = geneinfo.copy()
    else:
        df = pd.DataFrame(geneinfo)

    if df.shape[1] >= 6:
        try:
            int(df.iloc[0, 5])
        except (TypeError, ValueError):
            df.iloc[:, 5] = convertstrandinfo(df.iloc[:, 5].tolist())

    df["types"] = types
    df = df.dropna(subset=[df.columns[1]]).reset_index(drop=True)

    if colorby is not None:
        if colorbyrange is None:
            colorbyrange = (float(np.nanmin(colorby)), float(np.nanmax(colorby)))
        df["colors"] = maptocolors(colorby, colorbycol, num=100, rng=colorbyrange)
    else:
        df["colors"] = color

    if chrom is None or chromstart is None or chromend is None:
        chrom = df.iloc[0, 0]
        chromstart = float(np.nanmin(df.iloc[:, [1, 2]].values))
        chromend = float(np.nanmax(df.iloc[:, [1, 2]].values))
        chromstart -= abs(chromend - chromend) * 0.04
        chromend += abs(chromend - chromend) * 0.04

    bprange = abs(chromend - chromstart)
    wiggle = bprange * wigglefactor

    if df.shape[1] >= 4:
        names = df.iloc[:, 3].unique()
    else:
        names = np.arange(len(df))
    transcript_starts, transcript_stops, transcript_strands = [], [], []
    for n in names:
        sub = df[df.iloc[:, 3] == n] if df.shape[1] >= 4 else df.iloc[[0]]
        ts = float(np.nanmin(sub.iloc[:, [1, 2]].values))
        te = float(np.nanmax(sub.iloc[:, [1, 2]].values))
        transcript_starts.append(ts)
        transcript_stops.append(te)
        transcript_strands.append(int(sub.iloc[0, 5]) if df.shape[1] >= 6 else 1)
    transcript_info = pd.DataFrame({
        "names": names,
        "starts": transcript_starts,
        "stops": transcript_stops,
        "sizes": [te - ts for ts, te in zip(transcript_starts, transcript_stops)],
        "strand": transcript_strands,
    }).sort_values("sizes", ascending=False).reset_index(drop=True)

    if packrow:
        def _local_checkrow_g(data, alldata, maxrows, wiggle):
            thestart = min(float(data[0]), float(data[1])) - wiggle
            thestop = max(float(data[0]), float(data[1])) + wiggle
            for row in range(1, maxrows + 1):
                occ = alldata[alldata["plotrow"] == row]
                if occ.empty:
                    return row
                starts = occ["starts"].values
                stops = occ["stops"].values
                overlap = (
                    ((thestart >= starts) & (thestart <= stops)) |
                    ((thestop >= starts) & (thestop <= stops)) |
                    ((thestart <= starts) & (thestop >= stops))
                )
                if not overlap.any():
                    return row
            return None

        transcript_info["plotrow"] = 0
        for i in range(len(transcript_info)):
            transcript_info.at[i, "plotrow"] = _local_checkrow_g(
                [transcript_info.iloc[i]["starts"], transcript_info.iloc[i]["stops"]],
                transcript_info, maxrows, wiggle)
    else:
        transcript_info["plotrow"] = np.arange(len(transcript_info)) + 1

    transcript_info = transcript_info[transcript_info["plotrow"].notna()]
    transcript_info = transcript_info[
        ((transcript_info["starts"] > chromstart) & (transcript_info["starts"] < chromend)) |
        ((transcript_info["stops"] > chromstart) & (transcript_info["stops"] < chromend))
    ].reset_index(drop=True)

    if len(transcript_info) == 0:
        toprow = 1
    else:
        toprow = max(transcript_info["plotrow"])

    ax.set_xlim(chromstart, chromend)
    ax.set_ylim(0.5, toprow + 0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    if len(transcript_info) == 0:
        return [None, None]

    weight = {1: "normal", 2: "bold", 3: "italic", 4: "bold italic"}

    for i in range(len(transcript_info)):
        tname = transcript_info.iloc[i]["names"]
        sub = df[df.iloc[:, 3] == tname] if df.shape[1] >= 4 else df.iloc[[0]]
        col = sub["colors"].iloc[0] if isinstance(sub["colors"].iloc[0], str) else "navy"
        yval = transcript_info.iloc[i]["plotrow"]
        strand = int(transcript_info.iloc[i]["strand"])

        if labeltext:
            label_loc = (float(np.nanmax(sub.iloc[:, [1, 2]].values)) +
                         float(np.nanmin(sub.iloc[:, [1, 2]].values))) / 2.0
            if (labelat == "start" and strand == 1) or (labelat == "end" and strand == -1):
                label_loc = float(np.nanmin(sub.iloc[:, [1, 2]].values))
            if (labelat == "start" and strand == -1) or (labelat == "end" and strand == 1):
                label_loc = float(np.nanmax(sub.iloc[:, [1, 2]].values))
            adj = 1.0 if strand == 1 else 0.0
            ax.annotate(str(tname), xy=(label_loc, yval + labeloffset),
                        ha="right" if adj == 1.0 else "left",
                        va="center",
                        fontsize=10 * fontsize,
                        fontweight=weight.get(fonttype, "bold"),
                        annotation_clip=False)
            ax.annotate("", xy=(label_loc + strand * bprange * arrowlength / 4,
                                  yval + labeloffset),
                        xytext=(label_loc, yval + labeloffset),
                        arrowprops=dict(arrowstyle="->", lw=1, color=col))

        for typ in sub["types"].unique():
            sub_typ = sub[sub["types"] == typ]
            for _, ex in sub_typ.iterrows():
                h = bheight / 2.0 if typ == "utr" else bheight
                if plotgenetype == "box":
                    ax.add_patch(plt.Rectangle(
                        (ex.iloc[1], yval - h), ex.iloc[2] - ex.iloc[1], 2 * h,
                        facecolor=col, edgecolor=col, linewidth=0.5))
                elif plotgenetype == "arrow":
                    # R polygon: x=c(ex[2], ex[2]+offset, ex[3]+offset, ex[3], ex[3]+offset, ex[2]+offset, ex[2])
                    # R: offset = bprange * arrowlength * strand * -1
                    # When strand=+1: offset is NEGATIVE, so ex[2]+offset is LEFT of ex[2]
                    # When strand=-1: offset is POSITIVE, so ex[2]+offset is RIGHT of ex[2]
                    # This forms an "arrow head" pointing in the direction of strand.
                    offset = bprange * arrowlength * strand * -1
                    poly_x = [ex.iloc[1], ex.iloc[1] + offset, ex.iloc[2] + offset,
                              ex.iloc[2], ex.iloc[2] + offset, ex.iloc[1] + offset,
                              ex.iloc[1]]
                    poly_y = [yval, yval + h, yval + h, yval, yval - h, yval - h, yval]
                    # R's polygon() with the same x,y coords produces a flat-top
                    # arrow shape.  matplotlib's ax.fill() with a self-intersecting
                    # polygon can give weird results depending on fill rule.  Use
                    # Polygon with explicit fill rule for clean rendering.
                    from matplotlib.patches import Polygon as MplPolygon
                    arrow = MplPolygon(list(zip(poly_x, poly_y)),
                                          closed=True, facecolor=col, edgecolor=col,
                                          linewidth=0.5)
                    ax.add_patch(arrow)

    return [colorbyrange, colorbycol]


# ============================================================================
# plotManhattan
# ============================================================================
def plotManhattan(ax, bedfile, pvalues=None, chrom=None, chromstart=None,
                   chromend=None, genome=None,
                   col=SushiColors(5), space: float = 0.01,
                   ymax: float = 1.04, **plot_kwargs):
    """Plot a Manhattan-style -log10(pvalue) scatter."""
    if isinstance(bedfile, pd.DataFrame):
        df = bedfile.copy()
    else:
        df = pd.DataFrame(bedfile)

    if pvalues is None:
        if "pvalue" in df.columns:
            pvalues = df["pvalue"].values
        elif df.shape[1] >= 5:
            pvalues = df.iloc[:, 4].values
        else:
            pvalues = df.iloc[:, 1].values
    pvalues = np.asarray(pvalues, dtype=float)
    neglog_p = -np.log10(pvalues)

    if genome is not None:
        co = chromOffsets(genome, space)
        if isinstance(col, type(SushiColors(5))):
            colors = col(len(co))
        else:
            colors = list(col) * max(1, (len(co) // max(1, len(col)) + 1))
        df = df[df.iloc[:, 0].isin(co["chrom"].tolist())].reset_index(drop=True)
        if len(pvalues) > len(df):
            pvalues = pvalues[:len(df)]
        co_number = 0
        co_numbers = np.zeros(len(df), dtype=int)
        for i in range(len(co)):
            co_number += 1
            if co_number > len(colors):
                co_number = 1
            mask = (df.iloc[:, 0] == co.iloc[i]["chrom"]).values
            # Cast to int (the position column is int64 in the original
            # data; adding the float chrom offset and assigning back to an
            # int column raises "Invalid value for dtype 'int64'" in pandas).
            df.loc[mask, df.columns[1]] = (df.loc[mask, df.columns[1]].astype(float)
                                              + float(co.iloc[i]["start"])).astype(int)
            co_numbers[mask] = co_number
        sizes = (genome.iloc[:, 1].astype(float).values
                  if isinstance(genome, pd.DataFrame)
                  else np.asarray([g[1] for g in genome]))
        cumsums = np.cumsum(sizes)
        spacer = cumsums[-1] * space
        ax.scatter(df.iloc[:, 1].values, neglog_p,
                   c=[colors[i - 1] for i in co_numbers],
                   s=20, edgecolors="none", **plot_kwargs)
        ax.set_xlim(0, float(co["stop"].max()) + spacer)
        ax.set_ylim(0, float(np.nanmax(neglog_p)) * ymax)
        # R Manhattan adds chr number labels at the center of each chr.
        # We use the cumulative chromOffsets start/stop to position them.
        for i, row in co.iterrows():
            center = (row["start"] + row["stop"]) / 2.0
            # Strip "chr" prefix to match R's "1", "2", ... labels
            label = str(row["chrom"]).replace("chr", "")
            ax.text(center, -0.02 * ax.get_ylim()[1], label,
                    ha="center", va="top", fontsize=8,
                    transform=ax.transData, clip_on=False)
        # R sets y-axis ticks at pretty() values (4 6 8 10 12 14), with
        # 0 explicitly excluded. Use a FixedLocator with steps of 2.
        from matplotlib.ticker import FixedLocator, FixedFormatter
        ymax_val = ax.get_ylim()[1]
        # R uses axis(..., at=pretty(...), labels=-1*pretty(...)) which excludes
        # 0 (the minimum is always skipped on a Manhattan axis). Match that by
        # starting at step (2) when ymax>=2, or at 1 when ymax<2.
        tick_step = 2
        first_tick = tick_step  # R skips 0, so first tick is at tick_step
        last_tick = int(ymax_val // tick_step * tick_step)
        if last_tick < tick_step:
            last_tick = int(ymax_val)
            if last_tick < 1:
                first_tick = 0
        tick_positions = list(range(first_tick, last_tick + 1, tick_step))
        ax.yaxis.set_major_locator(FixedLocator(tick_positions))
        ax.yaxis.set_major_formatter(FixedFormatter([str(t) for t in tick_positions]))
        ax.tick_params(axis="y", labelsize=8)
        ax.tick_params(axis="x", labelbottom=False)
        # Add "-log10(P)" axis label (R uses mtext; matplotlib uses ylabel)
        ax.set_ylabel("log10(P)", fontsize=10, fontweight="bold")
    else:
        if chrom is not None:
            df = df[df.iloc[:, 0] == chrom].reset_index(drop=True)
            # Truncate pvalues to match filtered df (was missing in v0.1.5;
            # fixed in v0.1.6: scatter raises "x and y must be the same size")
            if len(neglog_p) > len(df):
                neglog_p = neglog_p[:len(df)]
        ax.scatter(df.iloc[:, 1].values, neglog_p,
                   c=col if not callable(col) else col(len(df)),
                   s=20, edgecolors="none", **plot_kwargs)
        ax.set_xlim(chromstart, chromend)
        ax.set_ylim(0, float(np.nanmax(neglog_p)) * ymax)
    # Note: R Manhattan examples also add mtext("chromosome", side=1)
    # and mtext("log10(P)", side=2). These are vignette-level additions,
    # not core plotManhattan behavior.


__all__ = ["plotBedgraph", "plotBed", "plotBedpe", "plotHic", "plotGenes", "plotManhattan"]
