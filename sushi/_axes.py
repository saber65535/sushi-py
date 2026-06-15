"""
sushi._axes -- axis labels, color legends, zoom-region overlays.

All functions take a matplotlib Axes and modify it in place. They never
create new figures.

Functions:
    labelgenome(ax, chrom, chromstart, chromend, ...)
    addlegend(ax, range_val, palette, ...)
    zoombox(ax, zoomregion=None, ...)
    zoomsregion(ax, region, ...)
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple, Callable

import numpy as np
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from ._helpers import (
    SushiColors,
    maptocolors,
    chromOffsets,
)


# ---------------------------------------------------------------------------
# labelgenome
# ---------------------------------------------------------------------------
def labelgenome(ax, chrom: str, chromstart: int, chromend: int,
                genome=None, space: float = 0.01,
                scale: str = "bp", side: int = 1,
                scipen: int = 20, n: int = 5,
                chromfont: int = 2, chromadjust: float = 0.015,
                chromcex: float = 1.0, chromline: float = 0.5,
                scalefont: int = 2, scaleadjust: float = 0.985,
                scalecex: float = 1.0, scaleline: float = 0.5,
                line: float = 0.18,
                edgeblankfraction: float = 0.10,
                **axis_kwargs):
    """Draw genomic coordinate axis ticks on an Axes.

    R signature:
        labelgenome(chrom, chromstart, chromend, ...)

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to draw on (the parent plot's x-axis).
    chrom, chromstart, chromend : str, int, int
        The region. If ``genome`` is None (single-chrom mode) the chrom name
        is drawn at chromadjust and the scale label ("bp", "Kb", "Mb") at
        scaleadjust. If ``genome`` is provided (multi-chrom mode) the
        chromosome names are drawn as axis tick labels at chrom centers.
    scale : {"bp", "Kb", "Mb"}
        Unit shown on the scale label (single-chrom mode only).
    n : int
        Desired number of tick labels.
    edgeblankfraction : float
        Fraction of the x-range to leave blank on each side for chrom/scale
        labels in single-chrom mode.
    """
    xmin, xmax = ax.get_xlim()
    rng = abs(xmax - xmin)

    if genome is not None:
        co = chromOffsets(genome, space)
        chromcenters = (co["start"] + co["stop"]) / 2.0
        labels = [c.replace("chr", "") for c in co["chrom"]]
        # draw as tick labels
        ax.set_xticks(chromcenters)
        ax.set_xticklabels(labels, **axis_kwargs)
        # ticks: short marks between chroms
        ax.set_xticks(np.concatenate([co["stop"].values[:-1] + space * co["size"].sum() / 2,
                                       []]), minor=True)
        return

    # single-chrom mode
    inner_start = chromstart + edgeblankfraction * rng
    inner_end = chromend - edgeblankfraction * rng

    fact = {"bp": 1, "Kb": 1000, "Mb": 1000000}.get(scale, 1)
    inner_start_lbl = inner_start / fact
    inner_end_lbl = inner_end / fact

    labels = np.array(_pretty([inner_start_lbl, inner_end_lbl], n=n))[1:-1]
    tick_positions = np.concatenate([[xmin], labels * fact, [xmax]])

    # matplotlib's axis() handles scientific notation via ScalarFormatter; we
    # use a FuncFormatter that prints plain integers.
    from matplotlib.ticker import FuncFormatter, FixedLocator, FixedFormatter

    def _fmt(v, pos):
        try:
            return f"{v / fact:g}"
        except Exception:
            return ""

    ax.xaxis.set_major_locator(FixedLocator(tick_positions))
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt))
    ax.tick_params(axis="x", which="major", direction="out", length=3,
                   labelsize=8 * scalecex, pad=18 + 10 * line)  # bumped from 4 to leave room for chrom/scale labels

    # chrom name + scale name
    # matplotlib: tick_params already moved the tick labels below the axis.
    # Place chrom and scale labels just below the tick labels.
    # Place the chrom/scale labels just below the axis. Note that to
    # render outside the axes bbox the caller must pass
    # ``fig.subplots_adjust(bottom=0.18)`` (or similar) before savefig.
    # In R, labelgenome uses axis() with `line=0.18` to place tick labels,
    # then mtext places the chromosome name and scale unit with their own
    # `line=` arguments.  In matplotlib, the cleanest analogue is to
    # use figure.text() anchored in figure-fraction coords, with the
    # caller responsible for having already set up the bottom/top margins
    # via subplots_adjust.  We position the chrom name at the LEFT-bottom
    # corner and the scale unit at the RIGHT-bottom corner, with a small
    # offset below the axis tick labels.
    # In practice, since the figure has subplots_adjust already, we use
    # ax.annotate with xycoords='axes fraction' and y close to 0 but
    # shifted enough to be BELOW the tick labels (which sit at y=0).
    # We rely on annotation_clip=False so the text can extend below the
    # axes bounding box.
    if side == 1:
        # Place labels just below the axis (R default line=.5)
        y_label = -0.16
        y_scale = -0.16
    else:
        y_label = 1.05 + 0.025 * chromline
        y_scale = 1.05 + 0.025 * scaleline

    weight = ["normal", "bold", "italic", "bold italic"][min(chromfont, 4) - 1]
    weight_s = ["normal", "bold", "italic", "bold italic"][min(scalefont, 4) - 1]

    ax.annotate(chrom.replace("chr", ""),
                xy=(chromadjust, y_label), xycoords="axes fraction",
                ha="left", va="top",
                fontweight=weight, fontsize=8 * chromcex,
                annotation_clip=False)
    ax.annotate(scale,
                xy=(scaleadjust, y_scale), xycoords="axes fraction",
                ha="right", va="top",
                fontweight=weight_s, fontsize=8 * scalecex,
                annotation_clip=False)


def _pretty(x, n=5):
    """Mimic R's pretty() for axis ticks."""
    a, b = min(x), max(x)
    if a == b:
        return np.array([a])
    # numpy's linspace is close enough to R pretty()
    return np.linspace(a, b, n + 1)


# ---------------------------------------------------------------------------
# addlegend
# ---------------------------------------------------------------------------
def addlegend(ax, range_val=None, title: str = "",
              palette: Callable = None,
              labels_digits: int = 1,
              side: str = "right", labelside: str = "left",
              xoffset: float = 0.1, width: float = 0.05,
              bottominset: float = 0.025, topinset: float = 0.025,
              tick_num: int = 5, tick_length: float = 0.01,
              txt_font: int = 1, txt_cex: float = 0.75,
              title_offset: float = 0.05,
              title_font: int = 2, title_cex: float = 1.0):
    """Draw a color-scaled legend on the right/left of an Axes.

    R signature:
        addlegend(range, title, palette, ...)

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to draw the legend on.
    range_val : (float, float)
        (min, max) values the colors map to.
    palette : callable
        A palette function like ``SushiColors(7)``.
    """
    if palette is None:
        palette = SushiColors(7)

    # 99 stops between min and max
    scalecol = maptocolors(np.linspace(range_val[0], range_val[1], 99), palette,
                            num=100, rng=range_val)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xrange = abs(xmax - xmin)
    yrange = abs(ymax - ymin)

    if side == "right":
        xleft = xmax + xrange * xoffset
        xright = xleft + xrange * width
    else:  # "left"
        xright = xmin - xrange * xoffset
        xleft = xright - xrange * width

    rect_edges = np.linspace(ymin + yrange * bottominset,
                              ymax - yrange * topinset,
                              num=len(scalecol) + 1)
    rect_tops = rect_edges[1:]
    rect_bots = rect_edges[:-1]

    # draw the color bars
    for i, color in enumerate(scalecol):
        ax.add_patch(plt.Rectangle((xleft, rect_bots[i]),
                                    xright - xleft,
                                    rect_tops[i] - rect_bots[i],
                                    facecolor=color, edgecolor="none",
                                    clip_on=False, zorder=10))

    # border
    ax.add_patch(plt.Rectangle((xleft, rect_bots[0]),
                                xright - xleft,
                                rect_tops[-1] - rect_bots[0],
                                facecolor="none", edgecolor="black",
                                clip_on=False, zorder=11))

    # ticks + labels
    if labelside == "left":
        tick_left = xleft - tick_length * xrange
        tick_right = xleft
        txt_x = xleft - tick_length * xrange
        txt_adj = (1.2, 0.5)
        title_x = xleft - tick_length * xrange - title_offset * xrange
        title_rot = 90
    else:  # "right"
        tick_left = xright
        tick_right = xright + tick_length * xrange
        txt_x = xright + tick_length * xrange
        txt_adj = (-0.2, 0.5)
        title_x = xright + tick_length * xrange + title_offset * xrange
        title_rot = 270

    # tick label positions
    tick_yvals = np.linspace(rect_edges[0], rect_edges[-1], num=tick_num)
    tick_labels = np.round(np.linspace(range_val[0], range_val[1], num=tick_num),
                            labels_digits)

    for y, lab in zip(tick_yvals, tick_labels):
        ax.plot([tick_left, tick_right], [y, y], color="black",
                clip_on=False, zorder=11)
        ax.annotate(f"{lab:g}", xy=(txt_x, y), xycoords="data",
                    ha="right" if labelside == "left" else "left",
                    va="center",
                    fontsize=10 * txt_cex,
                    fontweight=["normal", "bold", "italic", "bold italic"][min(txt_font, 4) - 1],
                    annotation_clip=False)

    # title
    if title:
        mid_y = (rect_edges[0] + rect_edges[-1]) / 2.0
        ax.annotate(title, xy=(title_x, mid_y),
                    xycoords="data", ha="center", va="center",
                    rotation=title_rot,
                    fontweight=["normal", "bold", "italic", "bold italic"][min(title_font, 4) - 1],
                    fontsize=10 * title_cex,
                    annotation_clip=False)


# ---------------------------------------------------------------------------
# zoombox
# ---------------------------------------------------------------------------
def zoombox(ax, zoomregion=None,
            lty: str = "--", lwd: float = 1.0, col: str = "black",
            topextend: float = 2.0, passthrough: bool = False):
    """Draw a zoom indicator above the axes (R zoombox).

    Use on the *second* plot of a zoom-in (the in-zoomed panel).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The zoomed-in axes.
    zoomregion : (start, stop) or None
        If a region is provided and ``passthrough`` is False, the top of the
        box has a gap between the in-region and the panel edges.
    topextend : float
        Fraction of the y-range by which the vertical lines extend upward
        beyond the panel top.
    passthrough : bool
        If True the vertical lines extend downward beyond the bottom too.
    """
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    yrange = abs(ymax - ymin)
    extend_top = topextend * yrange

    if passthrough:
        ax.plot([xmin, xmin], [ymin - extend_top, extend_top], linestyle=lty,
                linewidth=lwd, color=col, clip_on=False)
        ax.plot([xmax, xmax], [ymin - extend_top, extend_top], linestyle=lty,
                linewidth=lwd, color=col, clip_on=False)
    else:
        ax.plot([xmin, xmin], [ymax, extend_top], linestyle=lty,
                linewidth=lwd, color=col, clip_on=False)
        ax.plot([xmax, xmax], [ymax, extend_top], linestyle=lty,
                linewidth=lwd, color=col, clip_on=False)
        if zoomregion is None:
            ax.plot([xmin, xmax], [ymax, ymax], linestyle=lty,
                    linewidth=lwd, color=col, clip_on=False)
        else:
            zmin, zmax = zoomregion
            ax.plot([xmin, zmin], [ymax, ymax], linestyle=lty,
                    linewidth=lwd, color=col, clip_on=False)
            ax.plot([zmax, xmax], [ymax, ymax], linestyle=lty,
                    linewidth=lwd, color=col, clip_on=False)


# ---------------------------------------------------------------------------
# zoomsregion
# ---------------------------------------------------------------------------
def zoomsregion(ax, region: Sequence[float],
                chrom: Optional[str] = None,
                genome=None, space: float = 0.01,
                padding: float = 0.005,
                col: str = "white", zoomborder: str = "black",
                lty: str = "--", lwd: float = 1.0,
                extend=0.0, wideextend: float = 0.1,
                offsets=(0, 0), highlight: bool = False):
    """Draw a trapezoid connector from a parent panel to a zoomed panel.

    Use on the *first* plot of a zoom-in (the wide-region panel).

    R signature:
        zoomsregion(region, ...)

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The wide-region axes.
    region : (start, stop)
        The region that the zoom-in covers. If too narrow (< 2*padding of
        the panel width) it is widened symmetrically.
    chrom : str or None
        Required when genome is provided (multi-chrom mode).
    genome : DataFrame or None
        If provided, region coordinates are translated through chromOffsets.
    extend : float or (float, float)
        Vertical extension above and below the panel as a fraction of yrange.
    """
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xrange = abs(xmax - xmin)
    yrange = abs(ymax - ymin)

    if isinstance(extend, (int, float)):
        ext_up = float(extend) * yrange
        ext_lo = float(extend) * yrange
    else:
        ext_up, ext_lo = extend[0] * yrange, extend[1] * yrange

    plot_left = xmin + xrange * offsets[0]
    plot_right = xmax - xrange * offsets[1]

    # adjust region to multi-chrom offset
    if genome is not None:
        co = chromOffsets(genome, space)
        offset_row = co.loc[co["chrom"] == chrom]
        if offset_row.empty:
            raise ValueError(f"chrom {chrom!r} not found in genome")
        off = float(offset_row["start"].iloc[0])
        region = [region[0] + off, region[1] + off]
        # padding scaled to total layout width
        if abs(region[1] - region[0]) < 2 * padding * float(co["stop"].max()):
            mid = (region[0] + region[1]) / 2.0
            region = [mid - padding * float(co["stop"].max()),
                      mid + padding * float(co["stop"].max())]

    # add padding in single-chrom mode
    if (abs(region[1] - region[0]) / xrange) < (2 * padding):
        mid = (region[0] + region[1]) / 2.0
        region = [mid - padding * xrange, mid + padding * xrange]

    if highlight:
        plot_right = region[1]
        plot_left = region[0]
        min_y = ymin - ext_lo
        low_y = ymin - ext_lo
    else:
        min_y = ymin - ext_lo - (wideextend * yrange)
        low_y = ymin - ext_lo - (5 * yrange)

    # build the polygon (matching R's vertex order)
    poly_x = [
        region[0], region[0], region[1], region[1],
        plot_right, plot_right, plot_left, plot_left,
        region[0],
    ]
    poly_y = [
        ymin - ext_lo, ymax + ext_up, ymax + ext_up, ymin - ext_lo,
        min_y, low_y, low_y, min_y, ymin - ext_lo,
    ]
    ax.add_patch(Polygon(list(zip(poly_x, poly_y)),
                          closed=True, facecolor=col, edgecolor=zoomborder,
                          linestyle=lty, linewidth=lwd, fill=(col != "white" and col is not None),
                          clip_on=False))


__all__ = ["labelgenome", "addlegend", "zoombox", "zoomsregion"]
