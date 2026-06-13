"""
Sushi_Python_Version -- core helpers (no matplotlib side effects)
================================================================

Ported from Sushi R/Bioconductor package (Phanstiel et al., 2015).
https://github.com/PhanstielLab/Sushi

This module contains pure helpers:
    convertstrandinfo, sortChrom, chromOffsets,
    SushiColors, maptocolors, maptolwd, opaque,
    labelplot, _checkrow, _downsample_bedgraph.

Conventions:
    * Chromosome names are Python strings.
    * All genomic coordinates are 1-based inclusive (BED convention).
    * Plotting helpers return (Axes, ...) so callers can layer labels,
      legends, zooms.
"""

from __future__ import annotations

from typing import Callable, Iterable, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# convertstrandinfo
# ---------------------------------------------------------------------------
def convertstrandinfo(strandvector):
    """Convert '+/-' strand info to 1/-1 numeric (R semantics)."""
    out = list(strandvector)
    if len(out) == 0:
        return out
    if isinstance(out[0], str):
        return [(-1 if s == "-" else 1) for s in out]
    return out


# ---------------------------------------------------------------------------
# sortChrom
# ---------------------------------------------------------------------------
def sortChrom(genome):
    """Sort chroms as chr1, chr2, ..., chrX, chrY, chrM (numeric first)."""
    if not isinstance(genome, pd.DataFrame):
        genome = pd.DataFrame(list(genome), columns=["chrom", "size"])

    df = genome.copy()
    df["_bare"] = df["chrom"].astype(str).str.replace("chr", "", regex=False)

    numeric = pd.to_numeric(df["_bare"], errors="coerce")
    is_num = numeric.notna()
    num_df = df.loc[is_num].assign(_n=numeric[is_num]).sort_values("_n")
    str_df = df.loc[~is_num].sort_values("_bare")

    out = pd.concat([num_df, str_df], ignore_index=True)
    out["chrom"] = "chr" + out["_bare"].astype(str)
    return out[["chrom", "size"]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# chromOffsets
# ---------------------------------------------------------------------------
def chromOffsets(genome, space: float = 0.01) -> pd.DataFrame:
    """Per-chrom start/stop offsets for multi-chrom Manhattan plots."""
    sorted_df = sortChrom(genome).reset_index(drop=True)
    sizes = sorted_df["size"].astype(float).values
    cumsums = np.cumsum(sizes)
    total = cumsums[-1]
    spacer = total * space
    n = len(cumsums)
    additional = np.arange(1, n + 1) * spacer
    startpos = np.concatenate([[0.0], cumsums[:-1]]) + additional
    stoppos = cumsums + np.arange(1, n + 1) * spacer
    return pd.DataFrame({
        "chrom": sorted_df["chrom"].values,
        "size": sizes,
        "start": startpos,
        "stop": stoppos,
    })


# ---------------------------------------------------------------------------
# SushiColors -- discrete palette factory
# ---------------------------------------------------------------------------
_PALETTES = {
    2: ["blue", "#E5001B"],
    3: ["blue", "#E5001B", "orange"],
    4: ["blue", "#5900E5", "#E5001B", "orange"],
    5: ["blue", "#5900E5", "#E5001B", "orange", "yellow"],
    6: ["blue", "#5900E5", "#E5001B", "orange", "yellow"],
    7: ["black", "blue", "#5900E5", "#E5001B", "orange", "yellow", "white"],
}


class _SushiPaletteFactory:
    """Callable f(n) -> list of n hex colors from the Sushi discrete palette."""

    def __init__(self, name: str, stops: Sequence[str]):
        self.name = name
        self.stops = list(stops)
        self._cache = {}

    def __call__(self, n: int) -> list:
        if n == len(self.stops):
            return list(self.stops)
        if n in self._cache:
            return list(self._cache[n])
        idx = np.linspace(0, len(self.stops) - 1, n)
        if np.all(idx == idx.astype(int)):
            out = [self.stops[int(i)] for i in idx]
        else:
            out = [self.stops[int(round(i))] for i in idx]
        self._cache[n] = out
        return list(out)

    def __repr__(self) -> str:
        return f"SushiColors({self.name!r})"


def SushiColors(palette: Union[int, str] = "fire") -> Optional[Callable]:
    """Return a Sushi palette function (callable taking N -> list of N hex colors).

    R signature:
        SushiColors(palette = "fire")

    Parameters
    ----------
    palette : int or str
        * int in {2,3,4,5,6,7}  -- discrete Sushi palette
        * "list" -- return None (caller inspects ``_PALETTES``)
        * other  -- 2-color default (blue, #E5001B)
    """
    if palette == "list":
        return None
    if palette in _PALETTES:
        return _SushiPaletteFactory(name=str(palette), stops=_PALETTES[palette])
    return _SushiPaletteFactory(name="2", stops=_PALETTES[2])


# ---------------------------------------------------------------------------
# maptocolors
# ---------------------------------------------------------------------------
def maptocolors(vec, col: Callable, num: int = 100,
                range: Optional[Sequence[float]] = None) -> list:
    """Map a numeric vector to a discrete color palette (R semantics)."""
    arr = np.asarray(vec, dtype=float)
    if range is not None:
        arr = np.where(arr < range[0], range[0], arr)
        arr = np.where(arr > range[1], range[1], arr)
        lo, hi = range
    else:
        lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
    if lo == hi:
        return [col(1)[0]] * len(arr)
    breaks = np.linspace(lo, hi, num=num)
    palette = col(num + 1)
    edges = np.concatenate([[-np.inf], breaks, [np.inf]])
    idx = np.digitize(arr, edges) - 1
    idx = np.clip(idx, 0, len(palette) - 1)
    return [palette[int(i)] for i in idx]


# ---------------------------------------------------------------------------
# maptolwd
# ---------------------------------------------------------------------------
def maptolwd(lwdby, range=(1, 5)):
    """Map a numeric vector to line widths in `range`."""
    arr = np.asarray(lwdby, dtype=float)
    arr = arr - np.nanmin(arr)
    denom = np.nanmax(arr)
    if denom == 0:
        return [float(range[0])] * len(arr)
    frac = arr / denom
    widths = np.round(frac * abs(range[1] - range[0]), 1) + range[0]
    return widths.tolist()


# ---------------------------------------------------------------------------
# opaque
# ---------------------------------------------------------------------------
def opaque(color, transparency: float = 0.5):
    """Apply alpha to one or many colors (1 = opaque, 0 = transparent)."""
    import matplotlib.colors as mcolors

    def _one(c: str) -> str:
        rgba = mcolors.to_rgba(c)
        # mcolors.to_hex drops alpha; use to_rgba -> hex with keep_alpha=True
        return mcolors.to_hex((rgba[0], rgba[1], rgba[2], float(transparency)), keep_alpha=True)

    if isinstance(color, str):
        return _one(color)
    return [_one(c) for c in color]


# ---------------------------------------------------------------------------
# labelplot
# ---------------------------------------------------------------------------
def labelplot(ax, letter=None, title=None,
              letteradj=-0.05, titleadj=0.0,
              letterfont=2, titlefont=2,
              lettercex=1.2, titlecex=1.0,
              letterline=0.5, titleline=0.5,
              lettercol="black", titlecol="black"):
    """Add a paper letter + title to the top of an axes (R labelplot)."""
    weights = {1: "normal", 2: "bold", 3: "italic", 4: "bold italic"}
    if letter is not None:
        ax.text(letteradj, 1 + letterline / 10.0, letter,
                transform=ax.transAxes,
                ha="left", va="bottom",
                fontweight=weights.get(letterfont, "normal"),
                fontsize=10 * lettercex, color=lettercol)
    if title is not None:
        ax.text(titleadj, 1 + titleline / 10.0, title,
                transform=ax.transAxes,
                ha="left", va="bottom",
                fontweight=weights.get(titlefont, "normal"),
                fontsize=10 * titlecex, color=titlecol)


# ---------------------------------------------------------------------------
# _checkrow -- auto row packing
# ---------------------------------------------------------------------------
def _checkrow(data: Sequence, alldata: pd.DataFrame, maxrows: int,
              wiggle: float = 0, start_col: str = "start",
              stop_col: str = "stop") -> Optional[int]:
    """Return the lowest row index whose range does not overlap `data`."""
    thestart = float(min(data[0], data[1])) - wiggle
    thestop = float(max(data[0], data[1])) + wiggle
    for row in range(1, maxrows + 1):
        occupied = alldata.loc[alldata["plotrow"] == row]
        if occupied.empty:
            return row
        overlap = (
            (thestart >= occupied[start_col]) & (thestart <= occupied[stop_col])
        ) | (
            (thestop >= occupied[start_col]) & (thestop <= occupied[stop_col])
        ) | (
            (thestart <= occupied[start_col]) & (thestop >= occupied[stop_col])
        )
        if not overlap.any():
            return row
    return None


# ---------------------------------------------------------------------------
# _downsample_bedgraph
# ---------------------------------------------------------------------------
def _downsample_bedgraph(signaltrack: np.ndarray,
                         max_points: int = 8000) -> np.ndarray:
    """Halve signaltrack row-by-row until <= max_points (R while-loop)."""
    while signaltrack.shape[0] > max_points:
        if signaltrack.shape[0] % 2 != 0:
            signaltrack = signaltrack[:-1]
        starts = signaltrack[0::2, 0]
        stops = signaltrack[1::2, 1]
        values = (signaltrack[0::2, 2] + signaltrack[1::2, 2]) / 2.0
        signaltrack = np.stack([starts, stops, values], axis=1)
    return signaltrack


def _bedgraph_step_polygon(signaltrack: np.ndarray,
                            baseline: float = 0.0) -> tuple:
    """Build a closed step-function polygon from a bedGraph-like (n, 3) array.

    The polygon is the canonical "step plot" shape used by R Sushi's plotBedgraph:
    every bin spans its full [start, stop] at constant value, and adjacent
    bins share a vertical wall at the join. The polygon closes back to
    baseline on the left (start[0]-1) and right (stop[-1]+1).

    Returns
    -------
    (poly_x, poly_y) : np.ndarray
        1-D arrays of equal length, suitable for matplotlib ``ax.fill`` /
        ``ax.plot``. The outline follows the top edge then returns to
        baseline; ``ax.plot(poly_x, poly_y)`` draws the step outline.

    Why we DON'T use np.concatenate([starts, stops_reversed]) + matching
    reverse values: that pattern produces a single Z-folded polygon with
    a long diagonal from the last bin's top to the first bin's top, which
    fills the wrong area when adjacent bin values differ. The step polygon
    below avoids that diagonal by repeating each bin's start->stop edge.
    """
    if signaltrack.shape[0] == 0:
        return np.zeros(0), np.zeros(0)

    starts = signaltrack[:, 0].astype(float)
    stops = signaltrack[:, 1].astype(float)
    values = signaltrack[:, 2].astype(float)

    # Build top-edge as a step function:
    #   x: start[0], start[1], ..., start[n-1], stop[n-1]
    #   y: value[0], value[1], ..., value[n-1], value[n-1]
    top_x = np.concatenate([starts, [stops[-1]]])
    top_y = np.concatenate([values, [values[-1]]])

    # Close polygon: descend back to baseline in *reverse* order along
    # the same step function (so the right wall of bin k is on x=stop[k]).
    bottom_x = top_x[::-1]
    bottom_y = np.concatenate([top_y[::-1][1:], [top_y[::-1][0]]])
    # simpler: bottom_x = top_x[::-1], bottom_y = baseline * same length
    bottom_x = top_x[::-1]
    bottom_y = np.full_like(bottom_x, baseline)

    poly_x = np.concatenate([
        [starts[0] - 1.0],
        top_x,
        [stops[-1] + 1.0],
        [stops[-1] + 1.0],
        bottom_x,
        [starts[0] - 1.0],
    ])
    poly_y = np.concatenate([
        [baseline],
        top_y,
        [values[-1]],
        [baseline],
        bottom_y,
        [baseline],
    ])
    return poly_x, poly_y


__all__ = [
    "convertstrandinfo",
    "sortChrom",
    "chromOffsets",
    "SushiColors",
    "maptocolors",
    "maptolwd",
    "opaque",
    "labelplot",
    "_checkrow",
    "_downsample_bedgraph",
    "_bedgraph_step_polygon",
]
