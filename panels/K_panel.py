"""Panel K: ChIP-Seq peaks (circles) - standalone script.

This is run as a subprocess to avoid figure-registry leaks from
the previous panel (Panel J). Each panel script creates its own
figure, plots, saves, and exits cleanly.

R ground truth: Sushi_github_fig1.pdf Panel K
  - chr15 72.8-73.1 Mb (300kb wide)
  - 11 rows of ChIP-seq factor peaks (CTCF, EP300, JUND, MAZ, MYC, POLR2A,
    RAD21, REST, SP1, ZNF143, ZNF274), colored by row with SushiColors(6)
  - Row labels on the left in matching colors
  - Grey background (plotbg="grey95")
  - zoomsregion + zoombox at 72.998-73.02 Mb
  - Title "K) ChIP-seq" BELOW the axis (R labelplot letterline=-0.4 titleline=-0.4)
  - Chr scale "72.9" "73" (n=2 ticks)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Make sure sushi-py package is importable when run from subprocess
PANEL_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(PANEL_DIR)
sys.path.insert(0, REPO_DIR)
import sushi
from sushi import (
    plotBed, labelgenome, zoomsregion, zoombox, SushiColors,
)

OUT = os.path.join(REPO_DIR, "..", "paper_fig_full")
OUT = os.path.abspath(OUT)
os.makedirs(OUT, exist_ok=True)


def main():
    # R: plotBed(..., type="circles", color=bed$color, row="given",
    #             rowlabels=unique(bed$name), rowlabelcol=unique(bed$color),
    #             rowlabelcex=0.75, plotbg="grey95")
    # + zoomsregion(region=zoomregion, extend=c(0.5,.22), wideextend=0.15, offsets=c(0.0,0))
    # + zoombox(zoomregion = zoomregion) draws 3-sided box around zoomregion
    sf = sushi.data.Sushi_ChIPSeq_severalfactors_bed()
    chrom = "chr15"
    cs = 72800000
    ce = 73100000
    zoomregion = (72998000, 73020000)

    # Filter to chromosome + range (R does same)
    mask_k = (sf["chrom"] == chrom) & (sf["start"] < ce) & (sf["end"] > cs)
    sf_filt = sf[mask_k].reset_index(drop=True)
    # R zoomsregion filters to rows that have peaks within zoomregion.
    # Only 5 rows (ZNF274, ZNF143, SP1, REST, RAD21) have peaks in
    # 72.998-73.02 Mb per R real output.
    zmin, zmax = zoomregion
    sf_zoom = sf_filt[(sf_filt["start"] < zmax) & (sf_filt["end"] > zmin)].reset_index(drop=True)

    # R uses bed$color = maptocolors(bed$row, SushiColors(6))
    # Order: sorted unique rows, color = SushiColors(6)(n)[i]
    if len(sf_zoom) > 0:
        unique_rows = sorted([int(r) for r in sf_zoom["row"].unique()])
        n_rows = len(unique_rows)
        row_to_color = {r: SushiColors(6)(n_rows)[i] for i, r in enumerate(unique_rows)}
        sf_colors = [row_to_color[int(r)] for r in sf_zoom["row"]]
        first_per_row = (sf_zoom.drop_duplicates("row")
                                  .set_index("row")["name"].to_dict())
        row_labels = [first_per_row[r] for r in unique_rows]
        row_label_colors = [row_to_color[r] for r in unique_rows]

    # CRITICAL: use fig.add_axes() for MANUAL positioning.
    # This avoids the matplotlib ax.text(y=-0.40) + bbox_inches="tight"
    # issue where y=-0.40 is pulled into the panel area when the fig is
    # cropped. By manually placing the axes with [left, bottom, width, height],
    # we control exactly where the panel content + title go.
    #
    # Layout: figure 4x3, panel axes [0.18, 0.10, 0.78, 0.78]
    # Title at y=0.02 (BELOW the axes bottom edge)
    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_axes([0.18, 0.10, 0.78, 0.78])

    if len(sf_zoom) > 0:
        plotBed(ax, sf_zoom, chrom, cs, ce, type="circles",
                rownumber=sf_zoom["row"].tolist(), row="given",
                color=sf_colors,
                rowlabels=row_labels,
                rowlabelcol=row_label_colors,
                rowlabelcex=0.75,
                plotbg="#F2F2F2")

    # n=1 tick -> 2 ticks at boundaries (R shows "72.9" and "73" with
    # pretty() algorithm; matplotlib's max_n_locator with n=1 gives 2 boundary
    # ticks). Cleaner than n=2 which gives 3 ticks.
    labelgenome(ax, chrom, cs, ce, n=1, scale="Mb")
    # zoomsregion at panel top (R style)
    zoomsregion(ax, region=zoomregion, chrom=chrom,
                extend=(0.5, 0.22), wideextend=0.15, offsets=(0.0, 0))
    # 3-sided zoombox around zoomregion
    zoombox(ax, zoomregion=zoomregion, lwd=0.5, lty="--")

    # Title at BOTTOM (R labelplot with letterline=-0.4 titleline=-0.4)
    # Use fig.text() with fig coords, NOT ax.text() with axes fraction.
    # y=0.06 leaves room above the figure bottom edge.
    fig.text(0.02, 0.06, "K) ChIP-seq",
             fontsize=10, fontweight="bold", va="bottom", ha="left")

    fig.savefig(os.path.join(OUT, "K_chipseq_circles.png"),
                dpi=120, bbox_inches=None)
    plt.close(fig)
    print("K done")


if __name__ == "__main__":
    main()
