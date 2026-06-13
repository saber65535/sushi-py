"""Example 15: zoom-in figure with 1 wide + 2 narrow panels.
Output: /tmp/zoom.png
"""
import matplotlib.pyplot as plt
import sushi
from sushi import (
    plotBedgraph, labelgenome, zoombox, zoomsregion, SushiColors,
)

bg = sushi.data.Sushi_DNaseI_bedgraph()
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
fig.savefig("/tmp/zoom.png", dpi=120, bbox_inches="tight")
