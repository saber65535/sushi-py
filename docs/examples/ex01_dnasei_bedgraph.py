"""Example 1: plotBedgraph with gradient fill.
Output: /tmp/dnasei.png
"""
import matplotlib.pyplot as plt
import sushi
from sushi import plotBedgraph, labelgenome, SushiColors

bg = sushi.data.Sushi_DNaseI_bedgraph()

fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedgraph(ax, bg, "chr11", 1650000, 2350000, colorbycol=SushiColors(5))
labelgenome(ax, "chr11", 1650000, 2350000, n=4, scale="Mb", side=1)
fig.savefig("/tmp/dnasei.png", dpi=120, bbox_inches="tight")
