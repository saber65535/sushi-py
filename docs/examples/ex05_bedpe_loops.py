"""Example 5: BEDPE interactions as parabolic loops.
Output: /tmp/loops.png
"""
import matplotlib.pyplot as plt
import sushi
from sushi import plotBedpe, labelgenome, SushiColors
from matplotlib.lines import Line2D

bedpe = sushi.data.Sushi_5C_bedpe()
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedpe(ax, bedpe, "chr11", 1650000, 2350000,
          heights=bedpe["score"].tolist(), plottype="loops",
          colorby=bedpe["samplenumber"].tolist(),
          colorbycol=SushiColors(3))
labelgenome(ax, "chr11", 1650000, 2350000, n=3, scale="Mb", side=1)
handles = [Line2D([0], [0], marker="o", linestyle="",
                   markerfacecolor=c, markeredgecolor=c, label=l)
           for c, l in zip(SushiColors(3)(3),
                            ["K562", "HeLa", "GM12878"])]
ax.legend(handles=handles, loc="upper right", frameon=False)
fig.savefig("/tmp/loops.png", dpi=120, bbox_inches="tight")
