"""Example 6: BEDPE interactions as filled ribbons.
Output: /tmp/ribbons.png
"""
import matplotlib.pyplot as plt
import sushi
from sushi import plotBedpe, labelgenome, SushiColors

bedpe = sushi.data.Sushi_5C_bedpe()
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedpe(ax, bedpe, "chr11", 1650000, 2350000,
          heights=bedpe["score"].tolist(), plottype="ribbons",
          colorby=bedpe["samplenumber"].tolist(),
          colorbycol=SushiColors(3), border="black")
labelgenome(ax, "chr11", 1650000, 2350000, n=3, scale="Mb", side=1)
fig.savefig("/tmp/ribbons.png", dpi=120, bbox_inches="tight")
