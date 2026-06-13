"""Example 2: overlay two bedGraph tracks (CTCF + DNaseI).
Output: /tmp/overlay.png
"""
import matplotlib.pyplot as plt
import sushi
from sushi import plotBedgraph, labelgenome, SushiColors

fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotBedgraph(ax, sushi.data.Sushi_ChIPSeq_CTCF_bedgraph(),
             "chr11", 1955000, 1960000,
             transparency=0.5, color=SushiColors(2)(2)[0])
plotBedgraph(ax, sushi.data.Sushi_DNaseI_bedgraph(),
             "chr11", 1955000, 1960000,
             transparency=0.5, color=SushiColors(2)(2)[1],
             overlay=True, rescaleoverlay=True)
labelgenome(ax, "chr11", 1955000, 1960000, n=3, scale="Kb", side=1)
fig.savefig("/tmp/overlay.png", dpi=120, bbox_inches="tight")
