"""Example 3: flip a bedGraph track (mirror axis).
Output: /tmp/flip.png
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
             transparency=0.5, flip=True, color=SushiColors(2)(2)[1])
labelgenome(ax, "chr11", 1955000, 1960000, side=3, n=3, scale="Kb")
fig.savefig("/tmp/flip.png", dpi=120, bbox_inches="tight")
