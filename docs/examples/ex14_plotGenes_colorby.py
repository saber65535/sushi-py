"""Example 14: transcript models, box style, colored by log10(FPKM).
Output: /tmp/transcripts.png
"""
import matplotlib.pyplot as plt
import numpy as np
import sushi
from sushi import plotGenes, labelgenome, addlegend, SushiColors

trans = sushi.data.Sushi_transcripts_bed()
fig, ax = plt.subplots(figsize=(10, 5))
fig.subplots_adjust(bottom=0.22, right=0.82)
pg = plotGenes(ax, trans, "chr15", 72965000, 72990000,
               types=trans["type"].tolist(),
               colorby=np.log10(trans["score"].values + 0.001).tolist(),
               colorbycol=SushiColors(5), colorbyrange=(0, 1.0),
               labeltext=True, maxrows=50, bheight=0.4,
               plotgenetype="box")
labelgenome(ax, "chr15", 72965000, 72990000, n=3, scale="Mb", side=1)
addlegend(ax, range_val=pg[0], palette=pg[1], title="log10(FPKM)",
          side="right", bottominset=0.4, topinset=0, xoffset=-0.035,
          labelside="left", width=0.025, title_offset=0.055)
fig.savefig("/tmp/transcripts.png", dpi=120, bbox_inches="tight")
