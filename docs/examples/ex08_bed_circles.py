"""Example 8: BED peaks as one-circle-per-row (WashU style).
Output: /tmp/circles.png
"""
import matplotlib.pyplot as plt
import sushi
from sushi import plotBed, labelgenome, SushiColors

bed_full = sushi.data.Sushi_ChIPSeq_pol2_bed()
bed = bed_full[(bed_full["chrom"] == "chr11") &
               ((bed_full["start"] > 2281200) & (bed_full["start"] < 2282200) |
                (bed_full["end"] > 2281200) & (bed_full["end"] < 2282200))].reset_index(drop=True)
cb = bed["strand"].map({"+": 1, "-": -1}).tolist()

fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22, left=0.12)
plotBed(ax, bed, "chr11", 2281200, 2282200, type="circles",
        colorby=cb, colorbycol=SushiColors(2), row="auto", wiggle=0.001)
labelgenome(ax, "chr11", 2281200, 2282200, n=2, scale="Kb", side=1)
fig.savefig("/tmp/circles.png", dpi=120, bbox_inches="tight")
