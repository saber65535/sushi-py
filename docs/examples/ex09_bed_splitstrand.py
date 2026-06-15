"""Example 9: BED peaks split by strand (+ above axis, - below).
Output: /tmp/splitstrand.png
"""
import matplotlib.pyplot as plt
import sushi
from sushi import plotBed, labelgenome, SushiColors

bed_full = sushi.data.Sushi_ChIPSeq_pol2_bed()
bed = bed_full[(bed_full["chrom"] == "chr11") &
               ((bed_full["start"] > 2281200) & (bed_full["start"] < 2282200) |
                (bed_full["end"] > 2281200) & (bed_full["end"] < 2282200))].reset_index(drop=True)
# strand is stored numeric (-1.0/1.0) in the example data; convertstrandinfo
# returns -1/1 for both numeric and "+"/"-" string input.
cb = sushi.convertstrandinfo(bed["strand"].tolist())

fig, ax = plt.subplots(figsize=(10, 5))
fig.subplots_adjust(bottom=0.22, left=0.12)
plotBed(ax, bed, "chr11", 2281200, 2282200, colorby=cb,
        colorbycol=SushiColors(2), row="auto", splitstrand=True, wiggle=0.001)
labelgenome(ax, "chr11", 2281200, 2282200, n=2, scale="Kb", side=1)
fig.savefig("/tmp/splitstrand.png", dpi=120, bbox_inches="tight")
