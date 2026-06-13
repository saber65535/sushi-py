"""Example 10: multiple ChIP-seq factors on one axis, circles.
Output: /tmp/several_circles.png
"""
import matplotlib.pyplot as plt
import sushi
from sushi import plotBed, labelgenome, SushiColors, maptocolors

bed_full = sushi.data.Sushi_ChIPSeq_severalfactors_bed()
bed = bed_full[(bed_full["chrom"] == "chr15") &
               ((bed_full["start"] > 72800000) & (bed_full["start"] < 73100000) |
                (bed_full["end"] > 72800000) & (bed_full["end"] < 73100000))].reset_index(drop=True)
bed["color"] = maptocolors(bed["row"].tolist(), col=SushiColors(6))

fig, ax = plt.subplots(figsize=(10, 5))
fig.subplots_adjust(bottom=0.22)
plotBed(ax, bed, "chr15", 72800000, 73100000, type="circles",
        color=bed["color"].tolist(),
        rownumber=bed["row"].tolist(), row="given",
        plotbg="#F2F2F2",
        rowlabels=bed_full["name"].drop_duplicates().tolist()[:11],
        rowlabelcex=0.7)
labelgenome(ax, "chr15", 72800000, 73100000, n=3, scale="Mb", side=1)
fig.savefig("/tmp/several_circles.png", dpi=120, bbox_inches="tight")
