"""R Ex3 R ↔ Python 1:1 reproduction.

R real code (Sushi.Rnw line 245, Sushi-010.pdf):
    chrom = "chr11"; chromstart = 1650000; chromend = 2350000
    plotBedgraph(Sushi_DNaseI.bedgraph, chrom, chromstart, chromend, colorbycol=SushiColors(5))
    labelgenome(chrom, chromstart, chromend, n=4, scale="Mb")
    mtext("Read Depth", side=2, line=1.75, cex=1, font=2)
    axis(side=2, las=2, tcl=.2)
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
REPO_DIR = r"D:\奇奇怪怪的图\Sushi_Python_Version"
sys.path.insert(0, REPO_DIR)
import sushi
from sushi import plotBedgraph, labelgenome, SushiColors
OUT = r"D:\奇奇怪怪的图\Sushi_Python_Version\examples_R_real"
os.makedirs(OUT, exist_ok=True)
chrom = "chr11"; chromstart = 1650000; chromend = 2350000
bg = sushi.data.Sushi_DNaseI_bedgraph()
fig, ax = plt.subplots(figsize=(10, 5))
plotBedgraph(ax, bg, chrom, chromstart, chromend, colorbycol=SushiColors(5))
labelgenome(ax, chrom, chromstart, chromend, n=4, scale="Mb")
# R mtext("Read Depth", side=2, line=1.75, cex=1, font=2)
ax.set_ylabel("Read Depth", fontsize=10, fontweight="bold", labelpad=8)
# R axis(side=2, las=2, tcl=.2)
ax.tick_params(axis="y", labelsize=8)
out = os.path.join(OUT, "R_Ex3_dnasei_with_yaxis.png")
fig.savefig(out, dpi=100, bbox_inches="tight")
plt.close(fig)
print("saved:", out)
