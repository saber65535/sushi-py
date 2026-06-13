"""Example 12: Manhattan plot for GWAS p-values across all chroms.
Output: /tmp/manhattan.png
"""
import matplotlib.pyplot as plt
import sushi
from sushi import plotManhattan, labelgenome, SushiColors

gwas = sushi.data.Sushi_GWAS_bed()
genome = sushi.data.Sushi_hg18_genome()

fig, ax = plt.subplots(figsize=(12, 4))
fig.subplots_adjust(bottom=0.22)
plotManhattan(ax, gwas, pvalues=gwas.iloc[:, 4].values,
              genome=genome, col=SushiColors(6), space=0.01)
labelgenome(ax, "chr1", 0, 1, genome=genome, n=4, scale="Mb",
            edgeblankfraction=0.20, space=0.01)
fig.savefig("/tmp/manhattan.png", dpi=120, bbox_inches="tight")
