"""Example 13: gene model with arrow style + label.
Output: /tmp/genes.png
"""
import matplotlib.pyplot as plt
import sushi
from sushi import plotGenes, labelgenome

genes = sushi.data.Sushi_genes_bed()
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22)
plotGenes(ax, genes, "chr15", 72998000, 73020000,
          types=["exon"] * len(genes),
          maxrows=1, bheight=0.2, plotgenetype="arrow",
          bentline=False, labeloffset=0.4, fontsize=1.2,
          arrowlength=0.025, labeltext=True)
labelgenome(ax, "chr15", 72998000, 73020000, n=3, scale="Mb", side=1)
fig.savefig("/tmp/genes.png", dpi=120, bbox_inches="tight")
