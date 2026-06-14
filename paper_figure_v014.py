#!/usr/bin/env python3
"""sushi-py v0.1.4 equivalent of R PhanstielLab/Sushi vignettes/PaperFigure.R.

Generates a 14-panel layout that should match R's Figure_1.pdf 1:1
for visual diffing.
"""
import sys, os
sys.path.insert(0, r"C:\Users\Qianli\Desktop\Sushi_Python_Version")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sushi
from sushi import (plotBedgraph, plotBed, plotBedpe, plotHic, plotGenes,
                   plotManhattan, SushiColors, labelgenome, addlegend,
                   zoombox, zoomsregion, labelplot)

OUT_DIR = r"C:\Users\Qianli\Desktop\paper_fig_py_v014"
os.makedirs(OUT_DIR, exist_ok=True)

def save(fig, name):
    fig.savefig(os.path.join(OUT_DIR, name + ".png"), dpi=120, bbox_inches="tight")
    plt.close(fig)

# === Panel A: GWAS Manhattan ===
gw = sushi.data.Sushi_GWAS_bed()
genome = sushi.data.Sushi_hg18_genome()
fig, ax = plt.subplots(figsize=(14, 3))
fig.subplots_adjust(bottom=0.18, left=0.10, right=0.95, top=0.92)
plotManhattan(ax, gw, pvalues=gw["pval.GC.DBP"].values, genome=genome, col=SushiColors(6), space=0.05)
ax.set_xlabel("chromosome", fontsize=10, fontweight="bold")
ax.text(0.01, 1.05, "A)", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
ax.text(0.04, 1.05, "GWAS", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
save(fig, "01_panel_A_gwas")

# === Panel B: HiC ===
hic = sushi.data.Sushi_HiC_matrix()
import pandas as pd
pos = hic["positions"].astype(int); m = hic["matrix"]
full = np.full((m.shape[0], m.shape[0]), np.nan)
for i in range(m.shape[0]):
    for j in range(i, m.shape[1]):
        full[i, j] = m[i, j]
hic_df = pd.DataFrame(full, index=pos, columns=pos)
fig, ax = plt.subplots(figsize=(10, 4))
fig.subplots_adjust(bottom=0.22, right=0.82)
phic = plotHic(ax, hic_df, "chr11", 500000, 5050000, max_y=20, zrange=(0, 28),
               palette=sushi.SushiColors(7), flip=False)
addlegend(ax, range_val=phic[0], palette=phic[1], title="score",
          side="right", bottominset=0.4, topinset=0, xoffset=-0.035,
          labelside="left", width=0.025, title_offset=0.035)
labelgenome(ax, "chr11", 500000, 5050000, n=4, scale="Mb", edgeblankfraction=0.20)
ax.text(0.01, 1.05, "B)", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
ax.text(0.04, 1.05, "HiC", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
save(fig, "02_panel_B_hic")

# === Panel C: 5C ===
bedpe_5c = sushi.data.Sushi_5C_bedpe().reset_index(drop=True)
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22, right=0.85, top=0.92)
plotBedpe(ax, bedpe_5c, "chr11", 1650000, 2350000,
          heights=bedpe_5c["score"].tolist(), offset=0, flip=False,
          bty='n', lwd=1, plottype="loops",
          colorby=bedpe_5c["samplenumber"].tolist(), colorbycol=SushiColors(3))
labelgenome(ax, "chr11", 1650000, 2350000, n=3, scale="Mb")
ax.set_ylabel("Z-score", fontsize=10, fontweight="bold")
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(color=SushiColors(3)(3)[i], label=lbl)
                    for i, lbl in enumerate(["K562","HeLa","GM12878"])],
          loc="upper right", frameon=False, fontsize=8)
ax.text(0.01, 1.05, "C)", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
ax.text(0.04, 1.05, "5C", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
save(fig, "03_panel_C_5c")

# === Panel D: ChIA-PET (Pol2) ===
chiapet = sushi.data.Sushi_ChIAPET_pol2_bedpe().reset_index(drop=True)
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22, right=0.82, top=0.92)
distances = np.abs(chiapet["start1"] - chiapet["start2"])
plotBedpe(ax, chiapet, "chr11", 1650000, 2350000,
          heights=chiapet["Zscore"].tolist(),
          flip=True, bty='n', lwd=1, plottype="lines",
          colorby=distances.tolist(), colorbycol=SushiColors(5))
labelgenome(ax, "chr11", 1650000, 2350000, n=4, scale="Mb")
ax.text(0.01, 1.05, "D)", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
ax.text(0.04, 1.05, "ChIA-PET (Pol2)", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
save(fig, "04_panel_D_chiapet")

# === Panel E: DNaseI (with colorbycol gradient) ===
dnase = sushi.data.Sushi_DNaseI_bedgraph()
fig, ax = plt.subplots(figsize=(10, 3))
fig.subplots_adjust(bottom=0.22, top=0.92)
plotBedgraph(ax, dnase, "chr11", 1650000, 2350000, colorbycol=SushiColors(5))
labelgenome(ax, "chr11", 1650000, 2350000, n=4, scale="Mb")
ax.set_ylabel("Read Depth", fontsize=10, fontweight="bold")
ax.text(0.01, 1.05, "E)", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
ax.text(0.04, 1.05, "DnaseI", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
save(fig, "05_panel_E_dnaseI")

# === Panel F: ChIP-Seq / ChIP-Exo (overlay) ===
ctcf_chipseq = sushi.data.Sushi_ChIPSeq_CTCF_bedgraph().reset_index(drop=True)
ctcf_chipexo = sushi.data.Sushi_ChIPExo_CTCF_bedgraph().reset_index(drop=True)
fig, ax = plt.subplots(figsize=(6, 3))
fig.subplots_adjust(bottom=0.22, top=0.92)
plotBedgraph(ax, ctcf_chipseq, "chr11", 1860000, 1861000,
             transparency=0.5, color=SushiColors(2)(2)[1], overlay=False, rescaleoverlay=False)
plotBedgraph(ax, ctcf_chipexo, "chr11", 1860000, 1861000,
             transparency=0.5, color=SushiColors(2)(2)[1], overlay=True, rescaleoverlay=True)
labelgenome(ax, "chr11", 1860000, 1861000, n=3, scale="Mb", edgeblankfraction=0.2)
ax.legend(handles=[mpatches.Patch(color=SushiColors(2)(2)[0], label="ChIP-seq (CTCF)", alpha=0.5),
                    mpatches.Patch(color=SushiColors(2)(2)[1], label="ChIP-exo (CTCF)", alpha=0.5)],
          loc="upper right", frameon=False, fontsize=8)
ax.text(0.01, 1.05, "F)", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
ax.text(0.04, 1.05, "ChIP-Seq / ChIP-Exo", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
save(fig, "06_panel_F_chipseq_chipexo")

# === Panel G: ChIP-Seq pile-up (splitstrand) ===
pol2 = sushi.data.Sushi_ChIPSeq_pol2_bed().reset_index(drop=True)
fig, ax = plt.subplots(figsize=(6, 3))
fig.subplots_adjust(bottom=0.22, top=0.92)
plotBed(ax, pol2, "chr11", 2281000, 2282404, type="region",
        colorby=pol2["strand"].tolist(), colorbycol=SushiColors(2),
        wiggle=0.001, height=0.25)
labelgenome(ax, "chr11", 2281000, 2282400, n=2, scale="Mb")
ax.legend(handles=[mpatches.Patch(color=SushiColors(2)(2)[0], label="reverse"),
                    mpatches.Patch(color=SushiColors(2)(2)[1], label="forward")],
          loc="upper right", frameon=False, fontsize=8)
ax.text(0.01, 1.05, "G)", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
ax.text(0.04, 1.05, "ChIP-Seq", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
save(fig, "07_panel_G_chipseq_pileup")

print("First 7 panels done. Check images in", OUT_DIR)
