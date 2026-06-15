"""R Ex4-Ex15 R ↔ Python 1:1 reproduction.

Each R example is taken directly from Sushi.Rnw (the R vignette source)
and reproduced in Python using the sushi-py port. The R rendered figures
are in Sushi_R_real_render/Sushi-XXX.pdf.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg

REPO_DIR = r"D:\奇奇怪怪的图\Sushi_Python_Version"
sys.path.insert(0, REPO_DIR)
import sushi
from sushi import (
    plotBedgraph, plotBed, plotBedpe, plotHic, plotGenes, plotManhattan,
    labelgenome, addlegend, SushiColors, maptocolors, opaque,
)

OUT = r"D:\奇奇怪怪的图\Sushi_Python_Version\examples_R_real"
os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    out = os.path.join(OUT, name + ".png")
    fig.savefig(out, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print("  saved:", out)


# === R Ex4 (Sushi-011, par mfrow=2,1 2 bedgraphs stacked) ===
def R_Ex4():
    ctcf = sushi.data.Sushi_ChIPSeq_CTCF_bedgraph()
    dnase = sushi.data.Sushi_DNaseI_bedgraph()
    # R uses par(mfrow=c(2,1)) which puts CTCF on top, DNaseI on bottom
    # Reproduce by making 2 subplots stacked
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.subplots_adjust(hspace=0.15, left=0.10, right=0.95, bottom=0.10, top=0.95)
    chrom = "chr11"; chromstart = 1955000; chromend = 1960000
    # CTCF on top
    plotBedgraph(axes[0], ctcf, chrom, chromstart, chromend,
                 transparency=0.50, color=SushiColors(2)(2)[0])
    axes[0].set_ylabel("Read Depth", fontsize=8, fontweight="bold")
    # DNaseI on bottom
    plotBedgraph(axes[1], dnase, chrom, chromstart, chromend,
                 transparency=0.50, color=SushiColors(2)(2)[1],
                 overlay=True, rescaleoverlay=True)
    axes[1].set_ylabel("Read Depth", fontsize=8, fontweight="bold")
    labelgenome(axes[1], chrom, chromstart, chromend, n=3, scale="Kb")
    save(fig, "R_Ex4_par_mfrow2_stacked_overlay")


# === R Ex5 (Sushi-013, + legend "topright" with 2 filled rects) ===
def R_Ex5():
    import matplotlib.patches as mpatches
    ctcf = sushi.data.Sushi_ChIPSeq_CTCF_bedgraph()
    dnase = sushi.data.Sushi_DNaseI_bedgraph()
    chrom = "chr11"; chromstart = 1955000; chromend = 1960000
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.subplots_adjust(left=0.10, right=0.95, bottom=0.15, top=0.95)
    plotBedgraph(ax, ctcf, chrom, chromstart, chromend,
                 transparency=0.50, color=SushiColors(2)(2)[0])
    plotBedgraph(ax, dnase, chrom, chromstart, chromend,
                 transparency=0.50, color=SushiColors(2)(2)[1],
                 overlay=True, rescaleoverlay=True)
    labelgenome(ax, chrom, chromstart, chromend, n=3, scale="Kb")
    # R legend: top-right, 2 filled rects (blue + red)
    ax.legend(handles=[
        mpatches.Patch(facecolor=SushiColors(2)(2)[0], label="DNaseI"),
        mpatches.Patch(facecolor=SushiColors(2)(2)[1], label="ChIP-seq (CTCF)"),
    ], loc="upper right", frameon=False, fontsize=9)
    save(fig, "R_Ex5_overlay_with_legend_topright")


# === R Ex8 (Sushi-020, plotHic SushiColors(7) with addlegend right) ===
def R_Ex8():
    hic = sushi.data.Sushi_HiC_matrix()
    chrom = "chr11"; chromstart = 500000; chromend = 5050000
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.subplots_adjust(left=0.10, right=0.85, bottom=0.15, top=0.95)
    # plotHic returns [min_z, max_z] and palette when zrange given
    ret = plotHic(ax, hic, chrom, chromstart, chromend,
                   max_y=20, zrange=(0, 28), palette=SushiColors(7))
    if isinstance(ret, tuple) and len(ret) == 2:
        zrange_out, palette = ret
    else:
        zrange_out, palette = (0, 28), SushiColors(7)
    labelgenome(ax, chrom, chromstart, chromend, n=4, scale="Mb",
                edgeblankfraction=0.20)
    # R addlegend: side="right", bottominset=0.4, topinset=0, xoffset=-.035
    addlegend(ax, zrange_out, palette=palette, title="score", side="right",
              bottominset=0.4, topinset=0, xoffset=-0.035,
              labelside="left", width=0.025, title_offset=0.035)
    save(fig, "R_Ex8_hic_SushiColors7_with_legend_right")


# === R Ex10 (Sushi-023, plotBedpe 5C loops + 3-color legend) ===
def R_Ex10():
    import matplotlib.lines as mlines
    bedpe = sushi.data.Sushi_5C_bedpe()
    chrom = "chr11"; chromstart = 1650000; chromend = 2350000
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.subplots_adjust(left=0.10, right=0.85, bottom=0.15, top=0.95)
    plotBedpe(ax, bedpe, chrom, chromstart, chromend,
              heights=bedpe["score"].tolist(),
              plottype="loops",
              colorby=bedpe["samplenumber"].tolist(),
              colorbycol=SushiColors(3))
    labelgenome(ax, chrom, chromstart, chromend, n=3, scale="Mb")
    # R legend: 3 colored dots (K562, HeLa, GM12878)
    handles = [
        mlines.Line2D([], [], color=SushiColors(3)(3)[i], marker="o",
                      linestyle="None", markersize=8,
                      label=lbl)
        for i, lbl in enumerate(["K562", "HeLa", "GM12878"])
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False,
              fontsize=9, title_fontsize=9)
    ax.set_ylabel("Z-score", fontsize=8, fontweight="bold")
    ax.tick_params(axis="y", labelsize=7)
    save(fig, "R_Ex10_5C_loops_with_3color_legend")


# === R Ex11 (Sushi-024, plotBedpe 5C lines flip + addlegend right) ===
def R_Ex11():
    bedpe = sushi.data.Sushi_5C_bedpe()
    chrom = "chr11"; chromstart = 1650000; chromend = 2350000
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.subplots_adjust(left=0.10, right=0.85, bottom=0.15, top=0.95)
    ret = plotBedpe(ax, bedpe, chrom, chromstart, chromend,
                    heights=bedpe["score"].tolist(),
                    plottype="lines", flip=True,
                    colorby=bedpe["score"].tolist(),
                    colorbycol=SushiColors(5))
    labelgenome(ax, chrom, chromstart, chromend, side=3, n=3, scale="Mb")
    if isinstance(ret, tuple) and len(ret) == 2:
        zrange, palette = ret
    else:
        zrange, palette = (-3, 3), SushiColors(5)
    addlegend(ax, zrange, palette=palette, title="Z-score", side="right",
              bottominset=0.05, topinset=0.05, xoffset=-0.035,
              labelside="right", width=0.025, title_offset=0.045)
    save(fig, "R_Ex11_5C_lines_flip_with_legend")


# === R Ex12 (Sushi-026, plotBed Pol2 region + reverse/forward legend) ===
def R_Ex12():
    import matplotlib.patches as mpatches
    pol2 = sushi.data.Sushi_ChIPSeq_pol2_bed()
    chrom = "chr11"; chromstart = 2281200; chromend = 2282200
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.subplots_adjust(left=0.10, right=0.95, bottom=0.15, top=0.95)
    plotBed(ax, pol2, chrom, chromstart, chromend,
            colorby=pol2["strand"].tolist(),
            colorbycol=SushiColors(2), row="auto", wiggle=0.001)
    labelgenome(ax, chrom, chromstart, chromend, n=2, scale="Kb")
    # R legend: 2 colored rects (reverse, forward)
    ax.legend(handles=[
        mpatches.Patch(facecolor=SushiColors(2)(2)[0], label="reverse"),
        mpatches.Patch(facecolor=SushiColors(2)(2)[0], label="forward"),
    ], loc="upper right", frameon=False, fontsize=8)
    save(fig, "R_Ex12_Pol2_region_with_strand_legend")


# === R Ex14 (Sushi-029, plotBed severalfactors circles + mtext) ===
def R_Ex14():
    sf = sushi.data.Sushi_ChIPSeq_severalfactors_bed()
    chrom = "chr15"; chromstart = 72800000; chromend = 73100000
    sf = sf.copy()
    # R: bed$color = maptocolors(bed$row, SushiColors(6))
    unique_rows = sorted(sf["row"].unique())
    row_to_color = {r: SushiColors(6)(len(unique_rows))[i] for i, r in enumerate(unique_rows)}
    sf["color"] = sf["row"].map(row_to_color)
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.92)
    # rowlabels must match the visible rows (those with peaks in range)
    mask = (sf["chrom"] == chrom) & (sf["start"] >= chromstart) & (sf["start"] < chromend)
    sf_visible = sf[mask]
    rowlabels = list(sf_visible["name"].drop_duplicates())
    # For each visible row, the color is row_to_color[row]
    # unique rows in display order (R preserves first appearance)
    seen_rows = []
    for r in sf_visible["row"]:
        if r not in seen_rows:
            seen_rows.append(r)
    rowlabel_colors = [row_to_color[r] for r in seen_rows]
    plotBed(ax, sf, chrom, chromstart, chromend,
            rownumber=sf["row"].tolist(), type="circles",
            color=sf["color"].tolist(), row="given",
            plotbg="#F2F2F2",
            rowlabels=rowlabels,
            rowlabelcol=rowlabel_colors,
            rowlabelcex=0.75)
    labelgenome(ax, chrom, chromstart, chromend, n=3, scale="Mb")
    # R mtext("ChIP-seq", side=3, adj=-0.065, line=0.5, font=2)
    fig.text(0.02, 0.96, "ChIP-seq", fontsize=12, fontweight="bold",
             ha="left", va="top")
    save(fig, "R_Ex14_severalfactors_circles_with_mtext")


# === R Ex15 (Sushi-032, plotBed severalfactors region + mtext) ===
def R_Ex15():
    sf = sushi.data.Sushi_ChIPSeq_severalfactors_bed()
    chrom = "chr15"; chromstart = 72800000; chromend = 73100000
    sf = sf.copy()
    unique_rows = sorted(sf["row"].unique())
    row_to_color = {r: SushiColors(6)(len(unique_rows))[i] for i, r in enumerate(unique_rows)}
    sf["color"] = sf["row"].map(row_to_color)
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.92)
    mask = (sf["chrom"] == chrom) & (sf["start"] >= chromstart) & (sf["start"] < chromend)
    sf_visible = sf[mask]
    rowlabels = list(sf_visible["name"].drop_duplicates())
    seen_rows = []
    for r in sf_visible["row"]:
        if r not in seen_rows:
            seen_rows.append(r)
    rowlabel_colors = [row_to_color[r] for r in seen_rows]
    plotBed(ax, sf, chrom, chromstart, chromend,
            rownumber=sf["row"].tolist(), type="region",
            color=sf["color"].tolist(), row="given",
            plotbg="#F2F2F2",
            rowlabels=rowlabels,
            rowlabelcol=rowlabel_colors,
            rowlabelcex=0.75)
    labelgenome(ax, chrom, chromstart, chromend, n=3, scale="Mb")
    fig.text(0.02, 0.96, "ChIP-seq", fontsize=12, fontweight="bold",
             ha="left", va="top")
    save(fig, "R_Ex15_severalfactors_region_with_mtext")


print("=== R Ex4-15 1:1 reproduction ===")
for fn in [R_Ex4, R_Ex5, R_Ex8, R_Ex10, R_Ex11, R_Ex12, R_Ex14, R_Ex15]:
    name = fn.__name__
    try:
        fn()
    except Exception as e:
        import traceback
        print(f"  {name} FAILED: {e}")
        traceback.print_exc()
print("=== Done ===")
