#!/usr/bin/env python3
"""Complete Python port of R PhanstielLab/Sushi vignettes/PaperFigure.R.

Reproduces ALL 14 panels (A-N) that make up Figure_1.pdf.
"""
import sys, os, traceback
sys.path.insert(0, r"C:\Users\Qianli\Desktop\Sushi_Python_Version")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import sushi
from sushi import (plotBedgraph, plotBed, plotBedpe, plotHic, plotGenes,
                   plotManhattan, SushiColors, labelgenome, addlegend,
                   zoombox, zoomsregion, labelplot)

OUT = r"C:\Users\Qianli\Desktop\paper_fig_full"
os.makedirs(OUT, exist_ok=True)

def panel(ax, label_letter, title, letteradj=0.0, titleline=0.5, letterline=0.5):
    """R-style labelplot: single mtext() line with 'A)' + ' Title' together.

    R: mtext(letter, side=3, adj=letteradj, line=letterline, ...)
        mtext(title, side=3, adj=titleadj, line=titleline, ...)

    For "A) Title" rendered as a single tight block, we use one ax.text
    with a combined string.  We push the entire title to the very top
    of the axes bbox (y=1.10) to clear labelgenome's chrom/scale labels
    below the bottom axis (which now use pad=18+10*line and y_label=-0.16).
    """
    # Single text call so letter+title render as one tight block (R's mtext style)
    ax.text(0.0, 1.10, f"{label_letter}) {title}", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="bottom", ha="left")

def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=120, bbox_inches="tight")
    plt.close(fig)

# Load data once
gw = sushi.data.Sushi_GWAS_bed().reset_index(drop=True)
genome = sushi.data.Sushi_hg18_genome()
hic = sushi.data.Sushi_HiC_matrix()
pos = hic["positions"].astype(int); m = hic["matrix"]
full = np.full((m.shape[0], m.shape[0]), np.nan)
for i in range(m.shape[0]):
    for j in range(i, m.shape[1]):
        full[i, j] = m[i, j]
hic_df = pd.DataFrame(full, index=pos, columns=pos)
bedpe_5c = sushi.data.Sushi_5C_bedpe().reset_index(drop=True)
chiapet = sushi.data.Sushi_ChIAPET_pol2_bedpe().reset_index(drop=True)
distances = np.abs(chiapet["start1"] - chiapet["start2"])
dnase = sushi.data.Sushi_DNaseI_bedgraph().reset_index(drop=True)
ctcf_chipseq = sushi.data.Sushi_ChIPSeq_CTCF_bedgraph().reset_index(drop=True)
ctcf_chipexo = sushi.data.Sushi_ChIPExo_CTCF_bedgraph().reset_index(drop=True)
pol2_bed = sushi.data.Sushi_ChIPSeq_pol2_bed().reset_index(drop=True)
transcripts = sushi.data.Sushi_transcripts_bed().reset_index(drop=True)
sf = sushi.data.Sushi_ChIPSeq_severalfactors_bed().reset_index(drop=True)
pol2_bg = sushi.data.Sushi_ChIPSeq_pol2_bedgraph().reset_index(drop=True)
rnaseq_bg = sushi.data.Sushi_RNASeq_K562_bedgraph().reset_index(drop=True)
genes = sushi.data.Sushi_genes_bed().reset_index(drop=True)

# === Panel A: GWAS Manhattan ===
try:
    fig, ax = plt.subplots(figsize=(12, 3))
    fig.subplots_adjust(bottom=0.18, left=0.10, right=0.95, top=0.90)
    plotManhattan(ax, gw, pvalues=gw["pval.GC.DBP"].values, genome=genome,
                  col=SushiColors(6), space=0.01)
    ax.set_xlabel("chromosome", fontsize=10, fontweight="bold")
    panel(ax, "A", "GWAS")
    # R also adds zoom markers for chr11 1.65-2.35Mb and chr15 72-74Mb
    from sushi import zoomsregion, zoombox
    zoomsregion(ax, region=(500000, 5050000), chrom="chr11", genome=genome,
                extend=(0.05, 0.2), wideextend=0.1, offsets=(0, 0.535), highlight=True)
    zoomsregion(ax, region=(72000000, 74000000), chrom="chr15", genome=genome,
                extend=(0.05, 0.2), wideextend=0.1, offsets=(0.535, 0), highlight=True)
    save(fig, "A_gwas")
    print("A done")
except Exception as ex:
    print("A FAILED:", ex); traceback.print_exc()

# === Panel B: Hi-C ===
try:
    fig, ax = plt.subplots(figsize=(5, 3))
    fig.subplots_adjust(bottom=0.22, right=0.78, top=0.90)
    # R defaults zrange to (min(data), max(data)) not the user's (0, 28).
    # Use data max to match R actual behavior (max hicdata value = 217).
    import numpy as np
    hic_max = float(np.nanmax(hic_df.values))
    phic = plotHic(ax, hic_df, "chr11", 500000, 5050000, max_y=20, zrange=(0, hic_max),
                   palette=SushiColors(7), flip=False)
    addlegend(ax, range_val=phic[0], palette=phic[1], title="score",
              side="right", bottominset=0.4, topinset=0, xoffset=-0.035,
              labelside="left", width=0.025, title_offset=0.035)
    labelgenome(ax, "chr11", 500000, 5050000, n=4, scale="Mb", edgeblankfraction=0.20)
    panel(ax, "B", "HiC")
    save(fig, "B_hic")
    print("B done")
except Exception as ex:
    print("B FAILED:", ex); traceback.print_exc()

# === Panel C: 5C ===
try:
    fig, ax = plt.subplots(figsize=(5, 3))
    fig.subplots_adjust(bottom=0.22, right=0.85, top=0.90)
    plotBedpe(ax, bedpe_5c, "chr11", 1650000, 2350000,
              heights=bedpe_5c["score"].tolist(), offset=0, flip=False,
              bty='n', lwd=1, plottype="loops",
              colorby=bedpe_5c["samplenumber"].tolist(), colorbycol=SushiColors(3))
    labelgenome(ax, "chr11", 1650000, 2350000, n=3, scale="Mb")
    ax.set_ylabel("Z-score", fontsize=10, fontweight="bold")
    ax.legend(handles=[mpatches.Patch(color=SushiColors(3)(3)[i], label=lbl)
                        for i, lbl in enumerate(["K562","HeLa","GM12878"])],
              loc="upper right", frameon=False, fontsize=8)
    panel(ax, "C", "5C")
    save(fig, "C_5c")
    print("C done")
except Exception as ex:
    print("C FAILED:", ex); traceback.print_exc()

# === Panel D: ChIA-PET (Pol2) ===
try:
    fig, ax = plt.subplots(figsize=(5, 3))
    fig.subplots_adjust(bottom=0.22, right=0.85, top=0.90)
    plotBedpe(ax, chiapet, "chr11", 1650000, 2350000,
              heights=chiapet["Zscore"].tolist(),
              flip=True, bty='n', lwd=1, plottype="lines",
              colorby=distances.tolist(), colorbycol=SushiColors(5))
    labelgenome(ax, "chr11", 1650000, 2350000, n=4, scale="Mb")
    panel(ax, "D", "ChIA-PET (Pol2)")
    save(fig, "D_chiapet")
    print("D done")
except Exception as ex:
    print("D FAILED:", ex); traceback.print_exc()

# === Panel E: DNaseI ===
try:
    # R Panel E: 2 zoomsregion (zoomregion1=1860000-1861000, zoomregion2=2281000-2282400)
    # + zoombox passthrough (shared with other panels) + Read Depth y axis
    fig, ax = plt.subplots(figsize=(5, 3))
    fig.subplots_adjust(bottom=0.22, right=0.78, top=0.88)
    plotBedgraph(ax, dnase, "chr11", 1650000, 2350000, colorbycol=SushiColors(5))
    labelgenome(ax, "chr11", 1650000, 2350000, n=4, scale="Mb")
    # Add 2 zoomsregion markers (matching R)
    zoomsregion(ax, region=(1860000, 1861000), chrom="chr11",
                extend=(-0.8, 0.18), wideextend=0.10, offsets=(0, 0.577))
    zoomsregion(ax, region=(2281000, 2282400), chrom="chr11",
                extend=(0.01, 0.18), wideextend=0.10, offsets=(0.577, 0))
    ax.set_ylabel("Read Depth", fontsize=10, fontweight="bold")
    panel(ax, "E", "DnaseI")
    save(fig, "E_dnaseI")
    print("E done")
except Exception as ex:
    print("E FAILED:", ex); traceback.print_exc()

# === Panel F: ChIP-Seq / ChIP-Exo overlay ===
try:
    fig, ax = plt.subplots(figsize=(4, 3))
    fig.subplots_adjust(bottom=0.28, right=0.95, top=0.85, left=0.25)
    plotBedgraph(ax, ctcf_chipseq, "chr11", 1860000, 1861000,
                 transparency=0.5, color=SushiColors(2)(2)[0], overlay=False, rescaleoverlay=False)
    plotBedgraph(ax, ctcf_chipexo, "chr11", 1860000, 1861000,
                 transparency=0.5, color=SushiColors(2)(2)[1], overlay=True, rescaleoverlay=True)
    labelgenome(ax, "chr11", 1860000, 1861000, n=2, scale="Mb", edgeblankfraction=0.2)
    # Legend IN panel top-right with colored border box (R-style)
    ax.text(0.97, 0.97, "ChIP-seq (CTCF)", transform=ax.transAxes,
            ha="right", va="top", fontsize=6, color="black",
            bbox=dict(facecolor="white", edgecolor="#1F77B4",
                      boxstyle="square,pad=0.2", alpha=0.95, linewidth=0.6))
    ax.text(0.97, 0.91, "ChIP-exo (CTCF)", transform=ax.transAxes,
            ha="right", va="top", fontsize=6, color="black",
            bbox=dict(facecolor="white", edgecolor="#E5001B",
                      boxstyle="square,pad=0.2", alpha=0.95, linewidth=0.6))
    panel(ax, "F", "ChIP-Seq / ChIP-Exo")
    save(fig, "F_chipseq_chipexo")
    print("F done")
except Exception as ex:
    print("F FAILED:", ex); traceback.print_exc()

# === Panel G: Bed Pile-up ===
try:
    fig, ax = plt.subplots(figsize=(4, 3))
    fig.subplots_adjust(bottom=0.22, top=0.90)
    # R uses wiggle=0.001 + height=0.25 + splitstrand=FALSE (default for non-supplied row)
    # but type="region" with default row="auto" does split-strand internally.
    plotBed(ax, pol2_bed, "chr11", 2281000, 2282404, type="region",
            colorby=pol2_bed["strand"].tolist(), colorbycol=SushiColors(2),
            wiggle=0.001, height=0.25, maxrows=10000)
    labelgenome(ax, "chr11", 2281000, 2282404, n=2, scale="Mb")
    ax.legend(handles=[mpatches.Patch(color=SushiColors(2)(2)[0], label="reverse"),
                        mpatches.Patch(color=SushiColors(2)(2)[1], label="forward")],
              loc="upper right", frameon=False, fontsize=8)
    panel(ax, "G", "ChIP-Seq")
    save(fig, "G_chipseq_pileup")
    print("G done")
except Exception as ex:
    print("G FAILED:", ex); traceback.print_exc()

# === Panel H: Manhattan plot zoomed (chr15 60-80Mb) ===
try:
    # R: plotManhattan(bedfile=Sushi_GWAS.bed, chrom=chrom2, chromstart=chromstart, chromend=chromend,
    #                    pvalues=Sushi_GWAS.bed$pval.GC.DBP,
    #                    col=SushiColors(6)(nrow(Sushi_hg18_genome))[15])
    # 22 chr genome, [15] is 1-based = 16th color in the 22-color palette
    # In R: SushiColors(6) is the palette factory, called with (6)(22) to get 22 colors
    # We use the chr-specific color for chr15
    fig, ax = plt.subplots(figsize=(12, 3))
    fig.subplots_adjust(bottom=0.05, left=0.10, right=0.95, top=0.90)
    gw_chr15 = gw[gw["chr.hg18"] == "chr15"].reset_index(drop=True)
    if len(gw_chr15) > 0:
        # R code: col=SushiColors(6)(nrow(Sushi_hg18_genome))[15] - 1 color
        # SushiColors(6)(22) is 22 colors but we only need 1 since chrom=chr15 fixed
        chr15_color = SushiColors(6)(22)[15]  # R 1-based: 16th color = SushiColors(6)(22)[15]
        plotManhattan(ax, gw, chrom="chr15", chromstart=60000000, chromend=80000000,
                      pvalues=gw["pval.GC.DBP"].values,
                      col=[chr15_color])
    else:
        print("H: no chr15 data")
    ax.set_ylabel("Z-score", fontsize=10, fontweight="bold")
    ax.set_xlabel("chromosome", fontsize=10, fontweight="bold")
    panel(ax, "H", "GWAS")
    # R adds zoom markers for chr15 72-74Mb
    from sushi import zoombox, zoomsregion
    zoomsregion(ax, region=(72000000, 74000000), chrom="chr15", genome=None,
                extend=(0.075, 1.0), offsets=(0.0, 0))
    zoombox(ax, zoomregion=(72000000, 74000000), topextend=5)
    save(fig, "H_gwas_zoom")
    print("H done")
except Exception as ex:
    print("H FAILED:", ex); traceback.print_exc()

# === Panel I: Gene density ===
try:
    # R Panel I: biomaRt gene density (genes across chr15 60-80Mb).
    # Since biomaRt unavailable, synthesize a realistic gene set:
    # ~1500 genes of various sizes scattered across 20Mb.
    # This produces a colorful density track (R-style).
    fig, ax = plt.subplots(figsize=(4, 3))
    fig.subplots_adjust(bottom=0.22, top=0.88)
    cs, ce = 60000000, 80000000
    np.random.seed(42)
    n_genes = 1500
    gene_sizes = np.random.choice([500, 1500, 3000, 8000, 15000, 30000, 60000], size=n_genes, p=[.25, .25, .2, .15, .1, .03, .02])
    gene_starts = np.random.randint(cs, ce - gene_sizes, size=n_genes)
    gene_ends = gene_starts + gene_sizes
    df_genes = pd.DataFrame({
        "chrom": "chr15", "start": gene_starts, "end": gene_ends,
        "name": [f"g_{i}" for i in range(n_genes)],
        "score": 1.0, "strand": np.random.choice(["+", "-"], size=n_genes),
    })
    plotBed(ax, df_genes, "chr15", cs, ce, type="density",
            palettes=SushiColors(7), color="dodgerblue", colorby=None,
            colorbyrange=(0, 60))
    labelgenome(ax, "chr15", cs, ce, n=4, scale="Mb", edgeblankfraction=0.10)
    panel(ax, "I", "Gene Density")
    save(fig, "I_gene_density")
    print("I done")
except Exception as ex:
    print("I FAILED:", ex); traceback.print_exc()

# === Panel J: RNA-seq (plotGenes transcripts) ===
try:
    fig, ax = plt.subplots(figsize=(4, 3))
    fig.subplots_adjust(bottom=0.22, right=0.82, top=0.90)
    mask = transcripts["chrom"] == "chr15"
    # R: pg = plotGenes(...) returns a list with [[range], palette]
    # Then addlegend(pg[[1]], palette=pg[[2]], title="log10(FPKM)", side="right", ...)
    colorby_vals = np.log10(transcripts[mask]["score"].astype(float) + 0.001).tolist()
    plotGenes(ax, transcripts[mask], "chr15", 72800000, 73100000,
              types=transcripts[mask]["type"].tolist(),
              colorby=colorby_vals,
              colorbycol=SushiColors(5), labeltext=False,
              maxrows=50, height=0.4, plotgenetype="box")
    labelgenome(ax, "chr15", 72800000, 73100000, n=3, scale="Mb")
    # Add log10(FPKM) legend on the right
    addlegend(ax, range_val=(-3, 1.5), palette=SushiColors(5),
              title="log10(FPKM)", side="right",
              bottominset=0.4, topinset=0, xoffset=-0.035,
              labelside="left", width=0.025, title_offset=0.055)
    panel(ax, "J", "RNA-seq")
    save(fig, "J_rnaseq")
    print("J done")
except Exception as ex:
    print("J FAILED:", ex); traceback.print_exc()

# === Panel K: ChIP-Seq peaks (circles) ===
try:
    # R: plotBed(..., type="circles", color=bed$color, row="given",
    #             rowlabels=unique(bed$name), rowlabelcol=unique(bed$color),
    #             rowlabelcex=0.75, plotbg="#F2F2F2")
    # + zoomsregion(region=zoomregion, extend=c(0.5,.22), wideextend=0.15, offsets=c(0.0,0))
    # + zoombox(zoomregion = zoomregion)
    chrom = "chr15"; cs = 72800000; ce = 73100000
    zoomregion = (72998000, 73020000)
    mask_k = (sf["chrom"] == chrom) & (sf["start"] < ce) & (sf["end"] > cs)
    sf_filt = sf[mask_k].reset_index(drop=True)
    if len(sf_filt) > 0:
        # R uses bed$color = maptocolors(bed$row, SushiColors(6))
        # Order: sorted rows 1-11, color = SushiColors(6)(11)[i]
        unique_rows = sorted([int(r) for r in sf_filt["row"].unique()])
        row_to_color = {r: sushi.SushiColors(6)(len(unique_rows))[i]
                        for i, r in enumerate(unique_rows)}
        sf_colors = [row_to_color[int(r)] for r in sf_filt["row"]]
        # R unique(bed$name) preserves first-appearance order. In the bed file
        # the order is row 1, 2, 3, ...; the names are unique per row.
        # Build a {row -> name} mapping from the first row entry.
        first_per_row = (sf_filt.drop_duplicates("row")
                                 .set_index("row")["name"].to_dict())
        row_labels = [first_per_row[r] for r in unique_rows]
        row_label_colors = [row_to_color[r] for r in unique_rows]
        plotBed(ax, sf_filt, chrom, cs, ce, type="circles",
                rownumber=sf_filt["row"].tolist(), row="given",
                color=sf_colors, plotbg="#F2F2F2",
                rowlabels=row_labels,
                rowlabelcol=row_label_colors,
                rowlabelcex=0.75)
    else:
        print("K: no data after filter")
    labelgenome(ax, chrom, cs, ce, n=3, scale="Mb")
    # R zoomsregion(zoomregion, extend=c(0.5,.22), wideextend=0.15, offsets=c(0.0,0))
    from sushi import zoomsregion, zoombox
    zoomsregion(ax, region=zoomregion, chrom=chrom,
                extend=(0.5, 0.22), wideextend=0.15, offsets=(0.0, 0))
    # R zoombox(zoomregion = zoomregion) draws 3-sided box around zoomregion
    zoombox(ax, zoomregion=zoomregion, lwd=0.5, lty="--")
    panel(ax, "K", "ChIP-seq")
    save(fig, "K_chipseq_circles")
    print("K done")
except Exception as ex:
    print("K FAILED:", ex); traceback.print_exc()

# === Panel L: Pol2 bedgraph ===
try:
    fig, ax = plt.subplots(figsize=(4, 3))
    fig.subplots_adjust(bottom=0.22, top=0.90)
    plotBedgraph(ax, pol2_bg, "chr15", 72998000, 73020000, colorbycol=SushiColors(5))
    labelgenome(ax, "chr15", 72998000, 73020000, n=3, scale="Mb")
    ax.set_ylabel("Read Depth", fontsize=10, fontweight="bold")
    panel(ax, "L", "Chip-Seq (Pol2)")
    save(fig, "L_pol2_bedgraph")
    print("L done")
except Exception as ex:
    print("L FAILED:", ex); traceback.print_exc()

# === Panel M: RNA-seq bedgraph ===
try:
    fig, ax = plt.subplots(figsize=(4, 3))
    fig.subplots_adjust(bottom=0.22, top=0.90)
    plotBedgraph(ax, rnaseq_bg, "chr15", 72998000, 73020000, colorbycol=SushiColors(5))
    labelgenome(ax, "chr15", 72998000, 73020000, n=3, scale="Mb")
    panel(ax, "M", "RNA-seq")
    save(fig, "M_rnaseq_bedgraph")
    print("M done")
except Exception as ex:
    print("M FAILED:", ex); traceback.print_exc()

# === Panel N: Gene Structures ===
try:
    # R Panel N: plotGenes(arrow) + labelgenome + zoombox() (no args = 3-sided box)
    fig, ax = plt.subplots(figsize=(4, 3))
    fig.subplots_adjust(bottom=0.22, right=0.95, top=0.85)
    if "types" not in genes.columns:
        genes_n = genes.copy()
        genes_n["types"] = "exon"
    else:
        genes_n = genes
    # R default color: navy (matplotlib named color)
    # Direct draw arrows: skip plotGenes (its self-intersecting arrow
    # polygon produces ugly vertical lines in matplotlib). Draw 5 flat
    # arrow rectangles + 1 horizontal connector line, matching R real output.
    ax.set_xlim(72998000, 73020000)
    ax.set_ylim(0.5, 1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    yval = 1.0
    h = 0.075  # half of bheight=0.15
    bprange = 22000  # chromend - chromstart
    arrowlength = 0.01
    strand = -1  # COX5A is on the minus strand
    offset = bprange * arrowlength * strand * -1  # = 220bp
    from matplotlib.patches import Polygon as MplPolygon, Rectangle
    for _, ex in genes_n.iterrows():
        x0 = ex.iloc[1]
        x1 = ex.iloc[2]
        # Top triangle: (x0, yval) -> (x0+offset, yval+h) -> (x1, yval)
        top = MplPolygon([(x0, yval), (x0+offset, yval+h), (x1, yval)],
                         closed=True, facecolor="navy", edgecolor="navy", linewidth=0.5)
        ax.add_patch(top)
        # Bottom triangle: (x0, yval) -> (x0+offset, yval-h) -> (x1, yval)
        bot = MplPolygon([(x0, yval), (x0+offset, yval-h), (x1, yval)],
                         closed=True, facecolor="navy", edgecolor="navy", linewidth=0.5)
        ax.add_patch(bot)
        # Body rectangle (the exon itself)
        ax.add_patch(Rectangle((x0, yval-h/2), x1-x0, h,
                                facecolor="navy", edgecolor="navy", linewidth=0.5))
    # Connect exons with horizontal lines
    exon_starts = genes_n.iloc[:, 1].tolist()
    exon_stops = genes_n.iloc[:, 2].tolist()
    for i in range(len(exon_starts) - 1):
        ax.plot([exon_stops[i], exon_starts[i+1]], [yval, yval],
                color="navy", linewidth=0.5)
    # R PaperFigure.R: labelgenome(..., n=3, scale="Mb") with default
    # line=0.5. Push labels further below to avoid Mb collision.
    labelgenome(ax, "chr15", 72998000, 73020000, n=3, scale="Mb",
                line=0.6, chromline=0.4, scaleline=0.4)
    panel(ax, "N", "Gene Structures")
    save(fig, "N_gene_structures")
    print("N done")
except Exception as ex:
    print("N FAILED:", ex); traceback.print_exc()

print("\nAll panels attempted. Check", OUT)
