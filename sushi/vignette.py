"""
Reproduces all examples from Sushi R vignette (inst/doc/Sushi.R) using
the Python port.

Run: python -m sushi.vignette
Or:  python /path/to/sushi/vignette.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sushi
from sushi import (
    SushiColors, maptocolors, opaque, maptolwd,
    plotBedgraph, plotBed, plotBedpe, plotHic, plotGenes, plotManhattan,
    labelgenome, addlegend, labelplot, zoombox, zoomsregion,
)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vignette_output")
os.makedirs(OUT, exist_ok=True)
results = []


def run_example(name, fn):
    try:
        fn()
        results.append((name, "OK"))
        print(f"  [OK]   {name}")
    except Exception as e:
        import traceback
        results.append((name, f"FAIL: {e}"))
        print(f"  [FAIL] {name}: {e}")
        traceback.print_exc()


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=100, bbox_inches="tight")
    plt.close(fig)


# ---- 15 example functions ----
def ex01_dnasei_bedgraph():
    fig, ax = plt.subplots(figsize=(10, 3)); fig.subplots_adjust(bottom=0.22)
    bg = sushi.data.Sushi_DNaseI_bedgraph()
    plotBedgraph(ax, bg, "chr11", 1650000, 2350000, colorbycol=SushiColors(5))
    labelgenome(ax, "chr11", 1650000, 2350000, n=4, scale="Mb", side=1)
    save(fig, "01_dnasei_bedgraph")


def ex02_overlay_dnasei_ctcf():
    fig, ax = plt.subplots(figsize=(10, 3)); fig.subplots_adjust(bottom=0.22)
    ctcf = sushi.data.Sushi_ChIPSeq_CTCF_bedgraph()
    dnase = sushi.data.Sushi_DNaseI_bedgraph()
    plotBedgraph(ax, ctcf, "chr11", 1955000, 1960000, transparency=0.5, color=SushiColors(2)(2)[0])
    plotBedgraph(ax, dnase, "chr11", 1955000, 1960000, transparency=0.5, color=SushiColors(2)(2)[1], overlay=True, rescaleoverlay=True)
    labelgenome(ax, "chr11", 1955000, 1960000, n=3, scale="Kb", side=1)
    save(fig, "02_overlay_dnasei_ctcf")


def ex03_flip():
    fig, ax = plt.subplots(figsize=(10, 3)); fig.subplots_adjust(bottom=0.22)
    ctcf = sushi.data.Sushi_ChIPSeq_CTCF_bedgraph()
    dnase = sushi.data.Sushi_DNaseI_bedgraph()
    plotBedgraph(ax, ctcf, "chr11", 1955000, 1960000, transparency=0.5, color=SushiColors(2)(2)[0])
    plotBedgraph(ax, dnase, "chr11", 1955000, 1960000, transparency=0.5, flip=True, color=SushiColors(2)(2)[1])
    labelgenome(ax, "chr11", 1955000, 1960000, side=3, n=3, scale="Kb")
    save(fig, "03_flip_dnasei_ctcf")


def ex04_hic():
    hic = sushi.data.Sushi_HiC_matrix()
    # R stores HiC as lower-triangular (114 rows, 113 cols). Build a square
    # (114, 114) matrix padded with NaN in the upper triangle so that
    # pandas + plotHic are happy.
    import numpy as _np, pandas as _pd
    pos = hic["positions"].astype(int)
    m = hic["matrix"]
    full = _np.full((m.shape[0], m.shape[0]), _np.nan)
    for i in range(m.shape[0]):
        for j in range(i, m.shape[1]):
            full[i, j] = m[i, j]
    hic_df = _pd.DataFrame(full, index=pos, columns=pos)
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(bottom=0.22, right=0.82)
    phic = plotHic(ax, hic_df, "chr11", 500000, 5050000, max_y=20,
                   zrange=(0, 28), palette=SushiColors(7))
    addlegend(ax, range_val=phic[0], palette=phic[1], title="score", side="right",
              bottominset=0.4, topinset=0, xoffset=-0.035,
              labelside="left", width=0.025, title_offset=0.035)
    labelgenome(ax, "chr11", 500000, 5050000, n=4, scale="Mb",
                edgeblankfraction=0.20, side=1)
    save(fig, "04_hic_chr11")

def ex05_bedpe_loops():
    bedpe = sushi.data.Sushi_5C_bedpe()
    fig, ax = plt.subplots(figsize=(10, 3)); fig.subplots_adjust(bottom=0.22)
    plotBedpe(ax, bedpe, "chr11", 1650000, 2350000, heights=bedpe["score"].tolist(), plottype="loops", colorby=bedpe["samplenumber"].tolist(), colorbycol=SushiColors(3))
    labelgenome(ax, "chr11", 1650000, 2350000, n=3, scale="Mb", side=1)
    save(fig, "05_bedpe_5C_loops")


def ex06_bedpe_ribbons():
    bedpe = sushi.data.Sushi_5C.bedpe() if False else sushi.data.Sushi_5C_bedpe()
    fig, ax = plt.subplots(figsize=(10, 3)); fig.subplots_adjust(bottom=0.22)
    plotBedpe(ax, bedpe, "chr11", 1650000, 2350000, heights=bedpe["score"].tolist(), plottype="ribbons", colorby=bedpe["samplenumber"].tolist(), colorbycol=SushiColors(3), border="black")
    labelgenome(ax, "chr11", 1650000, 2350000, n=3, scale="Mb", side=1)
    save(fig, "06_bedpe_5C_ribbons")


def _filter_bed_by_region(bed, chrom, chromstart, chromend):
    """Pre-filter a BED dataframe to the plotting region (matches plotBed's
    internal filter). Returns (filtered_bed, strand_as_int_column)."""
    import pandas as _pd
    chrom_str = bed.iloc[:, 0].astype(str)
    start_in = (bed.iloc[:, 1] > chromstart) & (bed.iloc[:, 1] < chromend)
    end_in = (bed.iloc[:, 2] > chromstart) & (bed.iloc[:, 2] < chromend)
    mask = (chrom_str == chrom) & (start_in | end_in)
    return bed[mask].reset_index(drop=True)


def ex07_bed_region():
    bed_full = sushi.data.Sushi_ChIPSeq_pol2_bed()
    bed = _filter_bed_by_region(bed_full, "chr11", 2281200, 2282200)
    cb = bed["strand"].map({"+": 1, "-": -1}).tolist() if len(bed) > 0 else []
    fig, ax = plt.subplots(figsize=(10, 3)); fig.subplots_adjust(bottom=0.22, left=0.12)
    plotBed(ax, bed, "chr11", 2281200, 2282200, colorby=cb, colorbycol=SushiColors(2), row="auto", wiggle=0.001)
    labelgenome(ax, "chr11", 2281200, 2282200, n=2, scale="Kb", side=1)
    save(fig, "07_bed_pol2_region")


def ex08_bed_circles():
    bed_full = sushi.data.Sushi_ChIPSeq_pol2_bed()
    bed = _filter_bed_by_region(bed_full, "chr11", 2281200, 2282200)
    cb = bed["strand"].map({"+": 1, "-": -1}).tolist() if len(bed) > 0 else []
    fig, ax = plt.subplots(figsize=(10, 3)); fig.subplots_adjust(bottom=0.22, left=0.12)
    plotBed(ax, bed, "chr11", 2281200, 2282200, type="circles", colorby=cb, colorbycol=SushiColors(2), row="auto", wiggle=0.001)
    labelgenome(ax, "chr11", 2281200, 2282200, n=2, scale="Kb", side=1)
    save(fig, "08_bed_pol2_circles")


def ex09_bed_splitstrand():
    bed_full = sushi.data.Sushi_ChIPSeq_pol2_bed()
    bed = _filter_bed_by_region(bed_full, "chr11", 2281200, 2282200)
    cb = bed["strand"].map({"+": 1, "-": -1}).tolist() if len(bed) > 0 else []
    fig, ax = plt.subplots(figsize=(10, 3)); fig.subplots_adjust(bottom=0.22, left=0.12)
    plotBed(ax, bed, "chr11", 2281200, 2282200, colorby=cb, colorbycol=SushiColors(2), row="auto", splitstrand=True, wiggle=0.001)
    labelgenome(ax, "chr11", 2281200, 2282200, n=2, scale="Kb", side=1)
    save(fig, "09_bed_pol2_splitstrand")


def ex10_severalfactors_circles():
    bed_full = sushi.data.Sushi_ChIPSeq_severalfactors_bed()
    bed = _filter_bed_by_region(bed_full, "chr15", 72800000, 73100000)
    bed["color"] = maptocolors(bed["row"].tolist(), col=SushiColors(6)) if len(bed) > 0 else []
    fig, ax = plt.subplots(figsize=(10, 5)); fig.subplots_adjust(bottom=0.22)
    plotBed(ax, bed, "chr15", 72800000, 73100000, type="circles",
            color=bed["color"].tolist() if len(bed) > 0 else "navy",
            rownumber=bed["row"].tolist() if len(bed) > 0 else None,
            row="given", plotbg="#F2F2F2",
            rowlabels=bed_full["name"].drop_duplicates().tolist()[:11],
            rowlabelcex=0.7)
    labelgenome(ax, "chr15", 72800000, 73100000, n=3, scale="Mb", side=1)
    save(fig, "10_severalfactors_circles")


def ex11_severalfactors_region():
    bed_full = sushi.data.Sushi_ChIPSeq_severalfactors_bed()
    bed = _filter_bed_by_region(bed_full, "chr15", 72800000, 73100000)
    bed["color"] = maptocolors(bed["row"].tolist(), col=SushiColors(6)) if len(bed) > 0 else []
    fig, ax = plt.subplots(figsize=(10, 5)); fig.subplots_adjust(bottom=0.22)
    plotBed(ax, bed, "chr15", 72800000, 73100000, type="region",
            color=bed["color"].tolist() if len(bed) > 0 else "navy",
            rownumber=bed["row"].tolist() if len(bed) > 0 else None,
            row="given", plotbg="#F2F2F2",
            rowlabels=bed_full["name"].drop_duplicates().tolist()[:11],
            rowlabelcex=0.7)
    labelgenome(ax, "chr15", 72800000, 73100000, n=3, scale="Mb", side=1)
    save(fig, "11_severalfactors_region")


def ex12_manhattan():
    gwas = sushi.data.Sushi_GWAS_bed()
    genome = sushi.data.Sushi_hg18_genome()
    fig, ax = plt.subplots(figsize=(12, 4)); fig.subplots_adjust(bottom=0.22)
    plotManhattan(ax, gwas, pvalues=gwas.iloc[:, 4].values, genome=genome, col=SushiColors(6), space=0.01)
    # R vignette: in genome-only mode, labelgenome takes chrom as NA and
    # uses genome's chromOffsets. We emulate by setting chrom to chr1 and
    # letting ax's xlim come from plotManhattan.
    labelgenome(ax, chrom="chr1", chromstart=0, chromend=1, genome=genome, n=4,
                scale="Mb", edgeblankfraction=0.20, space=0.01)
    save(fig, "12_manhattan_gwas")


def ex13_genes():
    genes = sushi.data.Sushi_genes_bed()
    fig, ax = plt.subplots(figsize=(10, 3)); fig.subplots_adjust(bottom=0.22)
    plotGenes(ax, genes, "chr15", 72998000, 73020000, types=["exon"] * len(genes), maxrows=1, bheight=0.2, plotgenetype="arrow", bentline=False, labeloffset=0.4, fontsize=1.2, arrowlength=0.025, labeltext=True)
    labelgenome(ax, "chr15", 72998000, 73020000, n=3, scale="Mb", side=1)
    save(fig, "13_genes_chr15")


def ex14_transcripts_colorby():
    trans = sushi.data.Sushi_transcripts_bed()
    fig, ax = plt.subplots(figsize=(10, 5)); fig.subplots_adjust(bottom=0.22)
    plotGenes(ax, trans, "chr15", 72965000, 72990000, types=trans["type"].tolist(), colorby=np.log10(trans["score"].values + 0.001).tolist(), colorbycol=SushiColors(5), colorbyrange=(0, 1.0), labeltext=True, maxrows=50, bheight=0.4, plotgenetype="box")
    labelgenome(ax, "chr15", 72965000, 72990000, n=3, scale="Mb", side=1)
    save(fig, "14_transcripts_colorby")


def ex15_zoom():
    bg = sushi.data.Sushi_DNaseI_bedgraph()
    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])
    ax_top = fig.add_subplot(gs[0, :])
    ax_b1 = fig.add_subplot(gs[1, 0])
    ax_b2 = fig.add_subplot(gs[1, 1])
    fig.subplots_adjust(bottom=0.10, hspace=0.5, left=0.10, right=0.95)
    chrom, cs, ce = "chr11", 1900000, 2350000
    plotBedgraph(ax_top, bg, chrom, cs, ce, colorbycol=SushiColors(5))
    labelgenome(ax_top, chrom, cs, ce, n=4, scale="Mb", side=1)
    zoomsregion(ax_top, [1955000, 1960000], extend=(0.01, 0.13), wideextend=0.05, offsets=(0, 0.580))
    zoomsregion(ax_top, [2279000, 2284000], extend=(0.01, 0.13), wideextend=0.05, offsets=(0.580, 0))
    zr1 = (1955000, 1960000)
    plotBedgraph(ax_b1, bg, chrom, zr1[0], zr1[1], colorbycol=SushiColors(5))
    labelgenome(ax_b1, chrom, zr1[0], zr1[1], n=4, scale="Kb", edgeblankfraction=0.20)
    zoombox(ax_b1, zoomregion=None, lwd=1, col="black")
    zr2 = (2279000, 2284000)
    plotBedgraph(ax_b2, bg, chrom, zr2[0], zr2[1], colorbycol=SushiColors(5))
    labelgenome(ax_b2, chrom, zr2[0], zr2[1], n=4, scale="Kb", edgeblankfraction=0.20)
    zoombox(ax_b2, zoomregion=None, lwd=1, col="black")
    save(fig, "15_zoom_in_out")


if __name__ == "__main__":
    examples = [
        ("01 DNaseI bedgraph (gradient)", ex01_dnasei_bedgraph),
        ("02 overlay DNaseI + CTCF", ex02_overlay_dnasei_ctcf),
        ("03 flip CTCF + DNaseI", ex03_flip),
        ("04 plotHic + addlegend", ex04_hic),
        ("05 plotBedpe 5C loops", ex05_bedpe_loops),
        ("06 plotBedpe 5C ribbons", ex06_bedpe_ribbons),
        ("07 plotBed Pol2 region", ex07_bed_region),
        ("08 plotBed Pol2 circles", ex08_bed_circles),
        ("09 plotBed Pol2 splitstrand", ex09_bed_splitstrand),
        ("10 plotBed several-factors circles", ex10_severalfactors_circles),
        ("11 plotBed several-factors region", ex11_severalfactors_region),
        ("12 plotManhattan GWAS", ex12_manhattan),
        ("13 plotGenes (genes.bed)", ex13_genes),
        ("14 plotGenes transcripts (colorby)", ex14_transcripts_colorby),
        ("15 zoombox + zoomsregion", ex15_zoom),
    ]
    for name, fn in examples:
        print(f"--- {name} ---")
        run_example(name, fn)
    print()
    print("=== Vignette complete: %d/%d examples passed ===" % (
        sum(1 for _, s in results if s == "OK"), len(results)))
    for n, s in results:
        print(f"  {s[:3]:3s}  {n}")
    print("Output figures in:", OUT)
