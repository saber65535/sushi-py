"""Example 4: Hi-C lower-triangle heatmap with color legend.
Output: /tmp/hic.png
"""
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
import sushi

hic = sushi.data.Sushi_HiC_matrix()
pos = hic["positions"].astype(int)
m = hic["matrix"]
full = np.full((m.shape[0], m.shape[0]), np.nan)
for i in range(m.shape[0]):
    for j in range(i, m.shape[1]):
        full[i, j] = m[i, j]
hic_df = pd.DataFrame(full, index=pos, columns=pos)

fig, ax = plt.subplots(figsize=(10, 6))
fig.subplots_adjust(bottom=0.22, right=0.82)
phic = sushi.plotHic(ax, hic_df, "chr11", 500000, 5050000,
                     max_y=20, zrange=(0, 28), palette=sushi.SushiColors(7))
sushi.addlegend(ax, range_val=phic[0], palette=phic[1], title="score",
                side="right", bottominset=0.4, topinset=0, xoffset=-0.035,
                labelside="left", width=0.025, title_offset=0.035)
sushi.labelgenome(ax, "chr11", 500000, 5050000, n=4, scale="Mb",
                  edgeblankfraction=0.20, side=1)
fig.savefig("/tmp/hic.png", dpi=120, bbox_inches="tight")
