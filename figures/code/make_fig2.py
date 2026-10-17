#!/usr/bin/env python3
"""
Figure 2 — fit and separation for the 41 classifiable audit records.

Portable: every path is either a command-line argument or relative to this file.
Plotted values are computed exactly as in the previous generator (three-signature log
distances from sigma_obs and the case fraction); only lettering, legend text, size and
export formats changed.

Inputs (defaults, relative to this file):
  inputs/repository_scale_reanalysis.csv   frozen audit results (repository: results/)
  inputs/sim_tau.txt                        adequacy threshold tau (simulation layer 1)
Outputs (in --outdir, default: this file's directory):
  fig2_margins.eps                          vector
  fig2_margins.tiff                         1200 dpi, LZW-compressed raster
  fig2_margins.png                          600 dpi raster (embedded in the manuscript DOCX)
  fig2_plotted_values.csv                   the exact x/y values drawn, for verification

  python make_fig2.py [--reanalysis PATH] [--tau PATH] [--outdir DIR]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--reanalysis", default=os.path.join(HERE, "inputs", "repository_scale_reanalysis.csv"))
ap.add_argument("--tau", default=os.path.join(HERE, "inputs", "sim_tau.txt"))
ap.add_argument("--outdir", default=HERE)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

# ---------------------------------------------------------------- values (unchanged computation)
d = pd.read_csv(a.reanalysis)
exp_lin = np.sqrt(d.prev * (1 - d.prev))
exp_lo = 1 / np.sqrt(d.prev * (1 - d.prev))
D = np.c_[np.abs(np.log(d.sigma_obs / exp_lin)), np.abs(np.log(d.sigma_obs / 1.0)),
          np.abs(np.log(d.sigma_obs / exp_lo))]
names = np.array(["linear_probability", "standardised", "log_odds"])
o = np.argsort(D, axis=1)
d["assigned_scale"] = names[o[:, 0]]
d["distance_to_assigned"] = np.take_along_axis(D, o[:, :1], 1).ravel()
sep = np.take_along_axis(D, o[:, 1:2], 1).ravel() - d.distance_to_assigned
d["relative_separation"] = sep / np.abs(0.5 * np.log(d.prev * (1 - d.prev)))
d["disc"] = d.metadata_discordant.astype(str).isin(["True", "1", "1.0"])
# Deterministic row order: records at the attainable ceiling tie to floating-point precision, so
# order by separation rounded to 10 decimals, then by accession (independent of input row order).
d["_key"] = d.relative_separation.round(10)
g = d.sort_values(["_key", "id"], kind="mergesort").reset_index(drop=True)
TAU = float(open(a.tau).read().split("=")[-1])
g[["id", "assigned_scale", "disc", "distance_to_assigned", "relative_separation"]].assign(
    y_position=g.index).to_csv(os.path.join(a.outdir, "fig2_plotted_values.csv"), index=False)

# ---------------------------------------------------------------- style: three sizes, floor 8 pt
BASE, MID, SMALL = 9, 8.5, 8
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": BASE, "axes.titlesize": BASE,
                     "axes.labelsize": BASE, "xtick.labelsize": SMALL, "ytick.labelsize": SMALL,
                     "legend.fontsize": MID, "axes.spines.top": False, "axes.spines.right": False,
                     "ps.fonttype": 42, "pdf.fonttype": 42})
GREY = "#6b6b6b"
MK = {"linear_probability": ("o", "#1f5fa8", "Linear probability"),
      "log_odds": ("^", "#b8860b", "Log odds")}

# 174 mm = 6.85 in (double-column width); height within a 234 mm page
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.85, 8.7), sharey=True, gridspec_kw={"wspace": 0.06})
for sc, (mkr, col, lab) in MK.items():
    for ax, xcol in [(ax1, "distance_to_assigned"), (ax2, "relative_separation")]:
        s = g[(g.assigned_scale == sc) & (~g.disc)]
        ax.scatter(s[xcol], s.index, marker=mkr, s=34, facecolors=col, edgecolors="white", lw=.5,
                   zorder=3, label=lab if ax is ax1 else None)
        s2 = g[(g.assigned_scale == sc) & (g.disc)]
        ax.scatter(s2[xcol], s2.index, marker=mkr, s=60, facecolors="none", edgecolors="#c1121f", lw=1.3,
                   zorder=4, label=("Declared SD;\nindicated linear probability"
                                    if (ax is ax1 and sc == "linear_probability") else None))
ax1.axvline(TAU, color=GREY, lw=.8, ls="--", zorder=1)
ax1.text(TAU, len(g) - 1, " τ (post hoc)", fontsize=SMALL, color=GREY, va="center", ha="left")
ax2.axvline(1.0, color=GREY, lw=.8, ls=":", zorder=1)
ax1.set_xlabel("Log distance to assigned signature")
ax2.set_xlabel("Separation ÷ attainable ceiling")
ax1.set_title("Fit to the assigned signature", loc="left")
ax2.set_title("Separation relative to its ceiling", loc="left")
ax1.set_yticks(range(len(g)))
ax1.set_yticklabels(g.id)
for ax in (ax1, ax2):
    ax.margins(y=0.015, x=0.10)
ax1.set_xticks([0.0, 0.1, 0.2])
ax2.set_xticks([0.7, 0.8, 0.9, 1.0])
h, l = ax1.get_legend_handles_labels()
ax2.legend(h, l, frameon=False, loc="center left", handletextpad=0.4, borderaxespad=0.2)
fig.subplots_adjust(left=0.235, right=0.985, top=0.955, bottom=0.065, wspace=0.06)

# ---------------------------------------------------------------- render-then-verify (geometric)
fig.canvas.draw()
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
small = [t.get_text() for t, _ in texts if t.get_fontsize() < 8]
over = [(x.get_text(), y.get_text()) for i, (x, bx) in enumerate(texts)
        for y, by in texts[i + 1:] if bx.overlaps(by)]
outside = [t.get_text() for t, b in texts
           if not (fig.bbox.contains(b.x0, b.y0) and fig.bbox.contains(b.x1, b.y1))]
assert not small, small
assert not over, over
assert not outside, outside

fig.savefig(os.path.join(a.outdir, "fig2_margins.eps"), format="eps")
fig.savefig(os.path.join(a.outdir, "fig2_margins.png"), dpi=600, facecolor="white")
tif = os.path.join(a.outdir, "fig2_margins.tiff")
fig.savefig(tif, dpi=1200, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
from PIL import Image  # flatten to RGB (no alpha channel) for production
Image.open(tif).convert("RGB").save(tif, compression="tiff_lzw", dpi=(1200, 1200))
print("records", len(g), "| discordant", int(g.disc.sum()), "| tau", TAU,
      "| min font", min(t.get_fontsize() for t, _ in texts))
