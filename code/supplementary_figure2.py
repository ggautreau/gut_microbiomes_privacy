#!/usr/bin/env python3
"""Supplementary Figure 2: the links between the samples of Nielsen 2014 and Li 2014, and why they are not used.

    python3 -I supplementary_figure2.py LINKS.json NIELSEN_LI.json ROWS_V4.json OUTPUT_PREFIX      # writes OUTPUT_PREFIX.png (600 dpi) and .pdf

LINKS.json is the output of nielsen_li_links.py, NIELSEN_LI.json that of nielsen_li_overlap.py (Supplementary Note 1), ROWS_V4.json the donors of the 7 cohorts kept (data/donors.json).
A, rank of the sample of the other cohort with the same subject identifier among the 291,881 runs of the HGM set, when a donor of Nielsen queries the index (35 donors) and when a Li sample queries it (45 samples).
B, similarity to the query of the nearest run of other people, of the donor's other sample in its own cohort and of the sample of the same identifier in the other cohort (35 donors).
C, Bray-Curtis dissimilarity between species profiles. D, share of donors whose other sample is first, as counted and with the identifiers merged. E, what outranks the other sample of the Nielsen donors.
Palette and type follow Supplementary Figure 4: red = cohorts left out, grey = the rest.
"""
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

LINKS, NL, ROWS, PREFIX = sys.argv[1:5]
L = json.load(open(LINKS)); nl = json.load(open(NL)); kept = json.load(open(ROWS))
INK, INK2, INK3, RULE, GUIDE = "#1B2227", "#3F4A52", "#6B777F", "#C9D0D4", "#E3E7EA"
RED, RED_D, RED_L, GREY = "#C8553D", "#8E2F1C", "#E7ADA0", "#9AA5AC"
import os
from matplotlib import font_manager
FONTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fonts")                  # Liberation Sans shipped with the package (licences in fonts/), so that the drawing does not depend on the machine
font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if f.name != "Liberation Sans"]
for fn in sorted(os.listdir(FONTS)):
    if fn.lower().endswith(".ttf"):
        font_manager.fontManager.addfont(os.path.join(FONTS, fn))
plt.rcParams.update({"font.family": ["Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42})
W, H = 7.2, 6.9
fig = plt.figure(figsize=(W, H), dpi=300)


def ax_in(x0, y0, w, h):
    ax = fig.add_axes((x0 / W, y0 / H, w / W, h / H))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK3)
    ax.tick_params(length=2.5, width=0.5, color=INK3, labelcolor=INK2, labelsize=6.3)
    return ax


def title(x, y, letter, text):
    fig.text(x / W, y / H, letter, fontsize=10, fontweight="bold", color=INK, ha="left", va="bottom")
    fig.text((x + 0.2) / W, (y + 0.015) / H, text, fontsize=6.8, color=INK, ha="left", va="bottom")


rng = np.random.default_rng(11)
# ------------------------------------------------------------------ A: rank of the sample of the same identifier in the other cohort
F, R = L["forward"], L["reverse"]


def classes(ranks, n):
    rk = [x for x in ranks if x is not None]
    return [sum(x == 1 for x in rk), sum(2 <= x <= 10 for x in rk), sum(11 <= x <= 100 for x in rk), n - sum(x <= 100 for x in rk)]


bars = [("Nielsen donor queries\n(35 donors, identifier in Li)", classes(F["ranks"], F["with_li_sample"]), F["with_li_sample"]),
        ("Li sample queries\n(45 samples, identifier in Nielsen)", classes(R["ranks"], R["queries"]), R["queries"])]
axa = ax_in(1.75, 5.15, 2.2, 1.0)
cols = [RED_D, RED, RED_L, RULE]; labs = ["Rank 1", "2 to 10", "11 to 100", "Beyond 100\nor not found"]
for i, (lab, c, n) in enumerate(bars):
    y = 1 - i; left = 0
    for k, v in enumerate(c):
        if v:
            axa.barh(y, 100 * v / n, left=left, height=0.58, color=cols[k], linewidth=0, zorder=2)
            axa.text(left + 50 * v / n, y, str(v), ha="center", va="center", fontsize=6.0, color="white" if k < 2 else INK2, zorder=3, clip_on=False)
        left += 100 * v / n
axa.set_xlim(0, 100); axa.set_ylim(-0.55, 1.55); axa.set_yticks([1, 0]); axa.set_yticklabels([b[0] for b in bars], fontsize=6.4, color=INK)
axa.set_xticks([0, 25, 50, 75, 100]); axa.set_xlabel("Share of the queries (%); rank among the 291,881 runs of the HGM set", fontsize=6.4, color=INK, labelpad=2)
axa.grid(axis="x", color=GUIDE, linewidth=0.45, zorder=0); axa.set_axisbelow(True)
for k, (c_, l_) in enumerate(zip(cols, labs)):
    x0 = 0.05 + [0.0, 0.62, 1.34, 2.14][k]
    fig.patches.append(Rectangle((x0 / W, 6.3 / H), 0.09 / W, 0.09 / H, transform=fig.transFigure, facecolor=c_, edgecolor="none"))
    fig.text((x0 + 0.12) / W, 6.3 / H, l_.replace("\n", " "), fontsize=5.6, color=INK2, ha="left", va="bottom")
title(0.05, 6.52, "A", "The other cohort's sample is the nearest run")

# ------------------------------------------------------------------ B: similarity to the query
axb = ax_in(4.8, 5.15, 2.2, 1.0)
bi = np.array(F["best_imp"], float); wp = np.array(F["within_p"], float); cp = np.array([np.nan if v is None else v for v in F["p"]], float)
for a, b in zip(wp, cp):
    if not np.isnan(b):
        axb.plot([1, 2], [a, b], color=RULE, lw=0.6, zorder=1)
for x, v, col in ((0, bi, GREY), (1, wp, INK2), (2, cp, RED)):
    axb.scatter(x + rng.uniform(-0.12, 0.12, len(v)), v, s=6, color=col, alpha=0.75, linewidths=0, zorder=3)
    m = np.nanmedian(v); axb.plot([x - 0.3, x + 0.3], [m, m], color=INK, lw=1.4, solid_capstyle="butt", zorder=4)
axb.set_xlim(-0.5, 2.5); axb.set_ylim(0, 0.72); axb.set_xticks([0, 1, 2]); axb.set_xticklabels(["Nearest run of\nother people", "Own other\nsample (Nielsen)", "Same identifier\n(Li)"], fontsize=5.8)
axb.set_ylabel("Similarity to the query", fontsize=6.6, color=INK, labelpad=2); axb.grid(axis="y", color=GUIDE, linewidth=0.45, zorder=0); axb.set_axisbelow(True)
title(4.15, 6.52, "B", "More similar than the donor's own other sample")

# ------------------------------------------------------------------ C: species profiles
axc = ax_in(1.95, 2.75, 1.9, 1.35)
G = L["species"]["bc_groups"]; wn = L["species"]["within_nielsen_bc"]
rows_c = [("Same identifier,\nNielsen and Li (80 pairs)", "dots", G["matched"], RED), ("Same patient, two\nNielsen samples (78 pairs)", "dots", G["within_same"], INK2),
          ("Other pairs,\nNielsen and Li", "box", G["other_cross_quantiles"], GREY), ("Other pairs,\nwithin Nielsen", "box", G["within_diff_quantiles"], GREY)]
for i, (lab, kind, v, col) in enumerate(rows_c):
    y = 3 - i
    if kind == "dots":
        axc.scatter(v, y + rng.uniform(-0.24, 0.24, len(v)), s=5, color=col, alpha=0.6, linewidths=0, zorder=3)
        axc.plot([np.median(v)] * 2, [y - 0.36, y + 0.36], color=INK, lw=1.4, solid_capstyle="butt", zorder=4)
    else:
        axc.plot([v[0], v[4]], [y, y], color=col, lw=0.9, zorder=2)
        axc.add_patch(Rectangle((v[1], y - 0.22), v[3] - v[1], 0.44, facecolor=col, alpha=0.55, edgecolor="none", zorder=2))
        axc.plot([v[2]] * 2, [y - 0.36, y + 0.36], color=INK, lw=1.4, solid_capstyle="butt", zorder=4)
axc.set_xlim(0, 1); axc.set_ylim(-0.6, 3.6); axc.set_yticks([3, 2, 1, 0]); axc.set_yticklabels([r[0] for r in rows_c], fontsize=5.9, color=INK)
axc.set_xlabel("Bray-Curtis dissimilarity of the species profiles", fontsize=6.6, color=INK, labelpad=2); axc.grid(axis="x", color=GUIDE, linewidth=0.45, zorder=0); axc.set_axisbelow(True)
title(0.05, 4.3, "C", "The species profiles agree, without the sketch")

# ------------------------------------------------------------------ D: what counting does to the rate of the two cohorts
axd = ax_in(4.8, 2.75, 2.2, 1.35)
rates = {}
for r in kept:
    rates.setdefault(r["study"], []).append(r["rank_any"] == 1)
lo, hi = min(100 * np.mean(v) for v in rates.values()), max(100 * np.mean(v) for v in rates.values())
axd.add_patch(Rectangle((-0.55, lo), 3.9, hi - lo, facecolor=GUIDE, edgecolor="none", zorder=0))
axd.text(3.3, (lo + hi) / 2, f"7 cohorts kept\n{lo:.0f} to {hi:.0f}%", ha="right", va="center", fontsize=5.8, color=INK3, zorder=1)
vals = [(0, 100 * nl["nielsen_rank1_before"] / nl["nielsen_donors"], f"{nl['nielsen_rank1_before']}/{nl['nielsen_donors']}", RED_L), (1, 100 * nl["nielsen_rank1_after"] / nl["nielsen_donors"], f"{nl['nielsen_rank1_after']}/{nl['nielsen_donors']}", RED_D),
        (2.2, 100 * nl["li_rank1_before"] / nl["li_donors"], f"{nl['li_rank1_before']}/{nl['li_donors']}", RED_L), (3.2, 100 * nl["li_rank1_after"] / nl["li_donors"], f"{nl['li_rank1_after']}/{nl['li_donors']}", RED_D)]
for x, v, t, c in vals:
    axd.bar(x, v, width=0.78, color=c, linewidth=0, zorder=2)
    axd.text(x, v + 2, t, ha="center", va="bottom", fontsize=6.0, color=INK2, zorder=3)
axd.set_xlim(-0.55, 3.75); axd.set_ylim(0, 112); axd.set_yticks([0, 25, 50, 75, 100]); axd.set_xticks([0, 1, 2.2, 3.2]); axd.set_xticklabels(["as\ncounted", "identifiers\nmerged", "as\ncounted", "identifiers\nmerged"], fontsize=5.8)
axd.tick_params(axis="x", length=0)
for x_, t_ in ((0.5, f"Nielsen ({nl['nielsen_donors']} donors)"), (2.7, f"Li ({nl['li_donors']} donors)")):
    axd.text(x_, -0.27, t_, transform=blended_transform_factory(axd.transData, axd.transAxes), ha="center", va="top", fontsize=6.2, color=INK)
axd.set_ylabel("Other sample first (%)", fontsize=6.6, color=INK, labelpad=2); axd.grid(axis="y", color=GUIDE, linewidth=0.45, zorder=0); axd.set_axisbelow(True)
title(4.15, 4.3, "D", "Counting the other cohort's sample as another person")

# ------------------------------------------------------------------ E: who outranks the other sample of a Nielsen donor
OB = L["outranked_by"]; n_nd = sum(OB.values())
segs = [("first", "Other sample first", INK3), ("li_same_identifier", "Sample of Li, same identifier", RED_D), ("li_other_identifier", "Sample of Li, other identifier", RED),
        ("nielsen_other_person", "Other patient of Nielsen", RED_L), ("other_run", "Other run", RULE), ("not_found", "Not in the 5,000 nearest runs", "white")]
axe = ax_in(1.0, 0.45, 6.0, 0.4)
left = 0
for k, lab, col in segs:
    v = OB.get(k, 0)
    axe.barh(0, v, left=left, height=0.8, color=col, edgecolor=INK3 if col == "white" else "none", linewidth=0.5, zorder=2)
    if v:
        axe.text(left + v / 2, 0, str(v), ha="center", va="center", fontsize=6.4, color="white" if col in (INK3, RED_D, RED) else INK2, zorder=3)
    left += v
axe.set_xlim(0, n_nd); axe.set_ylim(-0.6, 0.6); axe.set_yticks([]); axe.set_xticks([0, 20, 40, 60, n_nd]); axe.set_xlabel(f"Nielsen donors (n = {n_nd}), as counted (a sample of Li is another person)", fontsize=6.6, color=INK, labelpad=2)
axe.spines["left"].set_visible(False)
ys = 1.3
for i, (k, lab, col) in enumerate(segs):
    x0 = 1.0 + [0.0, 2.1, 4.35, 0.0, 2.1, 4.35][i]; y0 = ys - 0.0 if i < 3 else ys - 0.28
    fig.patches.append(Rectangle((x0 / W, y0 / H), 0.09 / W, 0.09 / H, transform=fig.transFigure, facecolor=col, edgecolor=INK3 if col == "white" else "none", linewidth=0.5))
    fig.text((x0 + 0.13) / W, y0 / H, f"{lab} ({OB.get(k, 0)})", fontsize=5.9, color=INK2, ha="left", va="bottom")
title(0.05, 1.85, "E", "What outranks the other sample of a Nielsen donor")
fig.savefig(PREFIX + ".png", dpi=600)
fig.savefig(PREFIX + ".pdf")
print("wrote", PREFIX + ".png and .pdf")
