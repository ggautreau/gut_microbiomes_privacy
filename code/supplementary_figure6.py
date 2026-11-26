#!/usr/bin/env python3
"""Supplementary Figure 6: 95% cohort-then-donor resampling intervals for all three extrapolation models
of Fig. 1B and adaptive bootstrap sample-size selection by doubling B.

    python3 -I supplementary_figure6.py ADULT_RANKS.csv EXTRAPOLATION_BOOTSTRAP.json OUTPUT_PREFIX [EVT_RESAMPLING.json]

Writes OUTPUT_PREFIX.png (600 dpi), OUTPUT_PREFIX.pdf and OUTPUT_PREFIX_values.json.
EVT_RESAMPLING.json  output of evt_resampling.py. When it is given, a fourth panel shows the extreme-value model of the similarity tail (the tail drawn in Fig. 1B), and the
                     two lower panels show its bounds as B doubles, as for the three other models. On each resample the distribution of the shapes of the tails is
                     estimated again from the donors drawn, and the resample averages the resulting probabilities of these donors.
"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np

CSV_IN, BOOT_JSON, PREFIX = sys.argv[1:4]
EVT_JSON = sys.argv[4] if len(sys.argv) > 4 else None
os.makedirs(os.path.dirname(os.path.abspath(PREFIX)), exist_ok=True)

FONTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fonts")
if os.path.isdir(FONTS):
    font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if f.name != "Liberation Sans"]
    for fn in sorted(os.listdir(FONTS)):
        if fn.lower().endswith(".ttf"):
            font_manager.fontManager.addfont(os.path.join(FONTS, fn))

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.4,
    "ytick.major.size": 2.4,
    "xtick.minor.size": 1.4,
    "xtick.minor.width": 0.45,
    "xtick.major.pad": 2.0,
    "ytick.major.pad": 2.0,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Liberation Sans",
    "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold",
    "mathtext.default": "regular",
})

INK, INK2, INK3, RULE, GUIDE = "#14171A", "#3A4047", "#5C646C", "#9AA2AA", "#DCE0E4"
C_SK, C_SK_TINT = "#D55E00", "#F9D9BF"
C_PY, C_PY_INK = "#009988", "#00665B"
C_LN, C_LN_INK = "#882255", "#882255"
C_MS, C_MS_INK = "#A469F7", "#5A35A8"
C_BL, C_BL_INK = "#2F5FA8", "#1D3F75"
C_EV, C_EV_INK = "#111111", "#111111"                                           # the extreme-value model, as in Fig. 1B

bt = json.load(open(BOOT_JSON))
meta = bt["meta"]
M = meta["M"]
WORLD = meta["world"]
FIT_FROM = meta["fit_from"]
X_MAX = 1e12
THR_SE_W = meta["se_over_width_threshold"]

meas = bt["measured"]
Ns = np.asarray(meas["N"], float)
R_meas = np.asarray(meas["est"], float)
R_lo = np.asarray(meas["lo"], float)
R_hi = np.asarray(meas["hi"], float)
X_EXT = np.asarray(bt["x_ext"], float)

MODEL_SPECS = [
    ("pycmb", "Pitman–Yor correctness", C_PY, C_PY_INK, 0.93),
    ("line", "Straight line in log N", C_LN, C_LN_INK, 1.00),
    ("miss", "Power law of the miss rate", C_MS, C_MS_INK, 1.07),
]
MODEL_SPECS.sort(key=lambda m: -bt["models"][m[0]]["at_N"][str(int(WORLD))]["est"])
MODELS = {m[0]: bt["models"][m[0]] for m in MODEL_SPECS}                        # the data of each model and the grid of pool sizes of its curves
GRIDS = {m[0]: X_EXT for m in MODEL_SPECS}
if EVT_JSON:                                                                    # the extreme-value model comes last, as in the legend of Fig. 1B
    er = json.load(open(EVT_JSON))
    assert er["meta"]["M"] == M and er["meta"]["world"] == WORLD and er["meta"]["se_over_width_threshold"] == THR_SE_W
    MODEL_SPECS.append(("evt", "Extreme-value model", C_EV, C_EV_INK, 1.14))
    MODELS["evt"] = er["models"]["evt"]
    GRIDS["evt"] = np.asarray(er["x_ext"], float)
NTOP = len(MODEL_SPECS)
FS_TITLE, FS_BOX = (7.0, 5.8) if NTOP == 3 else (6.5, 5.4)                       # four panels are narrower


def style_ax(ax, grid="y"):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(RULE)
    if grid:
        ax.grid(axis=grid, color=GUIDE, linewidth=0.45, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=6.0, colors=INK2)


def add_risk_zone(ax):
    """Blue re-identification threshold area (9% solid, gradient up to 33%), matching Fig. 1B."""
    a_max, a_33, a_20, y_top = 0.28, 0.05, 0.20, 33.0
    yy = np.linspace(y_top, 9.0, 200)
    img = np.zeros((len(yy), 1, 4))
    img[..., :3] = matplotlib.colors.to_rgb(C_BL)
    img[:, 0, 3] = np.where(
        yy > 20.0,
        a_33 + (a_20 - a_33) * (33.0 - yy) / 13.0,
        a_20 + (a_max - a_20) * (20.0 - yy) / 11.0,
    )
    ax.imshow(img, extent=(2, X_MAX, 9.0, y_top), origin="upper", aspect="auto", interpolation="bilinear", zorder=0.5)
    ax.axhspan(0, 9, facecolor=C_BL, alpha=a_max, lw=0, zorder=0.5)
    ax.axhline(9, color=C_BL, lw=0.75, zorder=2.5)
    ax.axhline(33, color=C_BL, lw=0.35, zorder=2.5)


fig = plt.figure(figsize=(7.2, 5.15), dpi=300)
gs_top = fig.add_gridspec(1, NTOP, left=0.075, right=0.985, top=0.905, bottom=0.495, wspace=0.08)
gs_bot = fig.add_gridspec(1, 2, left=0.075, right=0.985, top=0.355, bottom=0.095, wspace=0.27)

WK = str(int(WORLD))
letters_top = ["A", "B", "C", "D"]
letters_bot = ("D", "E") if NTOP == 3 else ("E", "F")
for idx, (key, title, col, ink, _) in enumerate(MODEL_SPECS):
    ax = fig.add_subplot(gs_top[0, idx])
    style_ax(ax, "y")
    ax.set_xscale("log")
    add_risk_zone(ax)

    for xv in (M, WORLD):
        ax.plot([xv, xv], [0, 102], color=INK3, lw=0.45, ls=(0, (1.5, 2.0)), zorder=2)

    # Measured curve + 95% band (paler below FIT_FROM, darker on [FIT_FROM, M])
    below = Ns <= FIT_FROM
    above = Ns >= FIT_FROM
    ax.fill_between(Ns, R_lo, R_hi, color=C_SK_TINT, alpha=0.55, lw=0, zorder=1.5)
    if key == "evt":                                                             # this model is not fitted to the measured curve: the whole curve is pale
        ax.plot(Ns, R_meas, color=C_SK, lw=1.2, alpha=0.35, zorder=4, solid_capstyle="round")
    else:
        ax.plot(Ns[below], R_meas[below], color=C_SK, lw=1.2, alpha=0.35, zorder=4, solid_capstyle="round")
        ax.plot(Ns[above], R_meas[above], color=C_SK, lw=1.4, alpha=0.90, zorder=4.2, solid_capstyle="round")

    mdata = MODELS[key]
    xg = GRIDS[key]
    pt = np.asarray(mdata["point_curve"], float)
    lo_c = np.asarray(mdata["lo_curve"], float)
    hi_c = np.asarray(mdata["hi_curve"], float)
    ax.fill_between(xg, lo_c, hi_c, color=col, alpha=0.20 if key != "evt" else 0.13, lw=0, zorder=2.8)
    ax.plot(xg, lo_c, color=col, lw=0.55, alpha=0.75, zorder=3.5)
    ax.plot(xg, hi_c, color=col, lw=0.55, alpha=0.75, zorder=3.5)
    ax.plot(xg, pt, color=col, lw=1.25, ls=(0, (3.2, 1.5)), zorder=5)

    wstat = mdata["at_N"][WK]
    ax.plot([WORLD, WORLD], [wstat["lo"], wstat["hi"]], color=ink, lw=1.1, zorder=6)
    ax.plot([WORLD * 0.65, WORLD * 1.55], [wstat["lo"], wstat["lo"]], color=ink, lw=0.9, zorder=6)
    ax.plot([WORLD * 0.65, WORLD * 1.55], [wstat["hi"], wstat["hi"]], color=ink, lw=0.9, zorder=6)
    ax.plot([WORLD], [wstat["est"]], "o", ms=3.5, color=ink, mec="white", mew=0.5, zorder=7)

    # Summary box cleanly placed in open space between 33% line and 70% lower CI bound
    box_txt = (
        f"At 8×10$^9$: {wstat['est']:.1f}%\n"
        f"95% interval: {wstat['lo']:.1f}–{wstat['hi']:.1f}%\n"
        f"B = {mdata['B']:,} resamples"
    )
    ax.text(
        1.5e6, 51.0, box_txt,
        fontsize=FS_BOX, color=INK, ha="center", va="center", linespacing=1.28, zorder=8,
        bbox=dict(facecolor="white", edgecolor=col, lw=0.55, alpha=0.92, boxstyle="round,pad=0.30"),
    )

    ax.set_xlim(2, X_MAX)
    ax.set_ylim(0, 102)
    ax.set_xticks([10.0 ** e for e in range(1, 13)])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xticklabels([rf"$10^{{{e}}}$" if e in (1, 4, 7, 10) else "" for e in range(1, 13)], fontsize=5.8)
    ax.set_yticks(range(0, 101, 20))
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
    if idx == 0:
        ax.set_yticklabels([str(y) for y in range(0, 101, 20)], fontsize=6.0)
        ax.set_ylabel("Donors recovered (%)", fontsize=6.8, color=INK, labelpad=3)
        ax.text(3.2, 10.6, "9% risk", fontsize=5.3, color=C_BL_INK, fontweight="bold", ha="left", va="bottom", zorder=6)
        ax.text(3.2, 21.0, "20%", fontsize=5.1, color=C_BL_INK, ha="left", va="bottom", zorder=6)
        ax.text(3.2, 34.0, "33%", fontsize=5.1, color=C_BL_INK, ha="left", va="bottom", zorder=6,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
    else:
        ax.set_yticklabels([])

    ax.set_title(title, fontsize=FS_TITLE, fontweight="bold", color=ink, pad=14)
    ax.text(M, 103.5, "292k HGMs", fontsize=5.1, color=INK3, ha="center", va="bottom", clip_on=False)
    ax.text(WORLD, 103.5, "8×10$^9$", fontsize=5.1, color=INK3, ha="center", va="bottom", clip_on=False)
    ax.text(0.0, 1.07, letters_top[idx], transform=ax.transAxes, fontsize=8.8, fontweight="bold", color=INK, ha="left", va="bottom", clip_on=False)

fig.text(0.53, 0.435, "Pool size, N", fontsize=6.8, color=INK, ha="center", va="center")

# ---- Panel D: Adaptive stopping criterion (max MC SE / width) vs B
axD = fig.add_subplot(gs_bot[0, 0])
style_ax(axD, "y")
axD.set_xscale("log", base=2)

all_Bs = sorted({step["B"] for mkey, _, _, _, _ in MODEL_SPECS for step in MODELS[mkey]["doubling_trace"]})
max_ratio = max(step["max_se_over_width"] for mkey, _, _, _, _ in MODEL_SPECS for step in MODELS[mkey]["doubling_trace"])

for key, title, col, ink, x_jit in MODEL_SPECS:
    trace = MODELS[key]["doubling_trace"]
    bs = np.asarray([s["B"] for s in trace], float) * x_jit
    rs = [s["max_se_over_width"] for s in trace]
    axD.plot(bs, rs, color=col, lw=1.15, marker="o", ms=3.0, label=title, zorder=4)
    # Highlight stopping B with open white ring + inner dot
    axD.plot([bs[-1]], [rs[-1]], marker="o", ms=6.2, mfc="white", mec=ink, mew=1.1, zorder=5)
    axD.plot([bs[-1]], [rs[-1]], marker="o", ms=2.4, mfc=ink, mec="none", zorder=6)

axD.axhline(THR_SE_W, color=INK3, lw=0.8, ls=(0, (3.0, 2.0)), zorder=2.5)
axD.text(all_Bs[0] * 1.02, THR_SE_W - 0.008, f"{THR_SE_W:.2f} stopping criterion", fontsize=5.6, color=INK2, ha="left", va="top", zorder=6)
axD.set_xticks(all_Bs)
axD.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
axD.set_xticklabels([f"{b:,}" for b in all_Bs], fontsize=5.7)
axD.set_xlim(all_Bs[0] * 0.80, all_Bs[-1] * 1.25)
axD.set_ylim(0, max(max_ratio * 1.15, THR_SE_W * 2.2))
axD.set_xlabel("Number of resamples, B", fontsize=6.8, color=INK, labelpad=2.5)
axD.set_ylabel("Max MC SE / interval width", fontsize=6.8, color=INK, labelpad=3)
axD.set_title("Precision of the interval bounds vs. B", fontsize=7.0, fontweight="bold", color=INK, pad=6)
axD.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none", framealpha=0.85, fontsize=5.7, labelcolor=INK2, handlelength=1.6, labelspacing=0.35)
axD.text(-0.08, 1.04, letters_bot[0], transform=axD.transAxes, fontsize=8.8, fontweight="bold", color=INK, ha="left", va="bottom", clip_on=False)

# ---- Panel E: Lower and upper 95% CI bounds at 8e9 vs B with +/-1 MC SE error bars
axE = fig.add_subplot(gs_bot[0, 1])
style_ax(axE, "y")
axE.set_xscale("log", base=2)

min_lo = min(s["lo"] - s["se_lo"] for mkey, _, _, _, _ in MODEL_SPECS for s in MODELS[mkey]["doubling_trace"])
y_min_E = max(0.0, 5.0 * math.floor((min_lo - 4.0) / 5.0))

for key, title, col, ink, x_jit in MODEL_SPECS:
    trace = MODELS[key]["doubling_trace"]
    bx = np.asarray([s["B"] for s in trace], float) * x_jit
    lo_v = [s["lo"] for s in trace]
    hi_v = [s["hi"] for s in trace]
    se_l = [s["se_lo"] for s in trace]
    se_h = [s["se_hi"] for s in trace]
    axE.errorbar(bx, lo_v, yerr=se_l, color=col, lw=1.05, marker="o", ms=2.8, elinewidth=0.7, capsize=1.6, zorder=4)
    axE.errorbar(bx, hi_v, yerr=se_h, color=col, lw=0.95, ls=(0, (2.6, 1.4)), marker="^", ms=2.8, elinewidth=0.7, capsize=1.6, zorder=4)

axE.set_xticks(all_Bs)
axE.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
axE.set_xticklabels([f"{b:,}" for b in all_Bs], fontsize=5.7)
axE.set_xlim(all_Bs[0] * 0.80, all_Bs[-1] * 1.25)
axE.set_ylim(y_min_E, 102)
axE.set_xlabel("Number of resamples, B", fontsize=6.8, color=INK, labelpad=2.5)
axE.set_ylabel("95% interval bounds at 8×10$^9$ (%)", fontsize=6.8, color=INK, labelpad=3)
axE.set_title("95% interval bounds at 8×10$^9$ vs. B", fontsize=7.0, fontweight="bold", color=INK, pad=6)

handles_E = [
    Line2D([], [], color=INK2, lw=0.95, ls=(0, (2.6, 1.4)), marker="^", ms=3.0, label="97.5th pct. (±1 MC SE)"),
    Line2D([], [], color=INK2, lw=1.05, marker="o", ms=3.0, label="2.5th pct. (±1 MC SE)"),
]
axE.legend(handles=handles_E, loc="lower right", frameon=True, facecolor="white", edgecolor="none", framealpha=0.88, fontsize=5.6, labelcolor=INK2, handlelength=1.8, labelspacing=0.35)
axE.text(-0.08, 1.04, letters_bot[1], transform=axE.transAxes, fontsize=8.8, fontweight="bold", color=INK, ha="left", va="bottom", clip_on=False)

fig.savefig(PREFIX + ".png", dpi=600)
fig.savefig(PREFIX + ".pdf")

values_out = {
    "meta": meta,
    "models": {
        key: {
            "B": bt["models"][key]["B"],
            "converged": bt["models"][key]["converged"],
            "n_failed": bt["models"][key]["n_failed"],
            "n_warnings": bt["models"][key]["n_warnings"],
            "at_N": bt["models"][key]["at_N"],
            "doubling_trace": bt["models"][key]["doubling_trace"],
        }
        for key in ("pycmb", "line", "miss")
    },
}
if EVT_JSON:                                                                    # the extreme-value model: its interval and its bounds as B doubles
    values_out["models"]["evt"] = {k_: MODELS["evt"][k_] for k_ in ("B", "converged", "at_N", "doubling_trace")}
with open(PREFIX + "_values.json", "w", encoding="utf-8") as fh:
    json.dump(values_out, fh, indent=1)
print("wrote", PREFIX + ".png, .pdf and _values.json")
