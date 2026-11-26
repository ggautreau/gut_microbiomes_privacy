#!/usr/bin/env python3
"""Supplementary Figure 5: the three models of panel B of Figure 1 fitted for five windows of the measured curve, and the extreme-value model of the similarity tail for its five tails.

    python3 -I supplementary_figure5.py ADULT_RANKS.csv EXTRAPOLATION.json OUTPUT_PREFIX [EVT.json]      # writes OUTPUT_PREFIX.png (600 dpi), .pdf and _values.json (aggregates: share at 8e9, N9, errors of the thinning test, per fit)

ADULT_RANKS.csv  per-donor table of Figure 1 (rank_any_narrow); EXTRAPOLATION.json  output of extrapolation_tables.py (errors of the thinning test of the fifteen fits, key `holdout_windows`; the windows drawn in Figure 1B, key `holdout`).
Rows, the three models of Figure 1B in the order of its legend, the highest value at the population of the world first: the Pitman-Yor correctness model of Rocher et al. 2025 (PYC-MB), the straight line in log N and the power
law of the miss rate 1 - a N^b (Grother et al. 2019). Columns, the fit window, N in [lo, 291,881] for lo = 2, 10, 100, 1,000 and 10,000, least squares on the grid of pool sizes of Figure 1B (70 log-spaced values from 2 to 291,881 rounded to whole numbers, 68 distinct sizes).
Each panel shows the measured share of donors whose other sample is first in a pool of N runs drawn from the 291,881 runs of the HGM set (vermillion as in Figure 1, darker over
the window), the fit (dashed, in the colour of its curve in Figure 1B), its value at the population of the world and its root mean square error on thinned HGM sets at 100 times (in
parentheses). The same fit made on the upper bound of the measured share in small pools is computed for the values file (Supplementary Table 16) and not drawn. The panel of the fit drawn in Figure 1B is shaded. Shares are clipped to 0-100%; the vertical axis starts at 50%.
EVT.json  output of evt_model.py (aggregates/evt_summary.json). When it is given, a last row shows the extreme-value model of the similarity tail of Figure 1B, one panel per tail (the 100 to the 10,000 highest
similarities of the runs of other people to each query). This model is not fitted to the measured curve: it starts from the measured share at the size of the HGM set. Its panels show the value at the population
of the world without an error (the test of this model on thinned HGM sets is in Supplementary Table 18), and the values file does not change.
"""
import csv, json, math, os, sys
import numpy as np
from scipy.special import gammaln
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                  # python -I does not put the script's folder on the path
import pyc_mb as PY                                                             # Pitman-Yor correctness model of Rocher et al. 2025
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.offsetbox import AnnotationBbox, DrawingArea, HPacker, TextArea, VPacker

CSV_IN, EXT, PREFIX = sys.argv[1:4]
EVT = sys.argv[4] if len(sys.argv) > 4 else None
M, WORLD, X_MAX = 291_881, 8e9, 1e12
INK, INK2, INK3, RULE, GUIDE = "#14171A", "#3A4047", "#5C646C", "#9AA2AA", "#DCE0E4"
C_SK = "#D55E00"                                                                 # the data, as in Figure 1: vermillion
import os
from matplotlib import font_manager
FONTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fonts")                  # Liberation Sans shipped with the package (licences in fonts/), so that the drawing does not depend on the machine
font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if f.name != "Liberation Sans"]
for fn in sorted(os.listdir(FONTS)):
    if fn.lower().endswith(".ttf"):
        font_manager.fontManager.addfont(os.path.join(FONTS, fn))
plt.rcParams.update({"font.family": ["Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42, "mathtext.default": "regular"})

rows = list(csv.DictReader(open(CSV_IN)))
D = len(rows)
k = np.array([int(r["rank_any_narrow"]) - 1 if r["rank_any_narrow"] != "" else -1 for r in rows]); found = k >= 0
ext = json.load(open(EXT)); HW = ext["holdout_windows"]
FS = ("0.1", "0.03", "0.01")


def first_hyper(N):
    kk = k[found].astype(float); m = N - 1
    with np.errstate(all="ignore"):
        lg = gammaln(M - kk + 1) - gammaln(M - kk - m + 1) - (gammaln(M + 1) - gammaln(M - m + 1))
    return float(np.where(M - kk >= m, np.exp(lg), 0.0).sum() / D)


NS = np.unique(np.round(np.logspace(np.log10(2), np.log10(M), 70)).astype(int))      # the grid of Figure 1B
R1 = np.array([first_hyper(N) for N in NS])
WIN = (2, 10, 100, 1000, 10_000)
SP = ext["small_pools"]; cv = SP["curves"]                                      # the measured share and its upper bound on the grid (extrapolation_tables.py)
assert cv["N"] == [int(n_) for n_ in NS] and np.allclose(cv["measured"], R1, rtol=0, atol=1e-12), "the curves of extrapolation_tables.py are the measured curve of the ranks"
R_UP = np.array(cv["upper"])                                                    # upper bound: a donor whose other sample is not returned has only the runs of its list above that sample


def fit(key, lo, R=None):
    """The model fitted on the measured curve (or on the curve R) for N in [lo, M], as a function of N (share, not clipped)."""
    s = NS >= lo; x = NS[s].astype(float); y = (R1 if R is None else R)[s]
    if key == "line":
        c = np.polyfit(np.log2(x), y, 1); return lambda N: np.polyval(c, np.log2(np.asarray(N, float)))
    if key == "pycmb":
        f = PY.fit(x, y); return lambda N: np.atleast_1d(PY.kappa(np.asarray(N, float), f))
    c = np.polyfit(np.log(x), np.log(np.clip(1 - y, 1e-12, None)), 1)            # the miss rate 1 - R = a N^b, least squares on the logarithms
    return lambda N: 1 - np.exp(np.polyval(c, np.log(np.asarray(N, float))))


X_FINE = np.logspace(np.log10(M), 40, 6000)


def first_crossing(f1):
    """N at which the clipped curve falls to 9% (log-linear interpolation), or inf when it does not within 10^40."""
    y = 100 * np.clip(np.atleast_1d(f1(X_FINE)), 0, 1); y = np.where(np.isfinite(y), y, 100.0)
    i = np.where(y <= 9)[0]
    if len(i) == 0:
        return float("inf")
    j = i[0]
    if j == 0:
        return float(X_FINE[0])
    t = (y[j - 1] - 9) / (y[j - 1] - y[j])
    return float(10 ** (np.log10(X_FINE[j - 1]) + t * (np.log10(X_FINE[j]) - np.log10(X_FINE[j - 1]))))


rmse = lambda key, lo, f="0.01": math.hypot(HW[f]["error"][f"{key}|{lo}"], HW[f]["sd"][f"{key}|{lo}"])      # root mean square error over the thinnings, in points

# rows in the order of the legend of Figure 1B: name, key, colour of the curve, colour of the text, window drawn in Figure 1B
DRAWN_W = ext["holdout"]["0.01"]["windows"]                                     # the window of each fit drawn in Figure 1B: N >= 1,000 for the three models (Supplementary Note 3)
MODELS = [("Pitman–Yor\ncorrectness", "pycmb", "#009988", "#00665B", DRAWN_W["pycmb"]),
          ("Straight line\nin log N", "line", "#882255", "#882255", DRAWN_W["line"]),
          ("Power law of\nthe miss rate", "miss", "#A469F7", "#5A35A8", DRAWN_W["miss"])]
MODELS.sort(key=lambda m_: -float(np.clip(np.atleast_1d(fit(m_[1], m_[4])(WORLD))[0], 0, 1)))     # the highest value at 8e9 first, as in the legend of Figure 1B
assert all(m_[4] == 1000 for m_ in MODELS), "the three fits drawn in Figure 1B are those for N >= 1,000, where the measured share and its upper bound differ by at most 0.2 points"

EV = json.load(open(EVT)) if EVT else None                                      # the extreme-value model of the similarity tail (evt_model.py): a last row, one panel per tail
H3 = 1.54 + 1.22 * len(MODELS)                                                  # height with the rows of the fitted models alone: 1.22 in per row and room for a legend of three lines
ROW = (H3 - 0.563 - 0.98) / (len(MODELS) + 0.16 * (len(MODELS) - 1))            # height of one panel, in inches
GAP = 0.50                                                                      # between the rows of the fitted models and the last row, room for the headers of that row
fig = plt.figure(figsize=(7.2, H3 + (GAP + ROW if EV else 0)), dpi=300)
TOP, BOT = 1 - 0.563 / fig.get_figheight(), 0.98 / fig.get_figheight()          # the same margins in inches with and without the last row
BOT3 = BOT + ((GAP + ROW) / fig.get_figheight() if EV else 0)                   # bottom of the rows of the fitted models
gs = fig.add_gridspec(len(MODELS), 5, left=0.115, right=0.985, top=TOP, bottom=BOT3, wspace=0.10, hspace=0.16)


def frame(ax):
    """Axes, guide lines and dotted lines common to all the panels."""
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(RULE)
    ax.set_xscale("log")
    for yy in (60, 70, 80, 90, 100):
        ax.plot([2, X_MAX], [yy, yy], color=GUIDE, lw=0.4, zorder=0)
    for xv in (M, WORLD):
        ax.plot([xv, xv], [50, 102], color=INK3, lw=0.45, ls=(0, (1.5, 2.0)), zorder=2)


def ticks(ax, j, last):
    """Limits and ticks; the pool sizes are written under the last row and the shares at the left of the first column."""
    ax.set_xlim(2, X_MAX); ax.set_ylim(50, 102)
    ax.set_xticks([10.0 ** e for e in range(1, 13)]); ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xticklabels([r"$10^{%d}$" % e if e in (1, 4, 7, 10) else "" for e in range(1, 13)], fontsize=5.6)
    ax.set_yticks((50, 60, 70, 80, 90, 100)); ax.set_yticklabels([str(t) for t in (50, 60, 70, 80, 90, 100)], fontsize=5.6)
    ax.tick_params(length=2.0, width=0.5, colors=INK2, pad=1.4)
    if j:
        ax.set_yticklabels([])
    if not last:
        ax.set_xticklabels([])


def value(ax, v, text, ink, bold):
    """The value at the population of the world, under the cross."""
    ax.plot([WORLD], [v], "x", ms=4.2, mew=0.9, color=ink, zorder=6)
    row1 = TextArea(text, textprops=dict(fontsize=5.9, color=ink, fontweight="bold" if bold else "normal"))
    ax.add_artist(AnnotationBbox(row1, (WORLD / 1.6, v - 2.2), xycoords="data", box_alignment=(1, 1), frameon=True, pad=0.06,
                                 bboxprops=dict(facecolor="white", edgecolor="none", boxstyle="square,pad=0.0"), zorder=7))


def tint(ax, col, alpha=0.07):
    """The panel of a curve drawn in Figure 1B, tinted in the colour of that curve."""
    ax.set_facecolor(matplotlib.colors.to_rgba(col, alpha))
    for s_ in ("left", "bottom"):
        ax.spines[s_].set_color(col); ax.spines[s_].set_linewidth(1.0)


grid, models = {}, {}
for i, (name, key, col, ink, drawn) in enumerate(MODELS):
    for j, lo in enumerate(WIN):
        ax = fig.add_subplot(gs[i, j])
        f = fit(key, lo)
        v = 100 * float(np.clip(np.atleast_1d(f(WORLD))[0], 0, 1))
        fu = fit(key, lo, R_UP); vu = 100 * float(np.atleast_1d(fu(WORLD))[0])       # the same fit on the upper bound (not clipped: a miss rate above 1 gives a negative share)
        for c_, f_ in (("measured", f), ("upper", fu)):                              # the same fits as those of extrapolation_tables.py
            assert abs(100 * float(np.atleast_1d(f_(WORLD))[0]) - SP["at_8e9_windows"][c_][f"{key}|{lo}"]) < 1e-6, (key, lo, c_)
        grid[f"{key}|{lo}"] = dict(first_8e9=v, first_8e9_upper=vu, N9=first_crossing(f), holdout={fs: HW[fs]["error"][f"{key}|{lo}"] for fs in FS}, rmse_100=rmse(key, lo))
        if lo == drawn:
            models[key] = dict(grid[f"{key}|{lo}"], window=lo, residual_whole=100 * math.sqrt(float(np.mean((np.atleast_1d(f(NS.astype(float))) - R1) ** 2))))   # residual over the whole measured curve, in points
        frame(ax)
        inside = NS >= lo
        ax.plot(NS, 100 * R1, color=C_SK, lw=1.5, alpha=0.22, zorder=3, solid_capstyle="round")                       # the measured curve
        ax.plot(NS[inside], 100 * R1[inside], color=C_SK, lw=1.5, alpha=0.75, zorder=3.1, solid_capstyle="round")      # its part used by the fit
        X = np.logspace(np.log10(lo), np.log10(X_MAX), 400)
        ax.plot(X, 100 * np.clip(np.atleast_1d(f(X)), 0, 1), color=col, lw=0.9, ls=(0, (4.0, 1.8)), zorder=5)
        value(ax, v, f"{v:.1f}% ({rmse(key, lo):.1f})", ink, lo == drawn)                                            # with its error at 100 times
        ticks(ax, j, i == len(MODELS) - 1 and not EV)
        if lo == drawn:                                                                                               # the fit drawn in Figure 1B
            tint(ax, col)
        if i == 0:
            ax.set_title(f"Fitted for N ≥ {lo:,}", fontsize=7.0, color=INK, fontweight="bold", pad=17)
            ax.text(M, 102.5, "292k human gut\nmetagenomes", fontsize=5.0, color=INK3, ha="center", va="bottom", linespacing=1.05, clip_on=False)      # on top of their dotted lines, as in Figure 1B
            ax.text(WORLD, 102.5, "world\npopulation", fontsize=5.0, color=INK3, ha="center", va="bottom", linespacing=1.05, clip_on=False)
    fig.text(0.012, TOP - (i + 0.5) * (TOP - BOT3) / len(MODELS), name, rotation=90, fontsize=7.0, color=ink, fontweight="bold", ha="left", va="center", linespacing=1.15)
tails_out = []
if EV:                                                                          # the last row: the extreme-value model of the similarity tail, one panel per tail
    C_EV = "#111111"; NE = np.array(EV["N"]); DRAWN_T = EV["drawn"]
    assert EV["donors"] == D and len(EV["tails"]) == 5 and abs(NE[0] - M) < 1e-6 and NE[-1] >= X_MAX
    assert EV["curve"]["rate"] == EV["models"]["free"][str(DRAWN_T)]["curve"], "the tail drawn in Figure 1B"
    gs4 = fig.add_gridspec(1, 5, left=0.115, right=0.985, top=BOT + ROW / fig.get_figheight(), bottom=BOT, wspace=0.10)
    for j, r in enumerate(EV["tails"]):
        ax = fig.add_subplot(gs4[0, j]); m_ = EV["models"]["free"][str(r)]; c = np.array(m_["curve"]); v = m_["at_world"]
        assert abs(c[0] - 100 * R1[-1]) < 1e-6, "the model starts from the measured share at the size of the HGM set"
        assert abs(float(np.interp(np.log10(WORLD), np.log10(NE), c)) - v) < 1e-3, "the curve of the tail gives its value at the population of the world"
        frame(ax)
        ax.plot(NS, 100 * R1, color=C_SK, lw=1.5, alpha=0.22, zorder=3, solid_capstyle="round")                       # the measured curve: this model is not fitted to it
        ax.plot(NE, c, color=C_EV, lw=0.9, ls=(0, (4.0, 1.8)), zorder=5)
        value(ax, v, f"{v:.1f}%", C_EV, r == DRAWN_T)                                                                 # no error here: the test of this model is in Supplementary Table 18
        ticks(ax, j, True)
        if r == DRAWN_T:                                                                                              # the tail drawn in Figure 1B
            tint(ax, C_EV, 0.045)                                                                                     # black tints more than a colour
        ax.set_title(f"r = {r:,}\nhighest similarities", fontsize=7.0, color=INK, fontweight="bold", pad=4, linespacing=1.1)
        tails_out.append((r, v))
    fig.text(0.012, BOT + 0.5 * ROW / fig.get_figheight(), "Extreme-value\nmodel", rotation=90, fontsize=7.0, color=C_EV, fontweight="bold", ha="left", va="center", linespacing=1.15)
fig.text(0.062, (TOP + BOT) / 2, "Donors recovered (%)", rotation=90, fontsize=6.6, color=INK, ha="center", va="center")
fig.text(0.55, BOT - 0.27 / fig.get_figheight(), "Pool size, N", fontsize=6.8, color=INK, ha="center", va="center")
handles = [Line2D([], [], color=C_SK, lw=1.5, alpha=0.75, label="measured, other sample first, part used by the fit"),
           Line2D([], [], color=C_SK, lw=1.5, alpha=0.22, label="measured, not used by the fit"),
           Line2D([], [], color=INK3, lw=0.9, ls=(0, (4.0, 1.8)), label="model, one colour per model as in Fig. 1B"),
           Line2D([], [], color=INK3, ls="none", marker="x", ms=4.2, mew=0.9, label="value at 8×10$^9$ people (in parentheses, error at 100×)"),
           Patch(facecolor=matplotlib.colors.to_rgba(INK3, 0.10), edgecolor=INK3, lw=0.8, label="drawn in Fig. 1B"),
           Line2D([], [], color=INK3, lw=0.45, ls=(0, (1.5, 2.0)), label="291,881 HGMs, 8×10$^9$ people")]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.55, 0.0), ncol=2, frameon=False, fontsize=6.2, labelcolor=INK2, handlelength=2.2, columnspacing=1.6, labelspacing=0.45)
fig.savefig(PREFIX + ".png", dpi=600); fig.savefig(PREFIX + ".pdf")
json.dump(dict(meta=dict(donors=D, found=int(found.sum()), pi_found=100 * float(found.mean()), measured_first_at_M=100 * R1[-1], windows=list(WIN), grid_sizes=int(len(NS)), grid_values=70), models=models, grid=grid),
          open(PREFIX + "_values.json", "w"), indent=1)
print("wrote", PREFIX + ".png and .pdf")
for key_, v_ in grid.items():
    print(f"  {key_:16} first at 8e9 {v_['first_8e9']:5.1f}% | N9 {v_['N9']:.2e} | RMSE at 100x {v_['rmse_100']:.2f}")
for r_, v_ in tails_out:
    print(f"  extreme-value model, {r_:>6,} highest similarities: first at 8e9 {v_:5.1f}%")
for key_, v_ in models.items():
    print(f"  drawn {key_:10} N >= {v_['window']:>6}: residual over the whole curve {v_['residual_whole']:.3f} points")
