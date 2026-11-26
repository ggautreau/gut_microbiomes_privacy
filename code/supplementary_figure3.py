#!/usr/bin/env python3
"""Supplementary Figure 3: human DNA in the samples, and recovery. Vincent is included, to show why it was left out.

    python3 -I supplementary_figure3.py JOIN_TABLE.csv POOL_SUMMARY.json OUTPUT_PREFIX      # writes OUTPUT_PREFIX.png (600 dpi), .pdf and _values.json (aggregates only)

JOIN_TABLE.csv  one line per donor of the 558 (data/human_dna/donors_human_share.csv): `study`, `rank1` (the donor's other sample is the
                nearest of the 291,881 runs of the HGM set), `share_q` and `share_p` (% of the hashes of the sketch of each of the two samples that are among the hashes of the human reference, T2T-CHM13v2.0 or GRCh38.p14),
                `human_hashes_q` and `human_hashes_p` (their number), `stat_q` and `stat_p` (the human share of the same two runs in the taxonomic analysis that the NCBI runs on each run of the SRA, STAT:
                100 x human_hits / species_hits, the % of the species-level hits of STAT that are assigned to Homo sapiens; local/humancheck/all558/stat_per_run.json).
POOL_SUMMARY    aggregates/human_dna_pool_summary.json (the pool_final_summary.json of all558_pool.py): rank 1 in the pool of the 1,116 runs of the 558 donors, with all hashes, without the human hashes, and with these only.
A  share of human hashes in each sample, by cohort (a dot is a sample, a thin line joins the two samples of a donor, the bar is the median of the cohort)
B  share of donors recovered, by the human share of the more human of their two samples (Wilson 95% intervals)
C  share of donors whose other sample is the nearest of the 1,116 runs, with all hashes, without the human hashes and with the human hashes only, for the donors under and over 5% human
D  the human share of each of the 1,116 runs, by its human hashes and by the NCBI STAT analysis
E  share of donors recovered, by the STAT human share of the more human of their two samples (Wilson 95% intervals)
Palette and type follow Figure 1; grey = low in human DNA (cohorts kept, donors under 5%), red = human-rich (Vincent, donors with 5% or more).
"""
import csv, json, math, sys
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

JOIN, POOL, PREFIX = sys.argv[1:4]
INK, INK2, INK3, RULE, GUIDE = "#1B2227", "#3F4A52", "#6B777F", "#C9D0D4", "#E3E7EA"
C_RED, C_GREY = "#C8553D", "#6B777F"
import os
from matplotlib import font_manager
FONTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fonts")                  # Liberation Sans shipped with the package (licences in fonts/), so that the drawing does not depend on the machine
font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if f.name != "Liberation Sans"]
for fn in sorted(os.listdir(FONTS)):
    if fn.lower().endswith(".ttf"):
        font_manager.fontManager.addfont(os.path.join(FONTS, fn))
plt.rcParams.update({"font.family": ["Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42})
NAMES = {"VincentC_2016": "Vincent", "MehtaRS_2018": "Mehta", "CosteaPI_2017": "Costea", "HansenLBS_2018": "Hansen", "RaymondF_2016": "Raymond", "Bengtsson-PalmeJ_2015": "Bengtsson-Palme",
         "HMP_2019_ibdmdb": "HMP IBD", "Heitz-BuschartA_2016": "Heintz-Buschart"}
LEFT_OUT = {"VincentC_2016"}
FLOOR = 1e-3                                                                      # shares below 0.001% are drawn at 0.001%
W, H = 6.5, 7.6
H_OLD = 5.3                                                                       # height of the figure before panels D and E were added: A, B and C keep their size
yf = lambda y: (y * H_OLD + (H - H_OLD)) / H                                      # a bottom or top coordinate of the former figure, in the taller one
hf = lambda h: h * H_OLD / H


def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


def style(ax, grid=None):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK3)
    ax.tick_params(length=2.5, width=0.5, color=INK3, labelcolor=INK2, labelsize=6.5)
    if grid:
        ax.grid(axis=grid, color=GUIDE, linewidth=0.45, zorder=0); ax.set_axisbelow(True)


rows = list(csv.DictReader(open(JOIN)))
assert {r["study"] for r in rows} == set(NAMES)
by = {}
for r in rows:
    by.setdefault(r["study"], []).append((float(r["share_q"]), float(r["share_p"])))
allv = {c: np.array([x for p in v for x in p]) for c, v in by.items()}
order = sorted(by, key=lambda c: -np.median(allv[c]))                             # the most human cohort at the top
N = len(order)
fig = plt.figure(figsize=(W, H), dpi=300)
vals = dict(cohorts={})
# ------------------------------------------------------------------------------------------------ A: human hashes in each sample, by cohort
ax = fig.add_axes((0.185, yf(0.575), 0.565, hf(0.335))); style(ax, "x")
ax.set_xscale("log"); ax.set_xlim(0.7 * FLOOR, 130); ax.set_ylim(-0.6, N - 0.4)
rng = np.random.default_rng(3)
tr = blended_transform_factory(ax.transAxes, ax.transData)
for k, c in enumerate(order):
    y = N - 1 - k; col = C_RED if c in LEFT_OUT else C_GREY
    for a, b in by[c]:
        yy = y + rng.uniform(-0.30, 0.30); xa, xb = max(a, FLOOR), max(b, FLOOR)
        ax.plot([xa, xb], [yy, yy], color=col, lw=0.35, alpha=0.28, solid_capstyle="round", zorder=2)
        ax.scatter([xa, xb], [yy, yy], s=3.0, color=col, alpha=0.65, linewidths=0, zorder=3)
    v = allv[c]; more = np.array([max(a, b) for a, b in by[c]])                    # the more human of the two samples of a donor
    med = float(np.median(more)); q1, q3 = (float(x) for x in np.percentile(more, [25, 75]))
    ax.plot([max(med, FLOOR)] * 2, [y - 0.40, y + 0.40], color=INK, lw=1.6, solid_capstyle="butt", zorder=5)
    vals["cohorts"][NAMES[c]] = dict(donors=len(by[c]), samples=int(len(v)), median_more_human_sample_percent=med, q25_percent=q1, q75_percent=q3, max_percent=float(v.max()), median_of_all_samples_percent=float(np.median(v)),
                                     donors_5_percent_or_more=int((more >= 5).sum()), donors_50_percent_or_more=int((more >= 50).sum()))
    ax.text(1.085, y, f"{len(by[c])}", transform=tr, ha="right", va="center", fontsize=6.3, color=INK3)
    ax.text(1.265, y, ("<0.01" if med < 0.005 else (f"{med:.2f}" if med < 10 else f"{med:.0f}")), transform=tr, ha="right", va="center", fontsize=6.3, color=INK2)
    ax.text(1.43, y, f"{int((more >= 5).sum())}", transform=tr, ha="right", va="center", fontsize=6.3, color=INK2)
for x_, s_ in ((1.085, "Donors"), (1.265, "Median (%)"), (1.43, "Donors\n\u22655%")):
    ax.text(x_, N - 0.40, s_, transform=tr, ha="right", va="bottom", fontsize=5.8, color=INK3, clip_on=False, linespacing=1.1)
ax.axvline(5, color=INK3, lw=0.6, ls=(0, (1.5, 2.0)), zorder=1)
ax.set_yticks(range(N)); ax.set_yticklabels([NAMES[c] + ("  (left out)" if c in LEFT_OUT else "") for c in reversed(order)], fontsize=6.8, color=INK)
ax.set_xticks([FLOOR, 0.01, 0.1, 1, 10, 100]); ax.set_xticklabels(["≤0.001", "0.01", "0.1", "1", "10", "100"])
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_xlabel("Hashes of a sample that match a human genome (%)", fontsize=7, color=INK, labelpad=3)
ax.legend(handles=[Line2D([], [], marker="o", ls="-", lw=0.5, ms=2.6, color=C_GREY, alpha=0.9, label="cohort kept: a dot is a sample, a line joins the two samples of a donor"),
                   Line2D([], [], marker="o", ls="-", lw=0.5, ms=2.6, color=C_RED, alpha=0.9, label="Vincent, left out"),
                   Line2D([], [], color=INK, lw=1.6, label="median of the more human sample of each donor; dotted line, 5%")],
          loc="lower left", bbox_to_anchor=(-0.02, 1.0), ncol=1, frameon=False, fontsize=6.2, labelcolor=INK2, handlelength=2.0, borderaxespad=0.1, labelspacing=0.2)
# ------------------------------------------------------------------------------------------------ B: share recovered by the human share of the donor
more_all = np.array([max(float(r["share_q"]), float(r["share_p"])) for r in rows]); linked = np.array([r["rank1"] == "True" for r in rows])
bins = [(0, 0.5, "<0.5"), (0.5, 5, "0.5–5"), (5, 50, "5–50"), (50, 101, "≥50")]
ax = fig.add_axes((0.105, yf(0.105), 0.34, hf(0.29))); style(ax, "y")
vals["linked_by_human_share"] = []
for i, (lo, hi, lab) in enumerate(bins):
    m = (more_all >= lo) & (more_all < hi); k, n = int(linked[m].sum()), int(m.sum()); a, b = wilson(k, n)
    ax.bar(i, 100 * k / n, width=0.62, color=C_RED if lo >= 5 else C_GREY, linewidth=0, zorder=2)
    ax.plot([i, i], [a, b], color=INK, lw=0.8, solid_capstyle="butt", zorder=3)
    ax.text(i, b + 2.0, f"{100 * k / n:.0f}%", ha="center", va="bottom", fontsize=6.4, color=INK, zorder=4)
    vals["linked_by_human_share"].append(dict(label=lab, donors=n, linked=k, percent=100 * k / n, wilson_low=a, wilson_high=b))
ax.set_xticks(range(len(bins))); ax.set_xticklabels([f"{b[2]}\n(n = {v['donors']})" for b, v in zip(bins, vals["linked_by_human_share"])], fontsize=6.2)
ax.set_ylim(0, 112); ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xlabel("Human hashes in the more human sample of the donor (%)", fontsize=6.8, color=INK, labelpad=3)
ax.set_ylabel("Donors recovered (%)", fontsize=7, color=INK, labelpad=2)
# ------------------------------------------------------------------------------------------------ C: rank 1 in the pool with all hashes, without and with the human hashes only
P = json.load(open(POOL))["definitions"]
rowsC = [("all hashes", "All hashes"), ("exact human removed", "Human hashes\nremoved"), ("exact human only", "Human hashes\nonly")]
ax = fig.add_axes((0.675, yf(0.105), 0.285, hf(0.29))); style(ax, "x")
ypos = np.arange(len(rowsC))[::-1].astype(float); vals["pool_rank1_percent"] = {}
for y, (key, lab) in zip(ypos, rowsC):
    poor, rich = P[key]["top1_human_share_below_5"], P[key]["top1_human_share_5_or_more"]
    vals["pool_rank1_percent"][lab.replace("\n", " ")] = dict(under_5_percent_human=poor, five_percent_or_more=rich)
    ax.plot([rich["percent"], poor["percent"]], [y, y], color=RULE, lw=2.2, solid_capstyle="round", zorder=2)
    ax.scatter([poor["percent"]], [y], s=19, color=C_GREY, zorder=4, linewidths=0); ax.scatter([rich["percent"]], [y], s=19, color=C_RED, zorder=4, linewidths=0)
    if max(rich["percent"], poor["percent"]) < 15:                                  # both dots next to the axis: the numbers go above them
        ax.text(poor["percent"], y + 0.17, f"{poor['percent']:.0f}", fontsize=6.0, color=C_GREY, ha="center", va="bottom")
        ax.text(rich["percent"] + 1.5, y + 0.17, f"{rich['percent']:.0f}", fontsize=6.0, color=C_RED, ha="left", va="bottom")
        continue
    lo_, hi_ = sorted((rich["percent"], poor["percent"]))
    ax.text(rich["percent"] - 3.4 if rich["percent"] < poor["percent"] else rich["percent"] + 3.4, y, f"{rich['percent']:.0f}", fontsize=6.0, color=INK2, ha="right" if rich["percent"] < poor["percent"] else "left", va="center")
    ax.text(poor["percent"] + 3.4 if poor["percent"] >= rich["percent"] else poor["percent"] - 3.4, y, f"{poor['percent']:.0f}", fontsize=6.0, color=INK2, ha="left" if poor["percent"] >= rich["percent"] else "right", va="center")
ax.set_yticks(ypos); ax.set_yticklabels([r[1] for r in rowsC], fontsize=6.5, color=INK)
ax.set_xlim(-6, 112); ax.set_xticks([0, 25, 50, 75, 100]); ax.set_ylim(-0.6, len(rowsC) - 0.35)
ax.set_xlabel("Other sample nearest of the 1,116 runs (%)", fontsize=6.8, color=INK, labelpad=3)
nr, npo = P["all hashes"]["top1_human_share_5_or_more"]["n"], P["all hashes"]["top1_human_share_below_5"]["n"]
ax.legend(handles=[Line2D([], [], marker="o", ls="none", ms=3.8, color=C_GREY, label=f"under 5% human (n = {npo})"), Line2D([], [], marker="o", ls="none", ms=3.8, color=C_RED, label=f"5% or more (n = {nr})")],
          loc="lower right", bbox_to_anchor=(1.0, 0.22), frameon=False, fontsize=6.0, labelcolor=INK2, handletextpad=0.3, borderaxespad=0.0, labelspacing=0.3)
# ------------------------------------------------------------------------------------------------ D and E: the human share of the same runs in the NCBI STAT analysis (stat_per_run.json)
assert all(k in rows[0] for k in ("human_hashes_q", "human_hashes_p", "stat_q", "stat_p")), "donors_human_share.csv has no STAT columns"
oq = np.array([float(r["share_q"]) for r in rows]); op = np.array([float(r["share_p"]) for r in rows])
sq = np.array([float(r["stat_q"]) for r in rows]); sp = np.array([float(r["stat_p"]) for r in rows])
hq = np.array([int(r["human_hashes_q"]) for r in rows]); hp = np.array([int(r["human_hashes_p"]) for r in rows])
ours_r = np.concatenate([oq, op]); stat_r = np.concatenate([sq, sp]); hh_r = np.concatenate([hq, hp])
vin_r = np.array([r["study"] in LEFT_OUT for r in rows] * 2); coh_r = [r["study"] for r in rows] * 2
sc = dict(runs=int(len(ours_r)), spearman=float(spearmanr(ours_r, stat_r)[0]))
ge = hh_r >= 11                                                                  # a share that rests on 11 human hashes or more (below, counting noise dominates)
sc["spearman_11_or_more_human_hashes"] = dict(runs=int(ge.sum()), value=float(spearmanr(ours_r[ge], stat_r[ge])[0]))
sc["spearman_under_11_human_hashes"] = dict(runs=int((~ge).sum()), value=float(spearmanr(ours_r[~ge], stat_r[~ge])[0]))
sc["thresholds"] = {}
for t in (1, 5, 10, 50):
    A_, B_ = ours_r >= t, stat_r >= t
    tp, fp, fn, tn = int((A_ & B_).sum()), int((~A_ & B_).sum()), int((A_ & ~B_).sum()), int((~A_ & ~B_).sum())
    n_ = tp + fp + fn + tn; po = (tp + tn) / n_; pe = ((tp + fn) * (tp + fp) + (fp + tn) * (fn + tn)) / n_ ** 2
    sc["thresholds"][str(t)] = dict(both_high=tp, only_stat_high=fp, only_ours_high=fn, both_low=tn, kappa=float((po - pe) / (1 - pe)))
rich = ours_r >= 20; low = (ours_r >= 0.1) & (ours_r < 5)
sc["median_ratio_stat_over_ours_when_ours_20_or_more"] = dict(runs=int(rich.sum()), value=float(np.median(stat_r[rich] / ours_r[rich])))
sc["median_ratio_stat_over_ours_when_ours_0_1_to_5"] = dict(runs=int(low.sum()), value=float(np.median(stat_r[low] / ours_r[low])))
more_o = np.maximum(oq, op); more_s = np.maximum(sq, sp); vin_d = np.array([r["study"] in LEFT_OUT for r in rows])
sc["donors_5_percent_or_more"] = dict(vincent=dict(donors=int(vin_d.sum()), ours=int((more_o[vin_d] >= 5).sum()), stat=int((more_s[vin_d] >= 5).sum())),
                                      others=dict(donors=int((~vin_d).sum()), ours=int((more_o[~vin_d] >= 5).sum()), stat=int((more_s[~vin_d] >= 5).sum())))
sc["by_cohort"] = {}
for c in order:
    m = np.array([x == c for x in coh_r])
    sc["by_cohort"][NAMES[c]] = dict(runs=int(m.sum()), max_ours_percent=float(ours_r[m].max()), max_stat_percent=float(stat_r[m].max()), runs_5_percent_or_more_ours=int((ours_r[m] >= 5).sum()), runs_5_percent_or_more_stat=int((stat_r[m] >= 5).sum()))
# D: one dot per run
axD = fig.add_axes((0.105, 0.58 / H, 0.36, 1.45 / H)); style(axD)
xs, ys = np.maximum(ours_r, FLOOR), np.maximum(stat_r, FLOOR)
for m_, col_, z_ in ((~vin_r, C_GREY, 2), (vin_r, C_RED, 3)):
    axD.scatter(xs[m_], ys[m_], s=2.6, color=col_, alpha=0.42, linewidths=0, zorder=z_)
axD.plot([0.7 * FLOOR, 130], [0.7 * FLOOR, 130], color=INK3, lw=0.5, ls=(0, (1.5, 2.0)), zorder=1)
axD.axvline(5, color=RULE, lw=0.6, zorder=1); axD.axhline(5, color=RULE, lw=0.6, zorder=1)
axD.set_xscale("log"); axD.set_yscale("log"); axD.set_xlim(0.7 * FLOOR, 130); axD.set_ylim(0.7 * FLOOR, 130)
for axis_ in (axD.set_xticks, axD.set_yticks):
    axis_([FLOOR, 0.01, 0.1, 1, 10, 100])
axD.set_xticklabels(["≤0.001", "0.01", "0.1", "1", "10", "100"]); axD.set_yticklabels(["≤0.001", "0.01", "0.1", "1", "10", "100"])
axD.xaxis.set_minor_locator(matplotlib.ticker.NullLocator()); axD.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
axD.set_xlabel("Hashes of a run that match a human genome (%)", fontsize=6.8, color=INK, labelpad=3)
axD.set_ylabel("Human share by NCBI STAT (%)", fontsize=6.8, color=INK, labelpad=2)
same5 = sc["thresholds"]["5"]["both_high"] + sc["thresholds"]["5"]["both_low"]
axD.text(0.035, 0.965, f"Spearman {sc['spearman']:.2f}, {sc['runs']:,} runs\nsame side of 5%: {same5:,}", transform=axD.transAxes, ha="left", va="top", fontsize=6.2, color=INK2, linespacing=1.25)
axD.legend(handles=[Line2D([], [], marker="o", ls="none", ms=3.2, color=C_GREY, alpha=0.8, label="cohorts kept"), Line2D([], [], marker="o", ls="none", ms=3.2, color=C_RED, alpha=0.8, label="Vincent, left out")],
           loc="lower right", frameon=False, fontsize=6.0, labelcolor=INK2, handletextpad=0.3, borderaxespad=0.2, labelspacing=0.3)
# E: share of donors recovered by the STAT human share
axE = fig.add_axes((0.60, 0.58 / H, 0.36, 1.45 / H)); style(axE, "y")
vals_e = []
for i, (lo, hi, lab) in enumerate(bins):
    m = (more_s >= lo) & (more_s < hi); k, n = int(linked[m].sum()), int(m.sum()); a, b = wilson(k, n)
    axE.bar(i, 100 * k / n, width=0.62, color=C_RED if lo >= 5 else C_GREY, linewidth=0, zorder=2)
    axE.plot([i, i], [a, b], color=INK, lw=0.8, solid_capstyle="butt", zorder=3)
    axE.text(i, b + 2.0, f"{100 * k / n:.0f}%", ha="center", va="bottom", fontsize=6.4, color=INK, zorder=4)
    vals_e.append(dict(label=lab, donors=n, linked=k, percent=100 * k / n, wilson_low=a, wilson_high=b))
sc["linked_by_stat_human_share"] = vals_e
axE.set_xticks(range(len(bins))); axE.set_xticklabels([f"{b[2]}\n(n = {v['donors']})" for b, v in zip(bins, vals_e)], fontsize=6.2)
axE.set_ylim(0, 112); axE.set_yticks([0, 25, 50, 75, 100])
axE.set_xlabel("Human share by NCBI STAT in the more\nhuman sample of the donor (%)", fontsize=6.8, color=INK, labelpad=3, linespacing=1.2)
axE.set_ylabel("Donors recovered (%)", fontsize=7, color=INK, labelpad=2)
vals["stat_comparison"] = sc
for x_, y_, s_ in ((0.008, 0.992, "A"), (0.008, 0.447, "B"), (0.50, 0.447, "C")):
    fig.text(x_, yf(y_), s_, fontsize=9, fontweight="bold", color=INK, ha="left", va="top")
for x_, s_ in ((0.008, "D"), (0.50, "E")):
    fig.text(x_, (2.03 + 0.276) / H, s_, fontsize=9, fontweight="bold", color=INK, ha="left", va="top")
fig.savefig(PREFIX + ".png", dpi=600)
fig.savefig(PREFIX + ".pdf")
json.dump(vals, open(PREFIX + "_values.json", "w"), indent=1)
for n_, v in vals["cohorts"].items():
    print(f"{n_:16} donors {v['donors']:3d}  median {v['median_more_human_sample_percent']:7.3f}%  IQR {v['q25_percent']:.3f}-{v['q75_percent']:.3f}  max {v['max_percent']:5.1f}%  donors >= 5%: {v['donors_5_percent_or_more']:2d}  >= 50%: {v['donors_50_percent_or_more']:2d}")
print("recovered by human share:", [(d["label"], d["linked"], d["donors"]) for d in vals["linked_by_human_share"]])
print("pool rank 1:", {k: (round(v["under_5_percent_human"]["percent"], 1), round(v["five_percent_or_more"]["percent"], 1)) for k, v in vals["pool_rank1_percent"].items()})
print("STAT comparison:", f"Spearman {sc['spearman']:.3f}", "| 5% line", {k_: sc["thresholds"]["5"][k_] for k_ in ("both_high", "only_stat_high", "only_ours_high", "both_low")}, "| recovered by STAT share:", [(d["label"], d["linked"], d["donors"]) for d in sc["linked_by_stat_human_share"]])
print("wrote", PREFIX + ".png and .pdf")
