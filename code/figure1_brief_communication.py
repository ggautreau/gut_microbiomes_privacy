# Figure 1 of the Brief Communication. A = schematic of the analysis; B = recovery rate against the pool size, measured and extrapolated (Pitman-Yor correctness model, straight line in log N, power law of the miss rate and extreme-value model of the similarity tail, enlarged around the population of the world in a round zoom), with under it the table of the value at the population of the world of the three models against the window of their fit;
# C = share of donors with the other sample in the top X; D = recovery in each cohort; E = recovery against the time between the samples; F = recovery against the depth of the queried run; G = the donors of Raymond with and without a course of cephalosporins; H = similarity of the other sample and of the nearest run of other people when the donor is absent.
#!/usr/bin/env python3
"""Figure 1 of the Brief Communication, drawn from the per-donor tables (data/adult_ranks.csv and data/open_set_scores.csv).

    python3 -I figure1_brief_communication.py ARCHIVE_DIR OUTPUT_PREFIX VALUES_JSON          # ARCHIVE_DIR holds data/adult_ranks.csv, data/open_set_scores.csv, aggregates/raymond_antibiotic.json, aggregates/extrapolation_summary.json, aggregates/evt_summary.json and fonts/; VALUES_JSON = supplementary_figure5_values.json (run supplementary_figure5.py first)

A the analysis on the two samples of one donor (a drawing); B the recovery rate, the share of donors whose other sample is first in a pool of N runs, measured up to the whole HGM set of 291,881
runs and extrapolated by the Pitman-Yor correctness model of Rocher et al. 2025 (pyc_mb.py), a straight line in log N and a power law of the miss rate, all three fitted for N >= 1,000 (where the measured share and its upper bound differ by at most 0.2 points; Supplementary Note 3), and by the extreme-value model of the similarity tail (aggregates/evt_summary.json, written by evt_model.py; the tail selected by lowest thinning RMSE), with the interval of the straight line, a round zoom around the population of the world that holds the four values (the fits drawn are checked against VALUES_JSON and against the windows of aggregates/extrapolation_summary.json); C the share of donors whose other sample is among the top X runs of the HGM set; D recovery in each cohort, with the median size of the queried runs;
E recovery of the first sample alone against the days between the first and the last sample; F recovery against the bases of the queried run; G the donors of Raymond with and without a course of cephalosporins between their samples (E, F and G are drawn on one pair of axes: a bar has the same width in the three;
the bars of D to G are 95% Wilson score intervals of the share); H the similarity of the donor's other sample and of the nearest run of other people when the donor is absent (the same query with the donor's samples removed), the threshold at which as many absent persons are accepted as present donors are rejected, and the donors above it.

Input   data/adult_ranks.csv, one row per donor of the 7 cohorts kept (483), in the order of data/donors.json. Columns: study; age_class; q_lib (DNA);
        rank_any_narrow = 1 + the number of runs of the HGM set (291,881 human gut metagenomes of the Logan index, Methods) more similar to the query than the donor's other sample (empty = the other
        sample is not among the 5,000 nearest runs, which counts as not recovered); far_possible = 1 when the donor has another sample at least 30 days from the query sample; rank_far_narrow = the same rank when
        only those samples can be the other sample; days = days between the first and the last sample; rank_first_narrow = the rank when the first sample alone is the other sample (panel E); query_gbp = bases
        of the queried run in Gbp (panels D and F).
        aggregates/raymond_antibiotic.json (panel G): the donors of Raymond by the study_condition of their last stool sample (cephalosporins or control) and, for each group, those whose other sample is first (raymond_antibiotic.py).
        data/open_set_scores.csv (panel H), one row per donor in the same order: study; `partner` = the similarity of the donor's other sample (0 = beyond the 5,000 nearest runs);
        best_other = the similarity of the nearest run of other people among the 291,881 runs of the HGM set (make_open_set_scores.py).
Output  OUTPUT_PREFIX.png (600 dpi) and .pdf; OUTPUT_PREFIX_quadratic.json (the straight lines of panel B, their values and 95% intervals: Supplementary Note 3 and Table 15, and the numbers of the text) and OUTPUT_PREFIX_open_set.json (the numbers of panel H).
Random numbers: numpy default_rng(1); 100,000 resamples for the measured bands; first 800 resamples (selected by adaptive doubling in extrapolation_bootstrap.py) for the straight-line extrapolation.
"""
import csv, json, math, os, sys
import numpy as np
from scipy.special import gammaln
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerPatch
from matplotlib.transforms import Bbox
import matplotlib.patheffects as pe
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                  # python -I does not put the script's folder on the path
import pyc_mb as PY                                                             # Pitman-Yor correctness model of Rocher et al. 2025

HERE = os.path.abspath(sys.argv[1])                                   # the unpacked reproduction package (data/, fonts/)
PREFIX = sys.argv[2]
VALUES = sys.argv[3]                                                    # supplementary_figure5_values.json: the three models of Fig. 1B at the population of the world, as drawn
os.makedirs(os.path.dirname(os.path.abspath(PREFIX)), exist_ok=True)

FONTS = os.path.join(HERE, "fonts")                                     # Liberation Sans shipped with the package (licences in fonts/), so that the drawing does not depend on the machine
if os.path.isdir(FONTS):
    font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if f.name != "Liberation Sans"]
    for fn in sorted(os.listdir(FONTS)):
        if fn.lower().endswith(".ttf"):
            font_manager.fontManager.addfont(os.path.join(FONTS, fn))

M = 291_881                                                             # the HGM set: the runs of the index taken as human gut metagenomes (build_stock.py)
W, H, DPI = 7.2, 4.05, 600                                              # width of the canvas (inches), a height that H_AB + H_C replaces below, and the resolution (dpi)
CROP = 0.62                                                             # not used
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"], "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.6, "ytick.major.size": 2.6,
    "xtick.minor.size": 1.5, "xtick.minor.width": 0.45, "xtick.major.pad": 2.0, "ytick.major.pad": 2.0,
    "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic", "mathtext.bf": "Liberation Sans:bold", "mathtext.default": "regular",
})
INK, INK2, INK3, RULE, GUIDE = "#14171A", "#3A4047", "#5C646C", "#9AA2AA", "#DCE0E4"
C_SK, C_SK_TINT, C_SK_INK = "#D55E00", "#F9D9BF", "#8A3D00"             # k-mer sketches: vermillion, warm as it is the risk (the Okabe-Ito vermillion, safe for colour-blind readers against the blue of the thresholds); tint and dark ink
C_PY, C_PY_INK = "#009988", "#00665B"                                      # the Pitman-Yor extrapolation: teal (curve; diamond, value and name). The extrapolations, teal, wine, violet and black, were checked against the vermillion of the data and the blue of the thresholds for normal, deutan, protan and tritan vision
C_MS, C_MS_INK = "#A469F7", "#5A35A8"                                      # the power law of the miss rate: violet (curve; diamond, value and name)
C_EV, C_EV_INK = "#111111", "#111111"                                      # the extreme-value model of the similarity tail: black (curve; diamond, value and name), the one colour that stays apart from the five others for normal, deutan, protan and tritan vision (an olive falls on the vermillion for deutan vision)
C_LN, C_LN_TINT, C_LN_INK = "#882255", "#E9C9D8", "#882255"            # the straight line in log N: wine, dark against the teal so that their alternating dashes read; its band a pale rose, not a green zone
C_RED = "#C8553D"
C_BL, C_BL_INK = "#2F5FA8", "#1D3F75"                                      # the re-identification risk thresholds of panel B: cool, as the zone below them is the acceptable one
F_TICK, F_LAB, F_NOTE, F_NUM, F_LET = 6.4, 6.8, 6.2, 6.6, 9.0

with open(os.path.join(HERE, "data", "adult_ranks.csv"), newline="") as fh:
    rows = [r for r in csv.DictReader(fh) if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]
study = np.array([r["study"] for r in rows])
K = np.array([(int(r["rank_any_narrow"]) - 1) if r["rank_any_narrow"] != "" else M for r in rows], float)      # not found counts as a failure


def p_first(k, N):
    """Probability that a donor with k runs of the HGM set more similar than the other sample is first in a pool of N runs (the other sample and N - 1 runs drawn without replacement)."""
    n1 = N - 1
    ok = (M - k) >= n1
    out = np.zeros(len(k))
    out[ok] = np.exp(gammaln(M - k[ok] + 1) - gammaln(M - k[ok] - n1 + 1) - gammaln(M + 1) + gammaln(M - n1 + 1))
    return out


Ns = np.unique(np.round(np.logspace(np.log10(2), np.log10(M), 60)).astype(int))
R_all = np.array([p_first(K, N).mean() for N in Ns])
RESAMPLES = 100_000
RESAMPLES_LINE = 800                                                      # chosen by adaptive doubling in extrapolation_bootstrap.py (Supplementary Fig. 6)
rng = np.random.default_rng(1)
cohorts = sorted(set(study)); byc = {c: np.where(study == c)[0] for c in cohorts}
# cohorts resampled, then donors. RES_W[b, i] = share of donor i in resample b: a cohort of n donors drawn m times gives m n draws among its n donors.
m_c = rng.multinomial(len(cohorts), [1 / len(cohorts)] * len(cohorts), size=RESAMPLES)
RES_W = np.zeros((RESAMPLES, len(rows)))
for j, c in enumerate(cohorts):
    RES_W[:, byc[c]] = rng.multinomial(m_c[:, j] * len(byc[c]), [1 / len(byc[c])] * len(byc[c]))
RES_W /= RES_W.sum(axis=1, keepdims=True)
lo, hi = np.percentile(RES_W @ np.stack([p_first(K, N) for N in Ns], axis=1), [2.5, 97.5], axis=0)

# ---- the 10 nearest runs of the pool, and the extrapolation beyond the HGM set
from scipy.stats import hypergeom
k_all = np.array([(int(r["rank_any_narrow"]) - 1) if r["rank_any_narrow"] != "" else -1 for r in rows])            # runs of the HGM set above the other sample (-1: not among the 5,000 returned)
found_all = k_all >= 0
TOP = 10                                                                    # the other sample is among the TOP nearest runs of the pool
X_EXT = np.logspace(math.log10(M), 12, 400)                                 # the extrapolation starts at the last measured point
WORLD = 8e9


def top_inside(N):
    """Share of donors whose other sample is among the TOP nearest runs of a pool of N runs (N <= M + 1): the number of runs of the HGM set above the other sample among the N - 1 runs drawn without
    replacement is hypergeometric (count k in the HGM set); it must not exceed TOP - 1. A donor not found counts as a failure."""
    return float(hypergeom.cdf(TOP - 1, M, k_all[found_all], int(N) - 1).sum() / len(rows))


R_top = np.array([top_inside(N) for N in Ns]); R_top_M = R_top[-1]

# The three models, each fitted at the pool sizes N >= 1,000 of the grid of 70 log-spaced values from 2 to M rounded to whole numbers (68 distinct sizes, the grid of Supplementary Fig. 5): under 1,000 runs
# the measured share is only a lower bound, and from 1,000 runs on it is exact to 0.2 points (Supplementary Note 3). A straight line in log N fitted by least squares to each measured curve (the recovery
# rate and the share among the 10 nearest); shares clipped to 0-100%. Its interval repeats the fit on each of the resamples of cohorts, then donors, of the measured band. The Pitman-Yor correctness model
# and the power law of the miss rate are fitted to the same sizes (Supplementary Fig. 5).
NF = np.unique(np.round(np.logspace(math.log10(2), math.log10(M), 70)).astype(int))
FIT_FROM = 1_000
FIT_DEG = 1                                                                                   # a straight line in log N (Supplementary Note 3)
FIT_SEL = NF >= FIT_FROM
P1 = np.stack([p_first(K, N) for N in NF], axis=1)                                                                  # donors x pool sizes: first among the pool
P10 = np.stack([np.where(found_all, hypergeom.cdf(TOP - 1, M, np.where(found_all, k_all, 0), int(N) - 1), 0.0) for N in NF], axis=1)      # among the TOP nearest of the pool
lg2 = lambda x: np.log2(np.asarray(x, float))                                                            # the fits are in log N; the code takes the base 2, which only sets the unit of the coefficients (points per doubling of the pool)
quad_fit = lambda R: np.polyfit(lg2(NF[FIT_SEL]), R[FIT_SEL], FIT_DEG)
quad_raw = lambda c, N: np.polyval(c, lg2(N))
quad_at = lambda c, N: np.clip(quad_raw(c, N), 0, 1)
C1, C10 = quad_fit(P1.mean(0)), quad_fit(P10.mean(0))
Y1, Y10 = RES_W[:RESAMPLES_LINE] @ P1, RES_W[:RESAMPLES_LINE] @ P10
CB1, CB10 = np.polyfit(lg2(NF[FIT_SEL]), Y1[:, FIT_SEL].T, FIT_DEG), np.polyfit(lg2(NF[FIT_SEL]), Y10[:, FIT_SEL].T, FIT_DEG)      # the fit on each resample (coefficients x resamples)
at_all = lambda CB, N: 100 * np.clip(np.polyval(CB, lg2(N)), 0, 1)                                                  # the fits of the resamples at one pool size
band_ext = lambda CB: np.percentile(100 * np.clip(np.polyval(CB, lg2(X_EXT)[:, None]), 0, 1), [2.5, 97.5], axis=1)  # their band along X_EXT
Q1, Q10 = 100 * quad_at(C1, X_EXT), 100 * quad_at(C10, X_EXT)
Q1_draw = np.where(quad_raw(C1, X_EXT) > 0, Q1, np.nan)                                                            # stop the line where it reaches zero
SEL_PY = FIT_SEL                                                                              # the window of the Pitman-Yor model: the same sizes as the straight line
f_py = PY.fit(NF[SEL_PY].astype(float), P1.mean(0)[SEL_PY])                                                      # the Pitman-Yor correctness model (PYC-MB), two parameters, fitted for N >= 1,000 as in Supplementary Fig. 5
py_at = lambda N: 100 * np.clip(np.atleast_1d(PY.kappa(np.asarray(N, float), f_py)), 0, 1)
PY_draw = py_at(X_EXT); py_w = float(py_at(WORLD)[0])
SEL_PW = FIT_SEL                                                                              # the window of the power law of the miss rate: the same sizes as the straight line
c_ms = np.polyfit(np.log(NF[SEL_PW].astype(float)), np.log(np.clip(1 - P1.mean(0)[SEL_PW], 1e-12, None)), 1)      # the power law of the miss rate, 1 - R = a N^b, as in Supplementary Fig. 5
ms_at = lambda N: 100 * np.clip(1 - np.exp(np.polyval(c_ms, np.log(np.asarray(N, float)))), 0, 1)
MS_draw = ms_at(X_EXT); ms_w = float(ms_at(WORLD))
print(f"power law of the miss rate: first at 8e9 {ms_w:.1f}%, at the HGM set {MS_draw[0]:.2f}% (a = {math.exp(c_ms[1]):.4f}, b = {c_ms[0]:.4f})")
print(f"Pitman-Yor correctness model: first at 8e9 {py_w:.1f}%, at the HGM set {PY_draw[0]:.2f}% (h = {f_py.h:.1f} nats, gamma = {f_py.gamma:.4f})")
EVT = json.load(open(os.path.join(HERE, "aggregates", "evt_summary.json")))                                                 # the extreme-value model of the similarity tail (evt_model.py), for the tail drawn
ev_at = lambda N: np.interp(np.log10(np.asarray(N, float)), np.log10(EVT["curve"]["N"]), EVT["curve"]["rate"])            # its curve beyond the HGM set
EV_draw = ev_at(X_EXT); ev_w = float(EVT["models"]["free"][str(EVT["drawn"])]["at_world"])
assert abs(float(ev_at(WORLD)) - ev_w) < 0.05 and abs(EV_draw[0] - 100 * P1.mean(0)[-1]) < 1e-6, "the curve of the extreme-value model starts at the measured rate and gives its value at the population of the world"
print(f"extreme-value model of the similarity tail: first at 8e9 {ev_w:.1f}% (the {EVT['drawn']:,} highest similarities of each donor)")
lo_q1, hi_q1 = band_ext(CB1)
ci = lambda CB, N: [float(v) for v in np.percentile(at_all(CB, N), [2.5, 97.5])]
resid = lambda C, P: float(100 * (quad_at(C, M) - P.mean(0)[-1]))                                                  # fitted minus measured at N = M


def crossing(C, level):
    """Pool size above M at which the unclipped polynomial falls to `level` (share), or None."""
    roots = np.roots(list(C[:-1]) + [C[-1] - level]); roots = roots[np.isreal(roots)].real
    ok = [2.0 ** r for r in roots if 2.0 ** r > M]
    return float(min(ok)) if ok else None


Xc = np.logspace(math.log10(M), 14, 40000)
cross = np.where(quad_raw(C1, Xc) > quad_raw(C10, Xc))[0]                                   # beyond this pool size the recovery-rate line exceeds the 10-nearest one, which a share that includes it cannot do
N_cross = float(Xc[cross[0]]) if len(cross) else None
NTAB = (1e6, 1e7, 1e8, 1e9, 8e9, 1e10, 1e11)
quad_out = dict(fit_from=FIT_FROM, fit_points=int(FIT_SEL.sum()), resamples=RESAMPLES_LINE, coef_first=[float(v) for v in C1], coef_top10=[float(v) for v in C10],
                measured_first_at_M=100 * float(R_all[-1]), measured_top10_at_M=100 * float(R_top_M), residual_first_at_M=resid(C1, P1), residual_top10_at_M=resid(C10, P10),
                max_abs_residual_first=float(100 * np.abs(quad_at(C1, NF[FIT_SEL]) - P1.mean(0)[FIT_SEL]).max()), max_abs_residual_top10=float(100 * np.abs(quad_at(C10, NF[FIT_SEL]) - P10.mean(0)[FIT_SEL]).max()),
                N9_first=crossing(C1, 0.09), N0_first=crossing(C1, 0.0), N9_top10=crossing(C10, 0.09), N_cross=N_cross,
                first={str(int(N)): dict(est=float(100 * quad_at(C1, N)), lo=ci(CB1, N)[0], hi=ci(CB1, N)[1]) for N in NTAB},
                top10={str(int(N)): dict(est=float(100 * quad_at(C10, N)), lo=ci(CB10, N)[0], hi=ci(CB10, N)[1]) for N in NTAB})
json.dump(quad_out, open(PREFIX + "_quadratic.json", "w"), indent=1)
fmt_n = lambda v: "never" if v is None else f"{v:.2e}"
q1w, q10w = quad_out["first"][str(int(WORLD))]["est"], quad_out["top10"][str(int(WORLD))]["est"]
VAL = json.load(open(VALUES))["models"]
HOLD = json.load(open(os.path.join(HERE, "aggregates", "extrapolation_summary.json")))["holdout"]["0.01"]        # the test on thinned HGM sets at 100 times the thinned size (Supplementary Note 3)
rmse = lambda key: math.hypot(HOLD["error"][key], HOLD["sd"][key])
HW1 = json.load(open(os.path.join(HERE, "aggregates", "extrapolation_summary.json")))["holdout_windows"]["0.01"]      # the same test for the fifteen fits
rmse_w = lambda key, lo: math.hypot(HW1["error"][f"{key}|{lo}"], HW1["sd"][f"{key}|{lo}"])
assert all(abs(rmse_w(k_, w_) - rmse(k_)) < 1e-9 for k_, w_ in HOLD["windows"].items()), "the error of a fit drawn is that of its window"                                                  # rmse and rmse_w above: root mean square error over the thinnings, in points, as in Supplementary Note 3
assert HOLD["windows"] == {"pycmb": FIT_FROM, "line": FIT_FROM, "miss": FIT_FROM} and (SEL_PY == FIT_SEL).all() and (SEL_PW == FIT_SEL).all(), "the fits drawn are those tested on thinned HGM sets"
assert abs(q1w - VAL["line"]["first_8e9"]) < 1e-6 and all(VAL[k_]["window"] == w_ for k_, w_ in HOLD["windows"].items()) and abs(py_w - VAL["pycmb"]["first_8e9"]) < 1e-3 and abs(ms_w - VAL["miss"]["first_8e9"]) < 1e-6, "Fig. 1B and Supplementary Fig. 5 draw the same fits"
print(f"10 nearest of the pool: {100 * R_top[0]:.1f}% at N = {Ns[0]}, {100 * top_inside(10_000):.1f}% at 1e4, {100 * R_top_M:.1f}% at N = {Ns[-1]}")
print(f"straight line in log N on N >= {FIT_FROM} ({FIT_SEL.sum()} sizes): first at 8e9 {q1w:.1f}% ({quad_out['first'][str(int(WORLD))]['lo']:.1f}-{quad_out['first'][str(int(WORLD))]['hi']:.1f}), 10 nearest {q10w:.1f}% ({quad_out['top10'][str(int(WORLD))]['lo']:.1f}-{quad_out['top10'][str(int(WORLD))]['hi']:.1f}); "
      f"residual at M {quad_out['residual_first_at_M']:+.2f} / {quad_out['residual_top10_at_M']:+.2f} points; the two fits cross at {fmt_n(quad_out['N_cross'])}; first falls to 9% at {fmt_n(quad_out['N9_first'])}, to 0 at {fmt_n(quad_out['N0_first'])}")


def sci(x):
    e = int(math.floor(math.log10(x))); m = x / 10 ** e
    return f"{m:.1f}×10$^{{{e}}}$" if m < 9.95 else f"10$^{{{e + 1}}}$"


def style(a, grid="y"):
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        a.spines[sp].set_color(RULE)
    if grid:
        a.grid(axis=grid, color=GUIDE, linewidth=0.45, zorder=0)
    a.set_axisbelow(True)
    a.tick_params(axis="both", labelsize=F_TICK, colors=INK2)



# ======================================================================================================================
# Figure 1 of the Brief Communication: A, the analysis on one donor's pair of samples (schematic, drawn in inches);
# B, the recovery rate against the pool size.
# ======================================================================================================================
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle, Wedge
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
from matplotlib.transforms import blended_transform_factory


def HALF(cx, cy, r, z):                                                           # the other sample of the donor: a disc cut by a diagonal through its centre, half grey (upper left) and half orange (lower right), as two filled shapes without outline
    axA.add_patch(Circle((cx, cy), r, facecolor=C_SK, edgecolor="none", zorder=z, clip_on=False))
    axA.add_patch(Wedge((cx, cy), r + 0.0006, 45, 225, facecolor=INK3, edgecolor="none", zorder=z + 0.1, clip_on=False))

YB = 0.64                                                                              # bottom of the axes of the last row (inches): room for the counts, the rotated labels and the axis labels of E, F and G
H_AB, H_C = 4.35, YB + 1.60                                                  # heights (inches) of the block of A and B, with C and D under A, and of the last row below it (E to H)
H = H_AB + H_C
fig = plt.figure(figsize=(W, H), dpi=DPI)

# ------------------------------------------------------------------ panel A: the analysis
WA = 3.90
axA = fig.add_axes((0, H_C / H, WA / W, H_AB / H))
axA.set_xlim(-0.24, WA - 0.24); axA.set_ylim(0, H_AB); axA.axis("off")                 # the content is 0.24 in to the right: room for the left side
X0 = 0.30


def T(x, y, s, size=F_NOTE, color=INK2, ha="left", va="center", **kw):
    return axA.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, **kw)


# panel A, drawn compactly: the index on two lines
t = H_AB
days_all = np.array([float(r["days"]) for r in rows])
# The index as a rectangle with rounded corners and grey dots; the human gut metagenomes as a circle inside it (enlarged: to scale it would be 6% of the area), a crowd
# of discs: the runs of other people (white) and the other sample of the donor, half grey and half orange (is it first, recovered, or not?). Schematic:
# the positions and the numbers of dots and discs mean nothing.
A_UP = 0.29                                                                                 # the content of panel A 0.29 in higher: its two headings were removed
yc, R_INDEX, xc, xq = t - 0.875 + A_UP, 0.335, 2.08, 0.514                              # middle of the top row, half width of the right part of the index, its centre, the query sample
xs_, ys_ = xc - 0.2725, yc                                                             # centre of the HGM set: the centre of the index
HW_S, HH_S, RC_S = 0.43, 0.43, 0.43                                                    # its half width, half height and corner radius, all equal: a circle (the crowd and the grey dots use this rounded-box form)
xt = xc + R_INDEX + 0.09                                                               # the labels, to the right of the index
dots_rng = np.random.default_rng(7)                                                    # the dots only: the random numbers of the analysis are not touched
dots = []
IX_L, IX_R, IX_B, IX_T, IX_C = xc - R_INDEX - 0.545, xc + R_INDEX, yc - 0.515, yc + 0.515, 0.11          # the index: a rectangle with rounded corners
for _ in range(9000):
    p_ = (dots_rng.uniform(IX_L + 0.024, IX_R - 0.024), dots_rng.uniform(IX_B + 0.024, IX_T - 0.024))
    qx_, qy_ = max(IX_L + IX_C - p_[0], 0, p_[0] - (IX_R - IX_C)), max(IX_B + IX_C - p_[1], 0, p_[1] - (IX_T - IX_C))     # beyond the straight sides: in a corner
    if qx_ ** 2 + qy_ ** 2 <= (IX_C - 0.024) ** 2 and ((max(abs(p_[0] - xs_) - (HW_S - RC_S), 0) ** 2 + max(abs(p_[1] - ys_) - (HH_S - RC_S), 0) ** 2) ** 0.5 + min(max(abs(p_[0] - xs_) - (HW_S - RC_S), abs(p_[1] - ys_) - (HH_S - RC_S)), 0) >= RC_S + 0.022) and all((p_[0] - q_[0]) ** 2 + (p_[1] - q_[1]) ** 2 >= 0.045 ** 2 for q_ in dots):
        dots.append(p_)
axA.add_patch(FancyBboxPatch((IX_L, IX_B), IX_R - IX_L, IX_T - IX_B, boxstyle=f"round,pad=0,rounding_size={IX_C}", facecolor="white", edgecolor=INK3, lw=0.7, zorder=2))
for p_ in dots:
    axA.add_patch(Circle(p_, 0.0125, facecolor="#BBC1C7", edgecolor="none", zorder=3))
axA.add_patch(Circle((xs_, ys_), HW_S, facecolor=C_SK_TINT, edgecolor=C_SK, lw=0.6, zorder=3.2))
R_RANK = 0.054                                                                          # radius of a disc of the set, as the discs of the lists and the query sample
sc_rng = np.random.default_rng(1)                                              # a crowd of 40 discs in the HGM set: random starts, then pushed apart until none touches another
sc_P = np.column_stack([sc_rng.uniform(xs_ - HW_S + 0.05, xs_ + HW_S - 0.05, 40), sc_rng.uniform(ys_ - HH_S + 0.05, ys_ + HH_S - 0.05, 40)])
sc_P[0] = (xs_, ys_)                                                                   # one disc exactly at the centre of the set, kept there
sc_I = np.eye(40)
for sc_it in range(20000):
    sc_dv = sc_P[:, None, :] - sc_P[None, :, :]
    sc_dd = np.sqrt((sc_dv ** 2).sum(-1)) + sc_I
    sc_ov = np.clip(0.119 - sc_dd, 0, None) * (1 - sc_I)
    sc_P = sc_P + 0.25 * ((sc_ov / sc_dd)[:, :, None] * sc_dv).sum(1)
    sc_sx, sc_sy = np.where(sc_P[:, 0] >= xs_, 1.0, -1.0), np.where(sc_P[:, 1] >= ys_, 1.0, -1.0)             # back inside the set, a disc and 0.01 in of margin
    sc_dx, sc_dy = np.abs(sc_P[:, 0] - xs_) - (HW_S - RC_S), np.abs(sc_P[:, 1] - ys_) - (HH_S - RC_S)
    sc_qx, sc_qy = np.maximum(sc_dx, 0), np.maximum(sc_dy, 0)
    sc_nq = np.hypot(sc_qx, sc_qy)
    sc_k = np.where(sc_nq > RC_S - 0.064, (RC_S - 0.064) / np.where(sc_nq > 0, sc_nq, 1.0), 1.0)
    sc_P = np.column_stack([xs_ + sc_sx * np.where(sc_dx > 0, HW_S - RC_S + sc_qx * sc_k, np.abs(sc_P[:, 0] - xs_)), ys_ + sc_sy * np.where(sc_dy > 0, HH_S - RC_S + sc_qy * sc_k, np.abs(sc_P[:, 1] - ys_))])
    sc_P[0] = (xs_, ys_)
    if sc_ov.max() < 2e-4:
        break
assert sc_ov.max() < 2e-4, "the discs of the HGM set do not fit"
sc_pts = [tuple(q_) for q_ in sc_P]
sc_low = [q_ for q_ in sc_pts if ys_ - 0.25 <= q_[1] <= ys_ - 0.09]
sc_o = max(sc_low, key=lambda q_: q_[0])                                              # the other sample of the donor: the rightmost disc a little under the centre, near the height of its label
ranked = [sc_o] + [q_ for q_ in sc_pts if q_ is not sc_o]

for n_, (px_, py_) in enumerate(ranked, 1):
    if n_ == 1:
        HALF(px_, py_, R_RANK + 0.7 / 144, 5)                                         # the other sample of the donor, half grey, half orange
    else:
        axA.add_patch(Circle((px_, py_), R_RANK, facecolor="white", edgecolor=RULE, lw=0.7, zorder=5))      # a run of other people
axA.add_patch(Circle((xq, yc), 0.054, facecolor="white", edgecolor=C_SK, lw=1.0, zorder=3))
FA_SOLID = os.path.join(FONTS, "FontAwesome5Free-Solid-900.otf")                 # Font Awesome 5 Free, solid icons (SIL OFL 1.1; fonts/ of the package): the glyphs "female" and "male" of the left side, drawn as paths
T(xq, yc - 0.099, "query sample\n(last sample)", ha="center", va="top", linespacing=1.05)
axA.annotate("", (IX_L - 0.03, yc), (xq + 0.079, yc), arrowprops=dict(arrowstyle="-|>", lw=0.7, color=INK3, shrinkA=0, shrinkB=0, mutation_scale=6))
T((xq + 0.079 + IX_L - 0.03) / 2, yc + 0.05, "against the\nother public\nsamples", ha="center", va="bottom", color=INK3, linespacing=1.05)
T(xt, yc + 0.255, "Logan sketch index", color=INK, fontweight="bold")
T(xt, yc + 0.155, "4.8M metagenomes", color=INK2)                                       # rounded in the panel; the text gives 4,791,208
a_, y_stock, y_other = np.deg2rad(32), yc - 0.045, yc - 0.17
axA.plot([xs_ + HW_S, xt - 0.025], [ys_ + 0.03, y_stock], color=C_SK, lw=0.8, zorder=6, solid_capstyle="round")
T(xt, y_stock, "292k human gut metagenomes", color=C_SK_INK)                                  # rounded in the panel; the caption gives 291,881
axA.plot([ranked[0][0] + R_RANK * 0.95, xt - 0.025], [ranked[0][1] - 0.008, y_other], color=C_SK_INK, lw=0.8, zorder=6, solid_capstyle="round")
T(xt, y_other, "other sample of the donor", color=C_SK_INK, fontweight="bold")
AN_q = max([q_ for q_ in ranked[1:] if ys_ - 0.40 <= q_[1] <= ys_ - 0.20], key=lambda q_: q_[0])        # a white disc: a run of other people
y_people = yc - 0.29
axA.plot([AN_q[0] + R_RANK * 0.95, xt - 0.025], [AN_q[1] - 0.008, y_people], color=INK3, lw=0.8, zorder=6, solid_capstyle="round")
T(xt, y_people, "other people", color=INK2)
# Two cases, one chart each: the rank of a run (horizontal) and its similarity to the query (vertical), with a numbered disc on a stem for each of
# the five nearest runs, then the list goes on. Recovered: the other sample of the recovered donors, its median similarity (extrapolation_tables.py,
# key panel_a), then four runs of other people on a slight slope (illustrative: their medians, 0.190 to 0.181, are nearly flat). Not recovered: a
# schematic example, the other sample third.
pa_r = fig.canvas.get_renderer()
pa_wid = lambda a_: a_.get_window_extent(pa_r).width / fig.dpi                 # width of a text, in inches (the axes of panel A are drawn in inches)
pa_yh = t - 1.40 + A_UP
PA = json.load(open(os.path.join(HERE, "aggregates", "extrapolation_summary.json")))["panel_a"]
assert PA["recovered"] == int(np.sum(K == 0)) and len(PA["runs_of_other_people"]) == 5 and PA["other_sample"] > PA["runs_of_other_people"][0]
pa_vals = [[PA["other_sample"], 0.190, 0.178, 0.168, 0.160],                   # recovered: the median of the other sample, then a slight slope (illustrative)
           [0.182, 0.166, 0.150, 0.143, 0.138]]                       # not recovered: schematic, not data (the one donor at rank 3 has six similarities within 0.174-0.184)
pa_cases = [("Recovered", 1, C_SK_INK, 0, C_SK), ("Not recovered", 3, INK2, 2, INK3)]      # name, rank of the other sample, ink of the name, its position, its colour
pa_yb, pa_k, pa_cx, pa_dx, pa_end, pa_R = pa_yh - 0.81, 1.45, [0.61, 2.225], 0.158, 0.99, 0.054      # baseline, inches per unit of similarity, left end of each chart, step between two ranks, length of a baseline, radius of a disc (as in the index)
for pa_g in (0.2, 0.4):                                                        # light guides at 0.2 and 0.4, across the two charts
    axA.plot([pa_cx[0] - 0.02, pa_cx[1] + pa_end], [pa_yb + pa_g * pa_k] * 2, color=GUIDE, lw=0.45, zorder=1)
axA.plot([pa_cx[0] - 0.02] * 2, [pa_yb, pa_yb + 0.42 * pa_k], color=INK3, lw=0.6)
for pa_g in (0.0, 0.2, 0.4):
    axA.plot([pa_cx[0] - 0.045, pa_cx[0] - 0.02], [pa_yb + pa_g * pa_k] * 2, color=INK3, lw=0.6)
    T(pa_cx[0] - 0.07, pa_yb + pa_g * pa_k, f"{pa_g:.1f}" if pa_g else "0", ha="right", color=INK3)
axA.text(X0 + 0.01, pa_yb + 0.21 * pa_k, "Similarity", fontsize=F_NOTE, color=INK2, rotation=90, ha="center", va="center")
for (pa_name, pa_rank, pa_ink, pa_oth, pa_col), pa_c, pa_v6 in zip(pa_cases, pa_cx, pa_vals):
    pa_n = T(pa_c, pa_yh - 0.13, pa_name, color=pa_ink, fontweight="bold")
    T(pa_c + pa_wid(pa_n) + 0.05, pa_yh - 0.13, f"(rank {pa_rank})", color=INK3)
    axA.plot([pa_c - 0.02, pa_c + pa_end], [pa_yb, pa_yb], color=INK3, lw=0.6, solid_capstyle="butt")
    for pa_i, pa_v in enumerate(pa_v6):
        pa_x, pa_o = pa_c + 0.09 + pa_dx * pa_i, pa_i == pa_oth
        axA.plot([pa_x, pa_x], [pa_yb, pa_yb + pa_v * pa_k - pa_R], color=pa_col if pa_o else RULE, lw=1.0, solid_capstyle="butt", zorder=3)
        axA.add_patch(Circle((pa_x, pa_yb + pa_v * pa_k), pa_R, facecolor=pa_col if pa_o else "white", edgecolor=pa_col if pa_o else RULE, lw=0.7, zorder=5))
        axA.text(pa_x, pa_yb + pa_v * pa_k - 0.0015, str(pa_i + 1), fontsize=6.2, color="white" if pa_o else INK, fontweight="bold", ha="center", va="center", zorder=6)
    pa_xl = pa_c + 0.09 + pa_dx * 4                                          # the last lollipop
    for pa_d in range(3):                                                      # the list goes on, as far as the 5,000th run
        axA.add_patch(Circle((pa_xl + 0.10 + 0.03 * pa_d, pa_yb + 0.05), 0.010, facecolor=INK3, edgecolor="none"))


# ------------------------------------------------------------------ panel A, left side: the analysis is repeated for each of the 483 adult donors
from matplotlib.path import Path as pl_Path
assert int(np.sum(K >= 0)) == 483                                                  # the donors of the analysis


def pl_icon(glyph, cx, cy, h, color=INK2):                                          # an icon of Font Awesome 5 Free (solid), a vector path of height h inches centred on (cx, cy)
    tp_ = TextPath((0, 0), glyph, size=1.0, prop=font_manager.FontProperties(fname=FA_SOLID))
    e_ = tp_.get_extents()
    axA.add_patch(PathPatch(Affine2D().translate(-(e_.x0 + e_.width / 2), -(e_.y0 + e_.height / 2)).scale(h / e_.height).translate(cx, cy).transform_path(tp_),
                            facecolor=color, edgecolor="none", zorder=4, clip_on=False))


def pl_brace(x, y0, y1, w, color=INK3, lw=0.8):                                     # a curly brace pointing left: ends at x, tip at x - w, from y0 (bottom) to y1 (top)
    ym_, r_ = (y0 + y1) / 2, w / 2
    v_ = [(x, y1), (x - r_, y1), (x - r_, y1), (x - r_, y1 - r_), (x - r_, ym_ + r_), (x - r_, ym_), (x - r_, ym_), (x - w, ym_),
          (x - r_, ym_), (x - r_, ym_), (x - r_, ym_ - r_), (x - r_, y0 + r_), (x - r_, y0), (x - r_, y0), (x, y0)]
    c_ = [pl_Path.MOVETO] + [pl_Path.CURVE4] * 3 + [pl_Path.LINETO] + [pl_Path.CURVE4] * 6 + [pl_Path.LINETO] + [pl_Path.CURVE4] * 3
    axA.add_patch(PathPatch(pl_Path(v_, c_), facecolor="none", edgecolor=color, lw=lw, capstyle="round", joinstyle="round", zorder=4, clip_on=False))


pl_top, pl_bot = t - 0.07, t - 1.92                                                 # the left side spans the panel, from the top of the index to the baseline of the charts
pl_xr, pl_w = 0.25, 0.12                                                            # the right end of the brace and its depth
PL_F, PL_M = "\uf182", "\uf183"
PL_YD, PL_YN, PL_YA = 1.44, 1.64, 1.78                                              # under the top: the three dots, the number, the word
pl_gap_ir, pl_r, pl_lw, pl_gap_rr = 0.022, 0.043, 0.9, 0.012                        # icon to rounds, radius of the rounds, their line width (pt), space between the two rounds
pl_ro = pl_r + pl_lw / 144                                                          # outer radius of a round, in inches
pl_fw = TextPath((0, 0), PL_F, size=1.0, prop=font_manager.FontProperties(fname=FA_SOLID)).get_extents()
pl_fw = 0.14 * pl_fw.width / pl_fw.height                                           # width of the widest icon (the woman)
pl_W = pl_fw + pl_gap_ir + 4 * pl_ro + pl_gap_rr                                    # width of one row: the icon and the two rounds
PL_CX = (-0.20 + (pl_xr - pl_w) - 0.04) / 2                                         # the row is centred between the left margin and the tip of the brace
pl_x0 = PL_CX - pl_W / 2
pl_brace(pl_xr, pl_bot, pl_top, pl_w)
for pl_r_i in range(6):
    pl_y = pl_top - 0.24 - 0.20 * pl_r_i
    pl_icon(PL_F if pl_r_i % 2 == 0 else PL_M, pl_x0 + pl_fw / 2, pl_y, 0.14)
    pl_c1 = pl_x0 + pl_fw + pl_gap_ir + pl_ro                                       # centre of the first round
    pl_c2 = pl_c1 + 2 * pl_ro + pl_gap_rr                                           # centre of the second round
    HALF(pl_c1, pl_y, pl_ro, 3)          # the other sample, half grey, half orange
    axA.add_patch(Circle((pl_c2, pl_y), pl_r, facecolor="white", edgecolor=C_SK, lw=pl_lw, zorder=3, clip_on=False))        # the query sample: empty, second
for pl_d in (-0.030, 0.0, 0.030):                                                   # the series goes on
    axA.add_patch(Circle((PL_CX, pl_top - PL_YD + pl_d), 0.010, facecolor=INK3, edgecolor="none", clip_on=False))
axA.text(PL_CX, pl_top - PL_YN, "483", ha="center", va="center", fontsize=F_LAB, color=INK, fontweight="bold")
axA.text(PL_CX, pl_top - PL_YA, "adults", ha="center", va="center", fontsize=F_NOTE - 0.8, color=INK2)


# ------------------------------------------------------------------ panel A: the question (is the other sample of the donor first?) and its two answers
AR_yh = pa_yh - 0.13                                                                   # the row of the headers
AR_Q = (xs_, AR_yh)                                                                   # the question, between the two headers, under a stem from the index
axA.plot([AR_Q[0], AR_Q[0]], [IX_B, AR_Q[1] + 0.064], color=INK3, lw=0.8, solid_capstyle="butt", zorder=5)
for AR_x0, AR_x1, AR_xd, AR_col in ((AR_Q[0] - 0.064, 1.43, 1.62, C_SK), (AR_Q[0] + 0.064, 2.185, 1.995, INK3)):   # the two answers: recovered (orange) or not recovered (grey)
    axA.annotate("", (AR_x1, AR_yh), (AR_x0, AR_yh), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=AR_col, shrinkA=0, shrinkB=0, mutation_scale=6.5))
    axA.add_patch(Circle((AR_xd, AR_yh), 0.054, facecolor=AR_col, edgecolor=AR_col, lw=0.7, zorder=6))
HALF(AR_Q[0], AR_Q[1], 0.054 + 0.7 / 144, 6)                       # the question: the other sample, half grey, half orange, with a question mark
axA.text(AR_Q[0], AR_Q[1] - 0.002, "?", fontsize=6.8, color="white", fontweight="bold", ha="center", va="center", zorder=6.2)

# ------------------------------------------------------------------ panel B: the recovery rate against the pool size
XB0, WB, HB = 4.35, 2.65, 3.20                                                          # left edge, width and height of the axes of panel B (inches)
ax = fig.add_axes((XB0 / W, (H_C + H_AB - 0.42 - HB) / H, WB / W, HB / H))
style(ax, "y")
ax.set_xscale("log")
X_MAX = 1e12
ax.axvline(8e9, color=INK3, lw=0.6, ls=(0, (1.5, 2.0)), zorder=2)
ax.fill_between(Ns, 100 * lo, 100 * hi, color=C_SK_TINT, alpha=0.65, lw=0, zorder=1)
ax.fill_between(X_EXT, lo_q1, hi_q1, color=C_LN_TINT, alpha=0.45, lw=0, zorder=1)
ax.plot(Ns, 100 * R_all, color=C_SK, lw=1.7, zorder=5, solid_capstyle="round")
LW_EXT, DASH, GAP = 1.5, 2.6, 1.5                                                                # the extrapolations: dashes 1.5 pt thick, one colour per model
alt = lambda k_, m_: ((m_ - k_) % m_ * (DASH + GAP), (DASH, m_ * (DASH + GAP) - DASH))            # the k-th of m curves that overlap: the same short dashes, shifted so that the colours alternate with a gap between them
EV_HIGH = abs(ev_w - (py_w + q1w) / 2) < 0.8                                                     # the extreme-value model runs along the Pitman-Yor curve and the straight line (within the width of a line), or else under them, along the power law near the HGM set
if EV_HIGH:                                                                                      # three curves overlap: their three colours alternate, and the power law keeps the dash and the gap
    LS_PY, LS_EV, LS_LN, LS_MS = alt(0, 3), alt(1, 3), alt(2, 3), (0, (DASH, GAP))
else:                                                                                            # two pairs of curves overlap, the Pitman-Yor curve and the straight line, then the extreme-value model and the power law near the HGM set; the second pair is shifted by a quarter of the period from the first
    LS_PY, LS_LN, LS_MS, LS_EV = alt(0, 2), alt(1, 2), ((DASH + GAP) / 2, (DASH, DASH + 2 * GAP)), (3 * (DASH + GAP) / 2, (DASH, DASH + 2 * GAP))
LS_KEY = (0, (2.6, 5.6))                                                                         # the same dashes in the legend, starting with a dash
ZD = lambda k_, m_: ((m_ - k_) % m_ * DASH, (DASH, (m_ - 1) * DASH))                             # in the zoom the curves part: the dashes of the curves of a group alternate without a gap
LS_PY_Z, LS_EV_Z, LS_LN_Z, LS_MS_Z = (ZD(0, 3), ZD(1, 3), ZD(2, 3), (0, (DASH, DASH))) if EV_HIGH else (ZD(1, 2), ZD(1, 2), ZD(0, 2), ZD(0, 2))
ax.plot(X_EXT, Q1_draw, color=C_LN, lw=LW_EXT, ls=LS_LN, zorder=5)
ax.plot(X_EXT, PY_draw, color=C_PY, lw=LW_EXT, ls=LS_PY, zorder=6)                                                              # the Pitman-Yor extrapolation, without interval
ax.plot(X_EXT, MS_draw, color=C_MS, lw=LW_EXT, ls=LS_MS, zorder=5)                                                               # the power law of the miss rate, without interval
ax.plot(X_EXT, EV_draw, color=C_EV, lw=LW_EXT, ls=LS_EV, zorder=5)                                                               # the extreme-value model of the similarity tail, without interval
# Thresholds of re-identification risk, in blue: solid below 9% (advised for public release by the EMA and Health Canada; the CMS rule of 11 gives 1/11), and a fading area up to 33% (thresholds used mostly for trusted recipients, 20% the traditional one).
# The area starts as a very light blue at the thin line of 33%, deepens to the level of 20%, then to the level below 9%.
A_MAX, A_33, A_20, Y_TOP = 0.34, 0.06, 0.26, 33.0                              # opacity below 9%, at 33%, at 20%, top of the area (%): it starts at the thin line of 33%
_y = np.linspace(Y_TOP, 9, 400)                                                 # rows of the gradient, from the line of 33% down to 9%
_u = np.clip((Y_TOP - _y) / max(Y_TOP - 33, 1e-9), 0, 1); _soft = _u * _u * (3 - 2 * _u)    # smooth step for the rows above 33%: none, as the area starts at 33% (Y_TOP) without transition
_img = np.zeros((len(_y), 1, 4)); _img[..., :3] = matplotlib.colors.to_rgb(C_BL)
_img[:, 0, 3] = np.where(_y > 33, A_33 * _soft, np.where(_y > 20, A_33 + (A_20 - A_33) * (33 - _y) / 13, A_20 + (A_MAX - A_20) * (20 - _y) / 11))   # linear from 33% to 20%, then to 9%
ax.imshow(_img, extent=(2, 1e12, 9, Y_TOP), origin="upper", aspect="auto", interpolation="bilinear", zorder=0.6)
ax.axhspan(0, 9, facecolor=C_BL, alpha=A_MAX, lw=0, zorder=0.6)
ax.axhline(9, color=C_BL, lw=0.9, zorder=3)
ax.axhline(33, color=C_BL, lw=0.3, zorder=3)                                    # a thin line delimits the area at 33%
ax.axvline(M, color=INK3, lw=0.5, ls=(0, (1.5, 2.0)), zorder=2)
ax.plot([M], [100 * R_all[-1]], "o", ms=7.0, color=C_SK, mec="white", mew=0.7, zorder=8)                      # the measured share at the size of the HGM set, on its dotted line
ax.annotate(f"{100 * R_all[-1]:.1f}%", (M, 100 * R_all[-1]), xytext=(0, -7), textcoords="offset points", fontsize=F_NUM, color=C_SK_INK, fontweight="bold", ha="center", va="top", zorder=9, path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])      # under its dot, centred on the dotted line of the HGM set, with a white halo where the line crosses it
import matplotlib.patheffects as _pe
for N_, side_, off_, dx_ in ((10, -1, 4.0, 0), (100, -1, 4.0, 0), (1000, -1, 4.0, 0), (10_000, 1, 3.5, -1.2), (100_000, 1, 3.5, -1.7)):      # the measured share at five pool sizes, a dot on the curve and its value; the first three under it, the last two above it (clear of the value of the HGM set and of its dotted line)
    v_ = 100 * p_first(K, N_).mean()
    ax.plot([N_], [v_], "o", ms=3.6, color=C_SK, mec="white", mew=0.5, zorder=8)
    ax.annotate(f"{v_:.1f}%", (N_, v_), xytext=(dx_, off_ * side_), textcoords="offset points", fontsize=5.2, color=C_SK_INK, ha="center", va="top" if side_ < 0 else "bottom", zorder=9)
    print(f"pool N = {N_:,}: {v_:.2f}%")
ax.set_xlim(2, X_MAX); ax.set_ylim(0, 100)
ax.set_yticks(range(0, 101, 10)); ax.set_yticklabels([str(v) for v in range(0, 101, 10)])                  # a label every 10%; the blue 9% line carries its own label
ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(5)); ax.tick_params(axis="y", which="minor", length=1.6, width=0.45, color=RULE)      # and a tick every 5%
ax.set_xticks([10.0 ** k for k in range(1, 13)]); ax.set_xticklabels([r"$10^{%d}$" % k for k in range(1, 13)])         # every order of magnitude
ax.set_ylabel("Donors recovered (%)", fontsize=F_LAB, color=INK, labelpad=3)
ax.set_xlabel("Pool size, N", fontsize=F_LAB, color=INK, labelpad=2)
ax.text(M, 101.5, "human gut metagenomes,\n292k", fontsize=F_NOTE + 0.6, color=INK3, ha="center", va="bottom", clip_on=False)               # on top of its dotted line, as the world population
ax.text(8e9, 101.5, "world population,\n8×10$^9$", fontsize=F_NOTE + 0.6, color=INK3, ha="center", va="bottom", clip_on=False)
ax.text(2.6, 10.5, "9% (public release)", fontsize=F_NOTE, color=C_BL_INK, ha="left", va="bottom", fontweight="bold")          # above the line, clear of the 10% gridline
ax.text(2.6, 33.7, "33% (highly trusted data recipients)", fontsize=F_NOTE, color=C_BL_INK, ha="left", va="bottom", zorder=6, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))       # one line, on a white background, just above the thin line of 33% (it hides the dotted line of the HGM set behind it)
ax.text(2.6, 20.9, "20% (trusted recipients)", fontsize=F_NOTE, color=C_BL_INK, ha="left", va="bottom")
EXT_ = r"$\mathit{extrapolated}$: "; ETAL_ = r" $\mathit{et\ al}$.)"                               # the words of the legend set in italics
handles = [
    Line2D([], [], color=C_SK, lw=1.7, label=r"$\mathit{measured}$: donors recovered"),
    Line2D([], [], color=C_PY, lw=LW_EXT, ls=LS_KEY, label=EXT_ + "Pitman\u2013Yor\ncorrectness (Rocher" + ETAL_),
    Line2D([], [], color=C_LN, lw=LW_EXT, ls=LS_KEY, label=EXT_ + "straight line\nin log N (Baveja" + ETAL_),
    Line2D([], [], color=C_MS, lw=LW_EXT, ls=LS_KEY, label=EXT_ + "power law of\nmiss rate (Grother" + ETAL_),
    Line2D([], [], color=C_EV, lw=LW_EXT, ls=LS_KEY, label=EXT_ + "extreme-value\nmodel of the similarity tail"),
]
_ord = {C_PY: py_w, C_LN: q1w, C_MS: ms_w, C_EV: -1.0}; handles = handles[:1] + sorted(handles[1:], key=lambda h_: -_ord[h_.get_color()])       # the three models fitted to the measured curve by decreasing value at the population of the world, then the extreme-value model (their keys without the shift of the dashes)
handles.insert(2 + [h_.get_color() for h_ in handles[1:]].index(C_LN), Patch(facecolor=C_LN_TINT, alpha=0.45, edgecolor="none", label="95% interval of the\nstraight line"))      # the band of the straight line, right after it
handles.insert(1, Patch(facecolor=C_SK_TINT, alpha=0.65, edgecolor="none", label="95% interval of the\nmeasured curve"))      # the band of the measured curve, right after its line
for i_ in reversed([i_ for i_, h_ in enumerate(handles) if i_ > 1 and not isinstance(h_, Patch)]):                     # a little room before each curve that follows: a curve and its interval stay together
    handles.insert(i_, Line2D([], [], ls="none", label=" "))
HH_ = 0.3                                                                     # height of a row of the legend in font sizes: low, so that the empty rows between the groups are a few points high
SWATCH_ = lambda legend, orig_handle, xdescent, ydescent, width, height, fontsize: Rectangle((-xdescent, -ydescent + height / 2 - 0.35 * fontsize), width, 0.7 * fontsize)      # the swatch of a band keeps the height it had
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.015, 0.875), ncol=1, frameon=True, facecolor="white", edgecolor="none", framealpha=0.0, fontsize=F_NOTE + 0.2, handlelength=1.0, handleheight=HH_, handler_map={Patch: HandlerPatch(patch_func=SWATCH_)}, handletextpad=0.45, labelcolor=INK2, borderaxespad=0.0, borderpad=0.3, labelspacing=0.12)         # in the empty zone at the left, under the measured curve and clear of the intervals
for LT_t in ax.get_legend().get_texts():                                         # no background: the dotted line shows behind the legend, a white halo keeps the text clean
    LT_t.set_path_effects([pe.withStroke(linewidth=1.6, foreground="white")])
    if LT_t.get_text() == " ":                                                    # the empty entries that separate the groups: a few points high
        LT_t.set_fontsize(1.0)


# ------------------------------------------------------------------ zoom of panel B: the four extrapolations around the population of the world, in a round inset in the empty part of the panel,
# joined by two tangents to a small circle around their values at the population of the world; the values are written in the zoom only
z_x0 = math.log10(2); z_dx = math.log10(X_MAX) - z_x0
YB0 = H_C + H_AB - 0.42 - HB                                                    # bottom of panel B, inches from the bottom of the figure
z_in = lambda x_, y_: (XB0 + WB * (math.log10(x_) - z_x0) / z_dx, YB0 + HB * y_ / 100)
ZC, ZR = z_in(10 ** 8.9, 57.0), 0.55                                            # centre and radius of the zoom (inches)
z_mid = (min(py_w, q1w, ms_w, ev_w) + max(py_w, q1w, ms_w, ev_w)) / 2
ZY = (z_mid - 1.75, z_mid + 1.75)                                               # what it shows: 3.5 points centred on the four values, so that the curves that run together part,
Z_DEC = (ZY[1] - ZY[0]) / 2 * (HB / 100) / (WB / z_dx)                          # and as many inches of panel B across as up: the zoom enlarges both ways alike (half its width, in decades)
ZX = (WORLD / 10 ** Z_DEC, WORLD * 10 ** Z_DEC)
X_K = 2 * Z_DEC / 3.0                                                           # the distances across below, in decades, were set for a zoom 3 decades wide
Z_K = (ZY[1] - ZY[0]) / 6.2                                                     # the distances below, in points of the zoom, were set for a zoom 6.2 points high: scaled so that they stay the same on paper
assert ZY[0] + 0.5 * Z_K < min(py_w, q1w, ms_w, ev_w) and max(py_w, q1w, ms_w, ev_w) + 0.5 * Z_K < ZY[1], "the four values at the population of the world fall inside the zoom"
assert py_w >= q1w, "the Pitman-Yor value is written above its curve and the straight line's under its own"
axz = fig.add_axes(((ZC[0] - ZR) / W, (ZC[1] - ZR) / H, 2 * ZR / W, 2 * ZR / H))
axz.set_xscale("log"); axz.set_xlim(*ZX); axz.set_ylim(*ZY); axz.axis("off")
z_clip = Circle((0.5, 0.5), 0.5, transform=axz.transAxes, facecolor="white", edgecolor="none", zorder=0)
axz.add_patch(z_clip)
z_art = []
for yy_ in [y_ for y_ in range(80, 100) if ZY[0] + 0.3 * Z_K < y_ < ZY[1] - 0.3 * Z_K]:                     # a line every point
    z_art += axz.plot(ZX, [yy_, yy_], color=RULE, lw=0.4, zorder=1)
z_art += axz.plot([WORLD, WORLD], ZY, color=INK3, lw=0.5, ls=(0, (1.5, 2.0)), zorder=2)
for yv_, col_, ls_ in ((Q1_draw, C_LN, LS_LN_Z), (PY_draw, C_PY, LS_PY_Z), (MS_draw, C_MS, LS_MS_Z), (EV_draw, C_EV, LS_EV_Z)):
    z_art += axz.plot(X_EXT, yv_, color=col_, lw=1.8, ls=ls_, zorder=3)
for v_, ink_, al_ in ((q1w, C_LN_INK, 0.7), (py_w, C_PY_INK, 0.7), (ms_w, C_MS_INK, 0.7), (ev_w, C_EV_INK, 0.45)):
    z_art += axz.plot([WORLD], [v_], "D", ms=4.4, mfc=ink_, mec="white", mew=0.6, alpha=al_, zorder=4)      # translucent diamonds: the upper values overlap, and the black one, drawn last, lets the straight line show
for a_ in z_art:
    a_.set_clip_path(z_clip)
axz.add_patch(Circle((0.5, 0.5), 0.5, transform=axz.transAxes, facecolor="none", edgecolor=INK3, lw=0.7, zorder=6, clip_on=False))


def z_left(y_):
    """The pool size at the left edge of the round zoom, at the height y_ (%)."""
    fy = (y_ - ZY[0]) / (ZY[1] - ZY[0]) - 0.5; fx = 0.5 - math.sqrt(max(0.25 - fy * fy, 0.0))
    return 10 ** (math.log10(ZX[0]) + fx * (math.log10(ZX[1]) - math.log10(ZX[0])))


z_rend = fig.canvas.get_renderer(); z_inv = axz.transData.inverted(); z_ax = axz.transAxes.inverted()
z_box, z_curves = [], (lambda x_: 100 * quad_at(C1, x_), lambda x_: py_at(x_), ms_at, ev_at)


def z_fits(t_):
    """True when the text t_ stays inside the circle and meets neither a curve nor a text already placed (it is then kept)."""
    bb_ = t_.get_window_extent(z_rend); (xa_, ya_), (xb_, yb_) = z_inv.transform(bb_.get_points())
    if not all((u_ - 0.5) ** 2 + (v_ - 0.5) ** 2 < 0.25 for u_, v_ in z_ax.transform([(bb_.x0, bb_.y0), (bb_.x0, bb_.y1), (bb_.x1, bb_.y0), (bb_.x1, bb_.y1)])):
        return False
    xs_ = np.logspace(math.log10(xa_), math.log10(xb_), 40)
    if any(np.any((np.asarray(f_(xs_), float) > ya_ - 0.10 * Z_K) & (np.asarray(f_(xs_), float) < yb_ + 0.10 * Z_K)) for f_ in z_curves):
        return False
    if any(xa_ < c_ and a_ < xb_ and ya_ < d_ and b_ < yb_ for a_, b_, c_, d_ in z_box):
        return False
    z_box.append((xa_, ya_, xb_, yb_))
    return True


def z_place(label, options):
    """Writes the label at the first position (x, y, vertical alignment, text properties) where it fits; stops when none does."""
    for x_, y_, va_, kw_ in options:
        t_ = axz.text(x_, y_, label, ha="left", va=va_, **kw_)
        if z_fits(t_):
            return t_
        t_.remove()
    raise AssertionError(("no room in the zoom for", label))


z_lab = dict(fontsize=F_NUM, fontweight="bold", zorder=7, path_effects=[_pe.withStroke(linewidth=1.2, foreground="white")])
z_x = WORLD * 10 ** (0.13 * X_K)                                                        # the values start just right of the diamonds
def z_cands(f_, col_, side_):
    """Positions for the value of one curve: on the side side_ (+1 above, -1 under) first, near the population of the world first, then further along the curve."""
    out_ = []
    for d_ in (0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2):
        for sd_ in (side_, -side_):
            for x_ in (z_x,) + tuple(WORLD * 10 ** (o_ * X_K) for o_ in (-0.7, -0.9, -1.1, -1.25, 0.1, 0.35, 0.6)):
                out_.append((x_, float(np.asarray(f_(x_), float).ravel()[0]) + sd_ * d_ * Z_K, "bottom" if sd_ > 0 else "top", dict(z_lab, color=col_)))
    return out_


def z_bottom(t_):
    """The bottom of a text of the zoom, in points of the zoom."""
    return z_inv.transform(t_.get_window_extent(z_rend).get_points())[0][1]


# the values of the curves that run together (two decimals: the straight line and the extreme-value model differ in the second), in one column right of the
# population of the world, just under these curves and in their order from the top; then the power law, under its own curve
z_run = [(py_w, C_PY_INK, z_curves[1]), (q1w, C_LN_INK, z_curves[0])] + ([(ev_w, C_EV_INK, ev_at)] if EV_HIGH else [])
z_y = min(float(np.asarray(f_(x_), float).ravel()[0]) for _, _, f_ in z_run for x_ in np.logspace(math.log10(z_x), math.log10(z_x) + 0.8 * X_K, 20)) - 0.15 * Z_K
for v_, col_, _ in sorted(z_run, key=lambda t_: -t_[0]):
    z_y = z_bottom(z_place(f"{v_:.2f}%", [(z_x, z_y - 0.05 * i_ * Z_K, "top", dict(z_lab, color=col_)) for i_ in range(16)])) - 0.04 * Z_K
z_place(f"{ms_w:.2f}%", z_cands(ms_at, C_MS_INK, -1))
z_scale = 0
for yy_ in [y_ for y_ in range(80, 100) if ZY[0] + 0.6 * Z_K < y_ < ZY[1] - 0.6 * Z_K]:                            # the scale of the zoom, against the left edge of the circle: two lines labelled where there is room
    try:
        z_place(f"{yy_}%", [(max(z_left(yy_ + 0.1 * Z_K), z_left(yy_ + 0.5 * Z_K)) * 10 ** (0.06 * X_K), yy_ + 0.1 * Z_K, "bottom", dict(fontsize=5.0, color=INK3, zorder=5)),
                            (max(z_left(yy_ - 0.1 * Z_K), z_left(yy_ - 0.55 * Z_K)) * 10 ** (0.06 * X_K), yy_ - 0.1 * Z_K, "top", dict(fontsize=5.0, color=INK3, zorder=5))])
        z_scale += 1
    except AssertionError:
        pass
assert z_scale >= 1, "the zoom has no scale"
if not EV_HIGH:                                                                  # the extreme-value model, when it does not run with the two curves above: above its own curve
    z_place(f"{ev_w:.2f}%", z_cands(ev_at, C_EV_INK, +1))
z_lo, z_hi = min(py_w, q1w, ms_w, ev_w), max(py_w, q1w, ms_w, ev_w)
SC, SR = z_in(WORLD, (ZY[0] + ZY[1]) / 2), HB * (ZY[1] - ZY[0]) / 200          # the small circle: the part of panel B that the zoom enlarges, around the four values at the population of the world
fig.add_artist(Circle(SC, SR, transform=fig.dpi_scale_trans, facecolor="none", edgecolor=INK3, lw=0.6, zorder=20))
z_d = math.hypot(ZC[0] - SC[0], ZC[1] - SC[1]); z_u = ((ZC[0] - SC[0]) / z_d, (ZC[1] - SC[1]) / z_d); z_p = (-z_u[1], z_u[0])
z_c = (SR - ZR) / z_d; z_s = math.sqrt(1 - z_c * z_c)
for sg_ in (1, -1):                                                             # the two outer tangents of the two circles
    n_ = (z_c * z_u[0] + sg_ * z_s * z_p[0], z_c * z_u[1] + sg_ * z_s * z_p[1])
    fig.add_artist(Line2D([SC[0] + SR * n_[0], ZC[0] + ZR * n_[0]], [SC[1] + SR * n_[1], ZC[1] + ZR * n_[1]], transform=fig.dpi_scale_trans, color=INK3, lw=0.5, zorder=20))


# ------------------------------------------------------------------ panel H: the similarity of the other sample of every donor, and of the nearest run of other people when the donor is absent
# data/open_set_scores.csv: one line per donor, in the order of the analysis (`partner` = similarity of the other sample of the donor, 0 = beyond the 5,000 nearest runs; `best_other` = similarity of the nearest run of other people).
# The search of an absent person is the same query with the donor's samples removed, so its nearest run is the `best_other` run. The threshold sits between the two distributions, at the equal error rate: as many searches of an absent person are accepted as donors present are rejected.
with open(os.path.join(HERE, "data", "open_set_scores.csv"), newline="") as fh:
    g_rows = list(csv.DictReader(fh))
g_P = np.array([float(r["partner"]) for r in g_rows]); g_a = np.array([float(r["best_other"]) for r in g_rows]); g_n = len(g_rows)
g_cand = np.unique(np.concatenate([g_P[g_P > 0], g_a]))
g_far = np.array([np.mean(g_a >= c) for c in g_cand]); g_frr = np.array([np.mean(g_P < c) for c in g_cand])        # absent persons accepted; donors present rejected (a donor whose other sample is beyond the list counts as rejected)
g_t = float(g_cand[int(np.argmin(np.abs(g_far - g_frr)))])                                                         # equal error rate: the two shares are the same
g_edges = g_t + 0.01 * np.arange(-25, 71)
g_hs, _ = np.histogram(g_a, g_edges); g_hp, _ = np.histogram(g_P[g_P > 0], g_edges)
g_acc = g_edges[:-1] >= g_t - 1e-9
g_n_oth = int(np.sum(g_P >= g_t)); g_n_abs = int(np.sum(g_a >= g_t))
g_auc = float(np.mean(g_P[:, None] > g_a[None, :]) + 0.5 * np.mean(g_P[:, None] == g_a[None, :]))                # area under the ROC curve of the test: an other sample that is not returned (0) is below every nearest run of other people
g_ret = g_P[g_P > 0]; g_auc_ret = float(np.mean(g_ret[:, None] > g_a[None, :]) + 0.5 * np.mean(g_ret[:, None] == g_a[None, :]))   # the same over the donors whose other sample is returned
json.dump(dict(donors=g_n, other_sample_score_known=int((g_P > 0).sum()), other_sample_below_the_list=int((g_P == 0).sum()), threshold_equal_error=round(g_t, 5), error_rate_at_threshold=round(float(np.mean(g_a >= g_t)), 5),
               other_sample_at_or_above=g_n_oth, other_sample_at_or_above_and_top_hit=int(np.sum((g_P >= g_t) & (g_P >= g_a))), absent_person_at_or_above=g_n_abs,
               median_other_sample=round(float(np.median(g_P[g_P > 0])), 3), median_closest_other_person=round(float(np.median(g_a)), 3), auc=round(g_auc, 4), auc_when_the_other_sample_is_returned=round(g_auc_ret, 4)), open(PREFIX + "_open_set.json", "w"), indent=1)
C_STR, C_STR_TINT = "#7C868E", "#D3D8DC"
H_TOP0 = YB + 1.60 + 0.07                                                    # top of panel H (inches) if its axis went up to g_ymax
g_ymax = max(g_hs.max(), g_hp.max()) / 0.60                                  # the scale of panel H: the highest bar at 60% of that height
G_TOPC = 102                                                                 # the axis ends just above 100
G_SC = g_ymax / G_TOPC                                                       # axes fractions of the taller layout times G_SC keep their place in donors
H_TOP = 0.30 + (H_TOP0 - 0.30) / G_SC + 0.45                                        # the same inches per donor, plus the room that panel B leaves above panel H
axG = fig.add_axes((4.30 / W, 0.30 / H, 2.70 / W, (H_TOP - 0.30) / H))
style(axG, "y")
axG.stairs(g_hs, g_edges, fill=True, color=C_STR_TINT, lw=0, zorder=2); axG.stairs(np.where(g_acc, g_hs, 0), g_edges, fill=True, color=C_STR, alpha=0.75, lw=0, zorder=2); axG.stairs(g_hs, g_edges, color=C_STR, lw=0.8, zorder=3)
axG.stairs(g_hp, g_edges, fill=True, color=C_SK_TINT, alpha=0.8, lw=0, zorder=2); axG.stairs(np.where(g_acc, g_hp, 0), g_edges, fill=True, color=C_SK, alpha=0.6, lw=0, zorder=2); axG.stairs(g_hp, g_edges, color=C_SK, lw=0.8, zorder=3)
axG.set_xlim(0, 0.95); axG.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
axG.set_ylim(0, G_TOPC); axG.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(25))
axG.axvline(g_t, ymax=0.80 * G_SC, color=INK, lw=0.7, ls=(0, (3, 2)), zorder=4)
g_tr = blended_transform_factory(axG.transData, axG.transAxes)
axG.text(g_t, 0.815 * G_SC, f"{g_t:.4f}", fontsize=F_NOTE, color=INK, fontweight="bold", ha="center", va="bottom", transform=g_tr, zorder=5)
for g_col, g_above, g_y_below, g_y_above in ((C_SK_INK, g_n_oth, 0.655, 0.74), (INK3, g_n_abs, 0.74, 0.655)):   # share of each distribution below the threshold (left of the line) and above it (right); the larger share on top on each side
    axG.text(g_t - 0.03, g_y_below * G_SC, f"{100 * (g_n - g_above) / g_n:.1f}%", fontsize=F_NOTE + 0.6, color=g_col, fontweight="bold", ha="right", va="center", transform=g_tr, zorder=5)
    axG.text(g_t + 0.03, g_y_above * G_SC, f"{100 * g_above / g_n:.1f}%", fontsize=F_NOTE + 0.6, color=g_col, fontweight="bold", ha="left", va="center", transform=g_tr, zorder=5)
# inset: the ROC curve of the test over all thresholds (the dot is the threshold of the panel) and the area under it
axR = axG.inset_axes([0.61, 0.28 * G_SC, 0.36, 0.40 * G_SC])
axG.add_patch(matplotlib.patches.Rectangle((0.47, 0.13 * G_SC), 0.53, 0.60 * G_SC, transform=axG.transAxes, facecolor="white", edgecolor="none", zorder=1, clip_on=False))   # no grid line of the panel behind the inset and its labels
g_ts = np.concatenate([[np.inf], g_cand[::-1], [0.0]])
g_fpr = np.array([np.mean(g_a >= c) for c in g_ts]); g_tpr = np.array([np.mean(g_P >= c) for c in g_ts])
axR.plot([0, 1], [0, 1], color=RULE, lw=0.5, ls=(0, (2, 2)), zorder=1)
axR.plot(g_fpr, g_tpr, color=INK2, lw=0.9, zorder=3, solid_joinstyle="round")
axR.plot([np.mean(g_a >= g_t)], [np.mean(g_P >= g_t)], "o", ms=2.6, color=INK, zorder=4)
axR.set_xlim(-0.03, 1.03); axR.set_ylim(-0.03, 1.03); axR.set_xticks([0, 0.5, 1]); axR.set_yticks([0, 0.5, 1]); axR.set_xticklabels(["0", "50", "100"]); axR.set_yticklabels(["0", "50", "100"])
axR.tick_params(length=1.8, width=0.5, labelsize=F_NOTE + 0.4, colors=INK2, pad=1.4)
for sp_ in ("top", "right", "left", "bottom"):                               # a frame around the inset
    axR.spines[sp_].set_visible(True); axR.spines[sp_].set_color(INK3); axR.spines[sp_].set_linewidth(0.6)
axR.set_xlabel("Absent accepted (%)", fontsize=F_NOTE + 0.6, color=INK2, labelpad=1.4); axR.set_ylabel("Donors accepted (%)", fontsize=F_NOTE + 0.6, color=INK2, labelpad=1.4)
axR.text(0.95, 0.10, f"AUROC {g_auc:.3f}", fontsize=F_NOTE + 0.6, color=INK, fontweight="bold", ha="right", va="bottom", transform=axR.transAxes, zorder=5)
axR.text(float(np.mean(g_a >= g_t)) + 0.07, float(np.mean(g_P >= g_t)) - 0.04, f"{g_t:.4f}", fontsize=F_NOTE, color=INK, ha="left", va="top")
axG.legend(handles=[Patch(facecolor=C_SK_TINT, edgecolor=C_SK, lw=0.8, label=f"other sample of the donor ({int((g_P > 0).sum())} of {g_n})"),
                    Patch(facecolor=C_STR_TINT, edgecolor=C_STR, lw=0.8, label=f"nearest other people, donors absent ({g_n})")],
           loc="upper right", bbox_to_anchor=(0.99, 100 / G_TOPC - 0.015 * G_SC), fontsize=F_NOTE - 0.6, frameon=False, handlelength=1.1, handleheight=0.8, handletextpad=0.5, labelspacing=0.3, borderaxespad=0.0, labelcolor=INK2)
axG.set_xlabel("Similarity", fontsize=F_LAB, color=INK, labelpad=2)
axG.set_ylabel("Donors per 0.01 of similarity", fontsize=F_LAB, color=INK, labelpad=3)
# The searches whose nearest run of other people is more similar than every other sample (all from one cohort): a bracket under the label, and a thin arrow from the label to the bracket
g_big = g_a > g_P.max()
g_coh = {"CosteaPI_2017": "Costea"}
g_names = sorted(set(r["study"] for r, b in zip(g_rows, g_big) if b))
if g_big.any():                                                               # no such search once the rule of re-deposit has removed the Costea donors of Supplementary Note 1
    assert len(g_names) == 1 and g_names[0] in g_coh, g_names
    g_lo, g_hi = float(g_a[g_big].min()), float(g_a[g_big].max())
    g_by = g_ymax * 0.075                                                         # height of the bracket above the bars
    axG.plot([g_lo - 0.004, g_lo - 0.004, g_hi + 0.004, g_hi + 0.004], [g_by * 0.55, g_by, g_by, g_by * 0.55], color=INK3, lw=0.7, solid_capstyle="butt", zorder=5, clip_on=False)
    axG.annotate(f"{int(g_big.sum())} {g_coh[g_names[0]]} searches\n{g_lo:.2f} to {g_hi:.2f}", xy=((g_lo + g_hi) / 2, g_by + g_ymax * 0.012), xytext=((g_lo + g_hi) / 2, g_ymax * 0.135),                            # the label is centred above the bracket; the arrow goes from under the label to just above the bracket
                 fontsize=F_NOTE - 1.0, color=INK3, fontweight="bold", ha="center", va="bottom", linespacing=1.1, zorder=5, annotation_clip=False, bbox=dict(facecolor="white", edgecolor="none", pad=0.5),
                 arrowprops=dict(arrowstyle="-|>", color=INK3, lw=0.7, shrinkA=1.5, shrinkB=0, mutation_scale=5, relpos=(0.5, 0)))
fig.text(3.78 / W, (H_TOP + 0.03) / H, "H", fontsize=F_LET, fontweight="bold", color=INK, ha="left", va="top")


# ------------------------------------------------------------------ the middle row: C, share of donors whose other sample is in the top X; D, recovery in each cohort; the bottom row: E, recovery against the time between the first and the last sample
def wilson_bar(k, n):
    """The 95% Wilson score interval of the share k / n of n donors, in % (it exists for a share of 0 or 100%, and is asymmetric near them; the donors are taken as independent)."""
    p = k / n; z = 1.959964; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, 100 * (c - h)), min(100.0, 100 * (c + h))


days = np.array([float(r["days"]) for r in rows])
ok_first = np.array([r["rank_first_narrow"] == "1" for r in rows])                  # the first sample alone is first among the runs of the HGM set (not found counts as a failure)
weak = [c for c in cohorts if float(np.mean(K[byc[c]] == 0)) < 0.5]                  # the cohorts that recover less than half of their donors (none among the 7 cohorts kept)
is_weak = np.isin(study, weak)
NAMES = {"MehtaRS_2018": "Mehta", "CosteaPI_2017": "Costea", "NielsenHB_2014": "Nielsen", "VincentC_2016": "Vincent", "HansenLBS_2018": "Hansen", "PehrssonE_2016": "Pehrsson",
         "Bengtsson-PalmeJ_2015": "Bengtsson-Palme", "RaymondF_2016": "Raymond", "Heitz-BuschartA_2016": "Heintz-Buschart", "LiJ_2014": "Li", "HMP_2019_ibdmdb": "HMP IBD"}
BINS = [(0, 30, "\u226430"), (31, 90, "31\u201390"), (91, 180, "91\u2013180"), (181, 391, "181\u2013391")]                    # the last class takes the 3 donors at 366 to 391 days; no donor of the analysis is sampled further apart (Supplementary Table 12)

# ------------------------------------------------------------------ panel C: share of donors whose other sample is in the top X (cumulative distribution of the rank of the other sample among the runs of the HGM set)
XMAX_T = 5000                                                                      # the service returns the 5,000 nearest runs
rank_all = np.where(K < M, K + 1, M + 1).astype(int)                                # rank of the other sample among the runs of the HGM set (not found: beyond the 5,000 returned)
Xt = np.arange(1, XMAX_T + 1)
cum_share = lambda r: np.cumsum(np.bincount(np.minimum(r, XMAX_T + 1), minlength=XMAX_T + 2))[1:XMAX_T + 1] / len(r)      # share of donors with rank <= X, X = 1..5000
rec = 100 * cum_share(rank_all)
RANKS = np.unique(rank_all[rank_all <= XMAX_T])                                     # the share changes only at the ranks of the donors
lo_r, hi_r = np.percentile(100 * (RES_W @ (rank_all[:, None] <= RANKS[None, :]).astype(float)), [2.5, 97.5], axis=0)   # the resamples of the measured band
step_t = np.searchsorted(RANKS, Xt, side="right") - 1; lo_t, hi_t = lo_r[step_t], hi_r[step_t]
DY = 0.10 * (100 - math.floor(min(lo_t.min(), rec.min())))
D_LO = 70; D_CUT = 0.74 * (D_LO - 60) / (102 - 60)                         # panel D starts at 70% instead of 60%: the width saved (inches) goes to panel C
XC0 = 0.50; WC = 1.33 + D_CUT                                                                  # axes of panel C (inches)
H_CD = 1.72                                                                                   # height of the axes of C and D (inches), up to just under panel A
axC = fig.add_axes((XC0 / W, (H_C + 0.44) / H, WC / W, H_CD / H))
style(axC, "y")
axC.set_xscale("log")
axC.fill_between(Xt, lo_t, hi_t, step="post", color=C_SK_TINT, alpha=0.6, lw=0, zorder=1)
axC.plot(Xt, rec, drawstyle="steps-post", color=C_SK, lw=1.2, zorder=4)
for x_, col, ms_ in ((XMAX_T, INK3, 7.5), (1000, INK3, 7.5), (100, INK3, 7.5), (10, INK3, 7.5), (1, C_SK, 7.5)):             # big vermillion point at X = 1 (the first, as in panels A and B), big grey points at X = 10, 100, 1,000 and 5,000 (the whole list returned by the service)
    axC.plot([x_], [rec[x_ - 1]], "o", ms=ms_, color=col, mec="white", mew=0.7, zorder=6, clip_on=False)      # not clipped: the point of X = 5,000 sits at the edge of the axes
axC.text(1.7, rec[0] - DY, f"{rec[0]:.1f}%", fontsize=F_NOTE, color=C_SK_INK, fontweight="bold", ha="left", va="center")
axC.text(15.0, rec[9] - DY, f"{rec[9]:.1f}%", fontsize=F_NOTE, color=INK3, fontweight="bold", ha="left", va="center")
axC.text(100, rec[99] - 0.9 * DY, f"{rec[99]:.1f}%", fontsize=F_NOTE, color=INK3, fontweight="bold", ha="center", va="center")
axC.text(1000, rec[999] - 0.9 * DY, f"{rec[999]:.1f}%", fontsize=F_NOTE, color=INK3, fontweight="bold", ha="center", va="center")
axC.text(0.86 * XMAX_T, rec[XMAX_T - 1] + 0.95 * DY, f"{rec[XMAX_T - 1]:.1f}%", fontsize=F_NOTE, color=INK3, fontweight="bold", ha="right", va="center")      # above and left of the last point: the label of 1,000 is under its point
YLO =math.floor(min(lo_t.min(), rec.min()))                              # the y axis starts at the whole percent just under the lowest value drawn, band included: nothing is cut
axC.set_xlim(0.75, 6500); axC.set_ylim(YLO, 100)
axC.set_xticks([1, 10, 100, 1000, XMAX_T]); axC.set_xticklabels(["1", "10", "100", "1,000", "5,000"]); axC.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
from matplotlib.transforms import ScaledTranslation
for lab_, dx_ in zip(axC.get_xticklabels()[-2:], (-0.03, 0.03)):                       # the labels of 1,000 and 5,000 are 0.27 inch apart: each is moved 0.03 inch away from the other so that they do not touch
    lab_.set_transform(lab_.get_transform() + ScaledTranslation(dx_, 0, fig.dpi_scale_trans))
axC.set_yticks(list(range(5 * math.ceil(YLO / 5), 101, 5)))
axC.set_xlabel("Other sample in the top X", fontsize=F_LAB, color=INK, labelpad=2)
axC.set_ylabel("Donors recovered (%)", fontsize=F_LAB, color=INK, labelpad=3)
print(f"panel C (share of donors with the other sample in the top X): X = 1: {rec[0]:.1f}%; 2: {rec[1]:.1f}%; 3: {rec[2]:.1f}%; 5: {rec[4]:.1f}%; 10: {rec[9]:.1f}%; 100: {rec[99]:.1f}%; 172: {rec[171]:.1f}%; 1000: {rec[999]:.1f}%; 1721: {rec[1720]:.1f}%; 5000: {rec[-1]:.1f}%")

# ------------------------------------------------------------------ panel D: recovery in each cohort
gbp = np.array([float(r["query_gbp"]) for r in rows])
order = sorted(cohorts, key=lambda c: float(np.mean(K[byc[c]] == 0)))                       # the weakest at the bottom
XD0, WD = 2.72 + D_CUT, 0.74 - D_CUT                                           # the same scale per percent, the right edge where it was
G_X, N_X = 1 + 0.31 * 0.74 / WD, 1 + 0.59 * 0.74 / WD                          # the columns Gbp and n where they were
axD = fig.add_axes((XD0 / W, (H_C + 0.44) / H, WD / W, H_CD / H))
style(axD, None)
trD = blended_transform_factory(axD.transAxes, axD.transData)
for i, c in enumerate(order):
    idx = byc[c]; n = len(idx); k = int(np.sum(K[idx] == 0)); a, b = wilson_bar(k, n); col = C_RED if c in weak else C_SK
    assert a >= D_LO, (c, a)                                                          # every interval stays inside the axis
    axD.plot([a, b], [i, i], color=col, lw=0.9, solid_capstyle="butt", zorder=2)
    axD.plot([100 * k / n], [i], "o", ms=2.6 + 0.19 * math.sqrt(n), color=col, mec="white", mew=0.4, zorder=3)
    axD.text(G_X, i, f"{np.median(gbp[idx]):.2g}", transform=trD, ha="right", va="center", fontsize=F_NOTE, color=INK2)
    axD.text(N_X, i, str(n), transform=trD, ha="right", va="center", fontsize=F_NOTE, color=INK3)
axD.text(G_X, len(order) - 0.15, "Gbp", transform=trD, ha="right", va="bottom", fontsize=5.8, color=INK3, clip_on=False)
axD.text(N_X, len(order) - 0.15, "n", transform=trD, ha="right", va="bottom", fontsize=5.8, color=INK3, clip_on=False)
axD.set_yticks(range(len(order))); axD.set_yticklabels([NAMES[c] for c in order], fontsize=F_TICK, color=INK)
axD.set_xlim(D_LO, 102); axD.set_ylim(-0.7, len(order) - 0.3); axD.set_xticks(list(range(D_LO, 101, 10))); axD.grid(axis="x", color=GUIDE, linewidth=0.45, zorder=0)
axD.set_xlabel("Donors recovered (%)", fontsize=F_LAB, color=INK, labelpad=2)
axD.axvline(100 * float(np.mean(K == 0)), color=INK3, lw=0.5, ls=(0, (1.5, 2.0)), zorder=1)

# ------------------------------------------------------------------ panels E, F and G on one pair of axes: the same vertical scale and the same width of the bars
# E: recovery against the time between the first and the last sample; F: recovery against the bases of the queried run (Gbp); G: the donors of Raymond with and without a course of cephalosporins between their samples. One axes, so that a bar has the same width in the three.
BINS_F = [(0, 0.5, "<0.5"), (0.5, 1, "0.5\u20131"), (1, 2, "1\u20132"), (2, 4, "2\u20134"), (4, 1e9, "\u22654")]       # Gbp of the queried run (two classes of Supplementary Table 13 under 0.5 Gbp merged: 13 donors)
XEF0, WEF = 0.62, 3.05                                                                  # axes of E, F and G (inches); the letter of E is 0.55 to the left
X_F0 = len(BINS) + 0.8                                                                  # position of the first bar of F: 0.8 bar of gap after the last bar of E
X_R0 = X_F0 + len(BINS_F) - 1 + 1.8                                                     # the two bars of G: 0.8 bar of gap after the last bar of F
X_END = X_R0 + 1 + 0.6
UNIT = WEF / (X_END + 0.6)                                                              # inches between two bars
axE = fig.add_axes((XEF0 / W, YB / H, WEF / W, 1.60 / H))
style(axE, "y")
_mix = lambda c, a: tuple(a * v + (1 - a) for v in matplotlib.colors.to_rgb(c))      # the colour c over white at opacity a, opaque so that the gridlines do not show through
BAR_A = 0.70                                                                            # the bars of E, F and G are a softened vermillion with a full-colour outline
BW = 0.62                                                                               # width of a bar in units of the spacing between two bars: the same in E, F and G
trE = blended_transform_factory(axE.transData, axE.transAxes)
stats_c, stats_f = [], []
for j, (lo_d, hi_d, lab) in enumerate(BINS):
    in_bin = (days >= lo_d) & (days <= hi_d)
    n = int(in_bin.sum()); k = int(ok_first[in_bin].sum()); stats_c.append((lab, n, k))
    if n >= 1:                                                                          # every class is drawn, the small ones with their wide interval
        a, b = wilson_bar(k, n)
        axE.bar(j, 100 * k / n, width=BW, facecolor=_mix(C_SK, BAR_A), edgecolor=C_SK, linewidth=0.7, zorder=2)
        axE.plot([j, j], [a, b], color=INK, lw=0.8, solid_capstyle="butt", zorder=3)
    axE.text(j, -0.04, str(n), transform=trE, ha="center", va="top", fontsize=5.2, color=INK3)
for j, (lo_g, hi_g, lab) in enumerate(BINS_F):
    sel = (gbp >= lo_g) & (gbp < hi_g); n = int(sel.sum()); k = int(np.sum(K[sel] == 0)); stats_f.append((lab, n, k))
    a, b = wilson_bar(k, n)
    axE.bar(X_F0 + j, 100 * k / n, width=BW, facecolor=_mix(C_SK, BAR_A), edgecolor=C_SK, linewidth=0.7, zorder=2)
    axE.plot([X_F0 + j, X_F0 + j], [a, b], color=INK, lw=0.8, solid_capstyle="butt", zorder=3)
    axE.text(X_F0 + j, -0.04, str(n), transform=trE, ha="center", va="top", fontsize=5.2, color=INK3)
RAY = json.load(open(os.path.join(HERE, "aggregates", "raymond_antibiotic.json")))                  # panel G: raymond_antibiotic.py, the donors of Raymond by the study_condition of their last sample
RAY_ROWS = [("cephalosporins", "yes"), ("control", "no")]
assert sum(RAY["groups"][g_]["donors"] for g_, _ in RAY_ROWS) == int(np.sum(study == RAY["study"])) and sum(RAY["groups"][g_]["other_sample_first"] for g_, _ in RAY_ROWS) == int(np.sum(K[study == RAY["study"]] == 0)), "aggregates/raymond_antibiotic.json does not agree with data/adult_ranks.csv"
stats_r = []
for j, (g_, lab) in enumerate(RAY_ROWS):
    n = RAY["groups"][g_]["donors"]; k = RAY["groups"][g_]["other_sample_first"]; stats_r.append((lab, n, k))
    a, b = wilson_bar(k, n)
    axE.bar(X_R0 + j, 100 * k / n, width=BW, facecolor=_mix(C_SK, BAR_A), edgecolor=C_SK, linewidth=0.7, zorder=2)
    axE.plot([X_R0 + j, X_R0 + j], [a, b], color=INK, lw=0.8, solid_capstyle="butt", zorder=3)
    axE.text(X_R0 + j, -0.04, str(n), transform=trE, ha="center", va="top", fontsize=5.2, color=INK3)
axE.set_xticks(list(range(len(BINS))) + [X_F0 + j for j in range(len(BINS_F))] + [X_R0, X_R0 + 1])
axE.set_xticklabels([b[2] for b in BINS] + [b[2] for b in BINS_F] + [lab for _, lab in RAY_ROWS], fontsize=5.2, rotation=45, ha="right"); axE.tick_params(axis="x", pad=7.5)
axE.set_xlim(-0.6, X_END); axE.set_ylim(0, 108); axE.set_yticks(list(range(0, 101, 10))); axE.spines["left"].set_bounds(0, 100)
axE.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator([v_ for v_ in range(5, 101, 5)])); axE.tick_params(axis="y", which="minor", length=1.6, width=0.45, color=RULE)                         # a tick every 5, as the minor ticks of panel B (the positions of the labels are left to the major ticks)
axE.set_ylabel("Donors recovered (%)", fontsize=F_LAB, color=INK, labelpad=3)
x_in = lambda x: (XEF0 + (x + 0.6) * UNIT) / W                                         # figure fraction of a position on the axes of E, F and G
YL = YB - 0.373 - 0.035 + 0.012                                                        # top of the axis labels of E, F and G: 0.035 below the rotated labels (which end 0.373 below the axis)
fig.text(x_in((len(BINS) - 1) / 2), YL / H, "Days between the\nfirst and last sample", fontsize=F_LAB, color=INK, ha="center", va="top", linespacing=1.15)
fig.text(x_in(X_F0 + (len(BINS_F) - 1) / 2), YL / H, "Bases of the\nqueried run (Gbp)", fontsize=F_LAB, color=INK, ha="center", va="top", linespacing=1.15)
fig.text(x_in(X_R0 + 0.5), YL / H, "Cephalosporins\n(Raymond)", fontsize=F_LAB, color=INK, ha="center", va="top", linespacing=1.15)
print("panel G (other sample first, the donors of Raymond with a course of cephalosporins (yes) and without (no)):", "; ".join(f"{lab} {k}/{n}" for lab, n, k in stats_r))
print("panel E (first sample first, by days between the samples):", "; ".join(f"{lab} {k}/{n}" for lab, n, k in stats_c))
print("panel F (other sample first, by Gbp of the queried run):", "; ".join(f"{lab} {k}/{n}" for lab, n, k in stats_f), f"| width of a bar: {BW * UNIT:.3f} in in E and in F")


fig.text(0.008, 1 - 0.0534 / H, "A", fontsize=F_LET, fontweight="bold", color=INK, ha="left", va="top")
fig.text((XB0 - 0.40) / W, 1 - 0.0534 / H, "B", fontsize=F_LET, fontweight="bold", color=INK, ha="left", va="top")
fig.text(0.008, (H_C + 0.44 + H_CD + 0.19) / H, "C", fontsize=F_LET, fontweight="bold", color=INK, ha="left", va="top")
fig.text((1.80 + D_CUT) / W, (H_C + 0.44 + H_CD + 0.19) / H, "D", fontsize=F_LET, fontweight="bold", color=INK, ha="left", va="top")
fig.text((XEF0 - 0.55) / W, (H_C + 0.03) / H, "E", fontsize=F_LET, fontweight="bold", color=INK, ha="left", va="top")
fig.text(x_in(X_F0 - 0.5) - 0.14 / W, (H_C + 0.03) / H, "F", fontsize=F_LET, fontweight="bold", color=INK, ha="left", va="top", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))      # the white ground interrupts the gridline behind the letter
fig.text(x_in(X_R0 - 0.5) - 0.14 / W, (H_C + 0.03) / H, "G", fontsize=F_LET, fontweight="bold", color=INK, ha="left", va="top", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))

fig.savefig(PREFIX + ".png", dpi=DPI, pil_kwargs={"compress_level": 9})      # compression level 9, the same pixels: with the default level the compressed bytes of this figure happened to hold a string that the leak scan of the package flags
fig.savefig(PREFIX + ".pdf")
print("wrote", PREFIX + ".png and .pdf")
