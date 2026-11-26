#!/usr/bin/env python3
"""Thinning test of the three models of Figure 1B (Supplementary Note 3, Supplementary Table 16 and the errors in parentheses of the table of Fig. 1B), the small pools of Note 3 and the share
among the 10 nearest runs (aggregates only).

    python3 -I extrapolation_tables.py LISTS.npz OUT.json [B [WORKERS]]          # B thinnings at each fraction (400), in WORKERS processes (default all available cores minus two, at least 1); the result does not depend on WORKERS

LISTS.npz: for each donor of the analysis, `partner` (similarity of the donor's other sample, 0 = not among the 5,000 returned), `k_rows` (number of the 291,881 runs of the HGM set more similar than
it), `imp_nar` (similarities of the runs of other people in the HGM set, descending) and `study` (extract_lists.py).
  10 nearest of N : the other sample is among the 10 nearest runs of a pool of N <= M runs: hypergeometric on k (at most 9 runs of the HGM set above the other sample among the N - 1 runs drawn for the pool).
  small pools     : the measured share, and its upper bound when a donor whose other sample is not returned has only the runs of other people of its list above that sample.
                    The two curves on the grid of Fig. 1B (curves), and the three models at the population of the world fitted to each, for the windows drawn (at_8e9) and for the five
                    windows (at_8e9_windows); Supplementary Fig. 5 draws the fits to the measured curve, the column Upper of Supplementary Table 16 gives those to the upper bound.
  thinning test   : thin the runs of other people of each donor by Bernoulli(f) (M' = f M), keep the other samples; from the thinned lists alone predict the recovery rate measured at N = M
                    (10, 33 and 100 times M') with the three models of Fig. 1B, the Pitman-Yor correctness model (PYC-MB) of Rocher et al. 2025, the straight line in log N and the power law of
                    the miss rate 1 - a N^b, each fitted on the same share of the thinned set as [lo, M] is of the HGM set for lo = 2, 10, 100, 1,000 and 10,000 (key `holdout_windows`: Supplementary
                    Fig. 5 and Table 16). The fits drawn in Fig. 1B are the three models on [1,000, M], where the measured share and its upper bound for small pools differ
                    by at most 0.2 points (small pools above); key `holdout` repeats them. The test counts a donor whose other sample is not returned as a failure at every pool size, as the measured
                    curve does.
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import json, multiprocessing, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                   # python -I does not put the script's folder on the path
from scipy.special import gammaln
from scipy.stats import hypergeom
import pyc_mb as PY                                                              # Pitman-Yor correctness model of Rocher et al. 2025

LISTS, OUT = sys.argv[1], sys.argv[2]
B = int(sys.argv[3]) if len(sys.argv) > 3 else 400
AVAIL_CPUS = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
WORKERS = int(sys.argv[4]) if len(sys.argv) > 4 else max(1, AVAIL_CPUS - 2)
M = 291_881                                                                      # runs of the HGM set
d = np.load(LISTS, allow_pickle=True)
nar, partner, k, study, n_list = d["imp_nar"], d["partner"], d["k_rows"], d["study"], d["n_list"]
D = len(partner); found = partner > 0
out = dict(donors=int(D), found=int(found.sum()), k0=int((found & (k == 0)).sum()))

# ------------------------------------------------------------------ the 10 nearest runs of a pool of N <= M runs
TOP = 10


def top_inside(N):
    return float(hypergeom.cdf(TOP - 1, M, k[found], int(N) - 1).sum() / D)


meas10 = top_inside(M)
out["measured_top10"] = 100 * meas10
out["top10_inside"] = {str(N): 100 * top_inside(N) for N in (10, 100, 1000, 10_000, 100_000, M)}
cohorts = np.unique(study); members = {c: np.where(study == c)[0] for c in cohorts}
rng = np.random.default_rng(3); RES = 100_000; n_cohorts = len(cohorts)                           # cohorts resampled, then donors: a cohort drawn m times adds m binomials
m_c = rng.multinomial(n_cohorts, [1 / n_cohorts] * n_cohorts, size=RES); hits, tot = np.zeros(RES), np.zeros(RES)
for j, c in enumerate(cohorts):
    n_c = len(members[c]); hits += rng.binomial(m_c[:, j] * n_c, float(np.mean((found & (k <= TOP - 1))[members[c]]))); tot += m_c[:, j] * n_c
out["top10_ci_cohorts"] = [float(v) for v in np.percentile(100 * hits / tot, [2.5, 97.5])]

WIN = (2, 10, 100, 1000, 10_000)                                                   # the fit windows of the table of Fig. 1B: the same share of the thinned set as [lo, M] is of the HGM set (40 sizes from 2 for N >= 2, 30 otherwise)
DRAWN = dict(pycmb=1000, line=1000, miss=1000)                                     # the window of each fit drawn in Fig. 1B: N >= 1,000 for the three models, where the measured share is exact to 0.2 points

# ------------------------------------------------------------------ small pools: the measured share is a lower bound
# A donor whose other sample is beyond the list returned counts as a failure at every pool size. Such a donor has at least the runs of other
# people of its list (in the HGM set) above the other sample; with exactly that number the share is the upper bound below. The two meet for large pools.
s_list = np.array([len(nar[i]) for i in range(D)])


def first_share(N, kk):
    """Share of donors first in a pool of N; kk, runs of the HGM set above the other sample of each donor, -1 for a donor that fails at every pool size."""
    kv = kk[kk >= 0].astype(float); m = N - 1
    with np.errstate(all="ignore"):
        lg = gammaln(M - kv + 1) - gammaln(M - kv - m + 1) - (gammaln(M + 1) - gammaln(M - m + 1))
    return float(np.where(M - kv >= m, np.exp(lg), 0.0).sum() / D)


LIST_MAX = 5000                                                                  # the service returns at most 5,000 runs and stops earlier at a similarity of 0.05
# short and full: the runs of other people in the lists of the donors whose other sample is not returned, for the lists that stop at the similarity of 0.05 and for the lists of 5,000 runs
out["small_pools"] = dict(not_found=int((~found).sum()), list_runs=dict(min=int(s_list[~found].min()), median=float(np.median(s_list[~found])), max=int(s_list[~found].max()),
                                                                         short=sorted(int(v) for v in s_list[~found & (n_list < LIST_MAX)]), full=sorted(int(v) for v in s_list[~found & (n_list >= LIST_MAX)])),
                          share={str(N): dict(measured=100 * first_share(N, np.where(found, k, -1)), upper=100 * first_share(N, np.where(found, k, s_list))) for N in (2, 10, 100, 1000, 10_000, 100_000, M)})
pa_rec = found & (k == 0)                                                        # Fig. 1A: the other sample of the recovered donors and their five nearest runs of other people, medians
out["panel_a"] = dict(recovered=int(pa_rec.sum()), other_sample=float(np.median(partner[pa_rec])),
                      runs_of_other_people=[float(np.median([nar[i][j] for i in np.where(pa_rec)[0] if len(nar[i]) > j])) for j in range(5)])

# ------------------------------------------------------------------ thinning test
def plug(Ns, kp, Mst):
    kk = kp[found].astype(float); res = []
    for n in Ns:
        m = n - 1
        with np.errstate(all="ignore"):
            lg = gammaln(Mst - kk + 1) - gammaln(Mst - kk - m + 1) - (gammaln(Mst + 1) - gammaln(Mst - m + 1))
        res.append(np.where(Mst - kk >= m, np.exp(lg), 0.0).sum() / D)
    return np.array(res)


def fit_miss(Nfit, R):
    """Power law of the miss rate, 1 - R = a N^b (the miss rate is the share of donors whose other sample is not first), least squares on the logarithms over its window."""
    c = np.polyfit(np.log(np.asarray(Nfit, float)), np.log(np.clip(1 - np.asarray(R, float), 1e-12, None)), 1)
    return lambda N: float(1 - np.exp(np.polyval(c, np.log(float(N)))))


def fit_pycmb(Nfit, R):
    """PYC-MB of Rocher et al. 2025 (two parameters) on the window given; nan when the fit fails."""
    try:
        f = PY.fit(np.asarray(Nfit, float), R)
        return lambda N: float(PY.kappa(float(N), f))
    except Exception:
        return lambda N: float("nan")


def holdout_one(task):
    """One thinning at fraction f, seeded by (f, b) so that the result does not depend on the number of processes."""
    f, b = task
    rng = np.random.default_rng([20261022, int(round(1000 * f)), b])
    Mp = int(round(f * M)); target = M
    thin = [x[rng.random(len(x)) < f] for x in nar]
    kp = np.array([(int((t > partner[i]).sum()) if partner[i] > 0 else -1) for i, t in enumerate(thin)])
    pred = {}
    for lo in WIN:                                                                       # the three models of Fig. 1B fitted for N >= 2, 10, 100, 1,000 and 10,000 (Supplementary Fig. 5 and Table 16)
        g_ = np.unique(np.round(np.logspace(np.log10(2), np.log10(Mp + 1), 40)).astype(int)) if lo == 2 else np.unique(np.round(np.logspace(np.log10(max(2, Mp * lo / M)), np.log10(Mp + 1), 30)).astype(int))
        R_ = plug(g_, kp, Mp)
        pred[f"line|{lo}"] = float(np.polyval(np.polyfit(np.log2(g_), R_, 1), np.log2(target)))
        pred[f"pycmb|{lo}"] = fit_pycmb(g_, R_)(target)
        pred[f"miss|{lo}"] = fit_miss(g_, R_)(target)
    return pred


# the three models of Fig. 1B, with the windows drawn, fitted to the measured curve and to its upper bound for small pools (see "small pools" above)
NSB = np.unique(np.round(np.logspace(np.log10(2), np.log10(M), 70)).astype(int))          # the grid of Fig. 1B and of Supplementary Fig. 5
out["small_pools"]["at_8e9"] = {}                                                         # the three models of Fig. 1B at the population of the world, fitted to the measured curve and to its upper bound
for nm_, kk_ in (("measured", np.where(found, k, -1)), ("upper", np.where(found, k, s_list))):
    Rb = np.array([first_share(int(N), kk_) for N in NSB]); sel_ = {m_: NSB >= lo_ for m_, lo_ in DRAWN.items()}
    out["small_pools"]["at_8e9"][nm_] = dict(pycmb=100 * fit_pycmb(NSB[sel_["pycmb"]], Rb[sel_["pycmb"]])(8e9),
                                             line=100 * float(np.polyval(np.polyfit(np.log2(NSB[sel_["line"]]), Rb[sel_["line"]], 1), np.log2(8e9))),
                                             miss=100 * fit_miss(NSB[sel_["miss"]], Rb[sel_["miss"]])(8e9))
# the same for each of the five windows: how much a fit that starts in small pools depends on the donors whose other sample is not returned (not clipped: a miss rate above 1 gives a negative share)
out["small_pools"]["at_8e9_windows"] = {}; out["small_pools"]["curves"] = dict(N=[int(n_) for n_ in NSB])
for nm_, kk_ in (("measured", np.where(found, k, -1)), ("upper", np.where(found, k, s_list))):
    Rb = np.array([first_share(int(N), kk_) for N in NSB]); o_ = out["small_pools"]["at_8e9_windows"][nm_] = {}
    out["small_pools"]["curves"][nm_] = [float(v_) for v_ in Rb]                          # the share on the grid: the measured curve is drawn in Supplementary Fig. 5, its upper bound gives the column Upper of Supplementary Table 16 (aggregates only)
    for lo_ in WIN:
        s_ = NSB >= lo_
        o_[f"pycmb|{lo_}"] = 100 * fit_pycmb(NSB[s_], Rb[s_])(8e9)
        o_[f"line|{lo_}"] = 100 * float(np.polyval(np.polyfit(np.log2(NSB[s_]), Rb[s_], 1), np.log2(8e9)))
        o_[f"miss|{lo_}"] = 100 * fit_miss(NSB[s_], Rb[s_])(8e9)
out["holdout"] = {}
with multiprocessing.get_context("fork").Pool(max(1, min(WORKERS, B))) as pool:
    for f in (0.1, 0.03, 0.01):
        Mp = int(round(f * M)); ref = float(plug([M], k, M)[0])
        res = pool.map(holdout_one, [(f, b) for b in range(B)], chunksize=1)
        win_ = {m: np.array([r[m] for r in res], float) for m in res[0]}
        preds = {m_: win_[f"{m_}|{lo_}"] for m_, lo_ in DRAWN.items()}                  # the fits drawn in Fig. 1B: the same predictions as their window
        out.setdefault("holdout_windows", {})[str(f)] = dict(error={m: 100 * (float(np.nanmean(v)) - ref) for m, v in win_.items()}, sd={m: 100 * float(np.nanstd(v)) for m, v in win_.items()},
                                                          failed={m: int(np.isnan(v).sum()) for m, v in win_.items() if np.isnan(v).any()})
        out["holdout"][str(f)] = dict(Mprime=Mp, factor=float(M / Mp), measured=100 * ref, error={m: 100 * (float(np.nanmean(v)) - ref) for m, v in preds.items()}, sd={m: 100 * float(np.nanstd(v)) for m, v in preds.items()},
                                      failed={m: int(np.isnan(v).sum()) for m, v in preds.items() if np.isnan(v).any()}, thinnings=B, windows=dict(DRAWN))
json.dump(out, open(OUT, "w"), indent=1)
print("wrote", OUT)
