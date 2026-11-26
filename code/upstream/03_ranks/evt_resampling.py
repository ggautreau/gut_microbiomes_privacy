#!/usr/bin/env python3
"""Resampling band of the extreme-value model of the similarity tail, for the tail drawn in Fig. 1B (Supplementary Fig. 6). Aggregates only.

    python3 -I evt_resampling.py EVT_SIMILARITIES.npz EVT_SUMMARY.json OUT.json [WORKERS]

EVT_SIMILARITIES.npz  output of evt_similarities.py, as read by evt_model.py
EVT_SUMMARY.json      output of evt_model.py: the tail drawn in Fig. 1B, its curve and its intervals, which this script reproduces before it writes anything
OUT.json              for the tail drawn: the recovery rate of the model on the grid of pool sizes with its 95% interval at each size, and the bounds of the interval at the
                      population of the world on the first B resamples (B = 100, 200, 400 ...) with their Monte Carlo standard errors

The model is fitted once, on all the donors, by the functions of evt_model.py (cells, integrate_sigma, fit, stay, stay_cells and resample_values are those of that script).
A resample draws the cohorts with replacement, then the donors within each drawn cohort; the distribution Beta(a, b) of the shapes is estimated again from the donors drawn
(each weighted by the number of times it is drawn), which gives each donor a new posterior of its shape, and the probabilities of the donors drawn are averaged. The
posterior of the scale given each cell of the shape does not depend on a and b, so the survival of each donor given each cell is computed once.
The number of resamples follows the rule of extrapolation_bootstrap.py: B doubles from 100 until the Monte Carlo standard errors of the two bounds at the population of the
world (2,000 draws with replacement among the B resamples) are at most 5% of the width of the interval. The resamples are the first B of
the 100,000 that evt_model.py draws (the same generator, in the same order). This script stops if the rule does not select the number that evt_model.py uses, or if it does
not find the curve and the intervals of EVT_SUMMARY.json again. Numbers are written with 10 decimals."""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import json, multiprocessing, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                   # python -I does not put the script's folder on the path
from evt_model import fit, stay, stay_many, stay_cells, resample_values

M = 291_881
WORLD = 8e9
B_START = 100
N_SUB = 2_000
SE_OVER_W_THR = 0.05
NTAB = (1e6, 1e7, 1e8, 1e9, 8e9, 1e10, 1e11)

_RES_F = None
_RES_ALIVE = None


def _stay_at(v):
    return stay(_RES_F, _RES_ALIVE, v - M)


if __name__ == "__main__":
    SIM, SUMMARY, OUT = sys.argv[1:4]
    AVAIL_CPUS = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    WORKERS = int(sys.argv[4]) if len(sys.argv) > 4 else max(1, AVAIL_CPUS - 2)
    ev = json.load(open(SUMMARY))
    DRAWN = int(ev["drawn"]); USED = int(ev["resamples"]); POOL = int(ev["resamples_drawn"]); LO, HI = ev["ranges"]["free"]; N_GRID = np.array(ev["N"], float); ref = ev["models"]["free"][str(DRAWN)]
    d = np.load(SIM, allow_pickle=True)
    TOP, n_all, s_all, k_all, study = d["top"].astype(np.float64), d["n"].astype(np.float64), d["s"], d["k"], d["study"]
    D, K = TOP.shape
    alive_all = k_all == 0
    assert D == ev["donors"] and int(alive_all.sum()) == ev["recovered"] and ev["hgm_runs"] == M

    f = fit(TOP, n_all, s_all, DRAWN, LO, HI, workers=WORKERS)
    assert abs(f["a"] - ref["a"]) < 1e-9 and abs(f["b"] - ref["b"]) < 1e-9, "the distribution of the shapes is that of evt_model.py"
    P = stay_many(f, alive_all, N_GRID - M, workers=WORKERS)                         # donors x pool sizes
    _RES_F = f; _RES_ALIVE = alive_all
    with multiprocessing.get_context("fork").Pool(max(1, min(WORKERS, len(NTAB)))) as pool:
        p_tab = dict(zip(NTAB, pool.map(_stay_at, NTAB, chunksize=1)))
    point = 100 * P.mean(0)
    assert float(np.max(np.abs(point - np.array(ref["curve"])))) < 1e-9, "the curve is that of evt_model.py"
    assert abs(100 * float(p_tab[WORLD].mean()) - ref["at_world"]) < 1e-9 and abs(100 * float(p_tab[1e7].mean()) - ref["at_1e7"]) < 1e-9

    # ---------------------------------------------------------------- the resamples of evt_model.py: the cohorts with replacement, then the donors within each drawn cohort; the first USED give the intervals
    cohorts = np.unique(study); members = {c: np.where(study == c)[0] for c in cohorts}
    rng = np.random.default_rng(3)
    m_c = rng.multinomial(len(cohorts), [1 / len(cohorts)] * len(cohorts), size=POOL)
    W = np.zeros((POOL, D))
    for j, c in enumerate(cohorts):
        W[:, members[c]] = rng.multinomial(m_c[:, j] * len(members[c]), [1 / len(members[c])] * len(members[c]))
    W /= W.sum(1, keepdims=True)
    S_tab = np.concatenate([stay_cells(f, alive_all, [n_ - M], workers=WORKERS) for n_ in NTAB], 1)            # one pool size per call, as in evt_model.py
    S_grid = stay_cells(f, alive_all, N_GRID - M, workers=WORKERS)
    S_all = np.concatenate([S_tab, S_grid], 1)
    done = []                                                                       # recovery rates of the resamples computed so far (a and b estimated again on each)

    def values(b_):
        """Recovery rates (%) of the first b_ resamples at the pool sizes of NTAB, then of N_GRID."""
        have = sum(len(x) for x in done)
        if b_ > have:
            done.append(resample_values(f, S_all, W[have:b_], workers=WORKERS)[0])
        return np.concatenate(done)[:b_]

    V = values(USED)
    at_n = {}
    for j, n_ in enumerate(NTAB):
        lo_, hi_ = (float(v) for v in np.percentile(V[:, j], [2.5, 97.5]))
        if n_ in (1e7, WORLD):
            key = "interval_world" if n_ == WORLD else "interval_1e7"
            assert abs(lo_ - ref[key][0]) < 1e-9 and abs(hi_ - ref[key][1]) < 1e-9, f"the interval at {n_:g} is that of evt_model.py"
        at_n[str(int(n_))] = dict(est=100 * float(p_tab[n_].mean()), lo=lo_, hi=hi_, width=hi_ - lo_)
    band = np.percentile(V[:, len(NTAB):], [2.5, 97.5], axis=0).T                  # one pool size at a time, as evt_model.py computes an interval
    J_W = NTAB.index(WORLD)

    # ---------------------------------------------------------------- the number of resamples: the bounds on the first B, as in extrapolation_bootstrap.py
    trace, b_cur, step = [], B_START, 0
    while True:
        v = values(b_cur)[:, J_W]; lo_, hi_ = (float(x) for x in np.percentile(v, [2.5, 97.5]))
        idx = np.random.default_rng(20260 + step).integers(0, b_cur, size=(N_SUB, b_cur))
        qs = np.percentile(v[idx], [2.5, 97.5], axis=1)
        se_lo, se_hi = float(np.std(qs[0], ddof=1)), float(np.std(qs[1], ddof=1)); width = hi_ - lo_
        r_se = max(se_lo, se_hi) / width
        met = r_se <= SE_OVER_W_THR
        trace.append(dict(B=int(b_cur), lo=lo_, hi=hi_, width=width, se_lo=se_lo, se_hi=se_hi, max_se_over_width=r_se, passed=bool(met)))
        print(f"  first {b_cur:6,d} resamples: 95% interval [{lo_:5.2f}, {hi_:5.2f}], standard errors ({se_lo:.3f}, {se_hi:.3f}), largest over the width {r_se:.4f}{' -> criterion met' if met else ''}", flush=True)
        if met or b_cur * 2 > POOL:
            break
        b_cur *= 2; step += 1
    assert trace[-1]["passed"] and trace[-1]["B"] == USED, f"the rule selects {trace[-1]['B']} resamples and evt_model.py uses {USED}"
    w_ = at_n[str(int(WORLD))]
    assert abs(trace[-1]["lo"] - w_["lo"]) < 1e-9 and abs(trace[-1]["hi"] - w_["hi"]) < 1e-9, "the last step of the rule gives the interval of evt_model.py"

    rnd = lambda x: round(float(x), 10)


    def rounded(o):
        if isinstance(o, dict):
            return {k_: rounded(v_) for k_, v_ in o.items()}
        if isinstance(o, list):
            return [rounded(v_) for v_ in o]
        return rnd(o) if isinstance(o, float) else o


    out = dict(meta=dict(donors=int(D), recovered=int(alive_all.sum()), M=M, world=WORLD, tail=DRAWN, resamples=USED, resamples_drawn=POOL, seed=3, b_start=B_START, sub_draws=N_SUB,
                         se_over_width_threshold=SE_OVER_W_THR),
               x_ext=[float(v) for v in N_GRID],
               models=dict(evt=dict(B=USED, converged=True, point_curve=[float(v) for v in point],
                                    lo_curve=[float(v) for v in band[:, 0]], hi_curve=[float(v) for v in band[:, 1]], at_N=at_n, doubling_trace=trace)))
    json.dump(rounded(out), open(OUT, "w"), indent=1)
    print(f"tail of {DRAWN:,}: {w_['est']:.1f} % at 8e9, 95% interval {w_['lo']:.1f}-{w_['hi']:.1f} on the first {USED:,} of {POOL:,} resamples, the number that the rule selects")
    print("wrote", OUT)