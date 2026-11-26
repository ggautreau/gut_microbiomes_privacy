#!/usr/bin/env python3
"""Bootstrap confidence intervals for all three extrapolation models of Fig. 1B using adaptive doubling.

    python3 -I extrapolation_bootstrap.py ADULT_RANKS.csv OUTPUT.json [WORKERS]

For each of the three models fitted to the measured curve for N >= 1,000 (Pitman-Yor correctness PYC-MB,
straight line in log N, power law of the miss rate), runs the cohort-then-donor resampling scheme of
figure1_brief_communication.py (default_rng(1)) starting at B = 100 resamples and doubling B until
the Monte Carlo standard error of both the 2.5th and 97.5th percentile bounds at N = 8e9
(estimated from 2,000 sub-bootstrap draws) satisfies:
    max(SE_lo, SE_hi) / width <= 0.05.
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import csv
import json
import math
import multiprocessing as mp
import sys
import warnings
import numpy as np
from scipy.special import gammaln

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pyc_mb as PY

M = 291_881
WORLD = 8e9
FIT_FROM = 1_000
B_POOL = 100_000
B_START = 100
B_MAX_PY = 6_400
B_MAX_FAST = 51_200
N_SUB = 2_000
SE_OVER_W_THR = 0.05
NTAB = (1e6, 1e7, 1e8, 1e9, 8e9, 1e10, 1e11)

CSV_IN, JSON_OUT = sys.argv[1:3]
AVAIL_CPUS = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else max(1, AVAIL_CPUS - 2)

with open(CSV_IN, newline="") as fh:
    rows = [r for r in csv.DictReader(fh) if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]

D = len(rows)
study = np.array([r["study"] for r in rows])
K = np.array([(int(r["rank_any_narrow"]) - 1) if r["rank_any_narrow"] != "" else M for r in rows], float)


def p_first(k, N):
    """Probability that a donor with k runs more similar than the other sample is first in a pool of N runs."""
    n1 = N - 1
    ok = (M - k) >= n1
    out = np.zeros(len(k))
    out[ok] = np.exp(gammaln(M - k[ok] + 1) - gammaln(M - k[ok] - n1 + 1) - gammaln(M + 1) + gammaln(M - n1 + 1))
    return out


# The same grids as figure1_brief_communication.py
Ns = np.unique(np.round(np.logspace(np.log10(2), np.log10(M), 60)).astype(int))
NF = np.unique(np.round(np.logspace(math.log10(2), math.log10(M), 70)).astype(int))
FIT_SEL = NF >= FIT_FROM
X_FIT = NF[FIT_SEL].astype(float)
X_EXT = np.logspace(math.log10(M), 12, 200)
NTAB_ARR = np.asarray(NTAB, float)
EVAL_GRID = np.concatenate([X_EXT, NTAB_ARR])
N_EXT = len(X_EXT)

# Two-stage cohort-then-donor multinomial weights, identical RNG sequence to figure1_brief_communication.py
rng = np.random.default_rng(1)
cohorts = sorted(set(study))
byc = {c: np.where(study == c)[0] for c in cohorts}
m_c = rng.multinomial(len(cohorts), [1 / len(cohorts)] * len(cohorts), size=B_POOL)
RES_W = np.zeros((B_POOL, D))
for j, c in enumerate(cohorts):
    RES_W[:, byc[c]] = rng.multinomial(m_c[:, j] * len(byc[c]), [1 / len(byc[c])] * len(byc[c]))
RES_W /= RES_W.sum(axis=1, keepdims=True)

P_Ns = np.stack([p_first(K, N) for N in Ns], axis=1)
R_Ns = P_Ns.mean(axis=0)
lo_Ns, hi_Ns = np.percentile(RES_W @ P_Ns, [2.5, 97.5], axis=0)

P1 = np.stack([p_first(K, N) for N in NF], axis=1)
R_fit_pt = P1.mean(axis=0)[FIT_SEL]
Y1_fit = (RES_W @ P1)[:, FIT_SEL]


# Worker helpers for parallel PYC-MB fitting
_G_X_FIT = X_FIT
_G_EVAL_GRID = EVAL_GRID


def _py_worker(y_row):
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        try:
            f = PY.fit(_G_X_FIT, y_row)
            preds = 100.0 * np.clip(np.atleast_1d(PY.kappa(_G_EVAL_GRID, f)), 0.0, 1.0)
            ok = bool(np.all(np.isfinite(preds)))
        except Exception:
            preds = np.full(len(_G_EVAL_GRID), np.nan)
            ok = False
    return ok, len(ws), preds


def _fit_batch(key, y_batch, pool=None):
    """Fits batch of resampled curves (shape B_batch x len(X_FIT)) and returns (ok_mask, warn_count, preds_pct [B_batch x len(EVAL_GRID)])."""
    if key == "line":
        coef = np.polyfit(np.log2(X_FIT), y_batch.T, 1)
        log2_g = np.log2(EVAL_GRID)[:, None]
        preds = 100.0 * np.clip((coef[0][None, :] * log2_g + coef[1][None, :]).T, 0.0, 1.0)
        return np.all(np.isfinite(preds), axis=1), 0, preds
    if key == "miss":
        log_miss = np.log(np.clip(1.0 - y_batch.T, 1e-12, None))
        coef = np.polyfit(np.log(X_FIT), log_miss, 1)
        log_g = np.log(EVAL_GRID)[:, None]
        preds = 100.0 * np.clip(1.0 - np.exp(coef[0][None, :] * log_g + coef[1][None, :]).T, 0.0, 1.0)
        return np.all(np.isfinite(preds), axis=1), 0, preds
    # key == "pycmb"
    if pool is None:
        res = [_py_worker(row) for row in y_batch]
    else:
        res = pool.map(_py_worker, y_batch, chunksize=max(1, len(y_batch) // (max(1, WORKERS) * 4)))
    ok = np.array([r[0] for r in res], dtype=bool)
    n_warn = int(sum(r[1] for r in res))
    preds = np.stack([r[2] for r in res], axis=0)
    return ok, n_warn, preds


def _sub_boot_se(vals_at_world, step_idx):
    """Estimates 2.5/97.5 percentiles and their Monte Carlo SEs via N_SUB sub-bootstrap draws."""
    b_cur = len(vals_at_world)
    lo, hi = np.percentile(vals_at_world, [2.5, 97.5])
    sub_rng = np.random.default_rng(20260 + step_idx)
    idx = sub_rng.integers(0, b_cur, size=(N_SUB, b_cur))
    qs = np.percentile(vals_at_world[idx], [2.5, 97.5], axis=1)
    se_lo = float(np.std(qs[0], ddof=1))
    se_hi = float(np.std(qs[1], ddof=1))
    return float(lo), float(hi), se_lo, se_hi


def run_adaptive(key, b_max, pool=None):
    """Runs adaptive doubling bootstrap starting at B_START until SE and variance criteria pass."""
    idx_world = N_EXT + NTAB.index(WORLD)
    _, _, pt_preds = _fit_batch(key, R_fit_pt[None, :], pool=None)
    pt_preds = pt_preds[0]

    b_prev = 0
    b_cur = B_START
    preds_all = []
    ok_all = []
    n_warn_total = 0
    trace = []
    step_idx = 0
    converged = False

    while True:
        ok_b, n_w, pr_b = _fit_batch(key, Y1_fit[b_prev:b_cur], pool=pool)
        ok_all.append(ok_b)
        preds_all.append(pr_b)
        n_warn_total += n_w
        b_prev = b_cur

        ok_cat = np.concatenate(ok_all)
        pr_cat = np.concatenate(preds_all, axis=0)
        valid_pr = pr_cat[ok_cat]
        vals_w = valid_pr[:, idx_world]

        lo, hi, se_lo, se_hi = _sub_boot_se(vals_w, step_idx)
        width = hi - lo
        max_se_w = max(se_lo, se_hi) / width if width > 0 else float("inf")
        passed = max_se_w <= SE_OVER_W_THR

        trace.append({
            "B": int(b_cur),
            "lo": lo,
            "hi": hi,
            "width": float(width),
            "se_lo": se_lo,
            "se_hi": se_hi,
            "max_se_over_width": float(max_se_w),
            "passed": bool(passed),
        })
        print(f"  {key:5s} B={b_cur:5d}: 95% CI [{lo:5.2f}, {hi:5.2f}] (L={width:5.2f}), "
              f"SE=({se_lo:.3f}, {se_hi:.3f}), max_SE/L={max_se_w:.4f} -> {'PASS' if passed else 'double'}",
              flush=True)

        if passed:
            converged = True
            break
        if b_cur * 2 > b_max:
            print(f"  WARNING: {key} reached cap B_MAX={b_max} without meeting criterion", flush=True)
            break
        b_cur *= 2
        step_idx += 1

    ok_cat = np.concatenate(ok_all)
    pr_cat = np.concatenate(preds_all, axis=0)
    valid_pr = pr_cat[ok_cat]
    assert int((~ok_cat).sum()) == 0, f"{key}: unexpected failed fits: {int((~ok_cat).sum())}"

    q_curves = np.percentile(valid_pr[:, :N_EXT], [2.5, 50.0, 97.5], axis=0)
    at_n = {}
    for idx_n, n_val in enumerate(NTAB):
        col = valid_pr[:, N_EXT + idx_n]
        lo_n, med_n, hi_n = [float(x) for x in np.percentile(col, [2.5, 50.0, 97.5])]
        lo_s, hi_s, se_lo_n, se_hi_n = _sub_boot_se(col, 9000 + idx_n)
        assert abs(lo_n - lo_s) < 1e-12 and abs(hi_n - hi_s) < 1e-12
        at_n[str(int(n_val))] = {
            "est": float(pt_preds[N_EXT + idx_n]),
            "lo": lo_n,
            "med": med_n,
            "hi": hi_n,
            "width": float(hi_n - lo_n),
            "se_lo": se_lo_n,
            "se_hi": se_hi_n,
        }

    return {
        "B": int(b_cur),
        "converged": bool(converged),
        "n_failed": int((~ok_cat).sum()),
        "n_warnings": int(n_warn_total),
        "point_curve": [float(x) for x in pt_preds[:N_EXT]],
        "lo_curve": [float(x) for x in q_curves[0]],
        "med_curve": [float(x) for x in q_curves[1]],
        "hi_curve": [float(x) for x in q_curves[2]],
        "at_N": at_n,
        "doubling_trace": trace,
    }


def main():
    n_workers = max(1, WORKERS)
    print(f"Running adaptive doubling bootstrap for all 3 models (N >= {FIT_FROM:,}, {n_workers} workers for PYC-MB)...", flush=True)
    results = {}
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=n_workers) as pool:
        results["pycmb"] = run_adaptive("pycmb", B_MAX_PY, pool=pool)
    results["line"] = run_adaptive("line", B_MAX_FAST)
    results["miss"] = run_adaptive("miss", B_MAX_FAST)

    out = {
        "meta": {
            "donors": int(D),
            "M": int(M),
            "world": float(WORLD),
            "fit_from": int(FIT_FROM),
            "fit_points": int(FIT_SEL.sum()),
            "grid_sizes": int(len(NF)),
            "grid_values": 70,
            "pool_resamples": int(B_POOL),
            "sub_bootstrap_draws": int(N_SUB),
            "seed": 1,
            "se_over_width_threshold": float(SE_OVER_W_THR),
            "n_ext_points": int(N_EXT),
            "ntab": [float(x) for x in NTAB],
        },
        "measured": {
            "N": [int(n) for n in Ns],
            "est": [float(100.0 * r) for r in R_Ns],
            "lo": [float(100.0 * x) for x in lo_Ns],
            "hi": [float(100.0 * x) for x in hi_Ns],
        },
        "x_ext": [float(x) for x in X_EXT],
        "models": results,
    }
    os.makedirs(os.path.dirname(os.path.abspath(JSON_OUT)), exist_ok=True)
    with open(JSON_OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"Wrote {JSON_OUT}", flush=True)


if __name__ == "__main__":
    main()
