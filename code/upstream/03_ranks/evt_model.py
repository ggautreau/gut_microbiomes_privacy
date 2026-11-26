#!/usr/bin/env python3
"""Extreme-value model of the similarity tail: a fourth extrapolation of the recovery rate (Supplementary Note 3). Aggregates only.

    python3 -I evt_model.py EVT_SIMILARITIES.npz OUT.json [B] [WORKERS]

EVT_SIMILARITIES.npz  output of evt_similarities.py: for each of the 483 donors, the highest similarities of the runs of other people of the HGM set (`top`, decreasing),
                      their number (`n`), the similarity of the most similar other sample (`s`), the number of runs of other people above it (`k`), the cohort (`study`)
OUT.json              the model for five tails, its variants, its test on thinned HGM sets, the cohorts left out in turn and the curve of the tail drawn in Fig. 1B
B, WORKERS            thinnings of the test (default 400) and processes (default all available cores minus two, at least 1)

In a pool of N runs a donor is recovered when no run of other people has higher similarity to its query than its other sample. For a donor recovered in the HGM set
(M runs) the other sample lies above every run of other people, and the chance that one more run exceeds it is read on the upper tail of the similarities of the
runs of other people to this query. Above a threshold u_i, the (r+1)-th highest of the n_i similarities, the threshold exceedance probability zeta_i = P(X > u_i) has
a conjugate Beta(r + 1/2, n_i - r + 1/2) posterior under Jeffreys' Beta(1/2, 1/2) prior, and the excesses follow a generalised Pareto distribution with scale-invariant
reference prior p(sigma_i | xi_i) proportional to 1 / sigma_i on sigma_i > max(0, -xi_i * max_j y_{i,j}):
    P(X > u_i + y | zeta_i, sigma_i, xi_i) = zeta_i * (1 + xi_i y / sigma_i)_+^(-1/xi_i).
With lambda_i this probability at the similarity s_i of the other sample, a donor recovered in the HGM set stays recovered in a pool of N > M runs with
probability (1 - lambda_i)^(N - M); a donor that is not recovered in the HGM set is not recovered in a larger pool. R(N) is the mean over the donors (R(M) is the
measured rate). Each donor has its own shape, drawn from a distribution common to all the donors (hierarchical pooling), and its sign is free: xi_i ~ Beta(a, b) on
[-0.5, 0.5] (the interval of Martins and Stedinger 2000), with a and b estimated from all the donors by unpenalised exact maximum marginal likelihood (integrating out sigma_i directly over (-inf, inf) by adaptive Gauss-Kronrod quadrature
for each of 80 uniform shape cells); the survival probability of a donor is averaged jointly over the posterior of (zeta_i, sigma_i, xi_i) using direct adaptive
Gauss-Kronrod quadrature over sigma_i and 25-point Gauss-Legendre quadrature over the exact Beta quantile transform of zeta_i.
Intervals: the cohorts are drawn with replacement, then the donors within each drawn cohort; on each resample a and b are estimated again from the donors drawn (each
weighted by the number of times it is drawn), which gives each donor a new posterior of its shape, and the probabilities of the donors drawn are averaged. The posterior
of the scale given each cell of the shape does not depend on a and b, so the survival of each donor given each cell is computed once (stay_cells). Of the 100,000 resamples drawn, the first 400 are used, the number at which the rule of extrapolation_bootstrap.py is met for the tail drawn in Fig. 1B
(evt_resampling.py checks it).
Test on thinned HGM sets, as for the three other models: each run of other people is kept with probability f, the model is fitted on the thinned data alone,
with the same share of the runs in the tail, and predicts the recovery rate measured on the whole HGM set; a tail of fewer than ten similarities is not fitted.
The error of a prediction is split in two: what the model gives to the donors that are not recovered in the whole HGM set, and what it takes from those that are.
Fig. 1B draws the tail selected by lowest thinning RMSE across the tested distances (the 3,000 highest similarities)."""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import json, multiprocessing, sys
import numpy as np
from scipy.integrate import quad_vec
from scipy.optimize import minimize
from scipy.special import betainc, betaincinv, roots_legendre

M = 291_881
WORLD = 8e9
TAILS = (100, 300, 1000, 3000, 10_000)                                            # similarities above the threshold: 0.034 % to 3.4 % of the HGM set
VARIANT_TAILS = (300, 1000, 3000)                                                 # the tails of the two variants
RANGES = dict(free=(-0.5, 0.5), bounded=(-0.5, 0.0), unbounded=(0.0, 1.0))        # the shape: free in sign (the model); a bounded tail only; an unbounded tail only
THIN = (0.1, 0.03, 0.01)                                                          # the HGM set thinned 10, 33 and 100 times
MIN_EXCESSES = 10                                                                 # a thinned tail of fewer similarities is not fitted
DRAWN = 3000                                                                      # fallback tail drawn in Fig. 1B when B == 0; when B > 0, selected by lowest thinning RMSE
G = 80
_GL_X, _GL_W = roots_legendre(25)
GL_U = 0.5 * (_GL_X + 1.0)
GL_W = 0.5 * _GL_W
RESAMPLES = 100_000                                                               # resamples drawn
RESAMPLES_USED = 400                                                              # the first of them, used for the intervals: chosen by adaptive doubling in evt_resampling.py (Supplementary Fig. 6)

if __name__ == "__main__":
    SIM, OUT = sys.argv[1:3]
    B = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    AVAIL_CPUS = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    WORKERS = int(sys.argv[4]) if len(sys.argv) > 4 else max(1, AVAIL_CPUS - 2)
    d = np.load(SIM, allow_pickle=True)
    TOP, n_all, s_all, k_all, study = d["top"].astype(np.float64), d["n"].astype(np.float64), d["s"], d["k"], d["study"]
    D, K = TOP.shape
    alive_all = k_all == 0


def cells(lo, hi):
    """Uniform cells of the shape on [lo, hi]: edges and midpoints."""
    e = np.linspace(lo, hi, G + 1)
    return e, 0.5 * (e[:-1] + e[1:])


def _log_g_block(Yb, yb, x, smin, phi, buf):
    """Exact log-integrand g(phi) for phi = log(sigma - smin) under scale-invariant prior p(sigma) propto 1/sigma."""
    r = Yb.shape[1]
    phi_c = np.clip(phi, -500.0, 500.0)
    ep = np.exp(phi_c)
    sig = smin + ep
    if abs(x) < 1e-12:
        res = -(r + 1) * np.log(sig) + phi_c - Yb.sum(1) / sig
    elif x < 0:
        ax = -x
        axy = np.multiply(ax, Yb, out=buf)
        small = axy < (0.5 * sig)[:, None]
        lp = np.where(small, np.log1p(np.where(small, -axy / sig[:, None], 0.0)), np.log(np.maximum(ax * (yb[:, None] - Yb) + ep[:, None], 1e-300) / sig[:, None]))
        res = -(r + 1) * np.log(sig) + phi_c - (1.0 + 1.0 / x) * lp.sum(1)
    else:
        lp = np.log1p(np.multiply(x / sig[:, None], Yb, out=buf))
        res = -(r + 1) * np.log(sig) + phi_c - (1.0 + 1.0 / x) * lp.sum(1)
    return np.where(np.abs(phi) > 500.0, -np.inf, res)


def _F_H_block(Yb, yb, x, smin, phi, buf1, buf2):
    """Scaled gradient F (so dg/dphi = (ep/sig)*F, strictly decreasing) and negative second derivative dF/dphi = -(ep/sig)*H_in, plus -g''(phi) = H_phi."""
    r = Yb.shape[1]
    ep = np.exp(phi)
    sig = smin + ep
    if abs(x) < 1e-12:
        s = Yb.sum(1) / sig
        return -r + s, s, s
    if x < 0:
        ax = -x
        np.subtract(yb[:, None], Yb, out=buf1); np.multiply(buf1, ax, out=buf1); np.add(buf1, ep[:, None], out=buf1)
        np.multiply(Yb, ax, out=buf2); np.divide(buf2, buf1, out=buf2)
        c_xi = 1.0 / ax - 1.0
        wmax = smin / ep
        F = -r + c_xi * buf2.sum(1) + wmax
        np.add(buf2, 1.0, out=buf1); np.multiply(buf2, buf1, out=buf1)
        H_in = c_xi * buf1.sum(1) + wmax * (1.0 + wmax)
    else:
        np.multiply(Yb, x, out=buf2); np.add(buf2, sig[:, None], out=buf1); np.divide(buf2, buf1, out=buf2)
        c_xi = 1.0 + 1.0 / x
        F = -r + c_xi * buf2.sum(1)
        np.subtract(1.0, buf2, out=buf1); np.multiply(buf2, buf1, out=buf1)
        H_in = c_xi * buf1.sum(1)
    alpha = ep / sig
    return F, alpha * H_in, (alpha ** 2) * H_in


def _integrate_sigma_cell(Y, ymax, c0, x):
    """Exact scale marginalisation over sigma_i > max(0, -xi * max y) for one candidate shape `x`."""
    Dn, r = Y.shape
    ll_g = np.empty(Dn)
    phi_hat_g = np.empty(Dn)
    s_phi_g = np.empty(Dn)
    g0_g = np.empty(Dn)
    logZ_g = np.empty(Dn)
    step = max(1, 32_768 // r)                                                      # donors handled together, so that the arrays of a block stay in the cache of a processor (the results do not depend on it)
    for a in range(0, Dn, step):
        b = min(Dn, a + step)
        Yb = Y[a:b]; yb = ymax[a:b]; c0b = c0[a:b]; nb = b - a
        t1 = np.empty_like(Yb); t2 = np.empty_like(Yb)
        smin = (-x) * yb if x < 0 else np.zeros(nb)
        lo_b = c0b - 30.0
        hi_b = c0b + 30.0
        for _ in range(10):
            Fv = _F_H_block(Yb, yb, x, smin, lo_b, t1, t2)[0]
            need = Fv <= 0
            if not np.any(need):
                break
            lo_b = np.where(need, lo_b - 20.0, lo_b)
        for _ in range(10):
            Fv = _F_H_block(Yb, yb, x, smin, hi_b, t1, t2)[0]
            need = Fv >= 0
            if not np.any(need):
                break
            hi_b = np.where(need, hi_b + 20.0, hi_b)
        phi = np.clip(c0b, lo_b + 1e-6, hi_b - 1e-6)
        for _ in range(35):
            Fv, dF_neg, _ = _F_H_block(Yb, yb, x, smin, phi, t1, t2)
            lo_b = np.where(Fv > 0, phi, lo_b)
            hi_b = np.where(Fv <= 0, phi, hi_b)
            active = (np.abs(Fv) > 1e-9) & ((hi_b - lo_b) > 1e-10)
            if not np.any(active):
                break
            cand = phi + Fv / np.maximum(dF_neg, 1e-30)
            nxt = np.where((cand > lo_b) & (cand < hi_b), cand, 0.5 * (lo_b + hi_b))
            phi = np.where(active, nxt, phi)
        ph = phi
        _, _, H_hat = _F_H_block(Yb, yb, x, smin, ph, t1, t2)
        g0_b = _log_g_block(Yb, yb, x, smin, ph, t1)
        s_phi_b = 1.0 / np.sqrt(np.maximum(H_hat, 1e-30))
        Z_b, _ = quad_vec(
            lambda z_: np.exp(np.clip(_log_g_block(Yb, yb, x, smin, ph + s_phi_b * z_, t1) - g0_b, -700.0, 0.0)) * s_phi_b,
            -np.inf, np.inf, epsrel=1e-4, quadrature="gk15",
        )
        logZ_b = np.log(Z_b)
        ll_g[a:b] = g0_b + logZ_b
        phi_hat_g[a:b] = ph; s_phi_g[a:b] = s_phi_b; g0_g[a:b] = g0_b; logZ_g[a:b] = logZ_b
    return ll_g, phi_hat_g, s_phi_g, g0_g, logZ_g


_SIGMA_CTX = None


def _integrate_sigma_cell_worker(g):
    Y, xi, ymax, c0 = _SIGMA_CTX
    return _integrate_sigma_cell(Y, ymax, c0, float(xi[g]))


def integrate_sigma(Y, xi, workers=1):
    """Exact marginalisation over scale sigma_i > max(0, -xi * max y) under reference prior p(sigma_i) propto 1 / sigma_i."""
    global _SIGMA_CTX
    Gx = len(xi)
    ymax = Y.max(1)
    c0 = np.log(np.clip(Y.mean(1), 1e-12, None))
    if workers > 1:
        _SIGMA_CTX = (Y, xi, ymax, c0)
        with multiprocessing.get_context("fork").Pool(max(1, min(workers, Gx))) as cpool:
            parts = cpool.map(_integrate_sigma_cell_worker, range(Gx), chunksize=1)
        _SIGMA_CTX = None
    else:
        parts = [_integrate_sigma_cell(Y, ymax, c0, float(x)) for x in xi]
    ll, phi_hat, s_phi, g0, logZ = (np.column_stack([p[k] for p in parts]) for k in range(5))
    return dict(ll=ll, phi_hat=phi_hat, s_phi=s_phi, g0=g0, logZ=logZ, ymax=ymax, Y=Y)


def fit(top, n, s, r, lo, hi, prof=None, workers=1):
    """Tails of the donors. top: donors x (at least r + 1) highest similarities of the runs of other people, decreasing; n: their number; s: similarity of the other sample.
    prof: precomputed scale-marginalised representation from integrate_sigma, when already computed (those of a donor do not depend on the other donors).
    Returns the posterior over shape cells, exact Beta quantile nodes for zeta_i, and hyperparameters (a, b)."""
    u = top[:, r]; e, xi = cells(lo, hi)
    integ = integrate_sigma(np.clip(top[:, :r] - u[:, None], 1e-12, None), xi, workers=workers) if prof is None else prof
    ll = integ["ll"] - integ["ll"].max(1, keepdims=True)
    mass = lambda a, b: np.clip(np.diff(betainc(a, b, np.clip((e - lo) / (hi - lo), 0, 1))), 1e-300, None)      # mass of the Beta(a, b) on [lo, hi] in each cell

    def nlm(theta):                                                                  # minus the exact unpenalised marginal log-likelihood on log a, log b
        return -np.log((np.exp(ll) * mass(*np.exp(theta))[None, :]).sum(1)).sum()

    bounds = ((-12.0, 12.0), (-12.0, 12.0))
    best = min((minimize(nlm, x0, method="Nelder-Mead", bounds=bounds, options=dict(xatol=1e-4, fatol=1e-6, maxiter=2000)) for x0 in ((0.0, 0.0), (-1.0, 2.0), (1.0, 3.0), (3.0, 3.0))), key=lambda o: o.fun)
    a, b = np.exp(best.x)
    post = np.exp(ll) * mass(a, b)[None, :]; post /= post.sum(1, keepdims=True)
    az = np.full_like(n, r + 0.5, dtype=float)
    bz = np.maximum(n - r, 0.0) + 0.5
    zeta_nodes = betaincinv(az[:, None], bz[:, None], GL_U[None, :])
    return dict(post=post, xi=xi, a=float(a), b=float(b), xi_mean=(post * xi[None, :]).sum(1), u=u, s=s,
                zeta_nodes=zeta_nodes, integ=integ, L=np.exp(ll), e=e, lo=lo, hi=hi)


def _stay_cell(f, idx, extras_2d, g):
    """Contribution of shape cell g to survival probabilities for donors `idx` and extra pool sizes `extras_2d` (shape: len(idx) x M_ex)."""
    integ = f["integ"]
    x = float(f["xi"][g])
    Y = integ["Y"][idx]; ymax = integ["ymax"][idx]
    ph = integ["phi_hat"][idx, g]; sp = integ["s_phi"][idx, g]; g0 = integ["g0"][idx, g]; logZ = integ["logZ"][idx, g]
    post_g = f["post"][idx, g]; zn = f["zeta_nodes"][idx]
    m_ex = f["s"][idx] - f["u"][idx]; pos_m = m_ex > 0.0
    n_idx, m_cols = extras_2d.shape
    r = Y.shape[1]
    res = np.empty((n_idx, m_cols))
    step = max(1, 32_768 // r)
    for a in range(0, n_idx, step):
        b = min(n_idx, a + step)
        Yb = Y[a:b]; yb = ymax[a:b]; ph_b = ph[a:b]; sp_b = sp[a:b]; g0_b = g0[a:b]; logZ_b = logZ[a:b]
        post_b = post_g[a:b]; zn_b = zn[a:b]; m_b = m_ex[a:b]; pos_b = pos_m[a:b]; ex_b = extras_2d[a:b]
        smin_b = (-x) * yb if x < 0 else np.zeros(b - a)
        buf = np.empty_like(Yb)

        def integrand(z_):
            phi = ph_b + sp_b * z_
            w = np.exp(np.clip(_log_g_block(Yb, yb, x, smin_b, phi, buf) - g0_b, -700.0, 0.0) - logZ_b) * sp_b
            sig = smin_b + np.exp(np.clip(phi, -500.0, 500.0))
            if abs(x) < 1e-12:
                q = np.where(pos_b, np.exp(-m_b / sig), 1.0)
            else:
                base = 1.0 + x * m_b / sig
                q = np.where(pos_b, np.where(base > 0.0, np.power(np.clip(base, 1e-300, None), -1.0 / x), 0.0), 1.0)
            lam = np.clip(q[:, None] * zn_b, 0.0, 1.0 - 1e-16)
            eb = np.sum(GL_W[None, None, :] * np.exp(ex_b[:, :, None] * np.log1p(-lam[:, None, :])), axis=-1)
            return (post_b * w)[:, None] * eb

        val, _ = quad_vec(integrand, -np.inf, np.inf, epsrel=1e-3, epsabs=1e-5, quadrature="gk15")
        res[a:b] = val
    return res


_STAY_CTX = None


def _stay_cell_worker(g):
    f, idx, extras_2d = _STAY_CTX
    return _stay_cell(f, idx, extras_2d, g)


def stay(f, alive, extra, workers=1):
    """Probability that each donor is first once `extra` runs are added to the pool where `alive` says whether it is first."""
    global _STAY_CTX
    Dn = len(alive)
    ex_all = np.broadcast_to(np.asarray(extra, float), (Dn,))
    out = np.where(alive, 1.0, 0.0)
    idx = np.where(alive & (ex_all > 0.0))[0]
    if not len(idx):
        return out
    extras_2d = ex_all[idx, None]
    n_g = len(f["xi"])
    total = np.zeros((len(idx), 1))
    if workers > 1:
        _STAY_CTX = (f, idx, extras_2d)
        with multiprocessing.get_context("fork").Pool(max(1, min(workers, n_g))) as cpool:
            parts = cpool.map(_stay_cell_worker, range(n_g), chunksize=1)
        _STAY_CTX = None
        for part in parts:
            total += part
    else:
        for g in range(n_g):
            total += _stay_cell(f, idx, extras_2d, g)
    out[idx] = np.clip(total[:, 0], 0.0, 1.0)
    return out


def stay_many(f, alive, extras, workers=1):
    """Survival probabilities (donors x len(extras)) across a 1D grid of extra pool sizes `extras`, evaluating the scale log-likelihood once per adaptive node across all pool sizes."""
    global _STAY_CTX
    ex_grid = np.asarray(extras, float)
    Dn = len(alive)
    out = np.where(alive[:, None], 1.0, np.zeros((Dn, len(ex_grid))))
    pos_cols = np.where(ex_grid > 0.0)[0]
    idx = np.where(alive)[0]
    if not len(idx) or not len(pos_cols):
        return out
    extras_2d = np.broadcast_to(ex_grid[None, pos_cols], (len(idx), len(pos_cols)))
    n_g = len(f["xi"])
    if workers > 1:
        _STAY_CTX = (f, idx, extras_2d)
        with multiprocessing.get_context("fork").Pool(max(1, min(workers, n_g))) as cpool:
            parts = cpool.map(_stay_cell_worker, range(n_g), chunksize=1)
        _STAY_CTX = None
    else:
        parts = [_stay_cell(f, idx, extras_2d, g) for g in range(n_g)]
    total = np.sum(parts, axis=0)
    out[np.ix_(idx, pos_cols)] = np.clip(total, 0.0, 1.0)
    return out


def stay_cells(f, alive, extras, workers=1):
    """Probability that each donor is first once `extras[j]` runs are added to the pool, given each cell of the shape (donors x len(extras) x cells): the terms that
    stay() and stay_many() add up, before they are weighted by the posterior of the shape. They do not depend on a and b."""
    global _STAY_CTX
    ex_grid = np.atleast_1d(np.asarray(extras, float))
    Dn, n_g = len(alive), len(f["xi"])
    out = np.repeat(np.where(alive, 1.0, 0.0)[:, None, None], len(ex_grid), 1).repeat(n_g, 2)
    pos_cols = np.where(ex_grid > 0.0)[0]
    idx = np.where(alive)[0]
    if not len(idx) or not len(pos_cols):
        return out
    f1 = dict(f, post=np.ones_like(f["post"]))                                       # the posterior weight of each cell is applied later, once per resample
    extras_2d = np.broadcast_to(ex_grid[None, pos_cols], (len(idx), len(pos_cols)))
    if workers > 1:
        _STAY_CTX = (f1, idx, extras_2d)
        with multiprocessing.get_context("fork").Pool(max(1, min(workers, n_g))) as cpool:
            parts = cpool.map(_stay_cell_worker, range(n_g), chunksize=1)
        _STAY_CTX = None
    else:
        parts = [_stay_cell(f1, idx, extras_2d, g) for g in range(n_g)]
    out[np.ix_(idx, pos_cols, np.arange(n_g))] = np.clip(np.stack(parts, axis=-1), 0.0, 1.0)
    return out


def refit_shapes(f, w):
    """The Beta(a, b) of the shapes estimated again with donor weights w (the number of times a resample draws each donor, up to a factor), from the estimate on all the
    donors; returns the posterior of the cells of the shape of every donor (donors x cells) and (a, b)."""
    e, lo, hi, L = f["e"], f["lo"], f["hi"], f["L"]
    mass = lambda a, b: np.clip(np.diff(betainc(a, b, np.clip((e - lo) / (hi - lo), 0, 1))), 1e-300, None)
    keep = w > 0; Lk = L[keep]; wk = w[keep]

    def nlm(theta):
        return -(wk * np.log((Lk * mass(*np.exp(theta))[None, :]).sum(1))).sum()

    o = minimize(nlm, np.log([f["a"], f["b"]]), method="Nelder-Mead", bounds=((-12.0, 12.0), (-12.0, 12.0)), options=dict(xatol=1e-4, fatol=1e-6, maxiter=2000))
    a, b = np.exp(o.x)
    post = L * mass(a, b)[None, :]; post /= post.sum(1, keepdims=True)
    return post, float(a), float(b)


_RS_CTX = None


def _resample_worker(j):
    f, S, W = _RS_CTX
    post, a, b = refit_shapes(f, len(W[j]) * W[j])                                  # weights summing to the number of donors, as on the whole set
    return np.einsum("i,ig,ijg->j", W[j], post, S), a, b


def resample_values(f, S, W, workers=1):
    """Recovery rates (%) of the resamples W (rows: weights of the donors, summing to 1) at the pool sizes of S (output of stay_cells), with a and b estimated again on
    each resample; also returns the (a, b) of each resample."""
    global _RS_CTX
    _RS_CTX = (f, S, W)
    if workers > 1:
        with multiprocessing.get_context("fork").Pool(max(1, min(workers, len(W)))) as pool:
            res = pool.map(_resample_worker, range(len(W)), chunksize=max(1, len(W) // (4 * workers)))
    else:
        res = [_resample_worker(j) for j in range(len(W))]
    _RS_CTX = None
    return 100 * np.array([r_[0] for r_ in res]), np.array([[r_[1], r_[2]] for r_ in res])


def whole(arg):
    """The model on the whole HGM set, for one range of the shape and one tail."""
    name, r = arg
    return fit(TOP, n_all, s_all, r, *RANGES[name])


def thin_one(arg):
    """One thinned HGM set: each run of other people kept with probability fthin; the model fitted with the same share of the runs in the tail predicts the rate at M."""
    fthin, b = arg; g = np.random.default_rng(1000 * b + int(1000 * fthin))
    keep = g.random((D, K), dtype=np.float32) < fthin
    in_top = np.minimum(k_all, K)                                                     # runs of other people above the other sample that are among the highest K
    kt = np.array([keep[i, :in_top[i]].sum() for i in range(D)]) + g.binomial(np.maximum(k_all - K, 0), fthin)
    n_t = keep.sum(1) + g.binomial((n_all - K).astype(np.int64), fthin)
    order = np.argsort(~keep, axis=1, kind="stable")                                # the kept similarities first, in decreasing order
    res = {"no_fit": 100 * float((kt == 0).mean())}                                 # the share measured in the thinned set, taken as the prediction
    for r_full in TAILS:
        r = int(round(r_full * fthin))
        if r < MIN_EXCESSES:
            continue
        assert keep.sum(1).min() > r, "the thinned tail lies inside the similarities stored"
        ft = fit(np.take_along_axis(TOP, order[:, :r + 1], 1), n_t.astype(float), s_all, r, *RANGES["free"])
        p = stay(ft, kt == 0, M - n_t)
        res[r_full] = (100 * float(p.mean()), 100 * float(p[~alive_all].sum()) / D, 100 * float((p[alive_all] - 1).sum()) / D)      # the prediction; what it gives to the donors not recovered in the whole HGM set; what it takes from the donors recovered
    return res


if __name__ == "__main__":
    cohorts = np.unique(study); members = {c: np.where(study == c)[0] for c in cohorts}
    rng = np.random.default_rng(3)
    m_c = rng.multinomial(len(cohorts), [1 / len(cohorts)] * len(cohorts), size=RESAMPLES)      # the cohorts drawn with replacement, then the donors within each drawn cohort (as for the other models)
    W = np.zeros((RESAMPLES, D))
    for j, c in enumerate(cohorts):
        W[:, members[c]] = rng.multinomial(m_c[:, j] * len(members[c]), [1 / len(members[c])] * len(members[c]))
    W /= W.sum(1, keepdims=True)
    interval = lambda v: [float(x) for x in np.percentile(v, [2.5, 97.5])]
    N_GRID = np.logspace(np.log10(M), 12, 160)
    measured = 100 * float(alive_all.mean())

    out: dict = dict(donors=D, recovered=int(alive_all.sum()), hgm_runs=M, similarities_stored=K, tails=list(TAILS), tail_share={str(r): 100 * r / float(np.median(n_all)) for r in TAILS},
                     ranges={k_: list(v) for k_, v in RANGES.items()}, resamples=RESAMPLES_USED, resamples_drawn=RESAMPLES, min_excesses=MIN_EXCESSES, N=[float(v) for v in N_GRID], models={k_: {} for k_ in RANGES})
    tasks = [("free", r) for r in TAILS] + [(name, r) for name in ("bounded", "unbounded") for r in VARIANT_TAILS]
    ctx = multiprocessing.get_context("fork")
    fits = {(name, r): fit(TOP, n_all, s_all, r, *RANGES[name], workers=WORKERS) for name, r in tasks}
    for (name, r), f in fits.items():
        lo, hi = RANGES[name]; p_w = stay(f, alive_all, WORLD - M, workers=WORKERS); p_7 = stay(f, alive_all, 1e7 - M, workers=WORKERS)
        S_7, S_w = stay_cells(f, alive_all, [1e7 - M], workers=WORKERS), stay_cells(f, alive_all, [WORLD - M], workers=WORKERS)          # one pool size per call, as in evt_resampling.py
        v_rs, ab_rs = resample_values(f, np.concatenate([S_7, S_w], 1), W[:RESAMPLES_USED], workers=WORKERS)                             # a and b estimated again on each resample
        sm_rs = lo + (hi - lo) * ab_rs[:, 0] / ab_rs.sum(1)
        res: dict = dict(a=f["a"], b=f["b"], shape_mean=float(lo + (hi - lo) * f["a"] / (f["a"] + f["b"])), shape_quartiles=[float(q) for q in np.percentile(f["xi_mean"][alive_all], [25, 50, 75])],
                         shape_negative=float(np.mean(f["xi_mean"][alive_all] < 0)), at_1e7=100 * float(p_7.mean()), interval_1e7=interval(v_rs[:, 0]), at_world=100 * float(p_w.mean()), interval_world=interval(v_rs[:, 1]),
                         shape_mean_interval=interval(sm_rs), at_world_cells=100 * float(np.einsum("ig,ig->i", f["post"], S_w[:, 0]).mean()),
                         expected_recovered_at_M=float(stay(f, alive_all, M - 1, workers=WORKERS).sum()))                # of the donors recovered, those that stay first when the runs of the HGM set are drawn again from the fitted tails
        if name == "free":
            res["curve"] = (100 * stay_many(f, alive_all, N_GRID - M, workers=WORKERS).mean(0)).tolist()
            res["leave_one_out"] = {}
            for c in cohorts:                                                        # the model made again without each cohort in turn: the Beta is estimated on the donors kept
                keep = study != c; fc = fit(TOP[keep], n_all[keep], s_all[keep], r, lo, hi, prof={k_: v_[keep] for k_, v_ in f["integ"].items()})
                res["leave_one_out"][str(c)] = 100 * float(stay(fc, alive_all[keep], WORLD - M, workers=WORKERS).mean())
        out["models"][name][str(r)] = res
        print(f"{name:10s} tail {r:>6,d}: mean shape {res['shape_mean']:+.3f}; at 1e7 {res['at_1e7']:.1f} % ({res['interval_1e7'][0]:.0f}-{res['interval_1e7'][1]:.0f}); at 8e9 {res['at_world']:.1f} % "
              f"({res['interval_world'][0]:.0f}-{res['interval_world'][1]:.0f}); recovered donors that stay first when the HGM set is drawn again {res['expected_recovered_at_M']:.1f} of {int(alive_all.sum())}", flush=True)
    json.dump(out, open(OUT, "w"), indent=1)
    if B:
        out["thinning"] = {}; out["thinnings"] = B
        with ctx.Pool(max(1, min(WORKERS, B))) as pool:
            for fthin in THIN:
                runs = pool.map(thin_one, [(fthin, b) for b in range(B)], chunksize=1); o: dict = {}; out["thinning"][str(fthin)] = o
                nf = np.array([x["no_fit"] for x in runs]) - measured
                o["no_fit"] = dict(error=float(nf.mean()), sd=float(nf.std()), rmse=float(np.sqrt(np.mean(nf ** 2))), thinnings=len(nf))
                for r_full in TAILS:
                    v = np.array([x[r_full] for x in runs if r_full in x])
                    if len(v):
                        e = v[:, 0] - measured
                        o[str(r_full)] = dict(error=float(e.mean()), sd=float(e.std()), rmse=float(np.sqrt(np.mean(e ** 2))), given_to_not_recovered=float(v[:, 1].mean()), taken_from_recovered=float(v[:, 2].mean()),
                                              excesses=int(round(r_full * fthin)), thinnings=len(v))
                        print(f"thinned {1 / fthin:4.0f} times, tail {r_full:>6,d} ({int(round(r_full * fthin)):>5d} excesses): mean error {o[str(r_full)]['error']:+.2f} ({o[str(r_full)]['given_to_not_recovered']:+.2f} and {o[str(r_full)]['taken_from_recovered']:+.2f}), "
                              f"RMSE {o[str(r_full)]['rmse']:.2f}; no fit {o['no_fit']['error']:+.2f}, RMSE {o['no_fit']['rmse']:.2f}", flush=True)
                json.dump(out, open(OUT, "w"), indent=1)
    if B:
        best_by_f = {f: min((r_ for r_ in TAILS if str(r_) in out["thinning"][str(f)]), key=lambda r_: out["thinning"][str(f)][str(r_)]["rmse"]) for f in THIN}
        common = [r_ for r_ in TAILS if all(str(r_) in out["thinning"][str(f)] for f in THIN)]
        drawn = next(iter(best_by_f.values())) if len(set(best_by_f.values())) == 1 else min(common, key=lambda r_: sum(out["thinning"][str(f)][str(r_)]["rmse"] for f in THIN))
    else:
        drawn = DRAWN
    out["drawn"] = drawn
    out["curve"] = dict(N=out["N"], rate=out["models"]["free"][str(drawn)]["curve"])
    print(f"tail drawn: the {drawn:,} highest similarities ({out['tail_share'][str(drawn)]:.2f} % of the runs of other people), {out['models']['free'][str(drawn)]['at_world']:.1f} % at 8e9", flush=True)
    json.dump(out, open(OUT, "w"), indent=1)
    print("wrote", OUT)