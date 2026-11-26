"""pyc_mb.py -- Pitman-Yor correctness law (PYC) and its measurement-based fit (PYC-MB).  Validated module.

Reference: Rocher L., Hendrickx J.M., de Montjoye Y.-A. (2025), "A scaling law to model the effectiveness of
identification techniques", Nature Communications 16:347, doi:10.1038/s41467-024-55296-6 (open access, PMC11718298).
SI = 41467_2024_55296_MOESM1_ESM.pdf.  Below, "the paper" is this reference.  Numpy + scipy only (Python >= 3.9, scipy >= 1.7).

Public API:  kappa, fit, predict
    kappa(n, params, gamma=None, unit="nats")  expected fraction correctly identified in a closed-set gallery of n records
    fit(ns, kappas, weights=None, ...)         PYC-MB fit of (h, gamma) to t >= 2 measured points (n_i, kappa_i) -> PYCFit
    predict(params, n_target, unit="nats")     kappa at n' >> n with fitted (or given) parameters
  also: kappa_da (the law in (d, alpha)), hg_to_da / da_to_hg, loss, fit_all, fit_like_authors, fit_regimes,
        forecast_band (sensitivity of a forecast; NOT in the paper), PYCParams, PYCFit, bits_to_nats, psi_inv.

EQUATIONS (read in this work from the typeset PDF pages and the JATS/MathML of PMC11718298; numbers are the paper's)
  Eq. (1)  kappa(n) = sum_i pi_i [1 - (1 - pi_i)^n] / (n pi_i) = (1/n) sum_i [1 - (1 - pi_i)^n] = E[#distinct anonymity sets
           among n records]/n, with pi ~ PY(d, alpha) (Pitman-Yor), d in [0,1], alpha in [-d, +inf] (paper's text).
  Eq. (2)  h = E[H(pi)|d,alpha] = psi0(alpha+1) - psi0(1-d)  (NATS, digamma form);
           gamma = (psi0(1) - psi0(1-d)) / (psi0(alpha+1) - psi0(1-d)).  Inverse: u = 1-d = psi0^-1(psi0(1) - h gamma),
           v = 1+alpha = psi0^-1(psi0(1) + h - h gamma).
  Eq. (3)  AS PRINTED (page 3):  (v-1)/(n(u-1)) [ G(n-u+v) G(v-1) / (G(n-1+v) G(v-u)) - 1 ].  This is MINUS the correct value
           (kappa(1) = -1; |Eq.(3)_printed + Eq.(13)| < 1e-95 in 120-digit arithmetic): the denominator must read n(1-u).
           It is a typo of the printed formula; this module uses Eq. (13).
  Eq. (13) (Methods, p. 9; = SI Eq. S31)  E[kappa|d,alpha] = (1/(n d)) ( G(1+alpha)/G(d+alpha) * G(n+d+alpha)/G(n+alpha) - alpha ),
           d = 0 by its limit (alpha/n)(psi0(n+alpha) - psi0(alpha)).  Equivalently kappa(n) = E[K_n]/n with K_n the number of
           occupied tables of the (d, alpha) Chinese restaurant process, E[K_{m+1}] = E[K_m](1 + d/(alpha+m)) + alpha/(alpha+m)
           (checked to 1e-60 for d > 0, d = 0, d < 0, alpha < 0).  kappa(1) = 1 for every (h, gamma).
  Eq. (4)  (PYC-MB) R(h, gamma) = sum_{i=1}^t log(n_i) [kappa_hat(n_i) - E[kappa(n_i)|h,gamma]]^2, t >= 2 points n_i < n',
           minimised "using a derivative-free method (Nelder-Mead) and a logarithmic weighting"; the forecast is
           E[kappa(n')|h_hat, gamma_hat].  log = natural log (any base gives the same minimiser).  The n = 1 point has weight 0.
  Also read: SI S26 (stick-breaking beta_i ~ Beta(1-d, alpha+i d)), S28 (entropy), S31 (= Eq. 13, proof by the size-biased
  sample of S27), S2.4/S3.1 (uniform law), S3.2-S3.3 (50 log-spaced points, evaluation only if 0.01 < kappa_hat(n(t)) < 0.99).
  Not in the paper: start point, bounds, tolerances, restarts, units of h in the fit.  The authors' code (github.com/
  synthetic-society/dataless; model.py and extrapolate.py of commit 857a6b8d4a, 2024-11-14, read in full) uses start
  (h, gamma) = (12, 0.26), scipy's default Nelder-Mead options, loss mean(log(n)(kappa-model)^2) and NO bounds; the January 2026
  release (commit f2e4c645c3) restricts to h > 0, 0 <= gamma <= 1.  Units: h is in NATS here (the authors' code docstring says
  nats; "h = 41.54, gamma = 0.68 -> 75 % at 4e9" of the Results is reproduced in nats: 0.7438); the Discussion's "h = 40 bits:
  99 % (gamma = 0), 45 % (gamma = 1) at n = 7.53e9" is reproduced with h_nats = 40 ln 2 (0.9940, 0.4476).  The Introduction's
  worked example kappa(100) = 0.99, kappa(1k) = 0.80 -> "kappa ~ 0.21 for n = 10k" is reproduced (0.2145) only with gamma < 0
  (gamma in [0, 1] gives 0.346).  The authors' README outputs are reproduced to 5e-9 by fit_like_authors.

REGIMES of the discount d = 1 - u (all covered by the one analytic formula Eq. (13); h > 0 is exactly alpha + d > 0):
  gamma > 1        alpha < 0 (-d < alpha < 0, 0 < d < 1): heavy tail, a proper Pitman-Yor process (inside the paper's text domain)
  0 < gamma < 1    0 < d < 1, alpha > 0: heavy (power-law) tail, kappa ~ n^(d-1) for n >> alpha
  gamma = 0        d = 0: Dirichlet process, geometric tail (exact digamma limit)
  gamma < 0        d < 0: finite-support (Dirichlet-multinomial with m = alpha/|d| classes when m is an integer, an analytic
                   continuation otherwise; kappa ~ m/n beyond a knee).  The paper itself uses gamma = -1 ("finite uniform
                   distributions", its Fig. 4) although its text says d >= 0.  gamma -> -inf is the uniform law on e^h classes
                   (the entropy baseline ENT): kappa = (m/n)(1 - (1 - 1/m)^n), m = e^h; used exactly (error < 1e-15) when
                   psi0(1) - h gamma > 36.
  The map (h > 0, gamma real) -> (d < 1, alpha > -d) is onto the whole admissible Pitman-Yor domain.

FIT PROCEDURE (fit): the objective is exactly Eq. (4) (sum instead of the authors' mean: same minimiser; weights=None means
  log(n_i), an explicit weights array REPLACES the log weights).  The surface can be flat and the paper does not say how its
  single Nelder-Mead was stopped, so the global minimum is searched in (ln h, gamma): (i) profile scan: for each gamma of a
  43-node grid the loss is minimised over ln h (32-node scan + bounded Brent), the best separated minima of that profile (<= 4)
  are the starts; (ii) from each start and from the authors' start (12, 0.26): bounded trust-region least squares (accelerator
  only, same objective) then scipy Nelder-Mead (xatol 1e-8, fatol 1e-12 x loss scale) with up to 3 restarts from the incumbent;
  (iii) the lowest loss wins; a loss below 1e-16 (data reproduced to 1e-8) stops the search.  Default box h in [1e-6, 700] nats,
  gamma in [-50, 50] (the paper states none; the authors' 2024 code is unconstrained).  fit() warns when the optimum lies on a
  bound, when kappa_hat at the largest n is outside (0.01, 0.99) (SI S3.3 evaluates PYC-MB only there), and when exactly two
  points are given (exactly determined, possibly non-unique).  PYCFit.n_optima > 1 means a flat valley (non-identifiable).

KNOWN FAILURE CASES: (a) plateau curves (kappa_hat ~ 1 over the whole range) and curves much flatter
  than any Pitman-Yor law: the loss valley runs to h -> infinity (or gamma -> +/-inf) and the forecast is set by the box, not the
  data (flagged by bounds_hit); (b) t = 2 points: exactly determined, flat valley towards gamma -> -inf, the two inputs of the
  worked example are not attainable by ANY law (fitted 0.9775/0.8011); (c) the law decays to 0 as n -> inf for every (h, gamma)
  (no plateau); (d) forecasts are as uncertain as the data allow: the paper's own Table 1 RMSE for robust matching is 0.134
  (mu = 5 %) and 0.183 (mu = 1 %), and its theory covers exact matching only; (e) the paper evaluates 50 log-spaced points,
  few or irregular points change the weighting.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Optional, Sequence, Tuple, Union

import numpy as np
from scipy import optimize
from scipy.special import digamma, gammaln, zeta

__version__ = "validated-1.0"

__all__ = [
    "PYCParams", "PYCFit", "kappa", "kappa_da", "kappa_dalpha", "predict", "fit", "fit_all", "fit_like_authors",
    "fit_regimes", "forecast_band", "loss", "hg_to_da", "da_to_hg", "hg_to_dalpha", "dalpha_to_hg", "psi_inv",
    "bits_to_nats", "nats_to_bits", "REGIMES", "DEFAULT_START", "DEFAULT_H_BOUNDS", "DEFAULT_GAMMA_BOUNDS",
]

ArrayLike = Union[float, int, Sequence[float], np.ndarray]

LN2 = math.log(2.0)
PSI1 = float(digamma(1.0))                 # psi0(1) = -Euler's constant
DEFAULT_START = (12.0, 0.26)               # (h nats, gamma): INIT_STATE of the authors' code (not in the paper)
DEFAULT_H_BOUNDS = (1e-6, 700.0)           # nats; 700 = float64 range of psi0^-1 (exp(y) overflows beyond 709)
DEFAULT_GAMMA_BOUNDS = (-50.0, 50.0)
_Y_MAX = 700.0                             # largest argument of psi0^-1 evaluated in float64
_Y_UNIFORM = 36.0                          # psi0(1) - h gamma above this: uniform-law limit, error < 1e-15


def bits_to_nats(h_bits: ArrayLike):
    """Entropy in bits -> nats (the unit of the digamma formulas)."""
    return np.asarray(h_bits, dtype=float) * LN2


def nats_to_bits(h_nats: ArrayLike):
    return np.asarray(h_nats, dtype=float) / LN2


def _to_nats(h, unit):
    if unit == "nats":
        return h
    if unit == "bits":
        return h * LN2
    raise ValueError("unit must be 'nats' or 'bits'")


# ----------------------------------------------------------------------------------------------------------
# Asymptotic series for second differences of ln Gamma (DLMF 5.11.8):
#   lnG(z+a) - lnG(z) ~ a ln z + sum_{k>=2} c_k(a) z^(1-k),   c_k(a) = (-1)^k [B_k(a) - B_k(0)] / (k (k-1)),
#   B_k(a) - B_k(0) = sum_{j=1}^k C(k,j) B_{k-j} a^j   (no constant term, accurate for tiny a; B_1 = -1/2).
# ----------------------------------------------------------------------------------------------------------
_KMAX = 18
_Z0 = 20.0            # series used for z >= _Z0 and |a| <= z/8 (truncation error < 1e-16 relative)


def _bernoulli_numbers(m: int):
    B = [Fraction(0)] * (m + 1)
    B[0] = Fraction(1)
    for k in range(1, m + 1):
        B[k] = -sum(math.comb(k + 1, j) * B[j] for j in range(k)) / (k + 1)
    return B


def _series_table(kmax: int) -> np.ndarray:
    B = _bernoulli_numbers(kmax)
    tab = np.zeros((kmax + 1, kmax + 1))
    for k in range(2, kmax + 1):
        s = Fraction((-1) ** k, k * (k - 1))
        for j in range(1, k + 1):
            tab[k, j] = float(s * math.comb(k, j) * B[k - j])
    return tab


_TAB = _series_table(_KMAX)
_BERN = [float(b) for b in _bernoulli_numbers(24)]
_KK, _JJ = np.meshgrid(np.arange(_KMAX + 1), np.arange(_KMAX + 1), indexing="ij")
_MI = np.clip(_KK - _JJ, 0, None)


def _series_terms(a, z) -> np.ndarray:
    """T[k] = c_k(a) z^(1-k), k = 0.._KMAX (rows 0 and 1 are zero); a and z broadcast; shape (_KMAX+1, N)."""
    a, z = np.broadcast_arrays(np.atleast_1d(np.asarray(a, float)), np.atleast_1d(np.asarray(z, float)))
    eps = a / z
    iz = 1.0 / z
    ones = np.ones((1, z.size))
    pe = np.cumprod(np.concatenate([ones, np.tile(eps.ravel(), (_KMAX, 1))]), axis=0)
    pi = np.cumprod(np.concatenate([ones, np.tile(iz.ravel(), (_KMAX, 1))]), axis=0)
    prod = _TAB[:, :, None] * pe[None, :, :] * pi[_MI]
    return z.ravel() * prod.sum(axis=1)


def _lgd(z: np.ndarray, a: float) -> np.ndarray:
    """lnGamma(z + a) - lnGamma(z) for z > 0, z + a > 0 (z array, a scalar), accurate in every regime."""
    z = np.asarray(z, dtype=float)
    if a == 0.0:
        return np.zeros_like(z)
    out = np.empty_like(z)
    big = (z >= _Z0) & (abs(a) <= z / 8.0)
    if np.any(big):
        zb = z[big]
        out[big] = a * np.log(zb) + _series_terms(a, zb)[2:].sum(axis=0)
    rest = ~big
    if np.any(rest):
        zr = z[rest]
        if abs(a) > 2.0:
            out[rest] = gammaln(zr + a) - gammaln(zr)
        else:
            # shift z up to >= _Z0:  lnG(z+a)-lnG(z) = [lnG(z+N+a)-lnG(z+N)] - sum_{m<N} log1p(a/(z+m))
            nsh = np.maximum(0.0, np.ceil(_Z0 - zr))
            nmax = int(nsh.max()) if nsh.size else 0
            if nmax > 0:
                mm = np.arange(nmax)[:, None]
                acc = np.where(mm < nsh[None, :], np.log1p(a / (zr[None, :] + mm)), 0.0).sum(axis=0)
            else:
                acc = 0.0
            out[rest] = _lgd(zr + nsh, a) - acc
    return out


def _logratio(alpha: float, d: float, n: np.ndarray) -> np.ndarray:
    """S = ln[ (alpha+d)_n / (alpha)_n ] = lnG(alpha+d+n) - lnG(alpha+n) - lnG(alpha+d) + lnG(alpha)  (alpha > 0, d != 0,
    alpha + d > 0, n >= 0 array) without catastrophic cancellation:
      (i)  alpha >= Z0 and |d| <= alpha/8: series in d/alpha, any n;
      (ii) n <= min(alpha, alpha+d)/8 with min >= Z0: series in n/alpha (second difference symmetric in d <-> n);
      (iii) otherwise differences of lnG (no cancellation problem when S is not small)."""
    if alpha >= _Z0 and abs(d) <= alpha / 8.0:
        L = np.log1p(n / alpha)
        T = _series_terms(d, alpha)[2:, 0]
        ks = np.arange(2, _KMAX + 1)[:, None]
        return d * L + (T[:, None] * np.expm1(-(ks - 1) * L[None, :])).sum(axis=0)
    zmin = min(alpha, alpha + d)
    out = np.empty_like(n)
    m2 = (zmin >= _Z0) & (n <= zmin / 8.0)
    if np.any(m2):
        nn = n[m2]
        lg = np.log1p(d / alpha)
        Ta = _series_terms(nn, alpha)[2:]
        Tb = _series_terms(nn, alpha + d)[2:]
        out[m2] = nn * lg + (Tb - Ta).sum(axis=0)
    m3 = ~m2
    if np.any(m3):
        out[m3] = _lgd(alpha + n[m3], d) - _lgd(np.array([alpha]), d)[0]
    return out


def _psidiff(alpha: float, n: np.ndarray) -> np.ndarray:
    """psi0(alpha + n) - psi0(alpha) for alpha > 0, n >= 0 (relative accuracy also for alpha >> n)."""
    if alpha >= _Z0:
        L = np.log1p(n / alpha)
        out = L + n / (2.0 * alpha) / (alpha + n)
        ia2 = 1.0 / (alpha * alpha)
        p = ia2
        for k in range(1, 13):
            out = out - _BERN[2 * k] / (2.0 * k) * p * np.expm1(-2.0 * k * L)
            p = p * ia2
        return out
    return digamma(alpha + n) - digamma(alpha)


def _psi_inv_scalar(y: float) -> float:
    x = math.exp(y) + 0.5 if y >= -2.22 else -1.0 / (y - PSI1)
    for _ in range(80):
        xn = x - (float(digamma(x)) - y) / float(zeta(2.0, x))
        if not xn > 0.0:
            xn = 0.5 * x
        if abs(xn - x) <= 2e-15 * abs(xn):
            return xn
        x = xn
    return x


def psi_inv(y: ArrayLike):
    """Inverse digamma psi0^-1: R -> (0, inf) (psi0 is increasing and concave on (0, inf)).
    Minka (2003, App. C) starting point + safeguarded Newton iterations to full double precision; y <= 700."""
    y = np.asarray(y, dtype=float)
    if y.ndim == 0:
        yv = float(y)
        if yv > _Y_MAX:
            raise ValueError("psi_inv: argument > 700 (psi0^-1(y) ~ exp(y) overflows float64)")
        return _psi_inv_scalar(yv)
    y1 = y
    if np.any(y1 > _Y_MAX):
        raise ValueError("psi_inv: argument > 700 (psi0^-1(y) ~ exp(y) overflows float64)")
    return np.array([_psi_inv_scalar(float(v)) for v in y1.ravel()]).reshape(y1.shape)


def hg_to_da(h: float, gamma: float, unit: str = "nats") -> Tuple[float, float]:
    """(h, gamma) -> (d, alpha): d = 1 - u, alpha = v - 1 (Eq. 2/3).  Needs h > 0 and psi0(1) - h gamma, psi0(1) + h (1 - gamma)
    <= 700 (float64 range); d in (-inf, 1), alpha > -d."""
    h = float(_to_nats(float(h), unit))
    gamma = float(gamma)
    u = psi_inv(PSI1 - h * gamma)
    v = psi_inv(PSI1 + h * (1.0 - gamma))
    return 1.0 - u, v - 1.0


def da_to_hg(d: float, alpha: float, unit: str = "nats") -> Tuple[float, float]:
    """(d, alpha) -> (h, gamma), Eq. (2) (h in nats unless unit='bits')."""
    d = float(d)
    alpha = float(alpha)
    h = float(digamma(alpha + 1.0) - digamma(1.0 - d))
    gamma = float((PSI1 - digamma(1.0 - d)) / h)
    if unit == "bits":
        h /= LN2
    elif unit != "nats":
        raise ValueError("unit must be 'nats' or 'bits'")
    return h, gamma


hg_to_dalpha = hg_to_da
dalpha_to_hg = da_to_hg


def _kappa_raw(n: np.ndarray, d: float, alpha: float) -> np.ndarray:
    """Un-clipped evaluation of Eq. (13) for n >= 1 (array), d < 1, alpha > -d."""
    out = np.ones_like(n)
    msk = n > 1.0
    if not np.any(msk):
        return out
    nm = n[msk]
    if alpha > 0.0:
        if d == 0.0:
            ek = alpha * _psidiff(alpha, nm)                         # Dirichlet-process limit
        else:
            S = _logratio(alpha, d, nm)                              # E K_n = (alpha/d) [e^S - 1]
            if d > 0.0:
                ek = np.empty_like(S)
                sm = S <= 30.0
                ek[sm] = (alpha / d) * np.expm1(S[sm])
                if not np.all(sm):                                   # avoid overflow of e^S for tiny alpha
                    Sb = S[~sm]
                    ek[~sm] = np.exp(np.log(alpha / d) + Sb) * (-np.expm1(-Sb))
            else:
                ek = (alpha / d) * np.expm1(S)
    else:
        # alpha <= 0 (then d > -alpha >= 0): E K_n = [(alpha+d) R - alpha]/d, R = (z+d)_{n-1}/(z)_{n-1}, z = alpha+1
        z = alpha + 1.0
        R = np.exp(_lgd(z + nm - 1.0, d) - _lgd(np.array([z]), d)[0])
        ek = ((alpha + d) * R - alpha) / d
    out[msk] = ek / nm
    return out


def _check_n(n):
    n_arr = np.asarray(n, dtype=float)
    scalar = n_arr.ndim == 0
    n1 = np.atleast_1d(n_arr).ravel()
    if np.any(~(n1 >= 1.0)) or np.any(~np.isfinite(n1)):
        raise ValueError("gallery size n must be finite and >= 1")
    return n_arr, n1, scalar


def kappa_da(n: ArrayLike, d: float, alpha: float):
    """E[kappa(n) | d, alpha], Eq. (13) / SI S31, for n >= 1 (scalar or array).  Needs d < 1 and alpha + d > 0
    (d <= 0 is the finite-support regime, alpha < 0 the heavy-tail corner).  kappa(1) = 1."""
    d = float(d)
    alpha = float(alpha)
    if not (np.isfinite(d) and np.isfinite(alpha) and d < 1.0 and alpha + d > 0.0):
        raise ValueError("need d < 1 and alpha + d > 0, got d=%r alpha=%r" % (d, alpha))
    n_arr, n1, scalar = _check_n(n)
    k = _kappa_raw(n1, d, alpha)
    k = np.minimum(1.0, np.maximum(1.0 / n1, k))                    # exact bounds 1/n <= kappa <= 1
    return float(k[0]) if scalar else k.reshape(n_arr.shape)


kappa_dalpha = kappa_da


def _log1mexp(x: float) -> float:
    """log(1 - e^x) for x < 0."""
    return math.log(-math.expm1(x)) if x > -LN2 else math.log1p(-math.exp(x))


def _kappa_uniform(n: np.ndarray, h: float) -> np.ndarray:
    """Uniform law on m = e^h equiprobable classes: (m/n)(1 - (1-1/m)^n)  (S40; the gamma -> -inf limit)."""
    n = np.asarray(n, dtype=float)
    if h > _Y_MAX:
        return np.ones_like(n)                                      # m = e^h > 1e304 >> n: kappa = 1 to double precision
    lg = _log1mexp(-h)                                              # log(1 - 1/m)
    t = -np.expm1(n * lg)                                           # 1 - (1 - 1/m)^n
    k = math.exp(h) * t / n
    return np.minimum(1.0, np.maximum(1.0 / n, k))


def _kappa_hg(n1: np.ndarray, h: float, gamma: float) -> np.ndarray:
    """kappa(n | h, gamma) on a validated 1-d array n1 >= 1 (float64), h > 0."""
    yu = PSI1 - h * gamma
    yv = PSI1 + h * (1.0 - gamma)
    if yu > _Y_UNIFORM:
        return _kappa_uniform(n1, h)                                # |d| > e^36: uniform law, error < 1e-15
    if yv > _Y_MAX:
        return np.ones_like(n1)                                     # alpha > e^700 >> n: kappa = 1 to double precision
    u = psi_inv(yu)
    v = psi_inv(yv)
    d, a = 1.0 - u, v - 1.0
    k = _kappa_raw(n1, d, a)
    return np.minimum(1.0, np.maximum(1.0 / n1, k))


@dataclass(frozen=True)
class PYCParams:
    """Parameters of the PYC law: entropy h (nats, > 0) and tail complexity gamma (real)."""

    h: float
    gamma: float

    @property
    def d(self) -> float:
        return _dalpha_safe(self.h, self.gamma)[0]

    @property
    def alpha(self) -> float:
        return _dalpha_safe(self.h, self.gamma)[1]

    @property
    def h_bits(self) -> float:
        return self.h / LN2

    @property
    def regime(self) -> str:
        if PSI1 - self.h * self.gamma > _Y_UNIFORM:
            return "gamma -> -inf: uniform law on e^h classes (d -> -inf), kappa = (m/n)(1 - (1 - 1/m)^n)"
        d = self.d
        if abs(d) < 1e-12:
            return "gamma = 0: d = 0, Dirichlet process (geometric tail)"
        if d < 0:
            return "gamma < 0: d < 0, finite-support (Dirichlet-multinomial) regime, kappa ~ m/n beyond a knee"
        if self.alpha < 0:
            return "gamma > 1: alpha < 0 (-d < alpha < 0), heavy tail"
        return "0 < gamma <= 1: 0 < d < 1, alpha >= 0, heavy (power-law) tail kappa ~ n^(d-1)"

    def as_tuple(self) -> Tuple[float, float]:
        return (self.h, self.gamma)

    def __iter__(self):
        yield self.h
        yield self.gamma

    def kappa(self, n: ArrayLike):
        return kappa(n, self)

    @classmethod
    def from_hg(cls, h, gamma, unit="nats", **kw):
        return cls(h=float(_to_nats(float(h), unit)), gamma=float(gamma), **kw)

    @classmethod
    def from_da(cls, d, alpha, **kw):
        h, g = da_to_hg(d, alpha)
        return cls(h=float(h), gamma=float(g), **kw)


def _dalpha_safe(h, gamma):
    try:
        return hg_to_da(h, gamma)
    except ValueError:
        return (-math.inf, math.inf)


@dataclass(frozen=True)
class PYCFit(PYCParams):
    """Result of `fit`: the parameters plus diagnostics (aliases: nfev/n_evals, success/converged,
    at_bounds/bounds_hit)."""

    loss: float = float("nan")                  # R(h, gamma), Eq. (4), weights as used
    rmse: float = float("nan")                  # plain RMSE of the fitted points (n > 1)
    n_points: int = 0
    n_evals: int = 0
    converged: bool = False                     # Nelder-Mead convergence flag of the best run
    bounds_hit: Tuple[str, ...] = field(default_factory=tuple)   # subset of h_low, h_high, gamma_low, gamma_high
    n_optima: int = 1                           # distinct local optima found (flat valleys show up as several)
    start: Tuple[float, float] = (float("nan"), float("nan"))
    message: str = ""

    @property
    def nfev(self) -> int:
        return self.n_evals

    @property
    def success(self) -> bool:
        return self.converged

    @property
    def at_bounds(self) -> Tuple[str, ...]:
        return self.bounds_hit


def _as_params(params, gamma=None, unit="nats") -> Tuple[float, float]:
    if gamma is not None:
        return float(_to_nats(float(params), unit)), float(gamma)
    if isinstance(params, PYCParams):
        return params.h, params.gamma
    if isinstance(params, dict):
        return float(_to_nats(float(params["h"]), params.get("unit", unit))), float(params["gamma"])
    h, g = params
    return float(_to_nats(float(h), unit)), float(g)


def kappa(n: ArrayLike, params, gamma=None, unit: str = "nats"):
    """E[kappa(n) | h, gamma]: expected fraction correctly identified in a gallery of n >= 1 records (scalar or array).
    `params` = PYCParams / PYCFit / (h, gamma) / {'h':..,'gamma':..}, or kappa(n, h, gamma); h in `unit` ('nats' default,
    'bits' reads h * ln 2).  Eq. (13) in the (h, gamma) parametrisation of Eq. (2); kappa(1) = 1."""
    h, g = _as_params(params, gamma, unit)
    if not (h > 0.0 and np.isfinite(h) and np.isfinite(g)):
        raise ValueError("need h > 0 and finite gamma, got h=%r gamma=%r" % (h, g))
    n_arr, n1, scalar = _check_n(n)
    k = _kappa_hg(n1, h, g)
    return float(k[0]) if scalar else k.reshape(n_arr.shape)


def predict(params, n_target: ArrayLike, unit: str = "nats"):
    """Forecast E[kappa(n') | h_hat, gamma_hat] at gallery size(s) n' (PYC-MB prediction)."""
    return kappa(n_target, params, unit=unit)


# ----------------------------------------------------------------------------------------------------------
# PYC-MB fit
# ----------------------------------------------------------------------------------------------------------
def _prepare(ns, kappas, weights):
    ns = np.atleast_1d(np.asarray(ns, dtype=float)).ravel()
    ks = np.atleast_1d(np.asarray(kappas, dtype=float)).ravel()
    if ns.shape != ks.shape:
        raise ValueError("ns and kappas must have the same length")
    if not (np.all(np.isfinite(ns)) and np.all(np.isfinite(ks))):
        raise ValueError("non-finite input")
    if np.any(ns < 1.0):
        raise ValueError("gallery sizes must be >= 1")
    if np.any(ks < -1e-9) or np.any(ks > 1.0 + 1e-9):
        raise ValueError("kappas must lie in [0, 1]")
    ks = np.clip(ks, 0.0, 1.0)
    if weights is None:
        w = np.log(ns)                                      # Eq. (4): natural-log weighting
    else:
        w = np.atleast_1d(np.asarray(weights, dtype=float)).ravel()
        if w.shape != ns.shape or np.any(w < 0) or not np.all(np.isfinite(w)):
            raise ValueError("weights must be finite, non-negative, same length as ns")
    informative = (ns > 1.0) & (w > 0.0)                    # n = 1 carries no information: kappa(1) = 1 for every law
    if np.unique(ns[informative]).size < 2:
        raise ValueError("need t >= 2 measurements at distinct gallery sizes n > 1 with positive weight (two free parameters)")
    order = np.argsort(ns, kind="stable")
    return ns[order], ks[order], w[order]


def loss(params, ns, kappas, weights=None, unit: str = "nats") -> float:
    """R(h, gamma) = sum_i w_i [kappa_hat_i - E kappa(n_i | h, gamma)]^2, w_i = ln(n_i) unless given (Eq. 4)."""
    ns, ks, w = _prepare(ns, kappas, weights)
    r = ks - np.asarray(kappa(ns, params, unit=unit), float)
    return float(np.sum(w * r * r))


def _make_resid(nv, kv, wv):
    sw = np.sqrt(wv)

    def resid(h, g):
        try:
            with np.errstate(all="ignore"):
                r = sw * (kv - _kappa_hg(nv, h, g))
            return r if np.all(np.isfinite(r)) else np.full(nv.shape, 1e3)
        except (ValueError, FloatingPointError, OverflowError, ZeroDivisionError):
            return np.full(nv.shape, 1e3)
    return resid


_GAMMA_GRID = np.array([-50.0, -20.0, -10.0, -6.0, -4.0, -3.0, -2.2, -1.6, -1.2, -0.9, -0.6, -0.4, -0.25, -0.12, 0.0, 0.1,
                        0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.35, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0,
                        6.5, 8.0, 10.0, 13.0, 17.0, 25.0, 35.0, 50.0])


def _make_raw_loss(nv, kv, wv):
    def raw(h, g):
        try:
            with np.errstate(all="ignore"):
                r = kv - _kappa_hg(nv, h, g)
                val = float(np.dot(wv, r * r))
            return val if np.isfinite(val) else 1e30
        except (ValueError, FloatingPointError, OverflowError, ZeroDivisionError):
            return 1e30
    return raw


class _ExactFit(Exception):
    def __init__(self, x):
        super().__init__("exact fit")
        self.x = np.array(x, dtype=float)


_F_EXACT = 1e-16     # loss below which the data are reproduced to ~1e-8 in kappa: stop (noise-free curves)


def _profile_starts(raw, lo, hi, gnodes, nh, n_best):
    """Profile scan: for each gamma node minimise the loss over ln h (coarse scan + bounded Brent refinement), then return the
    `n_best` best separated minima of the profile P(gamma) as starting points (ln h, gamma).  The profile is much less
    fragile than the 2-d grid minima when the valley is narrow in h and long in gamma."""
    lh = np.linspace(lo[0], hi[0], nh)
    prof = []
    for gj in gnodes:
        vals = np.array([raw(math.exp(a), float(gj)) for a in lh])
        mins = [i for i in range(nh) if vals[i] <= (vals[i - 1] if i > 0 else np.inf) and vals[i] <= (vals[i + 1] if i < nh - 1 else np.inf)]
        mins = sorted(mins, key=lambda i: vals[i])[:2]
        if len(mins) == 2 and vals[mins[1]] > 1.5 * vals[mins[0]] + 1e-12:
            mins = mins[:1]
        for i in mins:
            a0, a1 = lh[max(i - 1, 0)], lh[min(i + 1, nh - 1)]
            best_a, best_v = lh[i], vals[i]
            if a1 > a0:
                r = optimize.minimize_scalar(lambda a: raw(math.exp(a), float(gj)), bounds=(a0, a1), method="bounded",
                                             options=dict(xatol=1e-6))
                if r.fun < best_v:
                    best_a, best_v = float(r.x), float(r.fun)
            prof.append((best_v, best_a, float(gj)))
    prof.sort(key=lambda t: (t[0], t[2]))
    picked = []
    for v, a, g in prof:
        if picked and (v > 2.0 * picked[0][0] + 1e-10 or picked[0][0] < _F_EXACT):   # cannot become the global optimum: skip
            break
        dup = False
        for pv, pa, pg in picked:
            same_floor = abs(v - pv) <= 1e-8 * max(pv, 1e-30) + 1e-18 and abs(a - pa) < 1e-3     # same flat valley floor
            near = abs(a - pa) <= 0.3 and abs(g - pg) <= 0.2 * max(1.0, abs(pg))
            if same_floor or near:
                dup = True
                break
        if not dup:
            picked.append((v, a, g))
        if len(picked) >= n_best:
            break
    return picked, (prof[0][0] if prof else np.inf)


def _search(raw, resid, hbnd, gbnd, starts, n_best, restarts, maxiter, xatol, fatol, grid, nh=32):
    """Global search of the loss in (ln h, gamma): profile scan -> best separated minima -> Nelder-Mead (+restarts) from them
    and from the user/authors' starts.  Returns a list of (loss, h, gamma, converged, message, start) sorted by loss and the
    number of loss evaluations."""
    hlo, hhi = hbnd
    glo, ghi = gbnd
    lo = np.array([math.log(hlo), glo])
    hi = np.array([math.log(hhi), ghi])
    nfev = [0]
    PEN = 1e4

    def obj(x):
        nfev[0] += 1
        xc = np.minimum(np.maximum(x, lo), hi)
        pen = float(np.sum((x - xc) ** 2))
        f = raw(math.exp(xc[0]), float(xc[1]))
        if f < _F_EXACT and pen == 0.0:
            raise _ExactFit(xc)
        return f + PEN * pen * (1.0 + f)

    cand = [(tuple(s), None) for s in starts]
    fscale = None
    if grid:
        gg = np.unique(np.concatenate([_GAMMA_GRID[(_GAMMA_GRID >= glo) & (_GAMMA_GRID <= ghi)], [glo, ghi]]))
        picked, fscale = _profile_starts(raw, lo, hi, gg, nh, n_best)
        nfev[0] += gg.size * (nh + 40)
        cand = [((math.exp(a), g), True) for _, a, g in picked] + cand
    if fatol is None:
        fatol = 1e-12 * max(fscale if fscale is not None else raw(*starts[0]), 1e-12)
    out = []
    for (h0, g0), from_grid in cand:
        x0 = np.array([math.log(min(max(h0, hlo), hhi)), min(max(g0, glo), ghi)])
        # accelerator: bounded Gauss-Newton/trust-region on the same objective (the minimiser of Eq. (4) does not depend on
        # the optimiser); Nelder-Mead below (the paper's method) then polishes from its result.
        try:
            xin = np.minimum(np.maximum(x0, lo + 1e-9), hi - 1e-9)
            lsq = optimize.least_squares(lambda x: resid(math.exp(x[0]), float(x[1])), xin, bounds=(lo, hi), method="trf",
                                         x_scale=np.array([1.0, max(1.0, abs(xin[1]))]), xtol=1e-14, ftol=1e-15, gtol=1e-15,
                                         max_nfev=200)
            nfev[0] += int(lsq.nfev) * 3
            if raw(math.exp(lsq.x[0]), float(lsq.x[1])) <= raw(math.exp(x0[0]), float(x0[1])):
                x0 = np.minimum(np.maximum(lsq.x, lo), hi)
        except Exception:
            pass
        if raw(math.exp(x0[0]), float(x0[1])) < _F_EXACT:
            out.append((raw(math.exp(x0[0]), float(x0[1])), math.exp(x0[0]), float(x0[1]), True, "exact fit (loss < %g)" % _F_EXACT, (h0, g0)))
            continue
        best_x, best_f, ok, msg = x0, None, False, ""
        for it in range(max(1, restarts if from_grid else min(restarts, 1)) + 1):
            if it == 0:
                sx, sg = (0.05, max(0.03, 0.02 * abs(x0[1])))
            else:
                sx, sg = 0.03, max(0.02, 0.01 * abs(best_x[1]))
            sim = np.array([best_x, best_x + [sx, 0.0], best_x + [0.0, sg]])
            try:
                res = optimize.minimize(obj, best_x, method="Nelder-Mead",
                                        options=dict(xatol=xatol, fatol=fatol, maxiter=maxiter, maxfev=2 * maxiter, initial_simplex=sim))
            except _ExactFit as e:
                best_x, best_f, ok, msg = e.x, raw(math.exp(e.x[0]), float(e.x[1])), True, "exact fit (loss < %g)" % _F_EXACT
                break
            ok, msg = bool(res.success), str(res.message)
            improved = best_f is None or res.fun < best_f - max(1e-12 * abs(best_f), 1e-17)
            if best_f is None or res.fun <= best_f:
                best_x, best_f = np.array(res.x, float), float(res.fun)
            if not improved:
                break
        xc = np.minimum(np.maximum(best_x, lo), hi)
        out.append((raw(math.exp(xc[0]), float(xc[1])), math.exp(xc[0]), float(xc[1]), ok, msg, (h0, g0)))
    out.sort(key=lambda t: t[0])
    return out, nfev[0]


def _bound_flags(h, g, hbnd, gbnd):
    flags = []
    tol_h = 1e-5 * (math.log(hbnd[1]) - math.log(hbnd[0]))
    tol_g = 1e-5 * (gbnd[1] - gbnd[0])
    if math.log(h) <= math.log(hbnd[0]) + tol_h:
        flags.append("h_low")
    if math.log(h) >= math.log(hbnd[1]) - tol_h:
        flags.append("h_high")
    if g <= gbnd[0] + tol_g and PSI1 - h * g <= _Y_UNIFORM:        # on the uniform-law plateau gamma_low is benign (exact limit)
        flags.append("gamma_low")
    if g >= gbnd[1] - tol_g:
        flags.append("gamma_high")
    return tuple(flags)


def fit_all(ns, kappas, weights=None, *, gamma_bounds=DEFAULT_GAMMA_BOUNDS, h_bounds=DEFAULT_H_BOUNDS, start=DEFAULT_START,
            starts=None, multistart: bool = True, grid=None, n_best: int = 4, restarts: int = 3, maxiter: int = 2000,
            xatol: float = 1e-8, fatol: Optional[float] = None, _warn: bool = False) -> list:
    """All distinct local optima found by the search of `fit`, best (lowest loss) first, as a list of PYCFit."""
    ns, ks, w = _prepare(ns, kappas, weights)
    if grid is not None:
        multistart = bool(grid)
    hb = (float(h_bounds[0]), float(h_bounds[1]))
    gb = (float(gamma_bounds[0]), float(gamma_bounds[1]))
    if not (0.0 < hb[0] < hb[1] and gb[0] < gb[1]):
        raise ValueError("bad bounds")
    sel = (ns > 1.0) & (w > 0.0)
    nv, kv, wv = ns[sel], ks[sel], w[sel]
    if _warn:
        k_last = float(kv[np.argmax(nv)])
        if not (0.01 < k_last < 0.99):
            warnings.warn("kappa_hat at the largest n is %.3g: the paper (SI S3.3) evaluates PYC-MB only where it lies in "
                          "(0.01, 0.99); plateau or vanishing curves leave (h, gamma) poorly determined" % k_last,
                          RuntimeWarning, stacklevel=3)
        if np.unique(nv).size == 2:
            warnings.warn("exactly two informative points for two parameters: the fit is exactly determined when it exists "
                          "and may be non-unique (flat valley); see fit_all and forecast_band", RuntimeWarning, stacklevel=3)
    raw = _make_raw_loss(nv, kv, wv)
    resid = _make_resid(nv, kv, wv)
    st = [tuple(float(x) for x in s) for s in (starts if starts is not None else [start])]
    found, nfev = _search(raw, resid, hb, gb, st, n_best if multistart else 0, restarts, maxiter, xatol, fatol, multistart)
    res = []
    for fv, h, g, ok, msg, s0 in found:
        if any(abs(math.log(h / q.h)) < 1e-3 and abs(g - q.gamma) < 1e-3 for q in res):
            continue
        pred = _kappa_hg(nv, h, g)
        rmse = float(np.sqrt(np.mean((kv - pred) ** 2)))
        res.append(PYCFit(h=float(h), gamma=float(g), loss=float(fv), rmse=rmse, n_points=int(nv.size), n_evals=int(nfev),
                          converged=bool(ok), bounds_hit=_bound_flags(h, g, hb, gb), n_optima=1, start=(float(s0[0]), float(s0[1])),
                          message=msg))
    # distinct optima whose loss is within 1e-6 relative of the best are the same valley seen from several starts
    nopt = len([q for q in res if q.loss <= res[0].loss * (1 + 1e-6) + 1e-14])
    return [PYCFit(**{**q.__dict__, "n_optima": nopt if i == 0 else 1}) for i, q in enumerate(res)]


def fit(ns, kappas, weights=None, *, gamma_bounds=DEFAULT_GAMMA_BOUNDS, h_bounds=DEFAULT_H_BOUNDS, start=DEFAULT_START,
        starts=None, multistart: bool = True, grid=None, n_best: int = 4, restarts: int = 3, maxiter: int = 2000,
        xatol: float = 1e-8, fatol: Optional[float] = None) -> PYCFit:
    """PYC-MB fit (Eq. 4): minimise  sum_i w_i [kappa_i - E kappa(n_i | h, gamma)]^2  over (h, gamma) and return a PYCFit.

    ns, kappas : the t >= 2 measured points (gallery size n_i >= 1, fraction correctly identified).
    weights    : None -> w_i = ln(n_i) (the paper; n = 1 has weight 0).  An explicit array REPLACES the log weights
                 (pass np.log(ns) * v to combine them with other weights v).
    gamma_bounds, h_bounds : box of the search, default gamma in [-50, 50], h in [1e-6, 700] nats (not in the paper; use
                 gamma_bounds=(0, 1) for the authors' 2026 code, (0, 50) for d >= 0).  fit() warns when the optimum is on a bound.
    start      : the authors' start (h, gamma) = (12, 0.26) (not in the paper); starts= a list of starts.
    multistart : also search a coarse grid and start Nelder-Mead from its `n_best` best separated cells (the objective is
                 unchanged; False = single start, paper style but with tight tolerances).
    restarts, maxiter, xatol, fatol : Nelder-Mead controls (scipy defaults would be 1e-4, 1e-4, 200 iterations per dimension).
    Raises ValueError for fewer than two informative points, kappa outside [0, 1], n < 1."""
    res = fit_all(ns, kappas, weights, gamma_bounds=gamma_bounds, h_bounds=h_bounds, start=start, starts=starts,
                  multistart=multistart, grid=grid, n_best=n_best, restarts=restarts, maxiter=maxiter, xatol=xatol,
                  fatol=fatol, _warn=True)
    best = res[0]
    if best.bounds_hit:
        warnings.warn("the optimum lies on the bound(s) %s of the search box (h in %s, gamma in %s): the loss is still "
                      "decreasing there, so the forecast depends on the box, not only on the data; see forecast_band"
                      % (",".join(best.bounds_hit), tuple(h_bounds), tuple(gamma_bounds)), RuntimeWarning, stacklevel=2)
    return best


def fit_like_authors(ns, kappas, weights=None, gamma_bounds=None) -> PYCFit:
    """The authors' published procedure, for comparison only: ONE scipy Nelder-Mead run with default options (xatol = fatol
    = 1e-4) in (h, gamma) from (12, 0.26) on mean_i(ln(n_i) (kappa_i - model_i)^2).  gamma_bounds=None: unconstrained as in
    their November-2024 code (only h > 0 enforced here); (0, 1): the 1e10 penalty of their 2026 release.  The kernel is this
    module's (their float formula loses accuracy for h >~ 15 at small gamma)."""
    ns, ks, w = _prepare(ns, kappas, weights)
    sel = (ns > 1.0) & (w > 0.0)
    nv, kv, wv = ns[sel], ks[sel], w[sel] / ns.size
    raw = _make_raw_loss(nv, kv, wv)
    nfev = [0]

    def f(x):
        nfev[0] += 1
        h, g = float(x[0]), float(x[1])
        if not (h > 0.0) or h > 700.0:
            return 1e10
        if gamma_bounds is not None and not (gamma_bounds[0] <= g <= gamma_bounds[1]):
            return 1e10
        v = raw(h, g)
        return v if v < 1e29 else 1e10

    res = optimize.minimize(f, np.array(DEFAULT_START, float), method="Nelder-Mead")
    h, g = float(res.x[0]), float(res.x[1])
    r = kv - _kappa_hg(nv, h, g) if h > 0 else np.full(nv.shape, np.nan)
    return PYCFit(h=h, gamma=g, loss=float(res.fun) * ns.size, rmse=float(np.sqrt(np.mean(r ** 2))), n_points=int(nv.size),
                  n_evals=int(res.nfev), converged=bool(res.success), bounds_hit=(), n_optima=1, start=DEFAULT_START,
                  message=str(res.message))


REGIMES = {
    "authors_2026_gamma_0_1": (0.0, 1.0),      # d in [0,1), alpha >= 0: the 2026 release of the authors' package
    "pitman_yor_d_ge_0": (0.0, 50.0),          # d >= 0, alpha > -d: the domain written in the paper's text
    "unconstrained_box": DEFAULT_GAMMA_BOUNDS,  # adds d < 0 (the paper uses gamma = -1 in its Fig. 4 and its worked example needs it)
}


def fit_regimes(ns, kappas, n_target=None, weights=None, regimes=None) -> dict:
    """Fit under several gamma ranges and return {label: (PYCFit, forecast at n_target or None)}: shows how much a forecast
    depends on the admitted range of the discount (worked example: 0.346 for gamma in [0,1], 0.2145 unconstrained)."""
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for label, gb in (regimes or REGIMES).items():
            r = fit(ns, kappas, weights, gamma_bounds=gb)
            out[label] = (r, None if n_target is None else kappa(n_target, r))
    return out


def forecast_band(ns, kappas, n_target, delta, weights=None, *, gamma_bounds=DEFAULT_GAMMA_BOUNDS,
                  h_bounds=DEFAULT_H_BOUNDS, best: Optional[PYCFit] = None, n_g: int = 60):
    """Sensitivity of a forecast (NOT in the paper): the smallest and largest kappa(n_target) over all (h, gamma) of the box
    whose loss, Eq. (4), is at most `delta` above its minimum.  For each gamma of a grid (n_g nodes + the fit's gamma
    neighbourhood) the feasible interval in ln h is found exactly (profile minimum by Brent, endpoints by bisection) and
    kappa(n_target) is evaluated at its ends; kappa is monotone in h at fixed gamma, so these are the extremes for that gamma.
    `delta` is in the units of Eq. (4); for independent errors of standard deviation s on each kappa_hat a heuristic
    one-sigma scale is delta ~ 2.3 s^2 mean(w).  Returns (kappa_min, kappa_max, best), arrays shaped like n_target.
    The extremes are over the gamma grid, so they slightly under-estimate the true range; the band is a profile-type
    sensitivity measure, not a confidence interval (errors of real curves are correlated)."""
    ns_, ks_, w = _prepare(ns, kappas, weights)
    sel = (ns_ > 1.0) & (w > 0.0)
    nv, kv, wv = ns_[sel], ks_[sel], w[sel]
    if best is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            best = fit(ns, kappas, weights, gamma_bounds=gamma_bounds, h_bounds=h_bounds)
    raw = _make_raw_loss(nv, kv, wv)
    thr = best.loss + float(delta)
    nt = np.atleast_1d(np.asarray(n_target, dtype=float))
    la, lb = math.log(h_bounds[0]), math.log(h_bounds[1])
    nodes = np.linspace(la, lb, 40)
    gs = np.concatenate([np.linspace(gamma_bounds[0], gamma_bounds[1], n_g), _GAMMA_GRID,
                         best.gamma + np.linspace(-1.0, 1.0, 61) * max(0.5, 0.5 * abs(best.gamma)), [best.gamma]])
    gs = np.unique(gs[(gs >= gamma_bounds[0]) & (gs <= gamma_bounds[1])])
    lo = np.full(nt.shape, np.inf)
    hi = np.full(nt.shape, -np.inf)

    def upd(a, g):
        nonlocal lo, hi
        k = _kappa_hg(nt, math.exp(a), float(g))
        lo = np.minimum(lo, k)
        hi = np.maximum(hi, k)

    def bisect(f, a_in, a_out):
        """a_in feasible (f <= thr), a_out infeasible: boundary between them."""
        for _ in range(32):
            mid = 0.5 * (a_in + a_out)
            if f(mid) <= thr:
                a_in = mid
            else:
                a_out = mid
        return a_in

    for g in gs:
        g = float(g)

        def f(a, g=g):
            return raw(math.exp(a), g)

        vals = np.array([f(a) for a in nodes])
        i = int(np.argmin(vals))
        a_star, p = float(nodes[i]), float(vals[i])
        r = optimize.minimize_scalar(f, bounds=(float(nodes[max(i - 1, 0)]), float(nodes[min(i + 1, nodes.size - 1)])),
                                    method="bounded", options=dict(xatol=1e-7))
        if r.fun < p:
            a_star, p = float(r.x), float(r.fun)
        if p > thr:
            continue
        upd(a_star, g)
        # left end: first infeasible node to the left of a_star
        left = [a for a in nodes[nodes < a_star][::-1]]
        a_in = a_star
        found = False
        for a in left:
            if f(a) <= thr:
                a_in = a
            else:
                upd(bisect(f, a_in, float(a)), g)
                found = True
                break
        if not found:
            upd(a_in, g)
        right = [a for a in nodes[nodes > a_star]]
        a_in = a_star
        found = False
        for a in right:
            if f(a) <= thr:
                a_in = a
            else:
                upd(bisect(f, a_in, float(a)), g)
                found = True
                break
        if not found:
            upd(a_in, g)
    kb = _kappa_hg(nt, best.h, best.gamma)
    lo = np.minimum(lo, kb)
    hi = np.maximum(hi, kb)
    return lo, hi, best
