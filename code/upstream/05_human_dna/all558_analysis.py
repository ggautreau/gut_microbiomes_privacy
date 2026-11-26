"""Is the similarity of the two samples of a donor carried by human DNA? Exact-hash analysis of every donor pair, with an optional layer
for the near-human 31-mers (human variants: aligned to T2T with at most 2 mismatches by all558_variants.py near).
For each pair: share of each sketch that is found in the human reference (T2T-CHM13v2.0 + GRCh38.p14), share of the shared hashes that is found in it,
the similarity recomputed after removing every hash that is found in it (exact, then exact + near), and the similarity of the human hashes alone.
usage: python3 -I all558_analysis.py WORKDIR SKDIR PAIRS.json OUT_PREFIX [NEAR.npz]
Writes OUT_PREFIX_table.csv (one line per donor pair) and OUT_PREFIX_summary.json (aggregates)."""
import csv, json, os, sys
import numpy as np

W, SK, PAIRS, OUT = sys.argv[1:5]
NEAR = sys.argv[5] if len(sys.argv) > 5 else None
pairs = json.load(open(PAIRS))
H = np.union1d(np.load(f"{W}/ref_t2t.npy"), np.load(f"{W}/ref_grch38.npy"))
N = np.load(NEAR)["h"] if NEAR else np.empty(0, np.uint64)
COMPLETE = set(np.load(NEAR)["complete_runs"].tolist()) if NEAR else set()      # runs whose non-exact-human hashes were all aligned
print(f"human reference hashes: union {len(H):,}; near-human hashes: {len(N):,}")


def load(run):
    z = np.load(f"{SK}/{run}.npz")
    return z["h"], int(z["nbases"])


cnt = lambda x, S: int(np.isin(x, S, assume_unique=True).sum())
rows, missing = [], 0
for p in pairs:
    if not (os.path.exists(f"{SK}/{p['q']}.npz") and os.path.exists(f"{SK}/{p['partner']}.npz")):
        missing += 1; continue
    A, bA = load(p["q"]); B, bB = load(p["partner"])
    inter = np.intersect1d(A, B, assume_unique=True)
    hA, hB, ih = cnt(A, H), cnt(B, H), cnt(inter, H)
    use_near = len(N) > 0 and p["q"] in COMPLETE and p["partner"] in COMPLETE
    nA_, nB_, ns = (cnt(A, N), cnt(B, N), cnt(inter, N)) if use_near else (np.nan, np.nan, np.nan)
    union = len(A) + len(B) - len(inter); union_h = hA + hB - ih; union_n = nA_ + nB_ - ns
    J = len(inter) / union
    Jn = (len(inter) - ih) / (union - union_h)
    Jn2 = (len(inter) - ih - ns) / (union - union_h - union_n) if use_near else float("nan")
    d = lambda a, b: a / b if (b and b == b) else float("nan")
    rows.append(dict(study=p["study"], subject=p["subject"], q=p["q"], partner=p["partner"], service=p["service"], nA=len(A), nB=len(B), gbA=bA / 1e9, gbB=bB / 1e9,
                     humanA=hA, humanB=hB, shared=len(inter), shared_human=ih, nearA=nA_, nearB=nB_, shared_near=ns,
                     J=J, J_nonhuman=Jn, J_nonhuman_near=Jn2, J_human_only=d(ih, union_h), J_near_only=d(ns, union_n), J_human_near_only=d(ih + ns, union_h + union_n)))
print(f"{len(rows)} pairs analysed, {missing} not yet available")
cols = list(rows[0])
with open(f"{OUT}_table.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)

a = lambda k: np.array([r[k] if r[k] is not None else np.nan for r in rows], float)
J, Jn, Jn2, Jh, Jv, Jhv = a("J"), a("J_nonhuman"), a("J_nonhuman_near"), a("J_human_only"), a("J_near_only"), a("J_human_near_only")
nA, nB, hA, hB, sh, ih, vA, vB, vs = a("nA"), a("nB"), a("humanA"), a("humanB"), a("shared"), a("shared_human"), a("nearA"), a("nearB"), a("shared_near")
share = 100 * np.maximum(hA / nA, hB / nB)                 # more human of the two sketches (%)
share_n = 100 * np.maximum(vA / nA, vB / nB)               # more near-human of the two sketches (%)
share_sh = 100 * ih / sh
d1, d2 = Jn - J, Jn2 - J
sv = a("service"); ok = ~np.isnan(sv)
q = lambda x, p: float(np.nanpercentile(x, p))
desc = lambda x: dict(median=q(x, 50), p5=q(x, 5), p95=q(x, 95), min=float(np.nanmin(x)), max=float(np.nanmax(x)))
strata = {"human share < 0.5%": share < 0.5, "0.5% to 5%": (share >= 0.5) & (share < 5), "5% or more": share >= 5}
S = dict(n_pairs=len(rows), similarity=desc(J),
         service_vs_sketch=dict(n=int(ok.sum()), pearson=float(np.corrcoef(J[ok], sv[ok])[0, 1]) if ok.sum() > 2 else None,
                                median_abs_diff=float(np.median(np.abs(J[ok] - sv[ok]))) if ok.any() else None, p95_abs_diff=q(np.abs(J[ok] - sv[ok]), 95) if ok.any() else None),
         human_share_of_sketch_percent=dict(**desc(share), n_over_0_1=int((share > 0.1).sum()), n_over_1=int((share > 1).sum()), n_over_5=int((share > 5).sum()), n_over_50=int((share > 50).sum())),
         human_share_of_shared_percent=desc(share_sh), human_hashes_shared=dict(median=q(ih, 50), max=int(ih.max()), pairs_with_none=int((ih == 0).sum())),
         delta_similarity_exact_human_removed=dict(**desc(d1), n_decrease=int((d1 < 0).sum()), n_increase=int((d1 > 0).sum()), n_abs_over_0_001=int((np.abs(d1) > 0.001).sum()),
                                                  n_abs_over_0_01=int((np.abs(d1) > 0.01).sum()), relative_max_percent=float(100 * np.nanmax(np.abs(d1) / J))),
         human_only_similarity=desc(Jh))
if len(N):
    S["near_human_share_of_sketch_percent"] = dict(**desc(share_n), n_over_1=int((share_n > 1).sum()))
    S["delta_similarity_exact_and_near_removed"] = dict(**desc(d2), n_decrease=int((d2 < 0).sum()), n_increase=int((d2 > 0).sum()), n_abs_over_0_01=int((np.abs(d2) > 0.01).sum()),
                                                       relative_max_percent=float(100 * np.nanmax(np.abs(d2) / J)), n_pairs=int((~np.isnan(d2)).sum()))
    S["near_only_similarity"] = desc(Jv)
for name, m in strata.items():
    if m.any():
        S.setdefault("by_human_share", {})[name] = dict(n=int(m.sum()), similarity_median=q(J[m], 50), similarity_without_human_median=q(Jn[m], 50),
                                                       delta_median=q(d1[m], 50), delta_min=float(d1[m].min()), delta_max=float(d1[m].max()),
                                                       human_only_similarity_median=q(Jh[m], 50), human_share_of_shared_median=q(share_sh[m], 50),
                                                       **({"near_share_median": q(share_n[m], 50), "delta_with_near_median": q(d2[m], 50), "delta_with_near_min": float(d2[m].min()),
                                                           "near_only_similarity_median": q(Jv[m], 50), "shared_near_median": q(vs[m], 50)} if len(N) else {}))
by = {}
for r, s in zip(rows, share):
    by.setdefault(r["study"], []).append(s)
S["by_cohort_human_share_percent"] = {k: dict(n=len(v), median=float(np.median(v)), p90=float(np.percentile(v, 90)), max=float(np.max(v)), n_over_5=int((np.array(v) >= 5).sum()))
                                      for k, v in sorted(by.items(), key=lambda kv: -len(kv[1]))}
json.dump(S, open(f"{OUT}_summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in S.items() if k not in ("by_cohort_human_share_percent",)}, indent=1))
print("by cohort (n, median % human, p90, max, n >= 5%):", {k: (v["n"], round(v["median"], 3), round(v["p90"], 2), round(v["max"], 2), v["n_over_5"]) for k, v in S["by_cohort_human_share_percent"].items()})
