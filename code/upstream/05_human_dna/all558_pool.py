"""Closed-set check inside the pool of the donors' runs: is the other sample of each donor still the nearest run of the donor's query when the human hashes are removed,
and what happens when only the human hashes are used? Similarity = Jaccard of the sketches (hash sets), as in the pair analysis.
Human hashes = hashes found exactly in T2T-CHM13v2.0 / GRCh38.p14; near-human hashes (optional) = 31-mers that align to T2T with at most 2 mismatches (human variants).
usage: python3 -I all558_pool.py WORKDIR SKDIR PAIRS.json OUT_PREFIX [NEAR.npz]
Writes OUT_PREFIX_summary.json (aggregates), OUT_PREFIX_ranks.csv (one line per donor) and OUT_PREFIX_unrelated.npz (similarities for the figure)."""
import csv, json, os, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata

W, SK, PAIRS, OUT = sys.argv[1:5]
NEAR = sys.argv[5] if len(sys.argv) > 5 else None
t0 = time.time()
pairs = [p for p in json.load(open(PAIRS)) if os.path.exists(f"{SK}/{p['q']}.npz") and os.path.exists(f"{SK}/{p['partner']}.npz")]
runs = sorted({r for p in pairs for r in (p["q"], p["partner"])})
idx = {r: i for i, r in enumerate(runs)}
H = np.union1d(np.load(f"{W}/ref_t2t.npy"), np.load(f"{W}/ref_grch38.npy"))
N = np.load(NEAR)["h"] if NEAR else np.empty(0, np.uint64)
hs = [np.load(f"{SK}/{r}.npz")["h"] for r in runs]
n = np.array([len(h) for h in hs], np.int64)
print(f"{len(pairs)} donors, pool of {len(runs)} runs, {n.sum() / 1e6:.0f} M hashes, {len(N):,} near-human hashes, {time.time() - t0:.0f} s", flush=True)

allh = np.concatenate(hs); del hs
uniq, inv = np.unique(allh, return_inverse=True); del allh
rows = np.repeat(np.arange(len(runs), dtype=np.int32), n)
d = np.bincount(inv, minlength=len(uniq))                       # in how many runs each hash occurs
hm = np.isin(uniq, H, assume_unique=True)                       # hashes found exactly in the human reference
nm = np.isin(uniq, N, assume_unique=True) if len(N) else np.zeros(len(uniq), bool)
nh = np.bincount(rows[hm[inv]], minlength=len(runs)).astype(np.int64)
nn = np.bincount(rows[nm[inv]], minlength=len(runs)).astype(np.int64)
ke = d[inv] >= 2                                                # a hash that occurs in one run only cannot be shared
kept = d >= 2
newcol = np.cumsum(kept, dtype=np.int64) - 1
X = sp.csr_matrix((np.ones(int(ke.sum()), np.float32), (rows[ke], newcol[inv[ke]].astype(np.int32))), shape=(len(runs), int(kept.sum())))
XH = X[:, np.nonzero(hm[kept])[0]]; XN = X[:, np.nonzero(nm[kept])[0]]
print(f"{len(uniq) / 1e6:.0f} M distinct hashes, {int(kept.sum()) / 1e6:.1f} M in two runs or more, {int(hm.sum())} human, {int(nm.sum())} near-human; matrix built {time.time() - t0:.0f} s", flush=True)

qi = np.array([idx[p["q"]] for p in pairs]); pj = np.array([idx[p["partner"]] for p in pairs]); arow = np.arange(len(qi))
I = (X[qi] @ X.T).toarray().astype(np.float64)                  # shared hashes, queries x runs
IH = (XH[qi] @ XH.T).toarray().astype(np.float64)               # shared exact-human hashes
IN = (XN[qi] @ XN.T).toarray().astype(np.float64) if len(N) else np.zeros_like(I)
print(f"products done {time.time() - t0:.0f} s", flush=True)
c = lambda v: (v[qi][:, None], v[None, :])
(nq, nj), (hq, hj), (vq, vj) = c(n), c(nh), c(nn)


def jac(i, a, b):
    den = a + b - i
    return np.divide(i, den, out=np.zeros_like(i), where=den > 0)


defs = {"all hashes": jac(I, nq, nj), "exact human removed": jac(I - IH, nq - hq, nj - hj), "exact human only": jac(IH, hq, hj)}
if len(N):
    defs.update({"exact and near human removed": jac(I - IH - IN, nq - hq - vq, nj - hj - vj), "exact and near human only": jac(IH + IN, hq + vq, hj + vj),
                 "near human only": jac(IN, vq, vj)})
study = np.empty(len(runs), object)
for p in pairs:
    study[idx[p["q"]]] = p["study"]; study[idx[p["partner"]]] = p["study"]
cols = np.arange(len(runs))[None, :]
unrel = (cols != qi[:, None]) & (cols != pj[:, None])           # neither the query itself nor the donor's other sample: a different person
same_c = unrel & (study[None, :] == study[qi][:, None]); diff_c = unrel & ~same_c
share = 100 * np.maximum(nh[qi] / n[qi], nh[pj] / n[pj])        # human share of the more human of the two samples (%)
rich = share >= 5
auc = lambda pos, neg: float((rankdata(np.concatenate([pos, neg]))[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))) if len(pos) and len(neg) else None
rs = np.random.default_rng(20261024)
pick = lambda v, m=30000: v[rs.permutation(len(v))[:m]]
res = dict(pool=dict(donors=len(pairs), runs=len(runs), donors_with_human_share_5_percent_or_more=int(rich.sum())), definitions={})
save, ranks = {}, {}
for name, M in defs.items():
    pv = M[arow, pj]
    row = dict(same_person=dict(median=float(np.median(pv)), p2_5=float(np.percentile(pv, 2.5)), p97_5=float(np.percentile(pv, 97.5))),
               unrelated_same_cohort=dict(median=float(np.median(M[same_c])), p97_5=float(np.percentile(M[same_c], 97.5)), n=int(same_c.sum())),
               unrelated_other_cohort=dict(median=float(np.median(M[diff_c])), p97_5=float(np.percentile(M[diff_c], 97.5)), n=int(diff_c.sum())))
    row["auc_vs_same_cohort"] = auc(pv, M[same_c]); row["auc_vs_other_cohort"] = auc(pv, M[diff_c])
    for tag, m in (("human_share_below_5", ~rich), ("human_share_5_or_more", rich)):
        if m.any():
            row["auc_vs_same_cohort_" + tag] = auc(pv[m], M[same_c & m[:, None]])
    Mr = M.copy(); Mr[arow, qi] = -1.0                          # a run is not its own neighbour
    ps = Mr[arow, pj]
    above = (Mr > ps[:, None]).sum(axis=1); ties = (Mr == ps[:, None]).sum(axis=1) - 1
    top1 = (above + ties) == 0; top10 = (above + ties) < 10      # ties count against the other sample
    row.update(top1=int(top1.sum()), top1_percent=float(100 * top1.mean()), top10=int(top10.sum()), top10_percent=float(100 * top10.mean()), median_rank_optimistic=float(np.median(above + 1)))
    for tag, m in (("human_share_below_5", ~rich), ("human_share_5_or_more", rich)):
        if m.any():
            row["top1_" + tag] = dict(n=int(m.sum()), top1=int(top1[m].sum()), percent=float(100 * top1[m].mean()))
    other = np.where(np.arange(M.shape[1])[None, :] == pj[:, None], -1.0, Mr).max(axis=1)
    row["margin_partner_minus_best_other_median"] = float(np.median(ps - other))
    res["definitions"][name] = row
    ranks[name] = (above + 1, ties)
    key = name.replace(" ", "_")
    save[f"partner__{key}"] = pv; save[f"same__{key}"] = pick(M[same_c]); save[f"diff__{key}"] = pick(M[diff_c])
    save[f"partner_rich__{key}"] = pv[rich]; save[f"partner_poor__{key}"] = pv[~rich]
    save[f"same_rich__{key}"] = pick(M[same_c & rich[:, None]]); save[f"same_poor__{key}"] = pick(M[same_c & ~rich[:, None]])
base = ranks["all hashes"][0]
res["rank_unchanged_when_removed"] = {k: dict(n=int((v[0] == base).sum()), of=len(pairs)) for k, v in ranks.items() if "removed" in k}
res["human_share_percent_of_donors"] = dict(median=float(np.median(share)), n_5_or_more=int(rich.sum()), cohorts_with_5_or_more={s: int(((study[qi] == s) & rich).sum()) for s in sorted(set(study[qi])) if ((study[qi] == s) & rich).any()})
np.savez(f"{OUT}_unrelated.npz", human_share=share, **save)
json.dump(res, open(f"{OUT}_summary.json", "w"), indent=1)
with open(f"{OUT}_ranks.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["q", "partner", "study", "human_share_percent"] + [f"{k}|{t}" for k in ranks for t in ("rank", "ties")])
    for k, p in enumerate(pairs):
        w.writerow([p["q"], p["partner"], p["study"], f"{share[k]:.4f}"] + [x for nme in ranks for x in (ranks[nme][0][k], ranks[nme][1][k])])
print(json.dumps({k: v for k, v in res.items() if k != "definitions"}, indent=1))
for name, row in res["definitions"].items():
    print(f"{name:30} top1 {row['top1']:4d}/{len(pairs)} ({row['top1_percent']:.1f}%)  top10 {row['top10_percent']:.1f}%  AUC vs same cohort {row['auc_vs_same_cohort']:.3f}  "
          f"same person median {row['same_person']['median']:.3f}  unrelated same cohort median {row['unrelated_same_cohort']['median']:.3f}")
print(f"done {time.time() - t0:.0f} s")
