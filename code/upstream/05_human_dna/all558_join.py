"""Does the human content of a donor's samples go with the result of the closed-set search in the HGM set of 291,881 runs?
Joins the pair table (human shares, similarities recomputed from the sketches) with the rank of the donor's other sample returned by the service (rows_v3.json,
rank_any_narrow = the rank-1 definition of the paper, here for the 558 donors of Supplementary Note 1). Aggregates only in the summary.
usage: python3 -I all558_join.py ANALYSIS_table.csv ROWS_V3.json OUT_PREFIX"""
import csv, json, sys
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

TABLE, ROWS, OUT = sys.argv[1:4]
R = {(r["study"], r["subject"]): r for r in json.load(open(ROWS))}
T = list(csv.DictReader(open(TABLE)))
f = lambda r, k: float(r[k]) if r[k] not in ("", "nan") else np.nan
D = []
for r in T:
    row = R[(r["study"], r["subject"])]
    rk = row["rank_any_narrow"]
    sq, sp = 100 * f(r, "humanA") / f(r, "nA"), 100 * f(r, "humanB") / f(r, "nB")
    D.append(dict(study=r["study"], q=r["q"], partner=r["partner"], rank=(rk if rk is not None else 10**6), rank1=rk == 1, share_q=sq, share_p=sp, share=max(sq, sp),
                  J=f(r, "J"), J_non=f(r, "J_nonhuman"), J_non_near=f(r, "J_nonhuman_near"), service=row["partner_any"], days=row["days"]))
a = lambda k: np.array([d[k] if d[k] is not None else np.nan for d in D], float)
share, rank1 = a("share"), np.array([d["rank1"] for d in D]); study = np.array([d["study"] for d in D])
S = dict(n=len(D), rank1=int(rank1.sum()), rank1_percent=float(100 * rank1.mean()))
bins = [(0, 0.5), (0.5, 5), (5, 50), (50, 101)]
S["by_human_share"] = {f"{lo}-{min(hi, 100)}%": dict(donors=int(m.sum()), rank1=int(rank1[m].sum()), percent=(float(100 * rank1[m].mean()) if m.any() else None),
                                                      cohorts={s: int((study[m] == s).sum()) for s in sorted(set(study[m]))})
                       for lo, hi in bins for m in [(share >= lo) & (share < hi)]}
S["by_cohort"] = {s: dict(donors=int((study == s).sum()), rank1=int(rank1[study == s].sum()), percent=float(100 * rank1[study == s].mean()),
                          human_share_median=float(np.median(share[study == s])), human_share_max=float(share[study == s].max())) for s in sorted(set(study))}
V = study == "VincentC_2016"
if V.any():
    sv = share[V]; rv = rank1[V]
    S["vincent"] = dict(donors=int(V.sum()), rank1=int(rv.sum()), human_share_median=float(np.median(sv)),
                        human_share_median_rank1=float(np.median(sv[rv])) if rv.any() else None, human_share_median_failures=float(np.median(sv[~rv])) if (~rv).any() else None)
    if rv.any() and (~rv).any():
        u = mannwhitneyu(sv[rv], sv[~rv]); S["vincent"]["mannwhitney_p_share_rank1_vs_failures"] = float(u.pvalue)
    svc = a("service")[V]; ok = ~np.isnan(svc)
    if ok.sum() > 3:
        c = spearmanr(sv[ok], svc[ok]); S["vincent"]["spearman_share_vs_service_similarity"] = dict(rho=float(c.statistic), p=float(c.pvalue), n=int(ok.sum()))
    sq, sp = a("share_q")[V], a("share_p")[V]
    c = spearmanr(sq, sp); S["vincent"]["spearman_share_query_vs_partner"] = dict(rho=float(c.statistic), p=float(c.pvalue))
    ratio = np.maximum(sq, sp) / np.maximum(np.minimum(sq, sp), 1e-3)
    S["vincent"]["pairs_with_a_share_ratio_over_2"] = int((ratio > 2).sum())
    S["vincent"]["pairs_one_sample_below_20_one_above_60"] = int(((np.minimum(sq, sp) < 20) & (np.maximum(sq, sp) > 60)).sum())
    terc = np.percentile(sv, [33.3, 66.7])
    S["vincent"]["rank1_by_share_tercile"] = {nm: dict(donors=int(m.sum()), rank1=int(rv[m].sum())) for nm, m in
                                              (("lowest third", sv <= terc[0]), ("middle third", (sv > terc[0]) & (sv <= terc[1])), ("highest third", sv > terc[1]))}
    S["vincent"]["similarity_median"] = dict(all_hashes=float(np.nanmedian(a("J")[V])), exact_human_removed=float(np.nanmedian(a("J_non")[V])),
                                             exact_and_near_removed=float(np.nanmedian(a("J_non_near")[V])) if not np.all(np.isnan(a("J_non_near")[V])) else None)
    S["vincent"]["days_between_samples_median"] = float(np.nanmedian(a("days")[V]))
json.dump(S, open(f"{OUT}_summary.json", "w"), indent=1)
with open(f"{OUT}_table.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(D[0])); w.writeheader(); w.writerows(D)
print(json.dumps(S, indent=1))
