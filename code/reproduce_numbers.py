"""Every number of main.tex that depends on the donor set, recomputed from the rows table (definitions copied from the build scripts and from the current text).
usage: python3 -I reproduce_numbers.py ROWS.json GBP.json EXTRAPOLATION.json FIG2_STATS.json QUADRATIC.json OUT.json
Run on the 483 donors of the analysis (data/donors.json) it gives the numbers of the main text and of the Methods."""
import json, math, sys
import numpy as np
from scipy.special import gammaln
from scipy.stats import hypergeom

ROWS, GBP, EXT, FIG2, QUAD, OUT = sys.argv[1:7]
rows = [r for r in json.load(open(ROWS)) if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]
gbp = json.load(open(GBP)); ext = json.load(open(EXT)); f2 = json.load(open(FIG2)); qd = json.load(open(QUAD))
M = 291_881
RESAMPLES = 100_000


def macro_boot(seed, by):
    """The mean of the cohort rates, on resamples of the cohorts and then of the donors inside each drawn cohort. A cohort drawn m times adds m binomials: only how many of its resampled donors are recovered matters."""
    rng = np.random.default_rng(seed); C = len(by)
    m = rng.multinomial(C, [1 / C] * C, size=RESAMPLES); tot = np.zeros(RESAMPLES)
    for j, v in enumerate(by.values()):
        tot += rng.binomial(m[:, j] * len(v), np.mean(v)) / len(v)
    return 100 * tot / C


N = len(rows)
pc = lambda k, n: 100 * k / n
res = dict(donors=N, cohorts=len({r["study"] for r in rows}))
nar = [r["rank_any_narrow"] for r in rows]
res["rank1_narrow"] = sum(1 for x in nar if x == 1); res["rank1_narrow_pct"] = pc(res["rank1_narrow"], N)
res["top10_narrow"] = sum(1 for x in nar if x is not None and x <= 10); res["top10_narrow_pct"] = pc(res["top10_narrow"], N)
# the same rank when only the donor's first sample counts as the other sample (the measure of Fig. 1E), against the most similar of all the other samples above
res["rank1_first_sample"] = sum(1 for r in rows if r["rank_v1_narrow"] == 1); res["rank1_first_sample_pct"] = pc(res["rank1_first_sample"], N)
res["rank1_broad"] = sum(1 for r in rows if r["rank_any"] == 1); res["rank1_broad_pct"] = pc(res["rank1_broad"], N)
res["absent"] = sum(1 for r in rows if r["rank_any"] is None); res["absent_pct"] = pc(res["absent"], N)
# pools of N runs (hypergeometric), as figure1_brief_communication.py
K = np.array([(x - 1) if x is not None else M for x in nar], float)
found = np.array([x is not None for x in nar]); kk = np.array([(x - 1) if x is not None else -1 for x in nar])


def p_first(k, n):
    n1 = n - 1; ok = (M - k) >= n1; out = np.zeros(len(k))
    out[ok] = np.exp(gammaln(M - k[ok] + 1) - gammaln(M - k[ok] - n1 + 1) - gammaln(M + 1) + gammaln(M - n1 + 1)); return out


top10 = lambda n: float(hypergeom.cdf(9, M, kk[found], int(n) - 1).sum() / N)
Ns = np.unique(np.round(np.logspace(np.log10(2), np.log10(M), 60)).astype(int)); n1e4 = int(Ns[np.argmin(np.abs(np.log10(Ns) - 4))])
res["pool"] = dict(first_N2=100 * float(p_first(K, 2).mean()), first_N1e2=100 * float(p_first(K, 100).mean()), first_N1e3=100 * float(p_first(K, 1000).mean()), first_N1e4=100 * float(p_first(K, n1e4).mean()), first_N10000=100 * float(p_first(K, 10_000).mean()), first_N1e5=100 * float(p_first(K, 100_000).mean()), first_NM=100 * float(p_first(K, M).mean()),
                   top10_N1e4=100 * top10(n1e4), top10_NM=100 * top10(M), N_near_1e4=n1e4, top10_NM_ci_cohorts=ext["top10_ci_cohorts"])
# cohorts: macro mean and interval (seed 1: cohorts, then donors, 100,000 resamples), rates on rank_any == 1 (all indexed runs)
by = {}
for r in rows:
    by.setdefault(r["study"], []).append(1 if r["rank_any"] == 1 else 0)
names = list(by); rates = {k: 100 * np.mean(v) for k, v in by.items()}
boot = macro_boot(1, by)
big = max(by, key=lambda k: len(by[k])); rest = [x for k, v in by.items() if k != big for x in v]
byn = {}
for r in rows:
    byn.setdefault(r["study"], []).append(1 if r["rank_any_narrow"] == 1 else 0)
restn = [x for k, v in byn.items() if k != big for x in v]
res["cohorts_detail"] = dict(rates_broad={k: (len(by[k]), 100 * np.mean(by[k])) for k in by}, rates_narrow={k: (len(byn[k]), 100 * np.mean(byn[k])) for k in byn},
                             macro_mean=float(np.mean(list(rates.values()))), macro_ci=[float(x) for x in np.percentile(boot, [2.5, 97.5])],
                             min=(min(rates, key=rates.get), min(rates.values())), max=(max(rates, key=rates.get), max(rates.values())),
                             largest=big, largest_n=len(by[big]), largest_share=pc(len(by[big]), N), without_largest_n=len(rest), without_largest_broad=pc(sum(rest), len(rest)), without_largest_narrow=pc(sum(restn), len(restn)))
# the same on the HGM set of 291,881 runs (rank_any_narrow), the measure of Fig. 1D and of the main text
byS = {}
for r in rows:
    byS.setdefault(r["study"], []).append(1 if r["rank_any_narrow"] == 1 else 0)
ratesS = {k: 100 * np.mean(v) for k, v in byS.items()}; bootS = macro_boot(1, byS)
restS = [x for k, v in byS.items() if k != big for x in v]
res["cohorts_detail_stock"] = dict(rates={k: (len(byS[k]), 100 * np.mean(byS[k])) for k in byS}, macro_mean=float(np.mean(list(ratesS.values()))), macro_ci=[float(x) for x in np.percentile(bootS, [2.5, 97.5])],
                                   min=(min(ratesS, key=ratesS.get), min(ratesS.values())), max=(max(ratesS, key=ratesS.get), max(ratesS.values())), without_largest_n=len(restS), without_largest=pc(sum(restS), len(restS)))
# 30-day subset (Methods, Time between samples)
far = [r for r in rows if r["far_possible"]]
res["far"] = dict(n=len(far), far_rank1_narrow=pc(sum(1 for r in far if r["rank_far_narrow"] == 1), len(far)), same_donors_any_partner=pc(sum(1 for r in far if r["rank_any_narrow"] == 1), len(far)),
                  all_donors=res["rank1_narrow_pct"])
days = np.array([float(r["days"]) for r in rows]); first = np.array([r["rank_v1_narrow"] == 1 for r in rows])
bins = [("<=30", 0, 30), ("31-90", 31, 90), ("91-180", 91, 180), ("181-391", 181, 391), (">391", 392, 1e9)]
res["first_sample_by_interval"] = {lab: (int(((days >= lo) & (days <= hi)).sum()), int(first[(days >= lo) & (days <= hi)].sum()), pc(first[(days >= lo) & (days <= hi)].sum(), max(1, ((days >= lo) & (days <= hi)).sum()))) for lab, lo, hi in bins}
inb = {}
for r_, d_ in zip(rows, days):                                            # which cohorts make the 181-391 day class of Fig. 1E
    if 181 <= d_ <= 391: inb[r_["study"]] = inb.get(r_["study"], 0) + 1
res["interval_181_391_by_cohort"] = dict(sorted(inb.items()))
res["days"] = dict(median=float(np.median(days)), q25=float(np.percentile(days, 25)), q75=float(np.percentile(days, 75)), max=float(days.max()))
# depth
g = np.array([gbp[f"{r['study']}|{r['subject']}"] for r in rows]); okb = np.array([r["rank_any"] == 1 for r in rows])
res["depth"] = dict(median_gbp=float(np.median(g)), below_half=int((g < 0.5).sum()), below_half_linked=int(okb[g < 0.5].sum()), above_half=int((g >= 0.5).sum()), above_half_linked=int(okb[g >= 0.5].sum()))
# CroCoDeEL (stats of supplementary_figure4.py) and the 30-day donors with an event with the same donor
A = f2["A"]
res["crocodeel"] = dict(none_all=A["none/all"], none_30d=A["none/30d"], same_all=A["same/all"], same_30d=A["same/30d"], other_all=A["other/all"], other_30d=A["other/30d"],
                        all_donors_broad_pct=res["rank1_broad_pct"], events=f2["n_events"], same=f2["n_same"], other=f2["n_other"])
# extrapolation
ho = ext["holdout"]
res["extrapolation"] = dict(world_first=qd["first"]["8000000000"], world_top10=qd["top10"]["8000000000"], line_overestimate_max=max(max(0.0, ho[f]["error"]["line"]) for f in ho),
                            line_errors={f: ho[f]["error"]["line"] for f in ho}, residual_first_at_M=qd["residual_first_at_M"], measured_first=qd["measured_first_at_M"], measured_top10=qd["measured_top10_at_M"])
json.dump(res, open(OUT, "w"), indent=1, default=float)


def show(d, ind=0):
    for k, v in d.items():
        if isinstance(v, dict): print(" " * ind + f"{k}:"); show(v, ind + 2)
        else: print(" " * ind + f"{k}: {v if not isinstance(v, float) else round(v, 2)}")


show(res)
