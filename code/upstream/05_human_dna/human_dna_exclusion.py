#!/usr/bin/env python3
"""The aggregates behind the exclusion of Vincent 2016 (Supplementary Note 1 and Supplementary Fig. 3): human DNA in the samples of each cohort and share of donors recovered.

    python3 -I human_dna_exclusion.py JOIN_FINAL_TABLE.csv POOL_FINAL_RANKS.csv JOIN_FINAL_SUMMARY.json POOL_FINAL_SUMMARY.json OUT.json

JOIN_FINAL_TABLE.csv     all558_join.py: one line per donor of the 558 (the 497 donors of the 7 cohorts kept and the 61 of Vincent), with `share_q` and `share_p`, the share (%) of the hashes of the sketch
                         of the queried run and of the donor's other sample that are found in a human genome (T2T-CHM13v2.0 or GRCh38.p14), and `rank1` (the other sample is the nearest of the 291,881 runs of the HGM set).
POOL_FINAL_RANKS.csv     all558_pool.py: rank of the donor's other sample among the 1,116 runs of the 558 donors, with all hashes and without those that are found in a human genome.
JOIN_FINAL_SUMMARY.json  all558_join.py: rank 1 by human share (the four strata).   POOL_FINAL_SUMMARY.json  all558_pool.py.
Only aggregates are written. The share of a donor is that of the more human of its two samples; 'discordant' = one sample under 20% human hashes and the other over 60%.
"""
import csv, json, statistics as S, sys

JOIN, RANKS, JSUM, PSUM, OUT = sys.argv[1:6]
J = list(csv.DictReader(open(JOIN))); R = list(csv.DictReader(open(RANKS))); js = json.load(open(JSUM)); ps = json.load(open(PSUM))
f = lambda r, k: float(r[k])
share = lambda r: max(f(r, "share_q"), f(r, "share_p"))
V = [r for r in J if r["study"] == "VincentC_2016"]; O = [r for r in J if r["study"] != "VincentC_2016"]
disc = [r for r in V if min(f(r, "share_q"), f(r, "share_p")) < 20 and max(f(r, "share_q"), f(r, "share_p")) > 60]
RV = [r for r in R if r["study"] == "VincentC_2016"]; RO = [r for r in R if r["study"] != "VincentC_2016"]
d = ps["definitions"]; ALL, REM = d["all hashes"], d["exact human removed"]
out = dict(
    vincent=dict(n=len(V), linked_stock=sum(r["rank1"] == "True" for r in V), share_median=S.median(share(r) for r in V), share_max=max(share(r) for r in V),
                 n_over5=sum(share(r) >= 5 for r in V), n_over50=sum(share(r) >= 50 for r in V), discordant=len(disc), discordant_linked=sum(r["rank1"] == "True" for r in disc),
                 pool_all=sum(r["all hashes|rank"] == "1" for r in RV), pool_removed=sum(r["exact human removed|rank"] == "1" for r in RV)),
    others=dict(n=len(O), share_median=S.median(share(r) for r in O), n_over5=sum(share(r) >= 5 for r in O), n_over5_linked=sum(share(r) >= 5 and r["rank1"] == "True" for r in O),
                linked_stock=sum(r["rank1"] == "True" for r in O)),
    strata=[dict(label=k.replace("-", " to ").replace("%", "").replace("0 to 0.5", "below 0.5").replace("50 to 100", "50 or more"), donors=v["donors"], linked=v["rank1"]) for k, v in js["by_human_share"].items()],
    pool=dict(rich_donors=ps["pool"]["donors_with_human_share_5_percent_or_more"], rich_all=ALL["top1_human_share_5_or_more"], rich_removed=REM["top1_human_share_5_or_more"],
              poor_all=ALL["top1_human_share_below_5"], poor_removed=REM["top1_human_share_below_5"], pool_runs=ps["pool"]["runs"], pool_donors=ps["pool"]["donors"],
              others_all=sum(r["all hashes|rank"] == "1" for r in RO), others_removed=sum(r["exact human removed|rank"] == "1" for r in RO), others_n=len(RO)))
json.dump(out, open(OUT, "w"), indent=1)
print("Vincent:", out["vincent"]); print("others:", out["others"])
