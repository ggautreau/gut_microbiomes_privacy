#!/usr/bin/env python3
"""Writes the per-donor table of Figure 1 (data/adult_ranks.csv) from the rows of the donors kept and the size of their queried runs.

    python3 -I make_adult_ranks.py ROWS_V4.json QUERY_GBP.json OUT.csv

ROWS_V4.json  rows_v4.json of exclude_donors.py (the 483 donors of the analysis);  QUERY_GBP.json  'study|subject' -> bases of the queried run in Gbp (ENA base_count of the run queried).
One row per donor, in the order of ROWS_V4.json. Columns: study; age_class; q_lib; rank_any_narrow (rank of the donor's other sample among the 291,881 runs of the HGM set, empty = not among
the 5,000 nearest runs); far_possible (1 = the donor has another sample at least 30 days from the query sample); rank_far_narrow (the same rank when only those samples can be the other sample); days (between the first and
the last sample); rank_first_narrow (the rank when the first sample alone is the other sample); query_gbp (rounded to 0.001 Gbp).
"""
import csv, json, sys

rows = json.load(open(sys.argv[1])); gbp = json.load(open(sys.argv[2]))
cell = lambda v: "" if v is None else v
with open(sys.argv[3], "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["study", "age_class", "q_lib", "rank_any_narrow", "far_possible", "rank_far_narrow", "days", "rank_first_narrow", "query_gbp"])
    for r in rows:
        w.writerow([r["study"], r["age_class"], r["q_lib"], cell(r["rank_any_narrow"]), 1 if r["far_possible"] else 0, cell(r["rank_far_narrow"]), "%g" % r["days"], cell(r["rank_v1_narrow"]), repr(round(gbp[f"{r['study']}|{r['subject']}"], 3))])
print(len(rows), "rows written")
