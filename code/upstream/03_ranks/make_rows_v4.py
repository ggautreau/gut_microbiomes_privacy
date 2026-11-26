#!/usr/bin/env python3
"""The 497 donors of the cohorts kept, before the two rules: the adult donors with a DNA query of the cohorts kept, from the rows of the 831 donors whose two samples are in the index.

    python3 -I make_rows_v4.py ROWS_V2_FULL.json OUT.json [--keep-vincent]

ROWS_V2_FULL.json  output of rows_with_partner.py on the final queries (831 donors whose first and last samples each have a run in the index; 12 cohorts).
A donor is kept when age_class is adult or adult_by_design (age of 18 or more in the curated metadata; Mehta 2018 is a cohort of adults by design) and the query is DNA (not a metatranscriptome),
and its cohort is not one of the four left out of the analysis (Supplementary Note 1): Pehrsson 2016, Nielsen 2014, Li 2014, Vincent 2016. The order of the rows is that of ROWS_V2_FULL.json.
Result: 497 donors in 7 cohorts. With --keep-vincent the 61 donors of Vincent are kept (558 donors in 8 cohorts, the table of the human-DNA check, all558_join.py).
"""
import json, sys

LEFT_OUT = {"PehrssonE_2016", "NielsenHB_2014", "LiJ_2014", "VincentC_2016"}
if "--keep-vincent" in sys.argv[3:]:
    LEFT_OUT.discard("VincentC_2016")
rows = json.load(open(sys.argv[1]))
kept = [r for r in rows if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA" and r["study"] not in LEFT_OUT]
json.dump(kept, open(sys.argv[2], "w"))
print(len(rows), "donors whose two samples are in the index;", len(kept), "kept in", len({r["study"] for r in kept}), "cohorts")
