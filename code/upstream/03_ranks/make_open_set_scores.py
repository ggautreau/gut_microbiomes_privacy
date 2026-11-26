#!/usr/bin/env python3
"""Writes data/open_set_scores.csv, the table of panel H of Figure 1: for each donor, the similarity of the donor's other sample and that of the nearest run of other people.

    python3 -I make_open_set_scores.py LISTS.npz OUT.csv          # LISTS.npz = lists_v4.npz of exclude_donors.py (for data/open_set_scores.csv), or lists_v4_497.npz of extract_lists.py (the 497 donors, for costea_same_people.py)

One row per donor, in the order of the analysis. Columns: study; `partner` (the similarity of the donor's other sample to the query, 0 = not among the 5,000 nearest runs);
`best_other` (the similarity of the nearest run of other people among the 291,881 runs of the HGM set). The search of a person absent from the index is the same query with the donor's samples removed:
its top hit is the `best_other` run.
"""
import csv, sys
import numpy as np

z = np.load(sys.argv[1], allow_pickle=True)
study, partner, impostors = z["study"], z["partner"], z["imp_nar"]                 # read each array once (every access to z[...] reads the whole array again)
n = len(study)
with open(sys.argv[2], "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["study", "partner", "best_other"])
    for i in range(n):
        w.writerow([study[i], f"{float(partner[i]):.6f}", f"{float(impostors[i][0]):.6f}"])
print(n, "rows written")
