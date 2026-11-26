#!/usr/bin/env python3
"""The 183 donors of the 4 cohorts left out (Pehrsson, Nielsen, Li, Vincent), for the samples studied of Supplementary Note 1, Supplementary Fig. 4 and Supplementary Table 8.

    python3 -I donors_left_out.py ROWS_V2_FULL.json OUT.json

ROWS_V2_FULL.json  the 831 donors whose two samples have a run in the index (rows_with_partner.py), with the query sample (v2s), the first sample (v1s) and the most similar other sample (best_any)
OUT.json           the donors of the 4 cohorts left out that are adults with a DNA query, in the format of the donors of the 7 cohorts kept, with only the fields that supplementary_figure4.py reads
                   (written next to the rows of the donors kept, since supplementary_figure4.py reads it from the folder of its ROWS.json; make_inputs.py writes data/donors_left_out.json from it)
"""
import json, sys
LEFT_OUT = ["PehrssonE_2016", "NielsenHB_2014", "LiJ_2014", "VincentC_2016"]
rows = [r for r in json.load(open(sys.argv[1])) if r["study"] in LEFT_OUT and r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]
out = [{k: r[k] for k in ("study", "subject", "v1s", "v2s", "best_any", "age_class", "q_lib")} for r in rows]
assert len(out) == 183, len(out)
json.dump(out, open(sys.argv[2], "w"))
print(len(out), "donors of the 4 cohorts left out")
