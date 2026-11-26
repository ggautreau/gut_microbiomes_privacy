#!/usr/bin/env python3
"""Raymond 2016: are the donors who took cephalosporins between their samples still recovered? The numbers of Fig. 1G.

    python3 -I raymond_antibiotic.py OUTPUT.json ROWS.json CMD_META.tsv

ROWS.json     rows_v4.json (exclude_donors.py): one row per donor with study, subject, rank_any_narrow (1 + the runs of the HGM set that are more similar than the donor's other sample; null = not among the
              5,000 nearest runs, which counts as a failure) and rank_v1_narrow (the same rank when the first sample alone is the other sample)
CMD_META.tsv  the harmonised sample metadata of curatedMetagenomicData (study_name, subject_id, sample_id, body_site, days_from_first_collection, study_condition)

In Raymond et al. (ISME J 2016) every donor has three stool samples, on days 0, 7 and 90. curatedMetagenomicData labels the samples taken after the start of the antibiotic course 'cephalosporins' and the
others 'control'. A donor is put in the group of the study_condition of the LAST stool sample (the queried one): 'cephalosporins' (a course of cephalosporins between the first and the last sample) or 'control'.
Counted by group: the donors, those whose other sample is first in the HGM set (the measure of panels D, F and G) and those whose first sample alone is first (the measure of panel E).
Output: counts and the pattern of the days of the samples, without names.
"""
import collections, csv, json, sys

OUT, ROWS, CMD = sys.argv[1:4]
STUDY = "RaymondF_2016"
rows = [r for r in json.load(open(ROWS)) if r["study"] == STUDY]
by = collections.defaultdict(list)
for r in csv.DictReader(open(CMD, encoding="utf-8"), delimiter="\t"):
    if r["study_name"] == STUDY and r["body_site"] == "stool":
        by[r["subject_id"]].append(r)
for v in by.values():
    v.sort(key=lambda r: float(r["days_from_first_collection"] or 0))

groups = collections.OrderedDict((g, {"donors": 0, "other_sample_first": 0, "first_sample_first": 0}) for g in ("cephalosporins", "control"))
days = collections.Counter()
for r in rows:
    s = by[r["subject"]]
    g = groups[s[-1]["study_condition"]]
    g["donors"] += 1
    g["other_sample_first"] += r["rank_any_narrow"] == 1
    g["first_sample_first"] += r["rank_v1_narrow"] == 1
    days[",".join(str(int(float(x["days_from_first_collection"]))) for x in s)] += 1
out = {"study": STUDY,
       "group": "study_condition of the last stool sample of the donor in curatedMetagenomicData: cephalosporins (a course between the first and the last sample) or control",
       "donors": len(rows), "groups": groups, "days_of_the_samples": dict(days)}
json.dump(out, open(OUT, "w"), indent=1)
print(f"wrote {OUT}:", "; ".join(f"{g}: {v['other_sample_first']} of {v['donors']} donors with the other sample first, {v['first_sample_first']} with the first sample alone" for g, v in groups.items()))
