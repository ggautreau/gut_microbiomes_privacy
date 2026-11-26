#!/usr/bin/env python3
"""Who each cohort is made of: the age label of the longitudinal cohorts of curatedMetagenomicData (Supplementary Note 1, Supplementary Tables 1 and 2).

    python3 -I cohort_age.py MATRICES_DIR OUT.csv

MATRICES_DIR  the 93 metadata tables <study>.metadata.tsv of curatedMetagenomicData (extract_all.R).
A longitudinal cohort is a study with at least 5 subjects that have two or more stool samples (30 of the 93 studies, 2,038 subjects). Per cohort the subjects counted are those and the age is that of their first
sample. curatedMetagenomicData records age in two columns and fills only one per study: `age` in years, and `infant_age` in days for the birth cohorts. Both are the string "NA" when absent, so a subject
counts as an infant (under 3 years) when either parses as a number and says so. A cohort is labelled only when at least 70% of its subjects carry an age: 'infants' when at least 80% of those are infants,
'adults' when at most 20% are, 'mixed' otherwise, 'not recorded' when the age is missing. Mehta 2018 records no age but all its 928 stool records are male: a cohort of adult men, labelled 'adults' by that field.
Cohorts are written by decreasing number of subjects. Aggregates only: no subject is named.
"""
import collections
import csv
import glob
import os
import sys

src, out = sys.argv[1:3]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):       # "NA", "", None
        return None


def first_samples(study):
    """The earliest stool sample of every subject who has at least two of them."""
    path = f"{src}/{study}.metadata.tsv"
    if not os.path.exists(path):
        return []
    rows = [r for r in csv.DictReader(open(path), delimiter="\t")
            if r.get("body_site") == "stool" and r.get("subject_id")]
    by = collections.defaultdict(list)
    for r in rows:
        by[r["subject_id"]].append(r)
    out_rows = []
    for rs in by.values():
        if len(rs) < 2:
            continue
        rs.sort(key=lambda r: num(r.get("days_from_first_collection")) or 0.0)
        if rs[0]["sample_id"] != rs[-1]["sample_id"]:
            out_rows.append(rs[0])
    return out_rows


COVERAGE = 0.7          # share of a cohort's subjects that must carry an age

# A cohort whose age field is empty but whose population is settled by another field of the same
# metadata. Only one so far: every one of the 928 stool records of Mehta 2018 carries gender "male",
# and the study follows adult men, so leaving it unlabelled would hide a fact the data state.
BY_OTHER_FIELD = {"MehtaRS_2018": ("adults", "all 928 records male, a cohort of adult men")}


def classify(rows):
    """share of subjects known to be infants, and the age label of the cohort.

    A cohort is only labelled when most of its subjects carry an age: the mother-infant
    studies record an age for the infants alone, and labelling those from the half that
    is recorded would call the whole cohort an infant cohort.
    """
    known = infant = 0
    for r in rows:
        years, days = num(r.get("age")), num(r.get("infant_age"))
        if years is None and days is None:
            continue
        known += 1
        if (years is not None and years < 3) or (days is not None and days <= 1095):
            infant += 1
    if not known or known < COVERAGE * len(rows):
        return known, round(100 * infant / known, 1) if known else 0.0, "not recorded"
    share = infant / known
    label = "infants" if share >= 0.8 else ("adults" if share <= 0.2 else "mixed")
    return known, round(100 * share, 1), label


MIN_SUBJECTS = 5
rows_out = []
for path in sorted(glob.glob(f"{src}/*.metadata.tsv")):
    study = os.path.basename(path)[:-len(".metadata.tsv")]
    subjects = first_samples(study)
    if len(subjects) < MIN_SUBJECTS:
        continue
    known, pct_infant, label = classify(subjects)
    note = ""
    if label == "not recorded" and study in BY_OTHER_FIELD:
        label, note = BY_OTHER_FIELD[study]
    rows_out.append([study, len(subjects), known, pct_infant, label, note])
rows_out.sort(key=lambda r: -r[1])                                      # the sort is stable: cohorts with the same number of subjects stay in alphabetical order
for r in rows_out:
    print(f"{r[0]:<26}{r[1]:>5} subjects  age known {r[2]:>4}  infants {r[3]:>5}%  -> {r[4]}" + (f"  ({r[5]})" if r[5] else ""))

with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["dataset", "subjects", "subjects_with_age", "pct_infant", "population", "from_other_field"])
    w.writerows(rows_out)
print("wrote", out, "-", len(rows_out), "cohorts")
