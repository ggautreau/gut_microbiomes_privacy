#!/usr/bin/env python3
"""HMP IBD: replaces the visit codes of the rows by the sample identifiers of curatedMetagenomicData, so that the rows join the CroCoDeEL tables and the sample metadata of the cohort.

    python3 -I hmp_curated_ids.py ROWS_V4.json ADULTES_DIR MATRICES_DIR OUT.json

ROWS_V4.json   output of make_rows_v4.py. For HMP IBD the samples of a donor are visit codes (build_hmp2_v2.py, from the ENA), because the curated metadata of this cohort hold no run accession.
ADULTES_DIR    runs_hmp2.json (run -> study, subject, visit), ibdmdb_runs_ena.tsv (read_count of each run), pairs_hmp2.json (the run of the first and last sample)
MATRICES_DIR   HMP_2019_ibdmdb.metadata.tsv (sample_id, subject_id, number_reads)

The link: a curated sample and a run of one subject are the same library when number_reads = 2 x read_count of the ENA. A visit with two curated samples is resolved by the run when the run is known (first and last sample);
a sample with no link keeps its visit code, which no curated table holds. Every other field and every other cohort are copied unchanged. extract_lists.py and every step that recognises the runs of a donor read ROWS_V4.json, with the visit codes.
"""
import collections, csv, json, sys

ROWS, ADU, MAT, OUT = sys.argv[1:5]
ST = "HMP_2019_ibdmdb"
by = collections.defaultdict(list)
for m in csv.DictReader(open(f"{MAT}/{ST}.metadata.tsv"), delimiter="\t"):
    by[(m["subject_id"], int(float(m["number_reads"])))].append(m["sample_id"])
ena = {e["run_accession"]: e for e in csv.DictReader(open(f"{ADU}/ibdmdb_runs_ena.tsv"), delimiter="\t")}
of_run = lambda subj, run: by.get((subj, 2 * int(ena[run]["read_count"])), []) if run in ena and ena[run].get("read_count") else []
visit = collections.defaultdict(set)
for run, (st, subj, vis) in json.load(open(f"{ADU}/runs_hmp2.json")).items():
    if st == ST:
        visit[(subj, vis)].update(of_run(subj, run))
pairs = {p["subject"]: p for p in json.load(open(f"{ADU}/pairs_hmp2.json")) if p["study"] == ST}
rows = json.load(open(ROWS)); n_out = 0; two = 0; left = []


def one(subj, vis, run=None):
    """The curated sample of a visit: the one of the run when the run is known and links, else the visit's only one (or the first of two by name), else the visit code."""
    global two
    c = sorted(visit.get((subj, vis), ()))
    if run:
        hit = [x for a in run.split(";") for x in of_run(subj, a) if x in c]
        if len(set(hit)) == 1:
            return hit[0]
    two += len(c) > 1
    return c[0] if c else vis


for r in rows:
    if r["study"] != ST:
        continue
    s = r["subject"]; p = pairs[s]
    r["v1s"], r["v2s"] = one(s, r["v1s"], p["v1_acc"]), one(s, r["v2s"], p["v2_acc"])
    for key in ("best_any", "best_far"):
        if r.get(key):
            r[key] = one(s, r[key])
    r["eligible"] = sorted({x for e in r["eligible"] for x in visit.get((s, e), ())})
    n_out += 1
    left += [k for k in ("v1s", "v2s", "best_any", "best_far") if r.get(k) and r[k] not in {x for v in by.values() for x in v}]
assert len(left) <= 1 and "v2s" not in left, left
json.dump(rows, open(OUT, "w"))
print(n_out, "donors of", ST, "; samples without a curated identifier:", left, "; samples resolved by name among two:", two)
