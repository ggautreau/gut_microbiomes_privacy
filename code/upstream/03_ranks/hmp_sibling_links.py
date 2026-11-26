#!/usr/bin/env python3
"""Runs of the two other projects of the HMP consortium that its public metadata table gives to a donor of HMP IBD (same person, other deposit).

    python3 -I hmp_sibling_links.py ROWS_V2_FULL.json HMP2_METADATA.csv ENA_PRJNA389280_RUNS.tsv ENA_STOCK_CANDIDATES.tsv OUT.json

ROWS_V2_FULL.json         the donors whose two samples have a run in the index (rows_with_partner.py): the subjects of HMP IBD
HMP2_METADATA.csv         the public metadata table of the consortium (hmp2_metadata_2018-08-20.csv, ibdmdb.org, downloads)
ENA_PRJNA389280_RUNS.tsv  run table of PRJNA389280 at the ENA (mehta_siblings.py fetch); ENA_STOCK_CANDIDATES.tsv the ENA records of the runs considered for the HGM set (build_stock.py)
OUT.json                  run accession -> [study, subject]: the run belongs to that donor; rows_with_partner.py and extract_lists.py (option --hmp) count it neither as
                          the donor's other sample nor as a run of another person
A run of PRJNA389280 is matched by its sample alias (the column Project of the table), a run of PRJNA306874 by its sample title SM-XXXXX (the column External ID without
its first letter); a run is kept only when the table gives it exactly one participant.
"""
import collections, csv, json, re, sys

ROWS, META, ENA389, CAND, OUT = sys.argv[1:6]
trows = list(csv.reader(open(META, encoding="utf-8", errors="replace")))
h = {c: i for i, c in enumerate(trows[0])}
proj, ext = collections.defaultdict(set), collections.defaultdict(set)
for r in trows[1:]:
    part, coll = r[h["Participant ID"]], r[h["site_sub_coll"]]
    proj[r[h["Project"]]].add((part, coll))
    m = re.match(r"^[A-Z](SM[0-9A-Z]{5})", r[h["External ID"]])
    if m:
        ext[m.group(1)].add((part, coll))
sib = {}
for r in csv.DictReader(open(ENA389), delimiter="\t"):
    sib[r["run_accession"]] = ("PRJNA389280", r["sample_alias"])
for r in csv.DictReader(open(CAND), delimiter="\t"):
    if r["study_accession"] == "PRJNA306874":
        sib[r["run_accession"]] = ("PRJNA306874", r["sample_title"])
    elif r["study_accession"] == "PRJNA389280" and r["run_accession"] not in sib:
        sib[r["run_accession"]] = ("PRJNA389280", r["sample_title"])
hmp_donors = {r["subject"] for r in json.load(open(ROWS)) if r["study"] == "HMP_2019_ibdmdb"}
links, stat = {}, collections.Counter()
for run, (pj, name) in sib.items():
    if pj == "PRJNA389280":
        cand = proj.get(name, set())
    else:
        m = re.match(r"^SM-([0-9A-Z]{5})$", name)
        cand = ext.get("SM" + m.group(1), set()) if m else set()
    parts = {a for a, _ in cand}
    if len(parts) != 1:
        stat["not in the table" if not parts else "several participants"] += 1
        continue
    a = next(iter(parts))
    if a in hmp_donors:
        links[run] = ["HMP_2019_ibdmdb", a]
        stat["given to a donor"] += 1
    else:
        stat["given to another participant"] += 1
json.dump(links, open(OUT, "w"))
print(json.dumps({"runs of the two projects": len(sib), "links": len(links), "donors": len({v[1] for v in links.values()}), **stat}))
