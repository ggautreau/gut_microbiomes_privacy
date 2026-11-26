#!/usr/bin/env python3
"""Runs of the Mehta cohort (PRJNA354235) that the curated metadata do not name, identified by the alias of their ENA sample.

In this project every stool sample was sequenced as a DNA and as an RNA library, and the two libraries sit in two BioSamples whose aliases carry the same subject code and the same sample name
(<subject>_<sample>.dna..., <subject>_<sample>.rna...). The curated metadata name one BioSample per sample, so the runs of the other one would count as runs of another person.

    python3 -I mehta_siblings.py fetch OUT.tsv [PROJECT]                # the run table of a project, PRJNA354235 by default (ENA portal; the request holds the project accession and nothing else)
    python3 -I mehta_siblings.py map DATA_DIR ENA.tsv OUT.json          # run -> [study, subject, sample] for the runs of the project that the curated BioSamples do not hold

DATA_DIR is the folder of rows_with_partner.py (curated_sample_biosample.json, matrices/MehtaRS_2018.metadata.tsv). A run is given the identity of the curated sample whose subject code and sample name
its alias carries; a run that names a sample outside the curated metadata stays unmapped (counted as another person). rows_with_partner.py, extract_lists.py and stock_sensitivity.py read OUT.json
from DATA_DIR/mehta_sibling_runs.json when it exists.
"""
import csv, json, re, sys, urllib.request

STUDY, PROJECT = "MehtaRS_2018", "PRJNA354235"
FIELDS = "run_accession,sample_accession,sample_alias,sample_title,library_name,library_strategy,library_source,read_count,base_count"
ALIAS = re.compile(r"^(\d+)_([A-Za-z0-9]+)\.(dna|rna)\.")
csv.field_size_limit(10**9)


def fetch(out, project=PROJECT):
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={project}&result=read_run&fields={FIELDS}&format=tsv&limit=0"
    open(out, "wb").write(urllib.request.urlopen(url, timeout=120).read())


def mapping(data, ena):
    cs = json.load(open(f"{data}/curated_sample_biosample.json"))
    subject = {r["sample_id"]: r["subject_id"] for r in csv.DictReader(open(f"{data}/matrices/{STUDY}.metadata.tsv"), delimiter="\t") if r["body_site"] == "stool" and r["subject_id"]}
    bs2s = {}
    for k, bs in cs.items():
        st, sid = k.split("|", 1)
        if st == STUDY and sid in subject:
            bs2s[bs] = (STUDY, subject[sid], sid)
    rows = list(csv.DictReader(open(ena), delimiter="\t", quoting=csv.QUOTE_NONE))
    key = lambda r: ALIAS.match(r["sample_alias"]).group(1, 2)
    ident = {}
    for r in rows:
        assert ALIAS.match(r["sample_alias"]), "alias outside the pattern"
        if r["sample_accession"] in bs2s:
            assert ident.setdefault(key(r), bs2s[r["sample_accession"]]) == bs2s[r["sample_accession"]], "two curated samples for one subject code and sample name"
    out = {r["run_accession"]: list(ident[key(r)]) for r in rows if r["sample_accession"] not in bs2s and key(r) in ident}
    return dict(sorted(out.items())), len(rows), sum(r["sample_accession"] not in bs2s for r in rows)


if __name__ == "__main__":
    if sys.argv[1] == "fetch":
        fetch(sys.argv[2], *sys.argv[3:4])
    else:
        out, n, n_out = mapping(sys.argv[2], sys.argv[3])
        json.dump(out, open(sys.argv[4], "w"), indent=0, sort_keys=True)
        print(f"{n} runs in the project, {n_out} outside the curated BioSamples, {len(out)} of them identified by their alias")
