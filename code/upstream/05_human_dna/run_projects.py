#!/usr/bin/env python3
"""BioProject (the ENA study accession) of every run of the human-DNA check, from the public ENA portal API: one request per 100 run accessions, and a request holds run accessions only.

    python3 -I run_projects.py PAIRS.json OUT.json

PAIRS.json  the 558 pairs written by all558_pairs.py (query run `q`, run of the other sample `partner`)
OUT.json    run accession -> study accession (PRJ...), sorted by run; every run of the pairs must be answered (the script stops otherwise)
"""
import json, sys, time, urllib.parse, urllib.request

PAIRS, OUT = sys.argv[1:3]
runs = sorted({r for p in json.load(open(PAIRS)) for r in (p["q"], p["partner"]) if r})
out = {}
for i in range(0, len(runs), 100):
    batch = runs[i:i + 100]
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode({"result": "read_run", "fields": "run_accession,study_accession", "format": "tsv", "limit": 0,
                                                                                   "query": " OR ".join(f'run_accession="{r}"' for r in batch)})
    txt = ""
    for att in range(4):
        try:
            txt = urllib.request.urlopen(url, timeout=120).read().decode(); break
        except Exception:
            time.sleep(3 * (att + 1))
    for l in txt.splitlines()[1:]:
        a = l.split("\t")
        if len(a) >= 2 and a[1].startswith("PRJ"):
            out[a[0]] = a[1]
    time.sleep(0.5)
missing = [r for r in runs if r not in out]
assert not missing, f"{len(missing)} runs without a BioProject, e.g. {missing[:3]}"
json.dump({r: out[r] for r in runs}, open(OUT, "w"), indent=0)
print(len(runs), "runs;", len(set(out.values())), "BioProjects")
