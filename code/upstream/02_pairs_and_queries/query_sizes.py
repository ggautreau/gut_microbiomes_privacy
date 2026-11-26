#!/usr/bin/env python3
"""Bases (Gbp) of the run queried for each donor, from the ENA portal API (public; a request holds run accessions only).

    python3 -I query_sizes.py RESULTS.jsonl OUT.json [CACHE.tsv ...]

RESULTS.jsonl  the answers of the service, one line per donor (study, subject, run = the run queried, results): the final queries, with the largest indexed run of each sample (merge_results.py).
CACHE.tsv      optional tables `run<TAB>read_count<TAB>base_count` already downloaded; only the runs missing from them are asked to the ENA, 100 per request
               (https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&fields=run_accession,read_count,base_count).
OUT.json       'study|subject' -> base_count / 1e9 (Gbp), the file query_gbp.json of the figures. The service sketches one run per query, so the size of that run, not of the sample, is what matters.
"""
import json, sys, time, urllib.parse, urllib.request

RES, OUT = sys.argv[1:3]
bases = {}
for path in sys.argv[3:]:
    for l in open(path):
        a = l.rstrip("\n").split("\t")
        if len(a) >= 3 and a[2].isdigit():
            bases[a[0]] = int(a[2])
queried = {}
for l in open(RES):
    try:
        q = json.loads(l)
    except Exception:
        continue
    if "results" in q:
        queried[(q["study"], q["subject"])] = q["run"]
missing = sorted({r for r in queried.values() if r not in bases})
for i in range(0, len(missing), 100):
    batch = missing[i:i + 100]
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode({"result": "read_run", "fields": "run_accession,read_count,base_count", "format": "tsv", "limit": 0,
                                                                                   "query": " OR ".join(f'run_accession="{r}"' for r in batch)})
    for att in range(4):
        try:
            txt = urllib.request.urlopen(url, timeout=120).read().decode(); break
        except Exception:
            time.sleep(3 * (att + 1)); txt = ""
    for l in txt.splitlines()[1:]:
        a = l.split("\t")
        if len(a) >= 3 and a[2].isdigit():
            bases[a[0]] = int(a[2])
    time.sleep(0.5)
still = [r for r in queried.values() if r not in bases]
assert not still, f"{len(still)} queried runs without a size, e.g. {still[:3]}"
json.dump({f"{st}|{sj}": bases[r] / 1e9 for (st, sj), r in queried.items()}, open(OUT, "w"))
print(len(queried), "donors;", len(missing), "runs asked to the ENA")
