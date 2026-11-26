"""Builds DATA_V2/results.jsonl: for every subject, the answer of the first queries (results.jsonl, then results_pass2.jsonl, first wins, as in
rows_with_partner.py), replaced by the result of the re-query with the largest indexed run when there is one. The results_pass2.jsonl written in DATA_V2 is empty."""
import json, os, sys
SP = os.environ.get("WORK", "local/work"); A = sys.argv[1]; D = sys.argv[2]            # A = folder of the first queries, D = folder written (adultes_data_v2)
res = {}
for fn in ("results.jsonl", "results_pass2.jsonl"):
    for l in open(f"{A}/{fn}"):
        try:
            q = json.loads(l)
        except Exception:
            continue
        if "results" in q and (q["study"], q["subject"]) not in res:
            res[(q["study"], q["subject"])] = q
new = {}
for l in open(f"{SP}/results_largest.jsonl"):
    try:
        q = json.loads(l)
    except Exception:
        continue
    if "results" in q:
        new[(q["study"], q["subject"])] = q
for k, q in new.items():
    res[k] = q
os.makedirs(D, exist_ok=True)
with open(f"{D}/results.jsonl", "w") as f:
    for q in res.values():
        f.write(json.dumps(q) + "\n")
open(f"{D}/results_pass2.jsonl", "w").close()
print(len(res), "subjects;", len(new), "with the largest-run query")
