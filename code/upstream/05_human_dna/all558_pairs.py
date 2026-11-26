"""The 558 donor pairs (query run, other sample of the same donor) with the similarity returned by the service, from the same tables as rows_with_partner.py.
usage: python3 -I all558_pairs.py LOCAL_DIR OUT.json"""
import csv, json, os, sys
SP, OUT = sys.argv[1], sys.argv[2]
ADU = SP + "/adultes_data_v2"; MAT = ADU + "/matrices"
runs = lambda s: [a for a in s.split(";") if a not in ("", "NA")]
P = json.load(open(f"{ADU}/pairs.json")) + [p for p in json.load(open(f"{ADU}/pairs_hmp2.json")) if p["study"] == "HMP_2019_ibdmdb"]
subset = {(p["study"], p["subject"]): p for p in P}
res = {}
for fn in ("results.jsonl", "results_pass2.jsonl"):
    for l in open(f"{ADU}/{fn}"):
        try: q = json.loads(l)
        except Exception: continue
        if "results" in q and (q["study"], q["subject"]) not in res: res[(q["study"], q["subject"])] = q
allruns = json.load(open(f"{ADU}/biosample_runs.json")); cs = json.load(open(f"{ADU}/curated_sample_biosample.json"))
meta_subj = {}
for study in {p["study"] for p in P}:
    for r in csv.DictReader(open(f"{MAT}/{study}.metadata.tsv"), delimiter="\t"):
        if r["body_site"] == "stool" and r["subject_id"]: meta_subj[(study, r["sample_id"])] = r["subject_id"]
bs2s = {}
for k, b in cs.items():
    st, sid = k.split("|", 1)
    if (st, sid) in meta_subj: bs2s[b] = (st, meta_subj[(st, sid)], sid)
acc = {run: bs2s[v["biosample"]] for run, v in allruns.items() if v["biosample"] in bs2s}
acc.update({run: tuple(v) for run, v in json.load(open(f"{ADU}/runs_hmp2.json")).items() if v[0] == "HMP_2019_ibdmdb"})
if os.path.exists(f"{ADU}/mehta_sibling_runs.json"):                                         # the other library of a Mehta sample sits in a BioSample that the curated metadata do not name: same mapping as rows_with_partner.py
    acc.update({run: tuple(v) for run, v in json.load(open(f"{ADU}/mehta_sibling_runs.json")).items() if run not in acc})
rows = {(r["study"], r["subject"]): r for r in json.load(open(f"{SP}/rows_v3.json"))}
keep = [k for k, r in rows.items() if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]
out = []
for key in keep:
    q = res[key]; study, subj = key; v2s = subset[key]["v2_sample"]
    same = [(sc, a) for a, sc in q["results"].items() if sc < 0.99 and acc.get(a) and acc[a][0] == study and acc[a][1] == subj and acc[a][2] != v2s]
    partner = max(same)[1] if same else (runs(subset[key]["v1_acc"]) or [None])[0]
    out.append(dict(study=study, subject=subj, q=q["run"], partner=partner, service=q["results"].get(partner), in_results=bool(same)))
json.dump(out, open(OUT, "w"))
ref = json.load(open(SP + "/humancheck/all558/all558_runs.json"))["pairs"]
assert [(o["study"], o["q"], o["partner"]) for o in out] == [tuple(x) for x in ref], "pair list differs from the one used for the download"
sv = [o["service"] for o in out if o["service"] is not None]
print(f"{len(out)} pairs, identical to the download list; service similarity known for {len(sv)} (other sample found among the 5,000 nearest runs), median {sorted(sv)[len(sv) // 2]:.3f}")
