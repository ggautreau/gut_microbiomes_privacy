#!/usr/bin/env python3
"""Why Nielsen and Li are left out: the two cohorts carry the same subject identifiers for 44 patients (Supplementary Note 1). Writes the aggregate numbers (no donor identifier) used by build_supplementary.py.

    python3 -I nielsen_li_overlap.py DATA ROWS.json OUT.json

DATA  folder of rows_with_partner.py (results.jsonl, biosample_runs.json, curated_sample_biosample.json, logan_gut_wgs_narrow_ids.txt.gz); the metadata tables are read from local/crocodeel_valles/matrices;
ROWS.json  the rows of all adult donors including Nielsen and Li (output of rows_with_partner.py).
Rank rule as in rows_with_partner.py: rank = 1 + the runs of the HGM set (the 291,881 runs taken as human gut metagenomes) that are not the donor's and are more similar than the donor's other sample.
'Merged': a sample of the paired cohort under the same subject identifier counts as the donor's own sample (a possible other sample) and not as a run of other people.
CroCoDeEL (crocodeel_valles/results) and the species tables give the contamination numbers of Nielsen.
"""
import collections, csv, gzip, json, os, sys
import numpy as np

DATA, ROWS, OUT = sys.argv[1:4]
CRO = "local/crocodeel_valles"
DUP = 0.99
NIEL, LI, LEC = "NielsenHB_2014", "LiJ_2014", "LeChatelierE_2013"
PAIR = {NIEL: LI, LI: NIEL}


def metadata(st):
    return list(csv.DictReader(open(f"{CRO}/matrices/{st}.metadata.tsv"), delimiter="\t"))


subjects = {st: {m["subject_id"] for m in metadata(st) if m["subject_id"] and m["body_site"] == "stool"} for st in (NIEL, LI, LEC)}
out = dict(shared_nielsen_li=len(subjects[NIEL] & subjects[LI]), shared_nielsen_lechatelier=len(subjects[NIEL] & subjects[LEC]), subjects_nielsen=len(subjects[NIEL]), subjects_li=len(subjects[LI]))

# identity of the runs of the two cohorts: run -> (study, subject, sample)
meta_subj = {(st, m["sample_id"]): m["subject_id"] for st in (NIEL, LI) for m in metadata(st) if m["body_site"] == "stool" and m["subject_id"]}
allruns = json.load(open(f"{DATA}/biosample_runs.json")); cs = json.load(open(f"{DATA}/curated_sample_biosample.json"))
bs2s = {}
for k, bs in cs.items():
    st, sid = k.split("|", 1)
    if (st, sid) in meta_subj:
        bs2s[bs] = (st, meta_subj[(st, sid)], sid)
acc = {run: bs2s[v["biosample"]] for run, v in allruns.items() if v["biosample"] in bs2s}
NARROW = set(gzip.open(f"{DATA}/logan_gut_wgs_narrow_ids.txt.gz", "rt").read().split("\n"))
res = {}
for fn in ("results.jsonl", "results_pass2.jsonl"):
    path = f"{DATA}/{fn}"
    if os.path.exists(path):
        for l in open(path):
            try:
                q = json.loads(l)
            except Exception:
                continue
            if "results" in q and (q["study"], q["subject"]) not in res:
                res[(q["study"], q["subject"])] = q
rows = {(r["study"], r["subject"]): r for r in json.load(open(ROWS)) if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"}

per = {NIEL: [], LI: []}
for (study, subj), r in rows.items():
    if study not in PAIR:
        continue
    q = res[(study, subj)]; v2s = r["v2s"]; other = PAIR[study]
    same, imp0, a0, imp1, a1 = [], [], [], [], []
    top = (0.0, None, None)
    for a, sc in q["results"].items():
        m = acc.get(a)
        if sc >= DUP:
            continue
        if m and m[0] == study and m[1] == subj:
            if m[2] != v2s:
                same.append((sc, m[0]))
            continue
        imp0.append(sc); a0.append(a)
        if sc > top[0]:
            top = (sc, m[0] if m else "unmapped", m[1] if m else None)
        if m and m[0] == other and m[1] == subj:
            same.append((sc, m[0]))
        else:
            imp1.append(sc); a1.append(a)
    i0, i1 = np.array(imp0), np.array(imp1)
    n0 = np.array([a in NARROW for a in a0], dtype=bool); n1 = np.array([a in NARROW for a in a1], dtype=bool)
    p_old = max([s for s, st in same if st == study], default=0.0); p_new = max([s for s, st in same], default=0.0)
    rk = lambda imp, nar, p: (int(((imp > p) & nar).sum()) + 1) if p > 0 else None
    per[study].append(dict(old=rk(i0, n0, p_old), new=rk(i1, n1, p_new), cross=any(st == other for s, st in same), top_study=top[1], top_same_id=(top[2] == subj)))
for st, tag in ((NIEL, "nielsen"), (LI, "li")):
    g = per[st]; n = len(g)
    out[f"{tag}_donors"] = n
    out[f"{tag}_linked_all"] = sum(r["rank_any"] == 1 for (stx, _), r in rows.items() if stx == st)          # rank 1 among all the runs of the index
    out[f"{tag}_with_same_id"] = sum(x["cross"] for x in g)
    out[f"{tag}_rank1_before"] = sum(x["old"] == 1 for x in g); out[f"{tag}_rank1_after"] = sum(x["new"] == 1 for x in g)
    out[f"{tag}_top10_before"] = sum(x["old"] is not None and x["old"] <= 10 for x in g); out[f"{tag}_top10_after"] = sum(x["new"] is not None and x["new"] <= 10 for x in g)
    fail = [x for x in g if x["old"] != 1 and x["old"] is not None]
    out[f"{tag}_outranked"] = len(fail)
    out[f"{tag}_best_other_cohort"] = sum(x["top_study"] == PAIR[st] for x in fail)
    out[f"{tag}_best_same_id"] = sum(x["top_study"] == PAIR[st] and x["top_same_id"] for x in fail)

# Nielsen: contamination between participants
md = {m["sample_id"]: m for m in metadata(NIEL)}; subj = {k: m["subject_id"] for k, m in md.items()}
ev = [(d["source"], d["target"], float(d["rate"])) for d in csv.DictReader([l for l in open(f"{CRO}/results/{NIEL}.conta.tsv") if not l.startswith("#")], delimiter="\t")]
oth = [e for e in ev if subj.get(e[0]) != subj.get(e[1])]
out.update(nielsen_events=len(ev), nielsen_events_other=len(oth), nielsen_median_rate_other=float(np.median([e[2] for e in oth])))
status = collections.Counter(); linked = collections.Counter()
for (study, subj_id), r in rows.items():
    if study != NIEL:
        continue
    Q = r["v2s"]; s_same = s_oth = False
    for s, t, _ in ev:
        if Q in (s, t):
            x = subj.get(t if s == Q else s)
            if x == subj_id: s_same = True
            elif x: s_oth = True
    k = "other" if s_oth else "same" if s_same else "none"
    status[k] += 1; linked[k] += int(r["rank_any"] == 1)
out.update(nielsen_status=dict(status), nielsen_status_linked=dict(linked))
lines = open(f"{CRO}/matrices/{NIEL}.species_relab.tsv").read().split("\n"); hdr = lines[0].split("\t")
tab = [l.split("\t") for l in lines[1:] if l.strip()]
X = np.array([[float(v) for v in r[1:]] for r in tab]).T; names = hdr[1:]
X = X / np.maximum(X.sum(1, keepdims=True), 1e-12)
sub = np.array([subj.get(x) for x in names]); n = len(names); iu = np.triu_indices(n, 1)
bc = np.array([np.abs(X[i] - X[j]).sum() / max((X[i] + X[j]).sum(), 1e-12) for i, j in zip(*iu)])
diff = sub[iu[0]] != sub[iu[1]]
out.update(nielsen_samples=n, nielsen_pairs_diff=int(diff.sum()), nielsen_pairs_diff_bc_below_010=int((diff & (bc < 0.10)).sum()), nielsen_median_bc_same=float(np.median(bc[~diff])), nielsen_median_bc_diff=float(np.median(bc[diff])))
json.dump(out, open(OUT, "w"), indent=1)
print("wrote", OUT); print(json.dumps(out, indent=0)[:1500])
