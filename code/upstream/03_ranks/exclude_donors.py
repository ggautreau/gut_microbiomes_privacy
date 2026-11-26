#!/usr/bin/env python3
"""The donors of the analysis: the 497 donors of the 7 cohorts kept (make_rows_v4.py) less the donors that two rules remove, applied in this order.

    python3 -I exclude_donors.py ROWS_497.json ROWS_497_CURATED.json LISTS_497.npz CROCODEEL_DIR OUT_DIR

  contamination  the samples studied of a donor are its query sample, its first sample and its most similar other sample (v2s, v1s and best_any of the curated rows). A sample
                 studied that takes part, as source or as target, in a CroCoDeEL event between two people (two different subjects of the curated metadata) at a rate above 5% is
                 removed. The donor leaves the analysis when the sample removed is its query sample or when no other stool sample of the donor remains in the curated metadata;
                 a donor that would keep other samples would need a new rank, and the script stops (it does not happen here).
  re-deposit     a donor whose nearest run of other people in the HGM set (LISTS: imp_nar[0], the runs given back to the donors of HMP IBD by hmp_sibling_links.py not
                 counted) is more similar to the query than the other sample of every donor still in the analysis: the searches of Costea whose nearest run of other people
                 carries the reads of the donor's own sample under another subject code (Supplementary Note 1). The script stops if one of them is not from Costea.
Writes OUT_DIR/rows_v4.json, rows_v4_curated.json and lists_v4.npz (the donors kept, in the order of ROWS_497) and donors_excluded.json (the curated rows of the donors removed,
in the same order, with the rule and the samples removed; make_inputs.py writes data/donors_excluded.json from it).
"""
import csv, json, sys
import numpy as np

ROWS, CUR, LISTS, CRO, OUT = sys.argv[1:6]
RATE = 0.05
rows, cur = json.load(open(ROWS)), json.load(open(CUR))
assert len(rows) == len(cur) == 497 and all(a["study"] == b["study"] and a["subject"] == b["subject"] for a, b in zip(rows, cur))
z = np.load(LISTS, allow_pickle=True)
L = {f: z[f] for f in z.files}
assert all((s, d) == (r["study"], r["subject"]) for s, d, r in zip(L["study"].tolist(), L["subject"].tolist(), rows)), "the lists must follow the order of the rows"
out, removed = {}, {}
for st in sorted({r["study"] for r in cur}):
    md = {m["sample_id"]: m for m in csv.DictReader(open(f"{CRO}/matrices/{st}.metadata.tsv"), delimiter="\t")}
    subj = {k: m["subject_id"] for k, m in md.items()}
    other = [(d["source"], d["target"], float(d["rate"])) for d in csv.DictReader([l for l in open(f"{CRO}/results/{st}.conta.tsv") if not l.startswith("#")], delimiter="\t")]
    other = [(s, t, rt) for s, t, rt in other if subj.get(s) != subj.get(t)]
    for i, r in enumerate(cur):
        if r["study"] != st:
            continue
        studied = {x for x in (r["v1s"], r["v2s"], r["best_any"]) if x}
        hit = sorted(x for x in studied if any(x in (s, t) and rt > RATE for s, t, rt in other))
        if not hit:
            continue
        left = [k for k in md if subj[k] == r["subject"] and k != r["v2s"] and k not in hit]
        if r["v2s"] not in hit and left:
            raise SystemExit(f"a donor of {st} keeps other samples once its contaminated samples are removed: it would need a new rank")
        out[i] = dict(rule="contamination", why="query sample" if r["v2s"] in hit else "no other sample left",
                      rates=sorted({round(rt, 4) for s, t, rt in other for x in hit if x in (s, t) and rt > RATE}))
        removed[i] = hit
left_idx = [i for i in range(len(rows)) if i not in out]
pmax = max(float(L["partner"][i]) for i in left_idx)
for i in left_idx:
    if float(L["imp_nar"][i][0]) > pmax:
        assert rows[i]["study"] == "CosteaPI_2017", rows[i]["study"]
        out[i] = dict(rule="re-deposit", nearest_other=round(float(L["imp_nar"][i][0]), 4))
keep = [i for i in range(len(rows)) if i not in out]
json.dump([rows[i] for i in keep], open(f"{OUT}/rows_v4.json", "w"))
json.dump([cur[i] for i in keep], open(f"{OUT}/rows_v4_curated.json", "w"))
np.savez_compressed(f"{OUT}/lists_v4.npz", **{f: (np.array([L[f][i] for i in keep], dtype=object) if f == "imp_nar" else L[f][keep]) for f in L})
json.dump([dict(cur[i], removed_samples=removed.get(i, []), **out[i]) for i in sorted(out)], open(f"{OUT}/donors_excluded.json", "w"), indent=1)
by = lambda rule: sum(1 for v in out.values() if v["rule"] == rule)
print(json.dumps(dict(candidates=len(rows), contamination=by("contamination"), re_deposit=by("re-deposit"), kept=len(keep), recovered=sum(rows[i]["rank_any_narrow"] == 1 for i in keep),
                      samples_removed=sum(len(v) for v in removed.values()))))
