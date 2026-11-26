"""HMP IBD pairs, version 3 (the rule of the other cohorts). The samples of a donor are its curated stool samples. curatedMetagenomicData lists no run for this cohort, so a visit counts as a
curated sample when a curated sample of the donor has number_reads = 2 x the read_count of one of its ENA runs (the link of hmp_curated_ids.py). The first and the last sample are the earliest and the latest of
these visits among those with such a run in the index; the query is the deepest of the runs of the last visit (as build_hmp2.py picks the deepest indexed run). A donor whose first and last visits do not
change keeps its pair. The earlier pairs (build_hmp2.py, kept by build_hmp2_v2.py) took the first and the last date among the runs whose ENA alias has no suffix, which left out some visits.
Reads the earlier pairs (DATA_DIR/pairs_hmp2.json) and writes OUT_DIR/pairs_hmp2.json, the same file with the entries of the adult donors of HMP IBD replaced, OUT_DIR/pairs_hmp2_v3_adults.json (the pairs of these adult donors alone) and OUT_DIR/requery_todo.json,
the donors whose query run changes (to be sent to the service again with requery_largest.py). No network.

    python3 -I build_hmp2_v3.py ADULT_DONORS_ROWS.json DATA_DIR ENA_TSV OUT_DIR
"""
import collections, csv, datetime, gzip, json, sys
ROWS, DATA, ENA, OUT = sys.argv[1:5]
idx = set(gzip.open(f"{DATA}/logan_index_ids.txt.gz", "rt").read().split("\n"))
runs = {k: v for k, v in json.load(open(f"{DATA}/runs_hmp2.json")).items() if v[0] == "HMP_2019_ibdmdb"}
ena = {e["run_accession"]: e for e in csv.DictReader(open(ENA), delimiter="\t")}
cur = collections.defaultdict(lambda: collections.defaultdict(list))
for m in csv.DictReader(open(f"{DATA}/matrices/HMP_2019_ibdmdb.metadata.tsv"), delimiter="\t"): cur[m["subject_id"]][int(float(m["number_reads"]))].append(m["sample_id"])
vis = collections.defaultdict(lambda: collections.defaultdict(list))                  # participant -> visit key -> [(run, date, reads)] of the indexed runs with a curated sample
for run, (st, part, key) in runs.items():
    e = ena.get(run)
    if run in idx and e and e.get("collection_date") and e.get("read_count") and cur[part].get(2 * int(e["read_count"])):
        vis[part][key].append((run, datetime.date.fromisoformat(e["collection_date"]).toordinal(), int(e["read_count"])))
old = {p["subject"]: p for p in json.load(open(f"{DATA}/pairs_hmp2.json")) if p["study"] == "HMP_2019_ibdmdb"}
adults = [r for r in json.load(open(ROWS)) if r["study"] == "HMP_2019_ibdmdb"]
new_pairs, todo, c = [], [], collections.Counter()
for r in adults:
    part = r["subject"]; V = vis[part]; vd = {k: min(d for _, d, _ in L) for k, L in V.items()}
    first, last = min(vd, key=vd.get), max(vd, key=vd.get); p = dict(old[part])
    if first == p["v1_sample"] and last == p["v2_sample"]: new_pairs.append(p); c["unchanged"] += 1; continue
    a = max(V[first], key=lambda x: x[2]); b = max(V[last], key=lambda x: x[2])
    c["query changed"] += last != p["v2_sample"]; c["first sample changed"] += first != p["v1_sample"]
    p.update(v1_acc=a[0], v2_acc=b[0], v1_sample=first, v2_sample=last, days=vd[last] - vd[first], v1_reads=str(a[2]), v2_reads=str(b[2]), n_samples=len(V)); new_pairs.append(p)
    if last != old[part]["v2_sample"]: todo.append({"study": "HMP_2019_ibdmdb", "subject": part, "new": b[0], "old": old[part]["v2_acc"], "old_visit": old[part]["v2_sample"], "new_visit": last, "reads": b[2]})
replaced = {p["subject"]: p for p in new_pairs}
json.dump([replaced.get(p["subject"], p) if p["study"] == "HMP_2019_ibdmdb" else p for p in json.load(open(f"{DATA}/pairs_hmp2.json"))], open(f"{OUT}/pairs_hmp2.json", "w"))
json.dump(new_pairs, open(f"{OUT}/pairs_hmp2_v3_adults.json", "w"), indent=0); json.dump(todo, open(f"{OUT}/requery_todo.json", "w"), indent=0)
print(dict(c), "| new queries:", len(todo), "| reads of the new queries (min, median):", min(t["reads"] for t in todo), sorted(t["reads"] for t in todo)[len(todo) // 2])
