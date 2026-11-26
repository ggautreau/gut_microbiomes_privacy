#!/usr/bin/env python3
"""Similarities of each donor's query to all the runs of the HGM set, computed from the hypervector database of Faure et al. (input of evt_model.py).

    python3 -I evt_similarities.py DB_DIR DATA_DIR ROWS_V4.json LISTS_V4.npz LINKS.json OUT.npz OUT.json

DB_DIR        big_db_reordered/ of hypervector_db.tar.gz (https://g-77cfc.ffdaa9.e229.data.globus.org/hypervector_db.tar.gz, about 10 Gb, 18 Gb unpacked):
              dimension.txt, vector_norms.txt (accession and norm), vectors.bin (one row of 2,048 values per accession). Not in the package.
DATA_DIR      the folder read by extract_lists.py (pairs.json, pairs_hmp2.json, results.jsonl, results_pass2.jsonl, logan_index_ids.txt.gz,
              logan_gut_wgs_narrow_ids.txt.gz = the HGM set, biosample_runs.json, curated_sample_biosample.json, runs_hmp2.json, mehta_sibling_runs.json, matrices/)
ROWS_V4.json  the 483 donors of the analysis; LISTS_V4.npz, the lists of extract_lists.py (their order, `k_rows`, `partner`); LINKS.json, hmp_sibling_links.json

The search service returns at most 5,000 runs for a query. Here the similarity of two runs is computed from their hypervectors a and b (divided by the square
root of the dimension) as the Jaccard index  a.b / (|a|^2 + |b|^2 - a.b), which is the similarity that the service returns (checked below on the pairs of the lists).
For each donor, with the identities and exclusions of extract_lists.py (runs of the donor, runs that public records give to the donor, re-deposits at 0.99 or more):
  s    similarity of its most similar other sample
  k    number of runs of other people of the HGM set above it, counted on the whole HGM set
  top  the TOP highest similarities of the runs of other people of the HGM set, in decreasing order
  n    number of runs of other people of the HGM set (present in the database)
OUT.npz holds these arrays in the order of ROWS_V4.json, with the cohort of each donor and whether its other sample was among the runs returned by the service;
OUT.json holds the checks (pairs compared with the service, largest difference, donors whose k equals that of the lists)."""
import csv, gzip, json, os, sys

_avail_cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
_default_threads = str(max(1, _avail_cpus - 2))
os.environ.setdefault("OPENBLAS_NUM_THREADS", _default_threads)
os.environ.setdefault("OMP_NUM_THREADS", _default_threads)
os.environ.setdefault("MKL_NUM_THREADS", _default_threads)

import numpy as np

DB, DATA, ROWS, LISTS, LINKS_F, OUT, OUTJ = sys.argv[1:8]
DUP = 0.99
TOP = 25_000                                                                      # enough for the test on thinned HGM sets of evt_model.py: at 100 times, about 250 of them are kept
M = 291_881

dim = int(open(f"{DB}/dimension.txt").read().split()[0])
names = [l.split()[0] for l in open(f"{DB}/vector_norms.txt") if l.strip()]
N = len(names); size = os.path.getsize(f"{DB}/vectors.bin")
assert size % (N * dim) == 0 and size // (N * dim) in (2, 4), (size, N, dim)
vec = np.memmap(f"{DB}/vectors.bin", dtype={2: np.int16, 4: np.int32}[size // (N * dim)], mode="r", shape=(N, dim))
row = {a: i for i, a in enumerate(names)}


def load(idx):
    """Rows of the database as float32, read in the order of the file."""
    idx = np.asarray(idx); order = np.argsort(idx); out = np.empty((len(idx), dim), np.float32)
    for s_ in range(0, len(idx), 20_000):
        o = order[s_:s_ + 20_000]; out[o] = vec[idx[o]]
    return out


def jaccard(d, a2, b2):
    return d / (a2 + b2 - d)


# ---------------------------------------------------------------- identities, as in extract_lists.py
runs = lambda s_: [a for a in s_.split(";") if a not in ("", "NA")]
P = json.load(open(f"{DATA}/pairs.json")) + [p for p in json.load(open(f"{DATA}/pairs_hmp2.json")) if p["study"] == "HMP_2019_ibdmdb"]
INDEX = set(gzip.open(f"{DATA}/logan_index_ids.txt.gz", "rt").read().split("\n"))
subset = {(p["study"], p["subject"]): p for p in P
          if runs(p["v1_acc"]) and runs(p["v2_acc"]) and any(a in INDEX for a in runs(p["v1_acc"])) and any(a in INDEX for a in runs(p["v2_acc"]))}
NARROW = [a for a in gzip.open(f"{DATA}/logan_gut_wgs_narrow_ids.txt.gz", "rt").read().split() if a]
assert len(NARROW) == M
res = {}
for fn in ("results.jsonl", "results_pass2.jsonl"):
    for l in open(f"{DATA}/{fn}"):
        try:
            q = json.loads(l)
        except Exception:
            continue
        if "results" in q and (q["study"], q["subject"]) not in res:
            res[(q["study"], q["subject"])] = q
allruns = json.load(open(f"{DATA}/biosample_runs.json")); cs = json.load(open(f"{DATA}/curated_sample_biosample.json"))
rows = json.load(open(ROWS)); keys = [(r["study"], r["subject"]) for r in rows]
meta_subj = {}
for study in {s_ for s_, _ in keys}:
    for r in csv.DictReader(open(f"{DATA}/matrices/{study}.metadata.tsv"), delimiter="\t"):
        if r["body_site"] == "stool" and r["subject_id"]:
            meta_subj[(study, r["sample_id"])] = r["subject_id"]
bs2s = {}
for k_, bs in cs.items():
    st, sid = k_.split("|", 1)
    if (st, sid) in meta_subj:
        bs2s[bs] = (st, meta_subj[(st, sid)], sid)
acc = {run: bs2s[v["biosample"]] for run, v in allruns.items() if v["biosample"] in bs2s}
acc.update({run: tuple(v) for run, v in json.load(open(f"{DATA}/runs_hmp2.json")).items() if v[0] == "HMP_2019_ibdmdb"})
if os.path.exists(f"{DATA}/mehta_sibling_runs.json"):
    acc.update({run: tuple(v) for run, v in json.load(open(f"{DATA}/mehta_sibling_runs.json")).items() if run not in acc})
LINKS = json.load(open(LINKS_F))
linked = {}
for a, who in LINKS.items():
    linked.setdefault(tuple(who), []).append(a)
own = {}
for run, m in acc.items():
    if run in row:
        own.setdefault((m[0], m[1]), []).append((run, m[2]))
lists = np.load(LISTS, allow_pickle=True)
assert [tuple(x) for x in zip(lists["study"], lists["subject"])] == keys, "the donors of the lists are those of the rows, in the same order"
k_pack, partner_pack = lists["k_rows"], lists["partner"]

# ---------------------------------------------------------------- the similarities
missing = [a for a in NARROW if a not in row]
nar_names = [a for a in NARROW if a in row]; pos = {a: j for j, a in enumerate(nar_names)}
H = load([row[a] for a in nar_names])
q_runs = [res[k_]["run"] for k_ in keys]
assert all(a in row for a in q_runs), "every queried run is in the database"
Q = load([row[a] for a in q_runs])
h2 = np.concatenate([np.einsum("ij,ij->i", c, c) for c in (H[s_:s_ + 40_000].astype(np.float64) for s_ in range(0, len(H), 40_000))]) / dim
q2 = np.einsum("ij,ij->i", Q.astype(np.float64), Q.astype(np.float64)) / dim
S = np.empty((len(Q), len(H)), np.float64)
for s_ in range(0, len(H), 40_000):
    d = (Q.astype(np.float64) @ H[s_:s_ + 40_000].astype(np.float64).T) / dim
    S[:, s_:s_ + 40_000] = jaccard(d, q2[:, None], h2[None, s_:s_ + 40_000])

# the similarity returned by the service, on every pair of the lists whose run is in the HGM set
diff = []; n_pairs = 0
for i, key in enumerate(keys):
    for a, sc in res[key]["results"].items():
        j = pos.get(a)
        if j is not None and sc < DUP:
            n_pairs += 1; diff.append(abs(S[i, j] - sc))
diff = np.array(diff)


def sim(a_run, b_run):
    va, vb = load([row[a_run]])[0].astype(np.float64), load([row[b_run]])[0].astype(np.float64); d = va @ vb / dim
    return float(jaccard(d, va @ va / dim, vb @ vb / dim))


top = np.zeros((len(keys), TOP), np.float32); n_imp = np.zeros(len(keys), np.int32); s_out = np.zeros(len(keys)); k_out = np.zeros(len(keys), np.int64)
same_k = 0; k_diff = []; partner_diff = 0.0
for i, key in enumerate(keys):
    study, subj = key; v2s = subset[key]["v2_sample"]; q_run = q_runs[i]; mine = own.get(key, [])
    same = [x for x in (sim(q_run, run) for run, sample in mine if sample != v2s and run != q_run) if x < DUP]
    partner = max(same, default=0.0)
    mask = S[i] < DUP
    for run, _ in mine:
        if run in pos:
            mask[pos[run]] = False
    for a in linked.get(key, []):
        if a in pos:
            mask[pos[a]] = False
    imp = S[i][mask]
    s_out[i] = partner; k_out[i] = int((imp > partner).sum()); n_imp[i] = len(imp)
    top[i] = -np.partition(-imp, TOP)[:TOP]; top[i] = -np.sort(-top[i])
    if partner_pack[i] > 0:
        partner_diff = max(partner_diff, abs(partner - partner_pack[i]))
        if k_out[i] == k_pack[i]:
            same_k += 1
        else:
            k_diff.append([int(k_pack[i]), int(k_out[i])])
returned = partner_pack > 0
np.savez_compressed(OUT + ".tmp.npz", top=top, n=n_imp, s=s_out, k=k_out, returned=returned, study=np.array([k_[0] for k_ in keys]))
os.replace(OUT + ".tmp.npz", OUT)                                                  # the file appears whole
checks = dict(database=dict(accessions=N, dimension=dim, bytes_per_value=size // (N * dim)), hgm_runs=M, hgm_runs_absent=len(missing), top=TOP,
              service=dict(pairs=n_pairs, max_abs_difference=float(diff.max()), median_abs_difference=float(np.median(diff))),
              returned=dict(donors=int(returned.sum()), same_k=same_k, k_lists_and_here=k_diff, max_abs_difference_other_sample=partner_diff),
              not_returned=dict(donors=int((~returned).sum()), k=sorted(int(v) for v in k_out[~returned]), similarity_other_sample=[float(v) for v in sorted(s_out[~returned])]),
              recovered=int((k_out == 0).sum()))
json.dump(checks, open(OUTJ, "w"), indent=1)
print(json.dumps({k_: v for k_, v in checks.items() if k_ != "not_returned"})[:600])
print("wrote", OUT, "and", OUTJ)
