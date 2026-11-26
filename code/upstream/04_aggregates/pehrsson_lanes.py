#!/usr/bin/env python3
"""Pehrsson et al. (PRJNA300541): samples of different participants with near-identical species profiles, and where they sit in the sequencing lanes.

    python3 -I pehrsson_lanes.py OUTPUT.json [ROWS.json [DATA_DIR]]

ROWS.json (optional, output of rows_with_partner.py) adds the 37 Pehrsson donors of the 680: how many have the query sample or the other sample in a near-duplicate pair, and how many of them are recovered.
DATA_DIR (optional, the DATA_DIR of rows_with_partner.py: results.jsonl, biosample_runs.json, curated_sample_biosample.json) adds what the searches of the index show for the adult donors whose queried sample is in one of
these pairs: for how many the nearest run, once the runs of the queried sample and the quasi-exact re-deposits are set aside as in rows_with_partner.py, is a run of the sample of the other participant of the pair,
and whether it is more similar than every other sample of the donor among the runs returned.

Reads, from crocodeel_valles/: the species table and the sample metadata of the cohort (curatedMetagenomicData), the CroCoDeEL result, and
pehrsson/pehrsson_run_alias.json (run accession -> run alias of the SRA, in the form <date>_ECP_Lane<n>/filtered_reads/Unpaired_<number>_filtered.txt).
The lane of a sample is the part of the alias of its unpaired run before '/filtered_reads'; its position in the lane is the rank of that run accession among
the samples of the lane in the table. Positions are inferred from accession order, not read from a barcode sheet. Only aggregates are written.
"""
import collections, csv, itertools, json, statistics as S, sys

CRO = "local/crocodeel_valles"
OUT = sys.argv[1]
ROWS = sys.argv[2] if len(sys.argv) > 2 else None
DATA = sys.argv[3] if len(sys.argv) > 3 else None
DUP = 0.99                                                                      # a run at this similarity or above is the query or a quasi-exact re-deposit (rows_with_partner.py)
md = {m["sample_id"]: m for m in csv.DictReader(open(f"{CRO}/matrices/PehrssonE_2016.metadata.tsv"), delimiter="\t")}
alias = json.load(open(f"{CRO}/pehrsson/pehrsson_run_alias.json"))
rows = [l.rstrip("\n").split("\t") for l in open(f"{CRO}/matrices/PehrssonE_2016.species_relab.tsv")]
samples = rows[0][1:]
X = {s: [float(r[1 + k]) for r in rows[1:]] for k, s in enumerate(samples)}
bc = lambda a, b: sum(abs(x - y) for x, y in zip(X[a], X[b])) / sum(x + y for x, y in zip(X[a], X[b]))
person = {s: md[s]["subject_id"] for s in samples}
assert len(samples) == 191 and set(samples) == set(md)

lane, acc = {}, {}
for s in samples:
    un = [r for r in md[s]["NCBI_accession"].split(";") if "/filtered_reads/Unpaired_" in alias.get(r, "")]
    if un:
        lane[s] = alias[un[0]].split("/")[0]; acc[s] = int(un[0][3:])
by_lane = collections.defaultdict(list)
for s in lane:
    by_lane[lane[s]].append(s)
pos = {}
for L, ss in by_lane.items():
    for i, s in enumerate(sorted(ss, key=lambda s: acc[s])):
        pos[s] = i
LANES = ["141118_ECP_Lane1", "141118_ECP_Lane2", "150102_ECP_Lane1", "150102_ECP_Lane2", "150102_ECP_Lane3", "150102_ECP_Lane4"]
RUN = "140805_ECP"                                                              # no lane in the file names

diff_all = []; same_all = []; near = []
for a, b in itertools.combinations(samples, 2):
    d = bc(a, b)
    (same_all if person[a] == person[b] else diff_all).append(d)
    if person[a] != person[b] and d < 0.10:
        near.append((a, b, d))
n_cross = sum(1 for a, b, _ in near if lane.get(a) != lane.get(b))

def off(a, b):
    return abs(pos[a] - pos[b])

per_lane = {}
for L in LANES + [RUN]:
    ss = by_lane[L]; at12 = [(a, b) for a, b in itertools.combinations(ss, 2) if person[a] != person[b] and off(a, b) == 12]
    per_lane[L] = dict(samples=len(ss), pairs12=len(at12), near12=sum(1 for a, b in at12 if bc(a, b) < 0.10),
                       near_all=sum(1 for a, b, _ in near if lane.get(a) == L and lane.get(b) == L))
offsets = collections.Counter(off(a, b) for a, b, _ in near if lane.get(a) in LANES and lane.get(a) == lane.get(b))
ctrl = [0, 0]
for L in LANES:
    for a, b in itertools.combinations(by_lane[L], 2):
        if person[a] != person[b] and off(a, b) in (1, 2, 11, 13):
            ctrl[0] += 1; ctrl[1] += bc(a, b) < 0.10
run_near = sorted(off(a, b) for a, b, _ in near if lane.get(a) == RUN and lane.get(b) == RUN)

donors = None
if ROWS:
    near_samples = {x for a, b, _ in near for x in (a, b)}
    R = [r for r in json.load(open(ROWS)) if r["study"] == "PehrssonE_2016" and r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]
    hit = [r for r in R if r["v2s"] in near_samples or r["best_any"] in near_samples]; rest = [r for r in R if r not in hit]
    donors = dict(n=len(R), linked=sum(r["rank_any"] == 1 for r in R), with_pair=len(hit), with_pair_linked=sum(r["rank_any"] == 1 for r in hit),
                  query_in_pair=sum(r["v2s"] in near_samples for r in R), query_in_pair_linked=sum(r["rank_any"] == 1 for r in R if r["v2s"] in near_samples),
                  without_pair=len(rest), without_pair_linked=sum(r["rank_any"] == 1 for r in rest))
search = None
if ROWS and DATA:
    # What the searches show: the nearest run of the donors whose queried sample is in a pair. Runs are assigned to samples through their BioSample, as in rows_with_partner.py.
    pair_of = {}
    for a, b, _ in near:
        assert a not in pair_of and b not in pair_of, "a sample is in two pairs: the other participant of the pair is not defined"
        pair_of[a] = b; pair_of[b] = a
    cs = json.load(open(f"{DATA}/curated_sample_biosample.json")); allruns = json.load(open(f"{DATA}/biosample_runs.json"))
    bs2s = {cs[f"PehrssonE_2016|{s}"]: s for s in samples if f"PehrssonE_2016|{s}" in cs}
    run_sample = {run: bs2s[v["biosample"]] for run, v in allruns.items() if v["biosample"] in bs2s}
    assert len(bs2s) == len(samples) and all(run_sample.get(a) == s for s in samples for a in md[s]["NCBI_accession"].split(";")), "a sample has no BioSample, or a curated run is not assigned to its sample"
    res = {}
    for fn in ("results.jsonl", "results_pass2.jsonl"):
        try:
            fh = open(f"{DATA}/{fn}")
        except FileNotFoundError:
            continue
        for l in fh:
            try:
                q = json.loads(l)
            except Exception:
                continue
            if "results" in q and q["study"] == "PehrssonE_2016" and q["subject"] not in res:
                res[q["subject"]] = q
    partner_first, above_own, sc_partner, sc_own = 0, 0, [], []
    Q = [r for r in R if r["v2s"] in near_samples]
    for r in Q:
        q = res[r["subject"]]; own, others = [], []
        assert run_sample.get(q["run"]) == r["v2s"], "the stored query is not a run of the queried sample of the row"
        for a, sc in q["results"].items():
            if sc >= DUP:
                assert run_sample.get(a) == r["v2s"], "a run at 0.99 or above is not a run of the queried sample: the Note says that only the runs of the queried sample are set aside"
                continue
            s = run_sample.get(a)
            if s == r["v2s"]:
                continue                                                         # another run of the queried sample: neither the donor's other sample nor another person
            (own if s is not None and person[s] == r["subject"] else others).append((sc, a))
        top_sc, top_run = max(others); best_own = max(own, default=(0.0, None))[0]
        assert len({run_sample.get(a) for sc, a in others if sc == top_sc}) == 1, "two samples tie for the nearest run"
        assert abs(top_sc - r["max_imp"]) < 1e-9 and abs(best_own - (r["partner_any"] or 0.0)) < 1e-9, "the lists are not read as rows_with_partner.py reads them"
        if run_sample.get(top_run) == pair_of[r["v2s"]]:
            partner_first += 1; sc_partner.append(top_sc)
            assert r["top_kind"] == "same_cohort_other_person", "rows_with_partner.py does not see this run as a person of the same cohort"
            if top_sc > best_own:
                above_own += 1
                if own:
                    sc_own.append(best_own)
    search = dict(query_in_pair=len(Q), nearest_is_pair_partner=partner_first, nearest_is_pair_partner_above_own=above_own, median_score_pair_partner=S.median(sc_partner) if sc_partner else None,
                  own_sample_returned=len(sc_own), median_score_best_own_sample=S.median(sc_own) if sc_own else None)
ev = [(d["source"], d["target"], float(d["rate"])) for d in csv.DictReader([l for l in open(f"{CRO}/results/PehrssonE_2016.conta.tsv") if not l.startswith("#")], delimiter="\t")]
ev_other = [(s, t, r) for s, t, r in ev if person[s] != person[t]]
out = dict(samples=len(samples), n_diff_pairs=len(diff_all), n_same_pairs=len(same_all), median_bc_diff=S.median(diff_all), median_bc_same=S.median(same_all),
           n_near=len(near), n_near_cross_lane=n_cross, n_near_in_lanes=sum(offsets.values()), near_offsets_in_lanes=dict(sorted(offsets.items())), near_run_offsets=run_near,
           control_offsets_1_2_11_13=ctrl, per_lane=per_lane,
           n_events=len(ev), n_events_other=len(ev_other), median_rate_events_other=S.median(r for _, _, r in ev_other),
           median_bc_events_other=S.median(bc(s, t) for s, t, _ in ev_other), donors=donors, search=search)
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps(out, indent=1))
