#!/usr/bin/env python3
"""How the cohorts of the analysis were chosen, from curatedMetagenomicData to the donors of the analysis, with the sex and the country of these donors (Supplementary Note 1 and Supplementary Tables 1 and 2).
Writes aggregate numbers only (no donor identifier) for build_supplementary.py.

    python3 -I cohort_selection.py ROWS_FINAL.json ROWS_ANALYSED.json OUT.json [EXCLUDED.json [ROWS_FINAL_WITHOUT_LINKS.json]]

ROWS_FINAL.json     rows of the 831 donors whose first and last samples are in the index (rows_v2_full.json of rows_with_partner.py: study, subject, age_class, q_lib)
ROWS_ANALYSED.json  rows of the donors of the analysis (rows_v4.json of exclude_donors.py, the 483 donors)
EXCLUDED.json       optional: donors_excluded.json of exclude_donors.py, the donors that the rules of contamination and re-deposit remove
ROWS_FINAL_WITHOUT_LINKS.json  optional: the rows of rows_with_partner.py without --hmp, to count the ranks that the runs given to the donors of HMP IBD change
Inputs read here: the 93 metadata tables and cohort_age.csv of crocodeel_valles/ (the species tables of curatedMetagenomicData, written on 1 October 2026), the pairs of adultes/ (first and last stool sample of each
donor of the adult cohorts) and the curated age category of all samples (cmd_meta.tsv, same release of the sample metadata).
Also read, at fixed paths: in local/adultes/ the index list (logan_index_ids.txt.gz), the BioSample tables (biosample_runs.json, curated_sample_biosample.json), the tables of the ENA (ena_PRJNA354235_runs.tsv, ibdmdb_runs_ena.tsv), the pairs and runs of HMP IBD (pairs_hmp2.json, runs_hmp2.json, runs_hmp2_dates.json) and the first answers (results.jsonl); the merged answers (local/adultes_data_v2/results.jsonl); and rows_v4_curated.json in the folder of ROWS_ANALYSED.
The longitudinal rule is that of cohort_age.py: a donor has at least two stool samples, a cohort is a study with at least five such donors.
"""
import collections, csv, datetime, glob, gzip, json, os, re, sys
import numpy as np

ROWS_FINAL, ROWS_ANALYSED, OUT = sys.argv[1:4]
EXCL = json.load(open(sys.argv[4])) if len(sys.argv) > 4 else None             # donors_excluded.json (exclude_donors.py): the donors that the rules of contamination and re-deposit remove
CRO = "local/crocodeel_valles"
A = "local/adultes"
CMD = "local/cmd_meta/cmd_meta.tsv"
MIN_DONORS = 5                                                                 # donors with two stool samples that a study needs to be a cohort (cohort_age.py, MIN_SUBJECTS)

# ---- the catalogue: every study of curatedMetagenomicData
studies = {}; stool_n = collections.defaultdict(collections.Counter)
for path in sorted(glob.glob(f"{CRO}/matrices/*.metadata.tsv")):
    st = os.path.basename(path)[:-len(".metadata.tsv")]
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    studies[st] = len(rows)
    for r in rows:
        if r["body_site"] == "stool" and r["subject_id"]:
            stool_n[st][r["subject_id"]] += 1
n_samples = sum(studies.values())
n_stool = sum(sum(c.values()) for c in stool_n.values())
n_stool_studies = sum(1 for st in studies if stool_n[st])
two = {st: sum(1 for k in c.values() if k >= 2) for st, c in stool_n.items()}
longit = {st: k for st, k in two.items() if k >= MIN_DONORS}
below = {st: k for st, k in two.items() if 0 < k < MIN_DONORS}
out: dict = dict(catalogue=dict(studies=len(studies), samples=n_samples, stool=n_stool, stool_studies=n_stool_studies),
           two_stool=dict(studies=sum(1 for k in two.values() if k), donors=sum(two.values())),
           longitudinal=dict(min_donors=MIN_DONORS, cohorts=len(longit), donors=sum(longit.values()), below_threshold_studies=len(below), below_threshold_donors=sum(below.values())))

# ---- the age label of each cohort (crocodeel_valles/reproduce/cohort_age.csv)
lab = {r["dataset"]: r for r in csv.DictReader(open(f"{CRO}/reproduce/cohort_age.csv"))}
assert set(lab) == set(longit) and all(int(lab[s]["subjects"]) == longit[s] for s in lab), "cohort_age.csv does not match the cohorts counted from the metadata tables"
groups = collections.defaultdict(list)
for s, r in lab.items():
    groups[r["population"]].append(s)
out["labels"] = {g: dict(cohorts=len(v), donors=sum(longit[s] for s in v)) for g, v in groups.items()}

# the cohorts that could not be labelled (fewer than 70% of the donors with an age): what the curated age category says about their stool samples
cat = collections.defaultdict(lambda: collections.defaultdict(set))
with open(CMD) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["body_site"] == "stool" and r["study_name"] in groups["not recorded"]:
            cat[r["study_name"]][r["age_category"]].add(r["subject_id"])
adult_cat = [s for s in groups["not recorded"] if set(cat[s]) == {"adult"}]
out["not_recorded_adult_category"] = dict(cohorts=len(adult_cat), donors=sum(longit[s] for s in adult_cat), names=sorted(adult_cat))

# ---- the 14 adult cohorts: from the donors with two stool samples to the donors of the analysis
pairs = json.load(open(f"{A}/pairs.json"))
n_pairs = collections.Counter(p["study"] for p in pairs)
final = json.load(open(ROWS_FINAL)); ana = json.load(open(ROWS_ANALYSED))
assert set(n_pairs) == set(groups["adults"]) and all(n_pairs[s] == longit[s] for s in n_pairs), "the pairs of pairs.json are not the donors of the 14 adult cohorts"
cohorts = {}
for st in groups["adults"]:
    f_ = [r for r in final if r["study"] == st]
    cohorts[st] = dict(donors=n_pairs[st], indexed=len(f_),
                       under18=sum(r["age_class"] == "under18_or_mixed" for r in f_), no_age=sum(r["age_class"] == "unknown" for r in f_), rna=sum(r["q_lib"] == "RNA" for r in f_),
                       adult_dna=sum(r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA" for r in f_),
                       analysed=sum(r["study"] == st for r in ana))
tot = {k: sum(c[k] for c in cohorts.values()) for k in next(iter(cohorts.values()))}
assert tot["indexed"] - tot["under18"] - tot["no_age"] - tot["rna"] == tot["adult_dna"], "the rules of age and library do not add up"
if EXCL is not None:                                                                  # the two rules that remove donors, between the adults with a DNA query of the cohorts kept and the analysis
    for st, c in cohorts.items():
        c["contamination"] = sum(1 for e in EXCL if e["study"] == st and e["rule"] == "contamination")
        c["re_deposit"] = sum(1 for e in EXCL if e["study"] == st and e["rule"] == "re-deposit")
        assert c["analysed"] == 0 or c["adult_dna"] - c["contamination"] - c["re_deposit"] == c["analysed"], st
    tot = {k: sum(c[k] for c in cohorts.values()) for k in next(iter(cohorts.values()))}
    out["excluded"] = dict(candidates=sum(c["adult_dna"] for c in cohorts.values() if c["analysed"]), contamination=sum(e["rule"] == "contamination" for e in EXCL),
                           contamination_query=sum(e["rule"] == "contamination" and e["why"] == "query sample" for e in EXCL),
                           contamination_no_other=sum(e["rule"] == "contamination" and e["why"] == "no other sample left" for e in EXCL),
                           samples_removed=sum(len(e["removed_samples"]) for e in EXCL), re_deposit=sum(e["rule"] == "re-deposit" for e in EXCL),
                           recovered_contamination=sum(e["rule"] == "contamination" and e["rank_any_narrow"] == 1 for e in EXCL),
                           recovered_re_deposit=sum(e["rule"] == "re-deposit" and e["rank_any_narrow"] == 1 for e in EXCL),
                           min_rate=min(x for e in EXCL for x in e.get("rates", [])), max_rate=max(x for e in EXCL for x in e.get("rates", [])),
                           by_cohort={st: [c["contamination"], c["re_deposit"]] for st, c in cohorts.items() if c["contamination"] or c["re_deposit"]})
out["adult_cohorts"] = cohorts
out["adult_totals"] = tot

# ---- sex and country of the donors of the analysis, from the curated metadata of their stool samples (Supplementary Table 2; aggregates only)
meta_ = collections.defaultdict(lambda: dict(gender=set(), country=set()))
for st in sorted({r["study"] for r in ana}):
    for r in csv.DictReader(open(f"{CRO}/matrices/{st}.metadata.tsv"), delimiter="\t"):
        if r["body_site"] == "stool" and r["subject_id"]:
            for k_ in ("gender", "country"):
                v_ = (r.get(k_) or "").strip()
                if v_ and v_ != "NA":
                    meta_[(st, r["subject_id"])][k_].add(v_)


def demo_(rows_):
    """Women, men, donors without a recorded sex, and countries (ISO 3166 alpha-3, as curatedMetagenomicData gives them) of a set of donors."""
    d_ = dict(donors=len(rows_), women=0, men=0, sex_not_recorded=0, countries=collections.Counter())
    for r in rows_:
        m_ = meta_[(r["study"], r["subject"])]
        assert len(m_["gender"]) <= 1 and m_["gender"] <= {"female", "male"}, ("one sex per donor", r["study"])
        d_["women" if m_["gender"] == {"female"} else "men" if m_["gender"] == {"male"} else "sex_not_recorded"] += 1
        assert len(m_["country"]) == 1, ("one country per donor", r["study"])
        d_["countries"][next(iter(m_["country"]))] += 1
    d_["countries"] = dict(sorted(d_["countries"].items(), key=lambda kv: (-kv[1], kv[0])))
    return d_


out["demographics"] = dict(by_cohort={st: demo_([r for r in ana if r["study"] == st]) for st in sorted({r["study"] for r in ana})}, total=demo_(ana))
assert out["demographics"]["total"]["donors"] == tot["analysed"] and all(out["demographics"]["by_cohort"][st]["donors"] == c["analysed"] for st, c in cohorts.items() if c["analysed"]), "the donors of the demographics are those of the analysis"
out["adult_by_recorded_age"] = sum(r["age_class"] == "adult" and r["q_lib"] == "DNA" for r in final)
out["adult_by_design_dna"] = sum(r["age_class"] == "adult_by_design" and r["q_lib"] == "DNA" for r in final)

# ---- what the text of Supplementary Note 1 adds about the closed set, the RNA queries and the interval between the two samples
out["ibd_rebuilt_pairs"] = len(json.load(open(f"{A}/pairs_hmp2.json")))                               # HMP IBD: pairs rebuilt from the ENA project (the curated metadata hold no run accession)
idx = set(gzip.open(f"{A}/logan_index_ids.txt.gz", "rt").read().split())
bs_runs = json.load(open(f"{A}/biosample_runs.json")); cs = json.load(open(f"{A}/curated_sample_biosample.json"))
runs_of = collections.defaultdict(list)
for run, v in bs_runs.items():
    runs_of[v["biosample"]].append(run)
acc = lambda t: [x for x in t.split(";") if x and x != "NA"]
extra = collections.Counter()                                                                          # donors with both samples in the index by BioSample but not by the curated accessions
for p_ in pairs:
    if p_["study"] == "HMP_2019_ibdmdb":
        continue
    cur = [any(x in idx for x in acc(p_[f"v{j}_acc"])) for j in (1, 2)]
    bio = [any(r in idx for r in runs_of.get(cs.get(f"{p_['study']}|{p_[f'v{j}_sample']}"), [])) for j in (1, 2)]
    if all(c or b_ for c, b_ in zip(cur, bio)) and not all(cur):
        extra[p_["study"]] += 1
out["biosample_only_in_index"] = dict(extra)
rna = [r for r in final if r["q_lib"] == "RNA"]; mdna = [r for r in final if r["study"] == "MehtaRS_2018" and r["q_lib"] == "DNA"]
out["rna_queries"] = dict(n=len(rna), linked=sum(r["rank_any"] == 1 for r in rna), mehta_dna_n=len(mdna), mehta_dna_linked=sum(r["rank_any"] == 1 for r in mdna))
dd = np.array([float(r["days"]) for r in ana])
out["analysed_interval"] = dict(n=len(dd), median=float(np.median(dd)), q25=float(np.percentile(dd, 25)), q75=float(np.percentile(dd, 75)), same_day=int((dd == 0).sum()), under_week=int((dd < 7).sum()),
                                same_day_cohorts=sorted({r["study"] for r in ana if float(r["days"]) == 0}))
# ---- the three precisions of Supplementary Notes 1 and 2: the donors outside the curated accessions that were queried anyway; Mehta, the second indexed run of the query sample and the RNA libraries left by the DNA rule, which tests the query only; HMP IBD, the visits left out of the pairs
queried = {(q["study"], q["subject"]) for q in map(json.loads, open(f"{A}/results.jsonl")) if "results" in q}
out["biosample_only_queried"] = {st: sum((p_["study"], p_["subject"]) in queried and not all(any(x in idx for x in acc(p_[f"v{j}_acc"])) for j in (1, 2)) and all(any(x in idx for x in acc(p_[f"v{j}_acc"])) or any(r in idx for r in runs_of.get(cs.get(f"{st}|{p_[f'v{j}_sample']}"), [])) for j in (1, 2)) for p_ in pairs if p_["study"] == st) for st in extra}
MERGED = "local/adultes_data_v2/results.jsonl"       # the answers that feed the HGM set (merge_results.py): the 851 stool queries
answered = {(q["study"], q["subject"]) for q in map(json.loads, open(MERGED)) if "results" in q}
final_keys = {(r["study"], r["subject"]) for r in json.load(open(ROWS_FINAL))}
assert final_keys <= answered
out["queried_outside_the_donors"] = dict(sorted(collections.Counter(st for st, _ in answered - final_keys).items()))                                  # donors left out of the 831 whose last sample was queried all the same (their lists only feed the HGM set)
mehta = [r for r in ana if r["study"] == "MehtaRS_2018"]
out["mehta_second_run"] = dict(n=len(mehta), with_second_indexed_run=sum(sum(x in idx for x in runs_of.get(cs.get(f"MehtaRS_2018|{r['v2s']}"), [])) >= 2 for r in mehta))
# The BioSample of the query sample of Mehta 2018 holds a second indexed run that the curated metadata do not list: is it smaller (fewer bases in the ENA) than the run that they list?
ena_mehta = {e["run_accession"]: int(e["base_count"]) for e in csv.DictReader(open(f"{A}/ena_PRJNA354235_runs.tsv"), delimiter="\t")}
pair_of = {(p_["study"], p_["subject"]): p_ for p_ in pairs}
second_smaller = 0
for r in mehta:
    indexed = [x for x in runs_of.get(cs.get(f"MehtaRS_2018|{r['v2s']}"), []) if x in idx]
    if len(indexed) >= 2:
        listed = [x for x in indexed if x in acc(pair_of[("MehtaRS_2018", r["subject"])]["v2_acc"])]
        unlisted = [x for x in indexed if x not in listed]
        assert listed and unlisted, "a donor of Mehta has several indexed runs in the BioSample of its query sample but none that the curated metadata do not list, or none that they list"
        second_smaller += max(ena_mehta[x] for x in unlisted) < max(ena_mehta[x] for x in listed)
out["mehta_second_run"]["second_run_smaller"] = second_smaller
rna_lib = lambda sample: any(".rna" in bs_runs[x].get("alias", "").lower() for x in runs_of.get(cs.get(f"MehtaRS_2018|{sample}"), []))
dna_lib = lambda sample: any(".rna" not in bs_runs[x].get("alias", "").lower() for x in runs_of.get(cs.get(f"MehtaRS_2018|{sample}"), []))
first_rna = [r for r in mehta if rna_lib(r["v1s"]) and not dna_lib(r["v1s"])]
far_rna = [r for r in mehta if r["far_possible"] and all(rna_lib(x) and not dna_lib(x) for x in r["eligible"])]
out["rna_other_samples"] = dict(first_sample=len(first_rna), first_sample_first_in_stock=sum(r["rank_v1_narrow"] == 1 for r in first_rna), only_partner_30d=len(far_rna), only_partner_30d_first_in_stock=sum(r["rank_far_narrow"] == 1 for r in far_rna), donors=len({r["subject"] for r in first_rna + far_rna}), same_donors={r["subject"] for r in first_rna} == {r["subject"] for r in far_rna})
ibd_ena = {e["run_accession"]: e for e in csv.DictReader(open(f"{A}/ibdmdb_runs_ena.tsv"), delimiter="\t")}; ibd_dates = json.load(open(f"{A}/runs_hmp2_dates.json")); visit_day, subject_days = {}, collections.defaultdict(list)
for run, (st_, sj, vis) in json.load(open(f"{A}/runs_hmp2.json")).items():
    if st_ == "HMP_2019_ibdmdb" and run in ibd_dates and "_MVX" not in ibd_ena.get(run, {}).get("sample_alias", ""):
        d_ = datetime.date.fromisoformat(ibd_dates[run]).toordinal(); visit_day[(sj, vis)] = d_; subject_days[sj].append(d_)
ibd = [r for r in ana if r["study"] == "HMP_2019_ibdmdb"]
out["ibd_pairs"] = dict(n=len(ibd), query_not_last=sum(visit_day.get((r["subject"], r["v2s"]), 1e9) < max(subject_days[r["subject"]]) for r in ibd), first_not_first=sum(visit_day.get((r["subject"], r["v1s"]), -1) > min(subject_days[r["subject"]]) for r in ibd))
# the same two counts by the curated days (days_from_first_collection of curatedMetagenomicData), which date also the 300 curated samples of HMP IBD that no run of the ENA project resolves (suffix _P; hmp_curated_ids.py)
cur_day, cur_subject, subject_cur_days = {}, {}, collections.defaultdict(list)
for m_ in csv.DictReader(open(f"{CRO}/matrices/HMP_2019_ibdmdb.metadata.tsv"), delimiter="\t"):
    if m_["body_site"] == "stool" and m_["subject_id"] and m_["days_from_first_collection"] not in ("", "NA"):
        cur_day[m_["sample_id"]] = float(m_["days_from_first_collection"]); cur_subject[m_["sample_id"]] = m_["subject_id"]; subject_cur_days[m_["subject_id"]].append(cur_day[m_["sample_id"]])
ibd_c = [r for r in json.load(open(os.path.join(os.path.dirname(os.path.abspath(ROWS_ANALYSED)), "rows_v4_curated.json"))) if r["study"] == "HMP_2019_ibdmdb"]
out["ibd_pairs"].update(curated_n=len(ibd_c), curated_query_not_last=sum(cur_day[r["v2s"]] < max(subject_cur_days[cur_subject[r["v2s"]]]) for r in ibd_c if r["v2s"] in cur_day),
                        curated_first_not_first=sum(cur_day[r["v1s"]] > min(subject_cur_days[cur_subject[r["v1s"]]]) for r in ibd_c if r["v1s"] in cur_day), curated_first_unresolved=sum(r["v1s"] not in cur_day for r in ibd_c),
                        curated_samples=len(cur_day), curated_samples_unresolved=sum(sid.endswith("_P") for sid in cur_day))
earliest_ = {sj_: min((d_, x_) for x_, d_ in cur_day.items() if cur_subject[x_] == sj_)[1] for sj_ in set(cur_subject.values())}               # the earliest curated sample of each participant
out["ibd_pairs"]["curated_first_not_first_earliest_unresolved"] = sum(cur_day[r["v1s"]] > min(subject_cur_days[cur_subject[r["v1s"]]]) and earliest_[cur_subject[r["v1s"]]].endswith("_P") for r in ibd_c if r["v1s"] in cur_day)
if len(sys.argv) > 5:                                                                  # rows_with_partner.py without --hmp: the ranks that the runs given back to the donors of HMP IBD changed
    nol = {(r["study"], r["subject"]): r for r in json.load(open(sys.argv[5]))}
    chg = lambda f: [r for r in ana if r[f] != nol[(r["study"], r["subject"])][f]]
    out["hmp_links_effect"] = dict(rank_any_narrow=len(chg("rank_any_narrow")), rank_v1_narrow=len(chg("rank_v1_narrow")), rank_far_narrow=len(chg("rank_far_narrow")),
                                  first_alone_now_first=sum(r["rank_v1_narrow"] == 1 for r in chg("rank_v1_narrow")), first_alone_was_first=sum(nol[(r["study"], r["subject"])]["rank_v1_narrow"] == 1 for r in chg("rank_v1_narrow")),
                                  studies=sorted({r["study"] for f in ("rank_any_narrow", "rank_v1_narrow", "rank_far_narrow") for r in chg(f)}))
# ---- the samples of each donor of the analysis along time (Supplementary Fig. 1, bottom): the number of donors that have each arrangement of the query sample, the first sample,
# the most similar other sample and the other samples, apart for the donors recovered and for the others. The stool samples of a donor are sorted as in build_pairs.py
# (by days from first collection, a missing day counting as day 0, the order of the metadata table breaking the ties). Counts only, no identifier.
sa_day = lambda m_: float(m_.get("days_from_first_collection") or 0) if re.fullmatch(r"-?[\d.]+", m_.get("days_from_first_collection") or "0") else 0.0
sa_rows = json.load(open(os.path.join(os.path.dirname(os.path.abspath(ROWS_ANALYSED)), "rows_v4_curated.json")))       # the donors of the analysis with the sample identifiers of the curated metadata
assert [(r["study"], r["subject"], r["rank_any_narrow"]) for r in sa_rows] == [(r["study"], r["subject"], r["rank_any_narrow"]) for r in ana], "rows_v4_curated.json does not hold the donors of the analysis"
SA_MAX = 4                                                                     # the arrangements are listed for the donors with at most four stool samples; the donors with more are counted together
sa_samples = collections.defaultdict(list)
for sa_st in sorted({r["study"] for r in sa_rows}):
    for m_ in csv.DictReader(open(f"{CRO}/matrices/{sa_st}.metadata.tsv"), delimiter="\t"):
        if m_["body_site"] == "stool" and m_["subject_id"]:
            sa_samples[(sa_st, m_["subject_id"])].append(m_)
sa_listed, sa_more, sa_ties = collections.Counter(), collections.defaultdict(list), 0
for r in sa_rows:
    sa_sorted = sorted(sa_samples[(r["study"], r["subject"])], key=sa_day)        # a stable sort, as in build_pairs.py
    sa_ids = [m_["sample_id"] for m_ in sa_sorted]
    assert len(sa_ids) == len(set(sa_ids)) >= 2 and r["v1s"] in sa_ids and r["v2s"] in sa_ids and r["v1s"] != r["v2s"] and r["best_any"] != r["v2s"] and (not r["best_any"] or r["best_any"] in sa_ids), "a sample of a donor is not in its metadata table"
    sa_rec = r["rank_any_narrow"] == 1
    if len(sa_ids) > SA_MAX:
        sa_more[(sa_rec, bool(r["best_any"]))].append((len(sa_ids), r["study"]))
        continue
    assert r["v1s"] == sa_ids[0] and r["v2s"] == sa_ids[-1], "the first and the query sample of a donor with at most four stool samples are its earliest and its last"
    sa_ties += len({sa_day(m_) for m_ in sa_sorted[1:-1]}) < len(sa_sorted) - 2   # two samples between the first and the query on the same day: their order is that of the metadata table
    sa_listed[(tuple("query" if x == r["v2s"] else ("first_and_most_similar" if x == r["best_any"] else "first" if r["best_any"] else "first_none_returned") if x == r["v1s"]
                     else "most_similar" if x == r["best_any"] else "other" for x in sa_ids), sa_rec)] += 1
assert sa_ties == 0, "two samples between the first and the query of a donor share a day: the arrangement would rest on the order of the metadata table"
out["sample_arrangements"] = dict(
    max_samples_listed=SA_MAX,
    listed=[dict(roles=list(k_), recovered=rec_, donors=n_) for (k_, rec_), n_ in sorted(sa_listed.items(), key=lambda kv: (not kv[0][1], len(kv[0][0]), -kv[1], kv[0][0]))],
    more=[dict(recovered=rec_, most_similar_returned=ret_, donors=len(v_), min_samples=min(x for x, _ in v_), max_samples=max(x for x, _ in v_), studies=sorted({s_ for _, s_ in v_}))
          for (rec_, ret_), v_ in sorted(sa_more.items(), key=lambda kv: (not kv[0][0], not kv[0][1]))])
assert sum(a_["donors"] for k_ in ("listed", "more") for a_ in out["sample_arrangements"][k_]) == len(ana), "every donor of the analysis is in one arrangement"
json.dump(out, open(OUT, "w"), indent=1)
print("wrote", OUT)
print(json.dumps({k: v for k, v in out.items() if k != "adult_cohorts"}, indent=0)[:2500])
for s, c in sorted(cohorts.items(), key=lambda kv: -kv[1]["donors"]):
    print(f"{s:24}", c)
