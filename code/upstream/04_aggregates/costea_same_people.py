#!/usr/bin/env python3
"""Supplementary Note 1 and Supplementary Table 9: the seven Costea searches whose nearest run of other people has a similarity above 0.7, and the deposits of the same reads under another subject code.

    python3 -I costea_same_people.py OUT.json ROWS.json OPEN_SET_SCORES.csv

ROWS.json             rows_v4_497.json (the 497 donors of the 7 cohorts kept, before the two rules, in the order of the analysis)
OPEN_SET_SCORES.csv   local/open_set_scores_497.csv (make_open_set_scores.py on lists_v4_497.npz: the two similarities of Fig. 1H for the 497 donors, in the same order)
Files of local/: the 5,000 nearest runs of every query (adultes_data_v2/results.jsonl), the HGM set (logan_gut_wgs_narrow_ids.txt.gz) and the curated metadata (cmd_meta.tsv).
Public records saved on 31 October 2026, read by accession only (costea_same_people/): the run lists of the ENA for four projects of EMBL (ena_PRJEB*.tsv), the BioSamples records of the queried and matched samples
(biosamples.json), the table of the four studies (ena_studies.tsv) and the counts of costea_reads_compare.py (reads_compare_first_tile.json) and of costea_trim_filter.py (trim_filter_first_tile.json, written on 8 November 2026). relabel.json names the one sample whose subject code we correct (see below). The output names the seven donors
by their subject identifier of curatedMetagenomicData and holds no run or sample accession of the deposits (only the accessions of the four projects).
"""
import collections, csv, datetime, gzip, json, math, re, sys

OUT, ROWS, OPEN_SET = sys.argv[1:4]
P = "local"
DATA = f"{P}/adultes_data_v2"
A = f"{P}/adultes"
CMD = f"{P}/cmd_meta/cmd_meta.tsv"
INV = f"{P}/costea_same_people"

STUDY = "CosteaPI_2017"
QUERY, RE_DEPOSIT, MATCH, THIRD = "PRJEB17632", "PRJEB24041", "PRJEB41102", "PRJEB50977"
OUTLIER = 0.7                                    # the seven stand apart: no other donor has a nearest run of other people above 0.5 (largest value 0.491)
DUP = 0.99                                       # quasi-exact re-deposits are removed from the lists at this similarity (Methods)
# one sample of the matched project is labelled with the subject code of another person, but its record (sex, age, family role) and its number continue the series of the donor's code
RELABEL = json.load(open(f"{INV}/relabel.json"))

rows = json.load(open(ROWS))
scores = list(csv.DictReader(open(OPEN_SET)))
assert len(rows) == len(scores) == 497
stock = set(gzip.open(f"{A}/logan_gut_wgs_narrow_ids.txt.gz", "rt").read().split("\n"))
ena = {}
for p in (QUERY, RE_DEPOSIT, MATCH, THIRD):
    for r in csv.DictReader(open(f"{INV}/ena_{p}.tsv"), delimiter="\t"):
        ena[r["run_accession"]] = r
studies = {r["study_accession"]: r for r in csv.DictReader(open(f"{INV}/ena_studies.tsv"), delimiter="\t")}
reads = {r["subject"]: r for r in json.load(open(f"{INV}/reads_compare_first_tile.json"))}
trim = {r["subject"]: r for r in json.load(open(f"{INV}/trim_filter_first_tile.json"))}                  # costea_trim_filter.py: the same first tile, trimmed or filtered?
biosample = {acc: {k: v[0]["text"] for k, v in d.get("characteristics", {}).items()} for acc, d in json.load(open(f"{INV}/biosamples.json")).items()}


def code(r):
    """Subject code of a run of the matched project (from the sample title) or of PRJEB50977 (from the sample name); None for the others."""
    if r["study_accession"] == MATCH and r["sample_alias"] in RELABEL:
        return RELABEL[r["sample_alias"]]
    m = None
    if r["study_accession"] == MATCH:
        m = re.search(r"subject (HD\.S\d+)", r["sample_title"])
    elif r["study_accession"] == THIRD:
        m = re.search(r"MyMicrobes_(HD\.S\d+)_ST_d\d+", r["sample_alias"])
    return m.group(1) if m else None


def day_of_third(r):
    m = re.search(r"_ST_d(\d+)$", r["sample_alias"])
    return int(m.group(1)) if m else None


meta = collections.defaultdict(list)
for r in csv.DictReader(open(CMD), delimiter="\t"):
    if r["study_name"] == STUDY:
        meta[r["subject_id"]].append(r)
res = {}
for line in open(f"{DATA}/results.jsonl"):
    q = json.loads(line)
    if q.get("study") == STUDY:
        res[q["subject"]] = q

# which donors, and what remains of the cohort
position = {(r["study"], r["subject"]): i for i, r in enumerate(rows)}
costea = [r for r in rows if r["study"] == STUDY]
seven = [r for r in costea if float(scores[position[(STUDY, r["subject"])]]["best_other"]) > OUTLIER]
assert len(seven) == 7 and all(r["subject"].startswith("DEU") for r in seven)
assert sum(float(s["best_other"]) > OUTLIER for s in scores) == 7, "no other donor of any cohort is above the cut"
countries = collections.Counter(re.sub(r"\d+", "", r["subject"]) for r in costea)
fail = [r for r in costea if r["rank_any_narrow"] != 1]
others = [r for r in fail if r not in seven]
assert all(r["subject"].startswith("KAZ") for r in others), "the failures left are all Kazakh donors"

donors = []
for r in seven:
    s = r["subject"]; i = position[(STUDY, s)]; q = res[s]; hits = q["results"]; qrec = ena[q["run"]]
    assert qrec["study_accession"] == QUERY
    partner = r["partner_any"] or 0.0
    assert abs(partner - float(scores[i]["partner"])) < 1e-9
    own_runs = [a for a, e in ena.items() if e["study_accession"] == QUERY and e["sample_alias"] == qrec["sample_alias"] and a != q["run"]]     # the other runs of the queried sample
    nr = sorted(((sc, a) for a, sc in hits.items() if sc < DUP and a in stock and a not in own_runs), reverse=True)        # the nearest runs of the HGM set: not the query, the other runs of its sample, nor a quasi-exact re-deposit; those above the donor's other sample are runs of other people
    top_sc, top = nr[0]
    assert abs(top_sc - float(scores[i]["best_other"])) < 1e-9 and ena[top]["study_accession"] == MATCH
    trec = ena[top]; cd = code(trec)
    # an exact re-deposit of the queried run in the project of the same centre
    exact = [a for a, sc in hits.items() if sc >= DUP and a != q["run"]]
    assert len(exact) == 1 and ena[exact[0]]["study_accession"] == RE_DEPOSIT and ena[exact[0]]["base_count"] == qrec["base_count"]
    # the other runs of the two samples
    match_runs = [a for a, e in ena.items() if e["study_accession"] == MATCH and e["sample_alias"] == trec["sample_alias"] and a != top]
    own_sc = [hits[a] for a in own_runs if a in hits]; match_sc = [hits[a] for a in match_runs if a in hits]
    # the reads, on the first tile of the two runs
    rd = reads[s]
    assert rd["same_tile"] and rd["contained"] == rd["common"]
    tf = trim[s]
    assert (tf["reads_q"], tf["reads_n"], tf["common"]) == (rd["reads_q"], rd["reads_n"], rd["common"]), "costea_trim_filter.py and costea_reads_compare.py did not read the same tile"
    assert tf["relation"]["queried_inside_matched"] == tf["common"] and tf["relation"]["ambiguous_position"] == 0       # every common read of the queried run is a strictly shorter part of its counterpart
    # sampling day: the Costea sample name is <person>-11-<day>-<n>; the sample names of PRJEB50977 give the code and the day
    mday = re.search(r"-11-(\d+)-\d+$", qrec["sample_alias"])
    assert mday
    day = int(mday.group(1))
    third = sorted({(e["sample_alias"], e["collection_date"]) for e in ena.values() if e["study_accession"] == THIRD and code(e) == cd and day_of_third(e) == day})
    assert len(third) == 1, (s, third)
    date_gap = abs((datetime.date.fromisoformat(third[0][1]) - datetime.date.fromisoformat(qrec["collection_date"])).days)
    # the days of the curated series against the days of the code in PRJEB50977
    curated_days = sorted({int(float(m["days_from_first_collection"])) for m in meta[s] if m["days_from_first_collection"] not in ("", "NA")})
    third_days = {day_of_third(e) for e in ena.values() if e["study_accession"] == THIRD and code(e) == cd}
    # records of the two samples
    bq, bt = biosample[qrec["sample_accession"]], biosample[trec["sample_accession"]]
    sex_q = {m["gender"] for m in meta[s] if m["gender"] not in ("", "NA")}; age_q = sorted({float(m["age"]) for m in meta[s] if m["age"] not in ("", "NA")})
    assert len(sex_q) == 1
    sex_t = bt.get("host sex", "")
    sex = "not recorded" if sex_t in ("", "not provided", "missing") else ("same" if sex_t == next(iter(sex_q)) else "different")
    age_at_query = age_q[0] + day / 365.25                                                          # curated age at the first sample, plus the time elapsed until the queried sample
    age_gap = abs(float(bt["host age"]) - age_at_query)
    # same centre, same site
    assert qrec["center_name"] == trec["center_name"]
    same_site = bq["geographic location (latitude)"] == bt["geographic location (latitude)"] and bq["geographic location (longitude)"] == bt["geographic location (longitude)"]
    same_year = qrec["collection_date"][:4] == trec["collection_date"][:4]
    # what stands above the donor's other sample in the HGM set, and what remains without the runs of the subject code and the new deposits of the reads of the donor (a sample of PRJEB24041 has the name of the Costea sample with a final underscore)
    above = [a for sc, a in nr if sc > partner]
    assert r["rank_any_narrow"] is None or r["rank_any_narrow"] == 1 + len(above)
    of_code = [a for a in above if a in ena and code(ena[a]) == cd]
    relabelled = [a for a in of_code if ena[a]["sample_alias"] in RELABEL]
    of_deposit = [a for a in above if a in ena and ena[a]["study_accession"] == RE_DEPOSIT and ena[a]["sample_alias"].split("-11-")[0] == qrec["sample_alias"].split("-11-")[0]]      # the person, as the first part of the Costea name
    assert not set(of_code) & set(of_deposit)
    rank_without = None if partner == 0 else 1 + len(above) - len(of_code) - len(of_deposit)
    donors.append(dict(subject=s, days=r["days"], partner=partner, rank_stock=r["rank_any_narrow"], closest_other=top_sc,
                       other_runs_of_queried_sample=[round(min(own_sc), 2), round(max(own_sc), 2)] if own_sc else None, other_runs_of_matched_sample=[round(min(match_sc), 2), round(max(match_sc), 2)],
                       reads_queried_first_tile=rd["reads_q"], reads_matched_first_tile=rd["reads_n"], reads_found_first_tile=rd["common"], reads_contained_first_tile=rd["contained"], median_length_queried=rd["median_len_q"], median_length_matched=rd["median_len_n"], trim_filter={k: v for k, v in tf.items() if k != "subject"},
                       sampling_day=day, third_deposit_same_day=True, date_gap_days=date_gap, same_year=same_year, same_site=same_site, sex=sex, age_gap_years=round(age_gap, 2),
                       curated_days_in_third_deposit=[sum(d in third_days for d in curated_days), len(curated_days)],
                       stock_runs_above_partner=None if partner == 0 else len(above), of_the_subject_code=None if partner == 0 else len(of_code), of_a_new_deposit=None if partner == 0 else len(of_deposit), of_a_relabelled_sample=None if partner == 0 else len(relabelled),
                       rank_without_the_code=rank_without))

# the donor whose other sample is not among the 5,000 nearest runs: the samples of its code, and those of another code
no_partner = [d for d in donors if d["partner"] == 0]
assert len(no_partner) == 1
subj = next(r["subject"] for r in seven if not r["partner_any"])
hits = res[subj]["results"]; floor = min(hits.values())
qrec = ena[res[subj]["run"]]
own_runs = {a for a, e in ena.items() if e["study_accession"] == QUERY and e["sample_alias"] == qrec["sample_alias"]}
cd = code(ena[max((a for a, sc in hits.items() if sc < DUP and a in stock and a not in own_runs), key=lambda a: hits[a])])
own_days = sorted({int(float(m["days_from_first_collection"])) for m in meta[subj] if m["days_from_first_collection"] not in ("", "NA")})
q_day = no_partner[0]["sampling_day"]
same_code_other_days = [e for e in ena.values() if e["study_accession"] == THIRD and code(e) == cd and day_of_third(e) != q_day]
assert all(a not in hits for a, e in ena.items() if e["study_accession"] == THIRD and code(e) == cd and day_of_third(e) != q_day), "the code's other days are not in the list"
other_code = collections.defaultdict(list)
for a, sc in hits.items():
    e = ena.get(a)
    if e and e["study_accession"] == THIRD and code(e) not in (None, cd) and day_of_third(e) in own_days and day_of_third(e) != q_day and e["collection_date"][:7] in {x["collection_date"][:7] for x in same_code_other_days}:
        other_code[code(e)].append((sc, day_of_third(e)))
best_code = max(other_code, key=lambda k: max(x[0] for x in other_code[k]))
other_days = sorted({d for _, d in other_code[best_code]})
other_code_scores = sorted(sc for sc, _ in other_code[best_code] if sc > 0.3)
same_day_other_code = [a for a, e in ena.items() if e["study_accession"] == THIRD and code(e) == best_code and day_of_third(e) == q_day]
assert all(a not in hits or hits[a] < 0.3 for a in same_day_other_code), "the same day of the other code does not resemble the query"
print("donor whose other sample is not returned:", no_partner[0]["subject"], "other code:", other_code[best_code], "same day of it:", [hits.get(a) for a in same_day_other_code], "floor", floor)

# Fig. 1H: the threshold at the equal error rate (as many searches of an absent person accepted as donors present rejected), as in figure1_brief_communication.py
best = [float(s["best_other"]) for s in scores]; partner_scores = [float(s["partner"]) for s in scores]
cand = sorted(set([x for x in partner_scores if x > 0] + best))
far = lambda c: sum(x >= c for x in best) / len(best); frr = lambda c: sum(x < c for x in partner_scores) / len(best)
t = min(cand, key=lambda c: abs(far(c) - frr(c)))
absent = sum(x >= t for x in best)
idx7 = {position[(STUDY, d["subject"])] for d in donors}
absent_without = sum(x >= t for i, x in enumerate(best) if i not in idx7)
assert min(x for i, x in enumerate(best) if i in idx7) > t and abs(far(t) - frr(t)) < 1.0 / len(best), "the seven are above the threshold, and the two error rates are equal"
regained = [d for d in donors if d["rank_without_the_code"] == 1]
others_partner = [r["partner_any"] or 0 for r in others]; others_best = [float(scores[position[(STUDY, r["subject"])]]["best_other"]) for r in others]
rest_of_donors = [x for i, x in enumerate(best) if i not in idx7]
out = dict(
    projects={p: dict(title=studies[p]["study_title"], centre=studies[p]["center_name"].split(";")[-1], first_public=studies[p]["first_public"]) for p in (QUERY, RE_DEPOSIT, MATCH, THIRD)},
    cohort=dict(donors=len(costea), german=countries["DEU"], kazakh=countries["KAZ"], linked_stock=sum(r["rank_any_narrow"] == 1 for r in costea), failures=len(fail), failures_german=len(fail) - len(others), failures_kazakh=len(others),
                failures_kazakh_max_partner=round(max(others_partner), 3), failures_kazakh_beyond_list=sum(not x for x in others_partner), failures_kazakh_max_closest_other=round(max(others_best), 3)),
    outlier_cut=OUTLIER, largest_closest_other_of_the_490=round(max(rest_of_donors), 3), largest_partner_of_all=round(max(partner_scores), 3),
    donors=donors,
    no_partner=dict(subject=no_partner[0]["subject"], floor_of_list=round(floor, 3), other_code_days=other_days, other_code_scores=[round(x, 2) for x in other_code_scores], own_days=own_days),
    relabelled=dict(samples=len(RELABEL), reason="labelled with the code of another person; record and number of the other code"),
    sensitivity=dict(regained=len(regained), costea_as_counted=[sum(r["rank_any_narrow"] == 1 for r in costea), len(costea)], costea_with_regained=[sum(r["rank_any_narrow"] == 1 for r in costea) + len(regained), len(costea)],
                     all_as_counted=[sum(r["rank_any_narrow"] == 1 for r in rows), len(rows)], all_with_regained=[sum(r["rank_any_narrow"] == 1 for r in rows) + len(regained), len(rows)]),
    open_set=dict(threshold=round(t, 4), absent_at_threshold=absent, absent_without_the_seven=[absent_without, len(best) - 7]),
)
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps({k: out[k] for k in ("cohort", "largest_closest_other_of_the_490", "largest_partner_of_all", "no_partner", "sensitivity", "open_set")}, indent=0))
for d in donors:
    print(d)
