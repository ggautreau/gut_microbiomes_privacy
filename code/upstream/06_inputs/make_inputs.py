#!/usr/bin/env python3
"""Writes the inputs of the figures and tables (data/ of the reproduction package) from the tables of the upstream steps.

    python3 -I make_inputs.py ROWS_V4.json QUERY_GBP.json CROCODEEL_DIR ADULTES_DIR NIELSEN_LI.conta.tsv JOIN_FINAL_TABLE.csv PAIRS_SERVICE.json RUN_PROJECTS.json STAT_PER_RUN.json OUT_DIR [DONORS_LEFT_OUT.json [DONORS_EXCLUDED.json]]

ROWS_V4.json        the 497 donors of the 7 cohorts kept, before the two rules (rows_with_partner.py, make_rows_v4.py, then the curated sample identifiers of HMP IBD: hmp_curated_ids.py), in the
                    order of the analysis; those of DONORS_EXCLUDED.json go to data/donors_excluded.json, the others to data/donors.json
QUERY_GBP.json      'study|subject' -> bases of the queried run in Gbp (ENA base_count)
CROCODEEL_DIR       crocodeel_valles/: matrices/<cohort>.metadata.tsv (curatedMetagenomicData) and results/<cohort>.conta.tsv (CroCoDeEL 1.3.1) of the 11 cohorts
ADULTES_DIR         adultes/: pairs.json, pairs_hmp2.json (run accessions of the two samples of each donor)
NIELSEN_LI.conta.tsv  CroCoDeEL on the species tables of Nielsen and Li merged (crocodeel_nielsen_li.py)
JOIN_FINAL_TABLE.csv  human-DNA check of the 558 donors (all558_join.py)
PAIRS_SERVICE.json  the 558 pairs of the human-DNA check: study, subject, query run `q`, run of the other sample `partner` (all558_pairs.py)
RUN_PROJECTS.json   run accession -> BioProject (run_projects.py)
STAT_PER_RUN.json   run accession -> human hashes of the sketch and human hits of the NCBI taxonomic analysis, STAT (local/humancheck/all558/stat_per_run.json; R. Faure)

Subjects and samples keep their identifiers of curatedMetagenomicData. The accessions of the two runs of each donor are in data/accessions.tsv.
What is kept: the cohort, the ranks, the similarities, the days, the cohort label of the age (adult or adult by design), the number of runs of the queried sample and the other samples at least 30 days from it.
What is left out: age, sex, country, the other run accessions, BioSamples, the header of the CroCoDeEL files (it records the user name and the host of the run), the species lists of the events.
DONORS_LEFT_OUT.json (optional) the 183 donors of the 4 cohorts left out (donors_left_out.py), with their query sample, first sample and most similar other sample; written as data/donors_left_out.json,
for the samples studied of Supplementary Note 1.
DONORS_EXCLUDED.json (optional) the donors that the rules of contamination and re-deposit remove (exclude_donors.py), written as data/donors_excluded.json with the rule.
"""
import csv, json, os, sys

ROWS, GBP, CRO, ADU, NLC, JOIN, PAIRS, PROJ, STAT, OUT = sys.argv[1:11]
DLO = sys.argv[11] if len(sys.argv) > 11 else None
DEX = sys.argv[12] if len(sys.argv) > 12 else None
dex = {(e["study"], e["subject"]): e for e in json.load(open(DEX))} if DEX else {}
COHORTS = ["MehtaRS_2018", "CosteaPI_2017", "HMP_2019_ibdmdb", "HansenLBS_2018", "Bengtsson-PalmeJ_2015", "RaymondF_2016", "Heitz-BuschartA_2016", "PehrssonE_2016", "NielsenHB_2014", "LiJ_2014", "VincentC_2016"]
GROUP = lambda st: "NL" if st in ("NielsenHB_2014", "LiJ_2014") else st                 # Nielsen and Li share their subject identifiers
runs = lambda s: [a for a in s.split(";") if a not in ("", "NA")]

known_subjects, known_samples = set(), set()                                        # (cohort or the pair Nielsen and Li, subject) and (cohort, sample) of the rows and of the CroCoDeEL tables


rows = json.load(open(ROWS)); gbp = json.load(open(GBP))
pairs = {(p["study"], p["subject"]): p for p in json.load(open(f"{ADU}/pairs.json")) + [q for q in json.load(open(f"{ADU}/pairs_hmp2.json")) if q["study"] == "HMP_2019_ibdmdb"]}
KEEP = ["study", "age_class", "q_lib", "rank_any", "rank_any_narrow", "rank_far_narrow", "rank_v1_narrow", "far_possible", "days", "partner_any", "max_imp", "top_kind"]
out_rows, out_gbp, out_ex = [], {}, []
for r in rows:                                                                       # the 497 donors, in the order of the analysis
    st = r["study"]
    new = {"study": st, "subject": r["subject"], "v1s": r["v1s"], "v2s": r["v2s"], "best_any": r["best_any"]}
    new["eligible"] = list(r["eligible"])                                           # the other samples at least 30 days from the query sample
    known_subjects.add((GROUP(st), r["subject"])); known_samples.update((st, x) for x in [r["v1s"], r["v2s"], r["best_any"]] + list(r["eligible"]))
    new.update({k: r[k] for k in KEEP if k != "study"})
    new["nruns"] = len(runs(pairs[(st, r["subject"])]["v2_acc"]))
    e = dex.get((st, r["subject"]))
    if e:                                                                          # removed by a rule: a file of its own
        new.update(rule=e["rule"], why=e.get("why"), removed_samples=list(e["removed_samples"]), rates=e.get("rates"), nearest_other=e.get("nearest_other"))
        out_ex.append(new)
        continue
    out_rows.append(new)
    out_gbp[f"{st}|{new['subject']}"] = round(gbp[f"{st}|{r['subject']}"], 3)                     # to 0.001 Gbp: the exact base count of a run is a public field of the ENA
os.makedirs(f"{OUT}/crocodeel", exist_ok=True); os.makedirs(f"{OUT}/human_dna", exist_ok=True)
json.dump(out_rows, open(f"{OUT}/donors.json", "w")); json.dump(out_gbp, open(f"{OUT}/query_gbp.json", "w"))
assert len(out_ex) == len(dex), "every donor removed by a rule must be one of the 497 donors"
if DEX:
    json.dump(out_ex, open(f"{OUT}/donors_excluded.json", "w"))

n_ev = {}
for st in COHORTS:
    with open(f"{CRO}/matrices/{st}.metadata.tsv", newline="") as fh:
        md = list(csv.DictReader(fh, delimiter="\t"))
    with open(f"{OUT}/crocodeel/{st}.metadata.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["sample_id", "subject_id", "days_from_first_collection"])
        for m in md:
            w.writerow([m["sample_id"], m["subject_id"], m["days_from_first_collection"]])
            known_samples.add((st, m["sample_id"])); known_subjects.add((GROUP(st), m["subject_id"]))
    ev = [l for l in open(f"{CRO}/results/{st}.conta.tsv") if not l.startswith("#")]
    rd = list(csv.DictReader(ev, delimiter="\t"))
    with open(f"{OUT}/crocodeel/{st}.conta.tsv", "w", newline="") as fh:
        fh.write("# CroCoDeEL 1.3.1, default settings\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["source", "target", "rate", "probability"])
        for d in rd:
            w.writerow([d["source"], d["target"], d["rate"], d["probability"]])
            known_samples.update({(st, d["source"]), (st, d["target"])})
    n_ev[st] = len(rd)

ev = [l for l in open(NLC) if not l.startswith("#")]
rd = list(csv.DictReader(ev, delimiter="\t"))
tag = {"N": "NielsenHB_2014", "L": "LiJ_2014"}
with open(f"{OUT}/crocodeel/nielsen_li_merged.conta.tsv", "w", newline="") as fh:
    fh.write("# CroCoDeEL 1.3.1, default settings, on the species tables of Nielsen and Li merged; sample names are N| or L| and the sample identifier\n")
    w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["source", "target", "rate", "probability"])
    for d in rd:
        a, b = d["source"], d["target"]
        assert a[1] == "|" and b[1] == "|"
        w.writerow([a, b, d["rate"], d["probability"]])
        known_samples.update({(tag[a[0]], a[2:]), (tag[b[0]], b[2:])})

with open(JOIN, newline="") as fh:                                                 # human-DNA check: the columns of supplementary_figure3.py
    jt = list(csv.DictReader(fh))
with open(f"{OUT}/human_dna/donors_human_share.csv", "w", newline="") as fh:
    stat = json.load(open(STAT))
    w = csv.writer(fh, lineterminator="\n"); w.writerow(["study", "rank1", "share_q", "share_p", "human_hashes_q", "human_hashes_p", "stat_q", "stat_p"])
    pc = lambda run: repr(100.0 * stat[run]["stat_human_hits"] / stat[run]["stat_species_hits"])                 # the human share of the NCBI STAT analysis, in percent
    for r in jt:
        w.writerow([r["study"], r["rank1"], r["share_q"], r["share_p"], stat[r["q"]]["human_hashes"], stat[r["partner"]]["human_hashes"], pc(r["q"]), pc(r["partner"])])
        assert abs(100.0 * stat[r["q"]]["human_hashes"] / stat[r["q"]]["hashes"] - float(r["share_q"])) < 1e-6 and abs(100.0 * stat[r["partner"]]["human_hashes"] / stat[r["partner"]]["hashes"] - float(r["share_p"])) < 1e-6
known_samples = {k for k in known_samples if k[1] not in (None, "", "NA")}                 # a donor whose other sample is beyond the 5,000 nearest runs names no such sample
# the accessions: one row per donor of the 558, in the order of the human-DNA check, with the subject of donors.json and the cohort; run accessions and BioProjects are public identifiers of the ENA
pairs_ = json.load(open(PAIRS)); proj = json.load(open(PROJ))
assert [(p["study"], p["q"], p["partner"]) for p in pairs_] == [(r["study"], r["q"], r["partner"]) for r in jt], "the pairs and the human-DNA table do not list the same donors in the same order"
analysed = {(r["study"], r["subject"]) for r in rows if (r["study"], r["subject"]) not in dex}        # the donors that the two rules remove are not in the analysis
assert analysed <= {(p["study"], p["subject"]) for p in pairs_}, "a donor of the analysis has no pair"
with open(f"{OUT}/accessions.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["subject", "cohort", "in_analysis", "query_run", "query_bioproject", "other_run", "other_bioproject", "other_in_5000"])
    for p in pairs_:
        w.writerow([p["subject"], p["study"], "yes" if (p["study"], p["subject"]) in analysed else "no", p["q"], proj[p["q"]], p["partner"], proj[p["partner"]], "yes" if p["in_results"] else "no"])
assert all((GROUP(p["study"]), p["subject"]) in known_subjects for p in pairs_), "a donor of the pairs is in no row and in no metadata table"
if DLO:                                                                              # the samples and subjects of these cohorts are those of the CroCoDeEL tables
    lo = json.load(open(DLO))
    out_lo = [{"study": r["study"], "subject": r["subject"], "v1s": r["v1s"], "v2s": r["v2s"], "best_any": r["best_any"], "age_class": r["age_class"], "q_lib": r["q_lib"]} for r in lo]
    assert all((GROUP(r["study"]), r["subject"]) in known_subjects and all(x in (None, "", "NA") or (r["study"], x) in known_samples for x in (r["v1s"], r["v2s"], r["best_any"])) for r in lo), "a sample or a subject of the donors left out is in no metadata table"
    json.dump(out_lo, open(f"{OUT}/donors_left_out.json", "w"))
print(len(out_rows), "donors;", len(pairs_), "rows of accessions;", len(known_subjects), "subjects;", len(known_samples), "samples;", "events:", n_ev, "; human-DNA table:", len(jt), "rows")
