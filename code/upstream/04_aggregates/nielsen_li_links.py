#!/usr/bin/env python3
"""The links between the samples of Nielsen 2014 and Li 2014 (Supplementary Note 1): the numbers behind Supplementary Fig. 2 and Supplementary Tables 6 and 7. Aggregates and anonymous arrays only: no subject identifier.

    python3 -I nielsen_li_links.py OUT.json [ROWS_V4.json NIELSEN_LI.json]      # the optional files give the sensitivity of the headline rate

Inputs (local): the rows of the 831 donors (rows_v2_full.json), the results of the queries of the donors, the BioSample tables and the list of the HGM set (adultes_data_v2/: results.jsonl, biosample_runs.json, curated_sample_biosample.json, logan_gut_wgs_narrow_ids.txt.gz), the answers of the reverse queries of the Li samples (nielsen_li_validation/reverse/answers.jsonl, made with the same
service and the same request as the first queries), the curated metadata and species tables of the two cohorts, and the run table of the ENA (public).
Definitions (as in nielsen_li_overlap.py and Supplementary Note 1). The two cohorts carry the same subject_id for 44 patients. For a query of one cohort, the OTHER SAMPLE ACROSS COHORTS is the most similar run of a stool
sample of the other cohort with the same subject_id; similarities of 0.99 or more (quasi-exact re-deposits) are ignored; the runs of every sample of that subject_id, in either cohort, are neither the other sample across cohorts nor runs of other people; a run of other people is any other run
of the HGM set (291,881 human gut metagenomes); rank = 1 + the runs of other people that are strictly more similar than the other sample across cohorts.
Species profiles: dissimilarity between the relative-abundance profiles of a Nielsen and a Li sample, Bray-Curtis (primary), Aitchison (Euclidean on CLR, pseudocount = half the smallest positive abundance of the two
tables) and 1 - Spearman; a sample is identified by the rank of the most similar sample of the same subject_id among all the samples of the other cohort.
"""
import collections, csv, gzip, json, re, statistics as S, sys
import numpy as np

OUT = sys.argv[1]
P = "local"
DATA = f"{P}/adultes_data_v2"
VAL = f"{P}/nielsen_li_validation"
CRO = "local/crocodeel_valles/matrices"
DUP = 0.99
ST = {"N": "NielsenHB_2014", "L": "LiJ_2014"}

meta = {t: [m for m in csv.DictReader(open(f"{CRO}/{st}.metadata.tsv"), delimiter="\t") if m["body_site"] == "stool" and m["subject_id"]] for t, st in ST.items()}
subj_of = {t: {m["sample_id"]: m["subject_id"] for m in meta[t]} for t in ST}
shared = sorted(set(subj_of["N"].values()) & set(subj_of["L"].values()))
bs_runs = json.load(open(f"{DATA}/biosample_runs.json")); cs = json.load(open(f"{DATA}/curated_sample_biosample.json"))
runs_of = collections.defaultdict(list)
for run, v in bs_runs.items():
    runs_of[v["biosample"]].append(run)
STOCK = set(gzip.open(f"{DATA}/logan_gut_wgs_narrow_ids.txt.gz", "rt").read().split())
who = {}                                                                        # run -> (cohort tag, subject_id, sample_id)
for t, st in ST.items():
    for sample, sj in subj_of[t].items():
        for r in runs_of.get(cs.get(f"{st}|{sample}"), []):
            who[r] = (t, sj, sample)

out: dict = dict(shared_subject_ids=len(shared))

# ------------------------------------------------------------------ the sketch index, Nielsen -> Li
res = {}
for l in open(f"{DATA}/results.jsonl"):
    try:
        q = json.loads(l)
    except Exception:
        continue
    if "results" in q:
        res[(q["study"], q["subject"])] = q
rows = [r for r in json.load(open(f"{P}/rows_v2_full.json")) if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]


def score_query(results, qrun, subject, dst):
    """(similarity p of the other sample across cohorts or None, rank or None, similarity of the nearest run of other people) for one query list."""
    cross, imp = [], []
    for a, sc in results.items():
        if a == qrun or sc >= DUP:
            continue
        w = who.get(a)
        if w and w[1] == subject:
            if w[0] == dst:
                cross.append(sc)
            continue
        if a in STOCK:
            imp.append(sc)
    p = max(cross, default=0.0)
    return (p if p > 0 else None), (1 + sum(1 for x in imp if x > p) if p > 0 else None), max(imp, default=0.0)


fwd = []
for r in rows:
    if r["study"] != ST["N"]:
        continue
    q = res[(r["study"], r["subject"])]
    has_li = r["subject"] in set(subj_of["L"].values())
    p, rk, bi = score_query(q["results"], q["run"], r["subject"], "L")
    fwd.append(dict(has_li=has_li, p=p, rank=rk, best_imp=bi, within_p=r["partner_any"], within_rank1=r["rank_any_narrow"] == 1))
f35 = [x for x in fwd if x["has_li"]]
out["forward"] = dict(donors=len(fwd), with_li_sample=len(f35), in_list=sum(x["p"] is not None for x in f35), rank1=sum(x["rank"] == 1 for x in f35), top10=sum(x["rank"] is not None and x["rank"] <= 10 for x in f35),
                      within_rank1=sum(x["within_rank1"] for x in f35), median_p=S.median(x["p"] for x in f35 if x["p"] is not None), median_within_p=S.median(x["within_p"] for x in f35),
                      median_best_impostor=S.median(x["best_imp"] for x in f35), ranks=sorted(x["rank"] for x in f35 if x["rank"] is not None),
                      p=[x["p"] for x in f35], within_p=[x["within_p"] for x in f35], best_imp=[x["best_imp"] for x in f35])
out["forward"]["li_closer_than_within"] = sum(1 for x in f35 if x["p"] is not None and x["p"] > x["within_p"])

# ------------------------------------------------------------------ the sketch index, Li -> Nielsen (new queries, one per Li sample of a shared patient)
rev = []
for l in open(f"{VAL}/reverse/answers.jsonl"):
    j = json.loads(l)
    if "results" not in j:
        continue
    w = who[j["run"]]
    p, rk, bi = score_query(j["results"], j["run"], w[1], "N")
    rev.append(dict(subject=w[1], p=p, rank=rk, best_imp=bi))
out["reverse"] = dict(queries=len(rev), patients=len({x["subject"] for x in rev}), in_list=sum(x["p"] is not None for x in rev), rank1=sum(x["rank"] == 1 for x in rev), top10=sum(x["rank"] is not None and x["rank"] <= 10 for x in rev),
                      median_p=S.median(x["p"] for x in rev if x["p"] is not None), median_best_impostor=S.median(x["best_imp"] for x in rev), ranks=sorted(x["rank"] for x in rev if x["rank"] is not None),
                      p=[x["p"] for x in rev], best_imp=[x["best_imp"] for x in rev])
donor_subjects = {r["subject"] for r in rows if r["study"] == ST["N"]}
f_rank = {r["subject"]: x["rank"] for r, x in zip([r for r in rows if r["study"] == ST["N"]], fwd)}
r_rank = collections.defaultdict(list)
for x in rev:
    r_rank[x["subject"]].append(x["rank"])
both = [s for s in donor_subjects if s in r_rank and f_rank.get(s) == 1 and min(v for v in r_rank[s] if v is not None) == 1]
out["reverse"]["mutual_rank1_among_donor_patients"] = len(both)
out["reverse"]["donor_patients_queried"] = sum(1 for s in donor_subjects if s in r_rank)
out["reverse"]["rank1_among_donor_patients"] = sum(1 for s in donor_subjects if s in r_rank and min(v for v in r_rank[s] if v is not None) == 1)

# ------------------------------------------------------------------ species profiles
def table(st):
    rws = [l.rstrip("\n").split("\t") for l in open(f"{CRO}/{st}.species_relab.tsv") if l.strip()]
    return rws[0][1:], {r[0]: np.array(r[1:], float) for r in rws[1:]}


names, tab = {}, {}
for t, st in ST.items():
    names[t], tab[t] = table(st)
species = sorted(set(tab["N"]) | set(tab["L"]))
X = {t: np.array([tab[t].get(sp, np.zeros(len(names[t]))) for sp in species]).T for t in ST}                    # samples x species
X = {t: X[t] / X[t].sum(1, keepdims=True) for t in ST}                                                          # profiles sum to 1
keepN = [i for i, s in enumerate(names["N"]) if s in subj_of["N"]]; keepL = [i for i, s in enumerate(names["L"]) if s in subj_of["L"]]
XN, XL = X["N"][keepN], X["L"][keepL]
sN = np.array([subj_of["N"][names["N"][i]] for i in keepN]); sL = np.array([subj_of["L"][names["L"][i]] for i in keepL])


def bc(A, B):
    return np.array([[np.abs(a - b).sum() / (a + b).sum() for b in B] for a in A])


def clr(A, pc):
    L = np.log(A + pc)
    return L - L.mean(1, keepdims=True)


def ranks_of(A):
    from scipy.stats import rankdata
    return np.array([rankdata(a) for a in A])


pc = 0.5 * min(XN[XN > 0].min(), XL[XL > 0].min())
CN, CL = clr(XN, pc), clr(XL, pc)
RN, RL = ranks_of(XN), ranks_of(XL)
DIST = {
    "Bray-Curtis": bc(XN, XL),
    "Aitchison": np.sqrt(((CN[:, None, :] - CL[None, :, :]) ** 2).sum(2)),
    "Spearman": 1 - np.corrcoef(np.vstack([RN, RL]))[: len(RN), len(RN):],
}
same = (sN[:, None] == sL[None, :])
shared_set = set(shared)
spec = {}
for m, Dm in DIST.items():
    # Li sample of a shared patient -> the 396 Nielsen samples
    r_l = []
    for j in range(len(sL)):
        if sL[j] not in shared_set:
            continue
        own = same[:, j]
        r_l.append(1 + int((Dm[~own, j] < Dm[own, j].min()).sum()))
    r_n = []
    for i in range(len(sN)):
        if sN[i] not in shared_set:
            continue
        own = same[i, :]
        r_n.append(1 + int((Dm[i, ~own] < Dm[i, own].min()).sum()))
    pos, neg = Dm[same], Dm[~same]
    from scipy.stats import rankdata
    allv = np.concatenate([pos, neg]); rk = rankdata(allv)
    auc = (rk[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))
    auc = 1 - auc                                                                           # a small distance means the same patient
    spec[m] = dict(li_to_nielsen_n=len(r_l), li_to_nielsen_rank1=sum(x == 1 for x in r_l), li_to_nielsen_top10=sum(x <= 10 for x in r_l), nielsen_to_li_n=len(r_n), nielsen_to_li_rank1=sum(x == 1 for x in r_n),
                   auc=float(auc), matched_pairs=int(same.sum()), median_matched=float(np.median(pos)), median_other_cross=float(np.median(neg)))
out["species"] = spec
# the point of comparison inside Nielsen with Bray-Curtis: two samples of one patient, and of two patients
BCN = bc(XN, XN); iu = np.triu_indices(len(sN), 1)
sameN = (sN[:, None] == sN[None, :])[iu]; dN = BCN[iu]
out["species"]["within_nielsen_bc"] = dict(same_patient_pairs=int(sameN.sum()), median_same_patient=float(np.median(dN[sameN])), median_different_patients=float(np.median(dN[~sameN])),
                                           q_same=[float(x) for x in np.percentile(dN[sameN], [5, 25, 50, 75, 95])], q_diff=[float(x) for x in np.percentile(dN[~sameN], [5, 25, 50, 75, 95])])
Dbc = DIST["Bray-Curtis"]
out["species"]["bc_groups"] = dict(matched=[float(x) for x in Dbc[same]], other_cross_quantiles=[float(x) for x in np.percentile(Dbc[~same], [5, 25, 50, 75, 95])],
                                   within_same=[float(x) for x in dN[sameN]], within_diff_quantiles=[float(x) for x in np.percentile(dN[~sameN], [5, 25, 50, 75, 95])])

# ------------------------------------------------------------------ the curated metadata of the shared patients (sex, age, BMI)
cmd = {}
with open(f"{P}/cmd_meta/cmd_meta.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["study_name"] in ST.values() and r["body_site"] == "stool":
            cmd.setdefault((r["study_name"], r["subject_id"]), []).append(r)
num = lambda x: float(x) if x not in ("", "NA", None) else None
conc = collections.Counter()
for sj in shared:
    n0 = min(cmd[(ST["N"], sj)], key=lambda r: (num(r["days_from_first_collection"]) or 0)); l0 = cmd[(ST["L"], sj)][0]
    conc["patients"] += 1
    conc["sex"] += n0["gender"] == l0["gender"]
    a, b = num(n0["age"]), num(l0["age"]); conc["age_le1"] += a is not None and b is not None and abs(a - b) <= 1; conc["age_equal"] += a is not None and a == b
    a, b = num(n0["BMI"]), num(l0["BMI"]); conc["bmi_known"] += a is not None and b is not None; conc["bmi_le05"] += a is not None and b is not None and abs(a - b) <= 0.5
    conc["bmi_equal"] += a is not None and b is not None and abs(a - b) < 1e-6
out["metadata_cmd"] = dict(conc)
point = lambda r: re.split(r"[-_]", r["sample_id"])[-1]                                      # the time point, as the last part of the sample name
two_one = sum(sorted(point(r) for r in cmd[(ST["N"], sj)]) == ["0", "2"] and [point(r) for r in cmd[(ST["L"], sj)]] == ["1"] for sj in shared)
out["time_points"] = dict(patients=len(shared), two_nielsen_one_li=two_one, others=len(shared) - two_one)

# ------------------------------------------------------------------ the identifiers are not a clean unit of person (curated metadata only)
byN = collections.defaultdict(list); byL = collections.defaultdict(list)
for (stn, sj), v in cmd.items():
    (byN if stn == ST["N"] else byL)[sj] = v
incons = [sj for sj, v in byN.items() if len(v) > 1 and (len({r["gender"] for r in v}) > 1 or (lambda ag: bool(ag) and max(ag) - min(ag) > 2)([num(r["age"]) for r in v if num(r["age"]) is not None]))]
rec = lambda r: (r["gender"], num(r["age"]), round(num(r["BMI"]), 6)) if num(r["age"]) is not None and num(r["BMI"]) is not None else None
firstN = {sj: rec(v[0]) for sj, v in byN.items()}; firstL = {sj: rec(v[0]) for sj, v in byL.items()}
groups = collections.defaultdict(set)
for sj, k in firstN.items():
    if k:
        groups[k].add(sj)
dup = [sorted(v) for v in groups.values() if len(v) > 1]
li_rec = collections.defaultdict(set)
for sj, k in firstL.items():
    if k:
        li_rec[k].add(sj)
other_code = [sj for sj, k in firstN.items() if k and sj not in byL and k in li_rec]
out["curation"] = dict(shared_with_inconsistent_records=sum(1 for sj in incons if sj in shared), nielsen_subjects_with_two_samples=sum(len(v) > 1 for v in byN.values()),
                       identical_record_groups=len(dup), identical_record_groups_two_series=sum(1 for g in dup if {x.split("_")[0] for x in g} == {"O2", "V1"}), identical_record_groups_healthy=sum(1 for g in dup if all(x.startswith("MH") for x in g)),
                       subjects_without_same_id_li_but_identical_record_under_another_code=len(other_code), of_which_in_series_O2_V1=sum(1 for x in other_code if x.startswith(("O2_", "V1_"))), of_which_V1_with_O2_counterpart=sum(1 for x in other_code if x.startswith("V1_") and all(y.startswith("O2_") for y in li_rec[firstN[x]])))
# ------------------------------------------------------------------ how the two cohorts were sequenced (ENA, public)
sample_of_run = {r: x[2] for r, x in who.items()}
tech = {}
for t, proj in (("N", "PRJEB1220"), ("L", "PRJEB5224")):
    rws = list(csv.DictReader(open(f"{VAL}/ena_{proj}.tsv"), delimiter="\t"))
    paired = [r for r in rws if r["library_layout"] == "PAIRED"]
    gb = collections.defaultdict(float)
    for r in rws:
        gb[sample_of_run.get(r["run_accession"], r["sample_accession"])] += int(r["base_count"]) / 1e9
    per_sample = list(gb.values())
    mate = [int(r["base_count"]) / int(r["read_count"]) / 2 for r in paired if int(r["read_count"]) > 0]
    tech[t] = dict(runs=len(rws), samples=len(per_sample), paired_runs=len(paired), models_paired=dict(collections.Counter(r["instrument_model"] for r in paired)), models_all=dict(collections.Counter(r["instrument_model"] for r in rws)),
                   mate_median=float(np.median(mate)), mate_p5_p95=[float(x) for x in np.percentile(mate, [5, 95])], gbp_median=float(np.median(per_sample)), gbp_iqr=[float(x) for x in np.percentile(per_sample, [25, 75])],
                   libraries=len({r["library_name"] for r in rws}), first_created=sorted({r["first_created"][:7] for r in rws}), centers=sorted({r["center_name"] for r in rws}), runs_in_stock=sum(1 for r in rws if r["run_accession"] in STOCK))
out["sequencing"] = tech
out["libraries_shared"] = len({r["library_name"] for r in csv.DictReader(open(f"{VAL}/ena_PRJEB1220.tsv"), delimiter="\t")} & {r["library_name"] for r in csv.DictReader(open(f"{VAL}/ena_PRJEB5224.tsv"), delimiter="\t")})
# ------------------------------------------------------------------ who outranks the other sample of a Nielsen donor (the first counting of Supplementary Note 1: a sample of Li is another person)
cat_ = collections.Counter()
for r in rows:
    if r["study"] != ST["N"]:
        continue
    q = res[(r["study"], r["subject"])]; sj = r["subject"]; same_s, imp0, top = [], [], (0.0, None, None)
    for a, sc in q["results"].items():
        if sc >= DUP:
            continue
        w = who.get(a)
        if w and w[0] == "N" and w[1] == sj:
            if w[2] != r["v2s"]:
                same_s.append(sc)
            continue
        imp0.append((sc, a))
        if sc > top[0]:
            top = (sc, w[0] if w else "unmapped", w[1] if w else None)
    p_old = max(same_s, default=0.0)
    if p_old == 0:
        cat_["not_found"] += 1; continue
    if 1 + sum(1 for sc, a in imp0 if a in STOCK and sc > p_old) == 1:
        cat_["first"] += 1; continue
    cat_["li_same_identifier" if top[1] == "L" and top[2] == sj else "li_other_identifier" if top[1] == "L" else "nielsen_other_person" if top[1] == "N" else "other_run"] += 1
out["outranked_by"] = dict(cat_)

# ------------------------------------------------------------------ the queried runs of the donors of the two cohorts (the largest indexed run of the last sample): layout and size from the ENA
ena_runs = {}
for proj in ("PRJEB1220", "PRJEB5224"):
    for r_ in csv.DictReader(open(f"{VAL}/ena_{proj}.tsv"), delimiter="\t"):
        ena_runs[r_["run_accession"]] = (r_["library_layout"], r_["instrument_model"], int(r_["base_count"]) / 1e9)
qr = {t: [ena_runs[res[(r["study"], r["subject"])]["run"]] for r in rows if r["study"] == ST[t]] for t in ST}
out["queried_runs"] = {t: dict(n=len(v), paired=sum(x[0] == "PAIRED" for x in v), models=dict(collections.Counter(x[1] for x in v)), gbp_median=float(np.median([x[2] for x in v])), below_half_gbp=sum(x[2] < 0.5 for x in v)) for t, v in qr.items()}
rows_kept = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None
if rows_kept:
    k_stock = sum(1 for r in rows_kept if r["rank_any_narrow"] == 1)
    nl_ = json.load(open(sys.argv[3]))
    asc = k_stock + nl_["nielsen_rank1_before"] + nl_["li_rank1_before"]; mrg = k_stock + nl_["nielsen_rank1_after"] + nl_["li_rank1_after"]; tot = len(rows_kept) + nl_["nielsen_donors"] + nl_["li_donors"]
    out["sensitivity"] = dict(kept_first=k_stock, kept=len(rows_kept), with_two_cohorts_as_counted=dict(first=asc, donors=tot), with_two_cohorts_merged=dict(first=mrg, donors=tot))
out["max_similarity_to_same_identifier_sample"] = max([x for x in out["forward"]["p"] if x is not None] + [x for x in out["reverse"]["p"] if x is not None])

json.dump(out, open(OUT, "w"), indent=1)


def show(d, ind=0):
    for k, v in d.items():
        if isinstance(v, dict):
            print(" " * ind + f"{k}:"); show(v, ind + 2)
        elif isinstance(v, list) and len(v) > 12:
            print(" " * ind + f"{k}: [{len(v)} values]")
        else:
            print(" " * ind + f"{k}: {round(v, 3) if isinstance(v, float) else v}")


show(out)
