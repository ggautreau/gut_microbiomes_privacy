#!/usr/bin/env python3
"""Writes supplementary.tex (Supplementary Notes 1 to 4 and Supplementary Tables 1 to 18, with the six Supplementary Figures drawn by their scripts) for the Brief Communication.

    python3 -I build_supplementary.py ROWS.json GBP.json OUT_DIR FIG2_STATS.json PEHRSSON.json EXTRAPOLATION.json FIG3_VALUES.json QUADRATIC.json NIELSEN_LI.json HUMAN_DNA.json SELECTION.json LINKS.json COSTEA.json HUMAN_GUT.json SENSITIVITY.json FIG3_STAT.json EVT.json EVT_SIMILARITIES.json BOOTSTRAP.json EVT_RESAMPLING.json EVT_REMOVED.json

ROWS.json     data/donors.json (one row per donor of the 7 cohorts kept: ranks, similarities, days, number of runs of the queried sample)
GBP.json      'study|subject' -> bases of the queried run in Gbp, from the ENA (read_run, base_count): the service sketches one run per query
OUT_DIR       supplementary.tex is written there; the title of the paper is read from OUT_DIR/title.txt (the \title{} of main.tex).
FIG2_STATS    stats of supplementary_figure4.py (Supplementary Note 1, Table 8); PEHRSSON.json: output of pehrsson_lanes.py (Supplementary Note 1, Table 5)
FIG3_VALUES   <prefix>_values.json of supplementary_figure5.py (Supplementary Fig. 5, Supplementary Note 3 and Table 16)
QUADRATIC     <prefix>_quadratic.json of figure1_brief_communication.py (the extrapolation of Fig. 1B: Supplementary Note 3 and Table 15)
NIELSEN_LI    output of nielsen_li_overlap.py (Supplementary Note 1: the two cohorts left out because they carry the same subject identifiers)
COSTEA        output of costea_same_people.py (Supplementary Note 1 and Table 9: the reads of seven Costea donors deposited under another subject code)
HUMAN_GUT     summary of build_stock.py (Supplementary Table 3); SENSITIVITY: output of stock_sensitivity.py (Supplementary Table 4)
FIG3_STAT     <prefix>_values.json of supplementary_figure3.py (Supplementary Note 1: the human share of the NCBI taxonomic analysis against ours, Supplementary Fig. 3D and E)
EVT           output of evt_model.py (Supplementary Note 3 and Tables 17 and 18: the extreme-value model of the similarity tail, the fourth model of Fig. 1B); EVT_SIMILARITIES: the checks of evt_similarities.py
EXTRAPOLATION output of extrapolation_tables.py (Supplementary Note 3 and Table 16: the small pools, the upper bound and the thinning test)
BOOTSTRAP     output of extrapolation_bootstrap.py (Supplementary Note 3 and Supplementary Fig. 6: adaptive doubling bootstrap intervals for all three models)
EVT_RESAMPLING  output of evt_resampling.py (Supplementary Fig. 6: the band of the extreme-value model and its bounds as the number of resamples doubles)
EVT_REMOVED   output of evt_removed_donors.py (last row of Supplementary Table 17: the extreme-value model with the donors that the two rules remove)
HUMAN_DNA     output of human_dna_exclusion.py (Supplementary Note 1: the human DNA in the samples of the cohorts kept and of Vincent)
SELECTION     output of cohort_selection.py (Supplementary Note 1, Tables 1 and 2: the selection of the cohorts and of the donors)
LINKS         output of nielsen_li_links.py (Supplementary Note 1, Supplementary Fig. 2, Tables 6 and 7)
Read beside the arguments: donors_excluded.json and open_set_scores.csv in the folder of ROWS.json, figure1_open_set.json in the folder of QUADRATIC.
Every number of the Notes and of the tables is computed here.
"""
import collections, csv, json, math, os, re, statistics as S, sys
import numpy as np
from scipy.stats import mannwhitneyu

ROWS, GBP, OVL = sys.argv[1], sys.argv[2], sys.argv[3]
FIG2 = sys.argv[4] if len(sys.argv) > 4 else None                               # stats of supplementary_figure4.py (optional)
PEH = sys.argv[5] if len(sys.argv) > 5 else None                                # pehrsson_lanes.py: the evidence for leaving Pehrsson out (Supplementary Note 1, Table 5)
EXT = sys.argv[6] if len(sys.argv) > 6 else None                                # extrapolation_tables.py: Supplementary Note 3 and Table 16
FIG3 = sys.argv[7] if len(sys.argv) > 7 else None                               # supplementary_figure5.py: predictions of each model at 8e9 (Table 16)
NL = sys.argv[9] if len(sys.argv) > 9 else None                                 # nielsen_li_overlap.py: Supplementary Note 1
VIN = sys.argv[10] if len(sys.argv) > 10 else None                              # aggregates of the human-DNA check behind the exclusion of Vincent (Supplementary Note 1)
LNK = sys.argv[12] if len(sys.argv) > 12 else None                              # nielsen_li_links.py: Supplementary Note 1, Fig. 2, Tables 6 and 7
SEL = sys.argv[11] if len(sys.argv) > 11 else None                              # cohort_selection.py: how the cohorts were chosen (Supplementary Note 1, Tables 1 and 2)
CST = sys.argv[13] if len(sys.argv) > 13 else None                              # costea_same_people.py: Supplementary Note 1 and Table 9
STK = sys.argv[14] if len(sys.argv) > 14 else None                              # build_stock.py: Supplementary Table 3
STS = sys.argv[15] if len(sys.argv) > 15 else None                              # stock_sensitivity.py: Supplementary Table 4
HDV = sys.argv[16] if len(sys.argv) > 16 else None                              # supplementary_figure3.py: the comparison with the human share of the NCBI taxonomic analysis (Supplementary Note 1, Fig. 3D and E)
EVT = sys.argv[17] if len(sys.argv) > 17 else None                              # evt_model.py: the extreme-value model of the similarity tail (Supplementary Note 3, Tables 17 and 18)
EVS = sys.argv[18] if len(sys.argv) > 18 else None                              # evt_similarities.py: the checks of the similarities computed from the hypervector database
BOOT = sys.argv[19] if len(sys.argv) > 19 else None                             # extrapolation_bootstrap.py: adaptive doubling bootstrap intervals for the three curve-fitted models (Supplementary Note 3, Fig. 6)
EVR = sys.argv[20] if len(sys.argv) > 20 else None                              # evt_resampling.py: the resampling band of the extreme-value model (Supplementary Fig. 6D to F)
EVX = sys.argv[21] if len(sys.argv) > 21 else None                              # evt_removed_donors.py: the extreme-value model with the donors that the two rules remove (last row of Supplementary Table 17)
INDEX_RUNS = "4,791,208"
QUAD = sys.argv[8] if len(sys.argv) > 8 else None                               # figure1_brief_communication.py: the straight line of Fig. 1B with its interval (Note 3, Table 15)

EXCLUDED = "PehrssonE_2016"                                                    # left out of the analysis (Supplementary Note 1)
SMALL = 0.5                                                                    # Gbp: below this a query is called small

rows = [r for r in json.load(open(ROWS)) if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]
_ex_path = os.path.join(os.path.dirname(os.path.abspath(ROWS)), "donors_excluded.json")
rows_ex = json.load(open(_ex_path)) if os.path.exists(_ex_path) else []           # the donors of the cohorts kept removed for contamination or re-deposit (exclude_donors.py)
gbp = json.load(open(GBP))
for r in rows:
    r["g"] = gbp[f"{r['study']}|{r['subject']}"]
    r["ok"] = r["rank_any"] == 1
    r["absent"] = r["rank_any"] is None
    r["out"] = (not r["ok"]) and (not r["absent"])

med = lambda L: S.median(L) if L else float("nan")
pct = lambda k, n: f"{100 * k / n:.0f}" if n else "--"
pct_rng = lambda kn: (lambda v: v[0] if v[0] == v[-1] else f"{v[0]}--{v[-1]}")(sorted((pct(k, n) for k, n in kn), key=float))   # the range of the shares of several classes, in whole percents
pct1 = lambda k, n: f"{100 * k / n:.1f}"
q = lambda x: f"{x:.2f}"
rp = lambda x: f"{round(100 * float(f'{x:.2f}')):.0f}"                                   # a rate (a fraction) as a whole percent, from its value to two decimals as the text had it
N = len(rows)
n_stock = sum(1 for r in rows if r["rank_any_narrow"] == 1)                      # the same donors as those recovered among all the runs of the index (the measure of Fig. 1)
ok = [r for r in rows if r["ok"]]; out = [r for r in rows if r["out"]]; ab = [r for r in rows if r["absent"]]
kinds = lambda L: collections.Counter(r["top_kind"] for r in L)
ko, kok, kab = kinds(out), kinds(ok), kinds(ab)
# Open-set test of Fig. 1H: the nearest run of other people of a donor is the run set against the other sample of an absent donor. figure1_open_set.json sits beside figure1_quadratic.json, open_set_scores.csv beside ROWS.
THR = json.load(open(os.path.join(os.path.dirname(QUAD), "figure1_open_set.json")))
OSS = list(csv.DictReader(open(os.path.join(os.path.dirname(ROWS), "open_set_scores.csv"))))
assert len(OSS) == len(rows) and all(abs(float(s["partner"]) - (r["partner_any"] or 0)) < 1e-4 and float(s["best_other"]) <= (r["max_imp"] or 0) + 1e-4 for s, r in zip(OSS, rows)), "open_set_scores.csv is not in the order of the donors"
accepted = [r for s, r in zip(OSS, rows) if float(s["best_other"]) >= THR["threshold_equal_error"]]
below = [r for s, r in zip(OSS, rows) if float(s["best_other"]) < THR["threshold_equal_error"]]
assert len(accepted) == THR["absent_person_at_or_above"], "the absent persons accepted do not match the numbers of Fig. 1H"
ka, kb = kinds(accepted), kinds(below)
mw_close = mannwhitneyu([r["max_imp"] for r in out], [r["max_imp"] for r in ok], alternative="two-sided", method="asymptotic")           # nearest run of other people, outranked against recovered donors (Note 2)
p_close, u_close = mw_close.pvalue, mw_close.statistic                             # the statistic U is that of the outranked donors
assert u_close == int(u_close), "the text prints U as a whole number"

# ---- Table 10: the three groups
def t1_row(name, L, k, sim_partner):
    n = len(L)
    miqr = lambda v: "{:.2f} ({:.2f}--{:.2f})".format(med(v), *np.percentile(v, [25, 75]))          # median (first quartile to third quartile)
    sp = miqr([r['partner_any'] for r in L]) if sim_partner else "--"
    return f"{name} & {n} & {pct1(n, N)} & {sp} & {miqr([r['max_imp'] for r in L])} & {k['same_cohort_other_person']} & {k['other_cohort']} & {k['unmapped']} \\\\"
T1 = "\n".join([t1_row("Recovered (rank 1)", ok, kok, True), t1_row("Rank 2 or more", out, ko, True), t1_row("Not in 5,000", ab, kab, False)])

# ---- Table 11: cohorts
by = collections.defaultdict(list)
for r in rows:
    by[r["study"]].append(r)
order = sorted(by, key=lambda s: sum(r["ok"] for r in by[s]) / len(by[s]))
gfmt = lambda x: f"{x:.2g}"
g2 = lambda x: (lambda t: t + ".0" if "." not in t and x < 10 else t)(f"{x:.2g}")                                  # two significant digits, with the final zero kept below 10 (4.0, not 4)
giqr = lambda v: "{} ({}--{})".format(g2(med(v)), *(g2(x) for x in np.percentile(v, [25, 75])))                       # Gbp: median (first quartile to third quartile)
diqr = lambda v: "{:.0f} ({:.0f}--{:.0f})".format(med(v), *np.percentile(v, [25, 75]))             # days: median (first quartile to third quartile)
T2 = []
for st in order:
    L = by[st]; k = sum(r["ok"] for r in L)
    T2.append(f"{st.replace('_', chr(92) + '_')} & {len(L)} & {k} ({pct(k, len(L))}) & {sum(r['absent'] for r in L)} & {diqr([r['days'] for r in L])} & {giqr([r['g'] for r in L])} & {'--' if st == 'HMP_2019_ibdmdb' else format(med([r['nruns'] for r in L]), '.0f')} \\\\")
T2.append(f"\\midrule\nAll cohorts & {N} & {len(ok)} ({pct(len(ok), N)}) & {len(ab)} & {diqr([r['days'] for r in rows])} & {giqr([r['g'] for r in rows])} & -- \\\\")
T2 = "\n".join(T2)
nfail = N - len(ok)
wk = order[0]                                                                  # the cohort with the lowest recovery rate
wk_n = len(by[wk]); wk_ok = sum(r["ok"] for r in by[wk])
rate_lo = min(sum(1 for r in by[s] if r["rank_any_narrow"] == 1) / len(by[s]) for s in by if s != wk)                 # the measure of Fig. 1D
rate_hi = max(sum(1 for r in by[s] if r["rank_any_narrow"] == 1) / len(by[s]) for s in by if s != wk)
assert wk_ok == sum(1 for r in by[wk] if r["rank_any_narrow"] == 1), "the weakest cohort is recovered as often on both measures"

# ---- Table 12: interval
bands = [("30 days or less", 0, 30), ("31--90", 30, 90), ("91--180", 90, 180), ("181--391", 180, 391), ("More than 391", 391, 1e9)]
T3 = []; first = {}                                                          # first: the first sample alone, first in the HGM set (Fig. 1E)
for lab, lo, hi in bands:
    L = [r for r in rows if lo < r["days"] <= hi] if lo else [r for r in rows if r["days"] <= hi]
    if not L and hi == 1e9:
        continue                                                                 # the open-ended class is shown only when it holds donors
    k = sum(r["ok"] for r in L)
    kf = sum(1 for r in L if r["rank_v1_narrow"] == 1)
    first[lab] = (kf, len(L))
    T3.append(f"{lab} & {len(L)} & {k} ({pct(k, len(L))}) & {sum(r['absent'] for r in L)} & {kf} ({pct(kf, len(L))}) \\\\")
T3 = "\n".join(T3)
over1 = [r for r in rows if r["days"] > 391]; over1_ok = sum(r["ok"] for r in over1)

# ---- Table 13: size of the queried run
cuts = [(0, 0.1, "Under 0.1"), (0.1, 0.5, "0.1--0.5"), (0.5, 1, "0.5--1"), (1, 2, "1--2"), (2, 4, "2--4"), (4, 1e9, "4 or more")]
T4 = []
for lo, hi, lab in cuts:
    g = [r for r in rows if lo <= r["g"] < hi]
    T4.append(f"{lab} & {len(g)} & {sum(r['ok'] for r in g)} ({pct(sum(r['ok'] for r in g), len(g))}) \\\\")
T4 = "\n".join(T4)
small = [r for r in rows if r["g"] < SMALL]; large = [r for r in rows if r["g"] >= SMALL]
ks, kl = sum(r["rank_any_narrow"] == 1 for r in small), sum(r["rank_any_narrow"] == 1 for r in large)                  # the measure of Fig. 1F
comp = collections.Counter(r["study"] for r in small)

# ---- Table 14: within cohort, below and above the cohort median
T5 = []; lo_k = lo_n = hi_k = hi_n = 0; half = {}; unsplit = []
for st in sorted(by):                                                          # every cohort: each is split at its own median
    g = by[st]
    m = med([r["g"] for r in g]); lo = [r for r in g if r["g"] <= m]; hi = [r for r in g if r["g"] > m]
    a, b = sum(r["ok"] for r in lo), sum(r["ok"] for r in hi)
    if not hi:                                                                 # all the queried runs of the cohort hold the same number of bases: it is not split, its donors are given across the two columns
        unsplit.append(st)
        T5.append(f"{st.replace('_', chr(92) + '_')} & {gfmt(m)} & \\multicolumn{{2}}{{c}}{{{a}/{len(lo)} ({pct(a, len(lo))})}} \\\\")
        continue
    lo_k += a; lo_n += len(lo); hi_k += b; hi_n += len(hi); half[st] = (pct(a, len(lo)), pct(b, len(hi)))
    T5.append(f"{st.replace('_', chr(92) + '_')} & {gfmt(m)} & {a}/{len(lo)} ({pct(a, len(lo))}) & {b}/{len(hi)} ({pct(b, len(hi))}) \\\\")
n5 = len(half)
assert unsplit == ["RaymondF_2016"], "the caption of Supplementary Table 14 says that RaymondF_2016 is the only cohort that is not split"
T5.append(f"\\midrule\nTotal, {n5} split cohorts & -- & {lo_k}/{lo_n} ({pct(lo_k, lo_n)}) & {hi_k}/{hi_n} ({pct(hi_k, hi_n)}) \\\\")
T5 = "\n".join(T5)
gq = sorted(r["g"] for r in rows)
iqr = (gq[int(.25 * (N - 1))], gq[int(.75 * (N - 1))])
b30 = [r for r in rows if r["days"] <= 30]; b30_ok = sum(r["ok"] for r in b30)
over1_peh = sum(1 for r in over1 if r["study"] == "PehrssonE_2016"); over1_rest = [r for r in over1 if r["study"] != "PehrssonE_2016"]; over1_rest_ok = sum(r["ok"] for r in over1_rest)
wo_peh = [r for r in rows if r["study"] != "PehrssonE_2016"]; wo_peh_ok = sum(r["ok"] for r in wo_peh)
peh = json.load(open(PEH)) if PEH else None
nl = json.load(open(NL)) if NL else None
nln = (nl['nielsen_donors'] + nl['li_donors']) if nl else 0; nl_linked = (nl['nielsen_linked_all'] + nl['li_linked_all']) if nl else 0
f2 = json.load(open(FIG2)) if FIG2 else None
vx = json.load(open(VIN)) if VIN else None
hdv = json.load(open(HDV)) if HDV else None
sc = json.load(open(SEL)) if SEL else None
st_ = json.load(open(STK)) if STK else None; se_ = json.load(open(STS)) if STS else None
M_TEX = format(st_["stock"], ",").replace(",", "{,}") if st_ else "291{,}881"          # the size of the HGM set, written for math mode
assert len(ok) == n_stock, "every run of the index competes: the donors recovered and the measure of Fig. 1 are the same"
lk = json.load(open(LNK)) if LNK else None
cs7 = json.load(open(CST)) if CST else None
costea_compare = (f"below the {min(d['closest_other'] for d in cs7['donors']):.2f} to {max(d['closest_other'] for d in cs7['donors']):.2f} found between a run and the same reads deposited again, less trimmed and filtered (below), and " if cs7 else "")

# Supplementary Note 1, "Independent estimate of the human share", and the caption of Supplementary Fig. 3D and E: every number is read from the values written by supplementary_figure3.py
stat_par = stat_cap = ""
if hdv:
    stc = hdv["stat_comparison"]; th5 = stc["thresholds"]["5"]; lst = stc["linked_by_stat_human_share"]; dvs = stc["donors_5_percent_or_more"]; bch = stc["by_cohort"]
    s11, s0 = stc["spearman_11_or_more_human_hashes"], stc["spearman_under_11_human_hashes"]
    r_lo, r_hi = stc["median_ratio_stat_over_ours_when_ours_0_1_to_5"], stc["median_ratio_stat_over_ours_when_ours_20_or_more"]
    runs5_o, runs5_s = th5["both_high"] + th5["only_ours_high"], th5["both_high"] + th5["only_stat_high"]
    assert all(bch[c]["runs_5_percent_or_more_ours"] == 0 and bch[c]["runs_5_percent_or_more_stat"] == 0 for c in ("Costea", "Raymond", "Heintz-Buschart")), "the sentence on the cohorts without a run at 5% would be wrong"
    assert sum(b["runs_5_percent_or_more_ours"] for b in bch.values()) == runs5_o and sum(b["runs_5_percent_or_more_stat"] for b in bch.values()) == runs5_s
    assert [d["label"] for d in lst] == ["<0.5", "0.5\u20135", "5\u201350", "\u226550"] and sum(d["donors"] for d in lst) == vx["pool"]["pool_donors"]
    stat_par = rf"""
\textbf{{Independent estimate of the human share.}} The NCBI runs a k-mer taxonomic analysis on every run of the SRA (STAT version 2.1.0\cite{{katz2021}}). We downloaded its result for the runs of the cohorts. The share of the species-level hits that STAT assigns to \textit{{Homo sapiens}} is a second human share of each run, computed on the run itself and not on our sketches. The two shares are correlated in rank (Spearman correlation {stc['spearman']:.2f} over the {stc['runs']:,} runs, {s11['value']:.2f} over the {s11['runs']} runs whose sketch holds 11 human hashes or more and {s0['value']:.2f} over the {s0['runs']} others, whose share rests on 10 human hashes or fewer; Supplementary Fig. 3D). They place {th5['both_high'] + th5['both_low']:,} of the {stc['runs']:,} runs on the same side of 5\%, {th5['both_high']} at 5\% or more by both and {th5['both_low']:,} below 5\% by both (Cohen's $\kappa$ {th5['kappa']:.2f}). Of the other {th5['only_ours_high'] + th5['only_stat_high']} runs, {th5['only_ours_high']} are at 5\% or more by the hashes only and {th5['only_stat_high']} by STAT only. STAT gives a lower share than the hashes when human DNA is scarce and a close one when it is abundant. Vincent holds {bch['Vincent']['runs_5_percent_or_more_ours']} of the {runs5_o} runs at 5\% or more by the hashes and {bch['Vincent']['runs_5_percent_or_more_stat']} of the {runs5_s} by STAT, and no run of Costea, Raymond or Heintz-Buschart reaches 5\% by either measure (largest share in Costea {bch['Costea']['max_ours_percent']:.1f}\% by the hashes and {bch['Costea']['max_stat_percent']:.1f}\% by STAT). By donor, the more human of the two samples is at 5\% or more for {dvs['vincent']['stat']} of the {dvs['vincent']['donors']} donors of Vincent by STAT and for {dvs['vincent']['ours']} by the hashes, and for {dvs['others']['stat']} of the {dvs['others']['donors']} donors of the cohorts kept by STAT and for {dvs['others']['ours']} by the hashes. The share recovered falls with the STAT human share as it does with ours, {pct(lst[0]['linked'], lst[0]['donors'])}\% below 0.5\% ({lst[0]['donors']} donors), {pct(lst[1]['linked'], lst[1]['donors'])}\% from 0.5 to 5\% ({lst[1]['donors']}), {pct(lst[2]['linked'], lst[2]['donors'])}\% from 5 to 50\% ({lst[2]['donors']}) and {pct(lst[3]['linked'], lst[3]['donors'])}\% at 50\% or more ({lst[3]['donors']}; Supplementary Fig. 3E). STAT is thus an independent indication, on another scale, that the human content is not an artefact of our references.
"""
    stat_cap = rf""" \textbf{{D}}, Share of the hashes of each of the {stc['runs']:,} runs of C that match a human genome, against the human share of the same run in the NCBI taxonomic analysis (STAT; Supplementary Note 1), for the {len(by)} cohorts kept (grey) and Vincent (red); the dotted line is the diagonal and the grey lines mark 5\%. \textbf{{E}}, As B, with the STAT human share."""

title = open(os.path.join(OVL, "title.txt"), encoding="utf-8").read().strip()                                                   # the \title{..} of main.tex
# Authors and affiliations under the title, as in main.tex (same order and wording); drop {AUTHORS} from the template for a double-anonymized submission
AUTHORS = (r"Guillaume Gautreau\textsuperscript{1}, Guillaume Kon Kam King\textsuperscript{1}, Nicolas Pons\textsuperscript{2}, Roland Faure\textsuperscript{3}" "\n\n"
           r"{\small\textsuperscript{1}Université Paris-Saclay, INRAE, MaIAGE, Jouy-en-Josas, France.\\" "\n"
           r"\textsuperscript{2}Université Paris-Saclay, INRAE, MetaGenoPolis, Jouy-en-Josas, France.\\" "\n"
           r"\textsuperscript{3}CNRS UMR3525, Institut Pasteur, Université Paris Cité, Paris, 75015, France.}")
n_out_same = ko["same_cohort_other_person"]
lo7 = min(d["closest_other"] for d in cs7["donors"]) if cs7 else None                                  # the nearest run of other people of the seven Costea donors of Supplementary Note 1
hi_out = [r for r in out if lo7 is not None and (r["max_imp"] or 0) >= lo7]
assert all(r["study"] == "CosteaPI_2017" for r in hi_out), "a donor outside Costea has a nearest run of other people as similar as the seven of Note 1"
hi_txt = (f"For {len(hi_out)} of the {ko['unmapped']} donors at rank 2 or more whose nearest run of other people lies outside the 12 cohorts, this run has a similarity of {lo7:.2f} to {max(r['max_imp'] for r in hi_out):.2f}. These {len(hi_out)} donors are among the seven Costea donors of Supplementary Note 1, whose reads sit in the archive under two subject codes." if hi_out else "")
band_wk = [r for r in by[wk] if 30 < r["days"] <= 90]; band_ot = [r for r in rows if r["study"] != wk and 30 < r["days"] <= 90]       # donors of the weakest cohort and of the others with 31--90 days between the samples
NAMES = {"PehrssonE_2016": "Pehrsson", "VincentC_2016": "Vincent", "NielsenHB_2014": "Nielsen", "LiJ_2014": "Li", "MehtaRS_2018": "Mehta", "CosteaPI_2017": "Costea", "HansenLBS_2018": "Hansen", "RaymondF_2016": "Raymond", "Bengtsson-PalmeJ_2015": "Bengtsson-Palme", "HMP_2019_ibdmdb": "HMP IBD", "Heitz-BuschartA_2016": "Heintz-Buschart", "IjazUZ_2017": "Ijaz"}
QO = ", ".join(f"{n} of {NAMES[st]}" for st, n in sc["queried_outside_the_donors"].items()) if sc else ""             # donors outside the 831 whose last sample was queried (cohort_selection.py)
multi = sorted(NAMES[st] for st in by if med([r["nruns"] for r in by[st]]) > 1)
multi_txt = ", ".join(multi[:-1]) + " and " + multi[-1]
assert sc is None or sc["mehta_second_run"]["second_run_smaller"] == sc["mehta_second_run"]["with_second_indexed_run"], "Supplementary Note 2 says that the second indexed run of the BioSample of the query sample is smaller than the listed run, for every such donor of Mehta"
note = rf"""
\section*{{Supplementary Note 2. Donors recovered and not recovered}}

Most donors are recovered, but not all. This Note looks at the donors that are not recovered, and at how the recovery rate varies between cohorts, with the time between the samples and with the depth of the query.

A donor is recovered when its other sample is at rank 1 (Methods). This is the case for {len(ok)} of {N} donors ({pct1(len(ok), N)}\%). The ranking is made among the 291,881 human gut metagenomes of the HGM set, which are part of the {INDEX_RUNS} runs of the index. The other {nfail} donors fall into two groups (Supplementary Fig. 1B and Supplementary Table 10). For {len(ab)} donors ({pct1(len(ab), N)}\%) the other sample is not among the runs that the Faure \textit{{et al.}} service returns. The service returns at most the 5,000 nearest runs, down to a similarity of 0.05, so these donors count as not recovered. For {len(out)} donors ({pct1(len(out), N)}\%) the other sample is returned, but at least one sample of another person has higher similarity to the query, placing the other sample at a median rank of {med([r['rank_any_narrow'] for r in out]):.0f}{(', including ' + str(sum(r['rank_any_narrow'] == 2 for r in out)) + ' for whom it is at rank 2') if any(r['rank_any_narrow'] == 2 for r in out) else ''}.

For the {len(out)} donors whose other sample is at rank 2 or more, the other sample is far less similar to the query than for the {len(ok)} recovered donors (median similarity {q(med([r['partner_any'] for r in out]))} compared with {q(med([r['partner_any'] for r in ok]))}). For each donor we also took the nearest run of other people among all the runs of the index. Its median similarity to the query is {q(med([r['max_imp'] for r in out]))} for the donors at rank 2 or more and {q(med([r['max_imp'] for r in ok]))} for the recovered donors (two-sided Mann--Whitney test, \textit{{U}} = {int(u_close):,}, \textit{{P}} = {p_close:.2f}). 

For the {len(ab)} donors whose other sample is not returned, the nearest run of other people has a median similarity of {q(med([r['max_imp'] for r in ab]))}. The nearest run of other people comes from the donor's own cohort for {n_out_same} of the {len(out)} donors at rank 2 or more ({pct(n_out_same, len(out))}\%), compared with {kok['same_cohort_other_person']} of the {len(ok)} recovered donors ({pct(kok['same_cohort_other_person'], len(ok))}\%). For {ko['unmapped']} of these {len(out)} donors it is a run outside the 12 adult cohorts whose runs we resolved (Supplementary Note 1). Such a run counts as another person, which can only lower the recovery rate. When the nearest run of other people comes from the donor's cohort, it can be a run of a relative or of a cohabitant. It can also be a run of the donor curated under a second identifier or deposited under the label of another sample.{' ' + hi_txt if hi_txt else ''}
In the open-set test of Fig. 1H (Methods), {len(accepted)} of the {N} searches for an absent person are accepted at its threshold. In {ka['same_cohort_other_person']} of these {len(accepted)} searches ({pct(ka['same_cohort_other_person'], len(accepted))}\%) the nearest run of other people comes from the donor's own cohort, compared with {kb['same_cohort_other_person']} of the {len(below)} searches that are not accepted ({pct(kb['same_cohort_other_person'], len(below))}\%).

The recovery rate differs between cohorts (Fig. 1D and Supplementary Table 11). {NAMES[wk]} has the lowest rate, with {wk_ok} of its {wk_n} donors recovered ({pct(wk_ok, wk_n)}\%), whereas {pct(int(rate_lo * 1000), 1000)}--{pct(int(rate_hi * 1000), 1000)}\% of the donors are recovered in each of the other {len(by) - 1} cohorts. To study the time between the samples, the first sample takes the place of the other sample (Fig. 1E and Supplementary Table 12). The first sample is at rank 1 in the HGM set for {pct_rng([first[c_] for c_ in ('30 days or less', '31--90', '91--180')])}\% of the donors whose first sample is up to 180 days from the query, and for {pct(first['181--391'][0], first['181--391'][1])}\% at 181--391 days.{'' if not over1 else ' Beyond 391 days ' + ('only one donor remains' if len(over1) == 1 else 'only ' + str(len(over1)) + ' donors remain') + '.'}

The recovery rate is lower only for queries of less than {SMALL} Gbp (Fig. 1F and Supplementary Tables 13 and 14). The Faure \textit{{et al.}} service sketches runs, so a query is one run. When the curated metadata list several runs for the query sample, as for all the query samples of {multi_txt} (Supplementary Table 11), the query is the largest indexed one. For {sc['mehta_second_run']['with_second_indexed_run']} of the {sc['mehta_second_run']['n']} donors of Mehta, the BioSample of the query sample holds a second indexed run, a smaller deposit of the same sample that the curated metadata do not list. This run counts neither as the donor's other sample nor as another person. Only {len(small)} queries hold less than {SMALL} Gbp, and {ks} of them are recovered ({pct(ks, len(small))}\%), compared with {kl} of the {len(large)} larger queries ({pct(kl, len(large))}\%; Fig. 1F and Supplementary Table 13).

Neither depth nor the interval between the samples explains the low rate of {NAMES[wk]}. The queries of {NAMES[wk]} hold a median of {gfmt(med([r['g'] for r in by[wk]]))} Gbp (Supplementary Table 11). Within {NAMES[wk]}, {sum(r['ok'] for r in by[wk] if r['g'] > med([x['g'] for x in by[wk]]))} of the {sum(r['g'] > med([x['g'] for x in by[wk]]) for r in by[wk])} donors whose query is above this median are recovered ({half[wk][1]}\%), compared with {sum(r['ok'] for r in by[wk] if r['g'] <= med([x['g'] for x in by[wk]]))} of the {sum(r['g'] <= med([x['g'] for x in by[wk]]) for r in by[wk])} donors whose query is at or below it ({half[wk][0]}\%; Supplementary Table 14). {('All its ' + str(wk_n) + ' donors have') if len(band_wk) == wk_n else ('Of its ' + str(wk_n) + ' donors, ' + str(len(band_wk)) + ' have')} 31--90 days between the two samples, and {pct(sum(r['ok'] for r in band_wk), len(band_wk))}\% of them are recovered, compared with {pct(sum(r['ok'] for r in band_ot), len(band_ot))}\% of the {len(band_ot)} donors of the other cohorts in that band."""
rows_far = [r for r in rows if r["far_possible"]]
note2 = fig2 = note3 = ""
if FIG2:
    A2 = f2["A"]
    assert lk, "Note 1 gives the number of shared identifiers: the twelfth argument (nielsen_li_links.py) is needed"
    xc = f2["nielsen_li_cross"]
    assert xc, "Supplementary Fig. 4C needs the run of CroCoDeEL on Nielsen and Li together (third argument of supplementary_figure4.py)"
    xw = xc["within"]
    r_ = lambda k: f"{A2[k][0]} of {A2[k][1]}"
    p_ = lambda k: pct1(A2[k][0], A2[k][1])
    note3 = ""
    if peh:
        pb = f2["by_cohort"]
        d = peh["donors"]; n_in = peh["n_near_in_lanes"]; n12 = peh["near_offsets_in_lanes"].get("12", 0); ctrl = peh["control_offsets_1_2_11_13"]; run = peh["near_run_offsets"]
        pl = peh["per_lane"]; lanes = [k for k in pl if k != "140805_ECP"]
        NUMW = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six", 7: "Seven"}
        assert sc, "Supplementary Note 1 describes how the cohorts were chosen: the eleventh argument (cohort_selection.py) is needed"
        cat_, lg, lb, ac, at = sc["catalogue"], sc["longitudinal"], sc["labels"], sc["adult_cohorts"], sc["adult_totals"]
        assert at["analysed"] == N and len(by) == sum(1 for c in ac.values() if c["analysed"]), "the donors of the analysis are not those of cohort_selection.py"
        n_left = sum(c["adult_dna"] for c in ac.values() if c["adult_dna"] and not c["analysed"]); k_left = sum(1 for c in ac.values() if c["adult_dna"] and not c["analysed"])
        n_idx = sum(1 for c in ac.values() if c["indexed"]); n_dna = sum(1 for c in ac.values() if c["adult_dna"])
        iv, rq, nr, ex_b = sc["analysed_interval"], sc["rna_queries"], sc["not_recorded_adult_category"], sc["biosample_only_in_index"]
        exs, hle, hl = sc["excluded"], sc["hmp_links_effect"], se_["hmp_links"]
        assert exs["candidates"] == at["adult_dna"] - n_left and exs["candidates"] - exs["contamination"] - exs["re_deposit"] == N == len(rows) and len(rows_ex) == exs["contamination"] + exs["re_deposit"], "the 497 donors of the cohorts kept, the donors removed by the two rules and the donors of the analysis must add up"
        seven = {d["subject"] for d in cs7["donors"]}
        seven_conta = [e for e in rows_ex if e["rule"] == "contamination" and e["subject"] in seven]; seven_redep = [e for e in rows_ex if e["rule"] == "re-deposit"]
        assert len(seven) == 7 and len(seven_conta) + len(seven_redep) == 7 and all(e["subject"] in seven for e in seven_redep) and hle["rank_any_narrow"] == 0 and hle["studies"] == ["HMP_2019_ibdmdb"], "the rule of re-deposit removes the Costea donors of Supplementary Note 1 that the rule of contamination leaves, and the runs given back change no rank of the other sample"
        assert iv["same_day_cohorts"] == ["MehtaRS_2018"] and list(ex_b) == ["MehtaRS_2018"], "the text of Note 1 names Mehta for the same-day pairs and for the donors outside the curated accessions"
        sel_note = rf"""We applied the rules below to the whole of curatedMetagenomicData, which holds {cat_['studies']} studies and {cat_['samples']:,} samples, {cat_['stool']:,} of them stool samples (Supplementary Fig. 1A and Supplementary Tables 1 and 2). They lead to {n_dna} cohorts, and we left {NUMW[k_left].lower()} of them out. A donor is a subject identifier of the curated metadata. We selected donors with at least two stool samples and considered their first and last samples, ordered by days since the first collection (with missing collection days assigned to day 0). A study with at least {NUMW[lg['min_donors']].lower()} such donors was considered a longitudinal cohort, yielding {lg['cohorts']} cohorts comprising {lg['donors']:,} donors. From here on, a donor of a cohort is a donor with at least two stool samples (but it can be more, Supplementary Fig. 1B).

We then labelled each cohort by the ages that the curated metadata record for its donors at their first sample. A donor under 3 years is an infant. A cohort is labelled when at least 70\% of its donors have an age. It is labelled as infants if at least 80\% of the donors with an age are infants, as adults if at most 20\% of them are, and as mixed otherwise. The label belongs to the cohort and not to its donors, so an adult cohort can hold minors, which the age rule below removes. Mehta records no age, but all its stool samples are from men and the study describes a cohort of adult men\cite{{mehta2018}}, so we counted it as adult. This gives {lb['infants']['cohorts']} infant cohorts ({lb['infants']['donors']} donors), {lb['adults']['cohorts']} adult cohorts ({lb['adults']['donors']}), {lb['mixed']['cohorts']} mixed cohort ({lb['mixed']['donors']}) and {lb['not recorded']['cohorts']} cohorts that could not be labelled for lack of ages ({lb['not recorded']['donors']}). We kept the adult cohorts, because the question is the linking of one adult, and did not query the infant cohorts. The {lb['not recorded']['cohorts']} cohorts without enough ages were not queried either. The curated age category marks all the stool samples of {nr['cohorts']} of them ({nr['donors']} donors) as adult. However, this category is filled even when no age is recorded, and these {nr['cohorts']} cohorts record no age, so we did not count them as adult.

A donor of the {lb['adults']['cohorts']} adult cohorts is analysed only if each of its two samples has a run in the index. This run must be among the runs listed in the curated metadata or, for HMP IBD, whose curated metadata list none, among the runs of its ENA project (below). This index rule holds for {at['indexed']} donors in {n_idx} cohorts, {ac['HMP_2019_ibdmdb']['indexed']} of them in HMP IBD. Another {ex_b['MehtaRS_2018']} donors of Mehta have an indexed run of each sample, but for at least one of the two samples this run is in the BioSample of the sample without being among its curated accessions. These donors are left out. {sum(sc['queried_outside_the_donors'].values())} donors that this rule leaves out were queried all the same, because the curated run of their last sample is indexed. Their lists of neighbours only complete the HGM set (Supplementary Table 3), and no rank of the analysis depends on them. Supplementary Table 4 gives the results with other sets of runs counted as other people. Two cohorts could not be tested. None of the reads of HMP 2012 ({ac['HMP_2012']['donors']} donors) is in the index, and the ENA lists them as controlled access (dbGaP). The curated samples of HMP 2019 T2D ({ac['HMP_2019_t2d']['donors']} donors) carry no run accession and its ENA project holds only amplicon runs. HMP IBD has no run accession in the curated metadata either, so we took the first and the last sample of each of its donors from its ENA project. This was possible for {sc['ibd_rebuilt_pairs']} of its {ac['HMP_2019_ibdmdb']['donors']} donors, those with at least two collection dates among the metagenomic runs of the project whose sample alias is a participant and visit code without suffix. For {ac['HMP_2019_ibdmdb']['indexed']} of these {sc['ibd_rebuilt_pairs']}, the first and the last of these dates each have one of these runs in the index. In the ENA project a sample is a visit of a participant. A visit counts as a curated sample when a curated sample of the donor has twice the read count of one of the ENA runs of the visit, because the curated metadata count both reads of a read pair where the ENA counts the read pair once. For the {sc['ibd_pairs']['n']} of these {ac['HMP_2019_ibdmdb']['indexed']} donors that stay in the analysis after the age rule below, the first and last samples are the earliest and the latest of the visits that have such a run in the index. The query is the deepest of these runs for the last visit. Fig. 1E uses these first and last samples. In the 30-day subset (Methods) the query is the same and the place of the other sample is taken by the most similar of the donor's visits at least 30 days from the query that have a run in the index.

Two donor rules select the donors of these cohorts. The first keeps the donors with a lowest recorded age of 18 years or more and those of Mehta, which records no age and was counted as adult (above). Of the {at['indexed']} donors it removes {at['under18']} under 18 and {at['no_age']} without an age. The {at['under18']} donors under 18 include the {ac['IjazUZ_2017']['under18']} donors of Ijaz that are in the index (all its donors are aged 7 to 14), so this cohort leaves the analysis. The second rule requires a DNA query. It removes the {at['rna']} donors of Mehta whose queried run is named as an RNA library (a library name ending in \_rna). Of these {rq['n']} donors, {rq['linked']} are recovered ({pct(rq['linked'], rq['n'])}\%; rank 1 among all the runs of the index), compared with {rq['mehta_dna_linked']} of the {rq['mehta_dna_n']} donors of Mehta with a DNA query ({pct(rq['mehta_dna_linked'], rq['mehta_dna_n'])}\%). The rule looks at the query only. For {sc['rna_other_samples']['donors']} donors of Mehta that stay in the analysis, the curated metadata name only RNA libraries for the first sample and for every other sample at least 30 days from the query. The DNA library of each of these samples is recognised as a run of the donor (Methods). The two donor rules leave {at['adult_dna']} donors in {n_dna} cohorts. Without the {NUMW[k_left].lower()} cohorts left out and the {exs['contamination'] + exs['re_deposit']} donors that the rules of contamination and re-deposit remove (below), the analysis holds {at['analysed']} donors in {len(by)} cohorts (Supplementary Table 2).

"""
        ps = peh.get("search")                                                     # pehrsson_lanes.py run with its DATA_DIR: the nearest run of the donors whose query sample is in a pair
        assert ps and ps["nearest_is_pair_partner_above_own"] == ps["nearest_is_pair_partner"] and peh["donors"]["query_in_pair"] == ps["query_in_pair"] and peh["donors"]["query_in_pair_linked"] == 0, "Supplementary Note 1 says that the run of the other participant is above every other sample of the donor and that none of these donors is recovered"
        note3 = rf"""
\section*{{Supplementary Note 1. The cohorts}}

The analysis requires donors with at least two stool samples, one to query the index and another to recover among the returned runs. We searched for such donors in curatedMetagenomicData, which compiles published metagenomic studies with curated metadata. This Note describes how we selected the donors for the analysis from this resource.

\subsection*{{Selection of the cohorts and donors}}

{sel_note}\subsection*{{Cohorts left out}}

Four of the {n_dna} cohorts that the rules lead to were left out of the analysis, for the reasons below. With these four cohorts included, the analysis would hold {N + nln + d['n'] + vx['vincent']['n']} donors in {len(by) + 4} cohorts, and {pct1(len(ok) + nl_linked + d['linked'] + vx['vincent']['linked_stock'], N + nln + d['n'] + vx['vincent']['n'])}\% of them would be recovered (rank 1 among all the runs of the index).

\subsubsection*{{Pehrsson}}

Its {peh['samples']} samples include pairs of samples that come from two different participants and have near-identical species profiles, so the sample of a participant is not always that participant's. The cohort was left out for these profiles.

In this cohort, {peh['n_near']} of the {peh['n_diff_pairs']:,} pairs of samples from two participants have a Bray-Curtis dissimilarity below 0.10 (e.g., samples 010202 and 020202, corresponding to SAMN04261148, participant 1\_2, 68-year-old female and SAMN04261163, participant 2\_2, 4-year-old female, both 7 days after the first collection, have a dissimilarity of 0.04). The median dissimilarity is {peh['median_bc_diff']:.2f} over all these pairs and {peh['median_bc_same']:.2f} between two samples of one participant. Of the {peh['n_events']} events that CroCoDeEL reports in this cohort (below), {peh['n_events_other']} are between two participants (median rate {rp(peh['median_rate_events_other'])}\%, median dissimilarity of their two samples {peh['median_bc_events_other']:.2f}), against {f2['n_other']} of the {f2['kept']['events']} events of all the samples of the {len(by)} cohorts kept. Of the {ps['query_in_pair']} adult donors of this cohort whose query sample is in one of these {peh['n_near']} pairs, {ps['nearest_is_pair_partner']} have as nearest run of the index, once the runs of the query sample itself are set aside, a run of the sample of the other participant of the pair (median similarity {ps['median_score_pair_partner']:.2f}). This run is more similar than every other sample of the donor among the runs returned. For {ps['own_sample_returned']} of these {ps['nearest_is_pair_partner']} donors at least one such sample is returned, and the most similar one has a similarity of {ps['median_score_best_own_sample']:.2f} at the median. None of the {ps['query_in_pair']} donors is recovered.

We located these pairs in the sequencing lanes (Supplementary Table 5). Each sample has an unpaired run whose name in the SRA carries a date and, for {sum(v['samples'] for k, v in pl.items() if k != '140805_ECP')} of the {peh['samples']} samples, a lane. None of the {peh['n_near']} pairs spans two lanes. When the samples of a lane are ordered by the accession of their unpaired run, {n12} of the {n_in} pairs that sit inside the {len(lanes)} named lanes join samples exactly 12 positions apart. Of the {sum(pl[k]['pairs12'] for k in lanes)} pairs of samples from two participants that are 12 positions apart in these lanes, {sum(pl[k]['near12'] for k in lanes)} are below 0.10, compared with {ctrl[1]} of the {ctrl[0]} pairs that are 1, 2, 11 or 13 positions apart. The other {len(run)} of the {peh['n_near']} pairs are inside the sequencing run dated 140805, which has no lane name. The Methods of Pehrsson \textit{{et al.}}\cite{{pehrsson2016}} report demultiplexing by barcode with no mismatches, with reads whose pair did not carry the same barcode kept as unpaired reads.

Of the {d['n']} adult donors of this cohort, {d['with_pair']} have the query sample or the other sample in one of the {peh['n_near']} pairs below 0.10; {d['with_pair_linked']} of them is recovered, compared with {d['without_pair_linked']} of the {d['without_pair']} others.

Although the cause cannot be clearly established from the public data, this pattern may result from a mixing of DNA or barcodes between two sets of 12 libraries of a lane during library preparation. We left this cohort out of the analysis as a precaution. The corresponding author of Pehrsson \textit{{et al.}} was informed of these observations on 17 November 2026 and is looking into the matter. The authors were able to reproduce the near-identical profiles and to confirm the issue (e-mail from the second author, 23 November 2026).

\subsubsection*{{Nielsen and Li}}

These two cohorts carry the same subject identifiers for {nl['shared_nielsen_li']} patients. The samples that the two cohorts hold under one identifier are different samples with different run accessions, two in Nielsen and a third in Li for {lk['time_points']['two_nielsen_one_li']} of the {nl['shared_nielsen_li']} identifiers (below). For {nl['nielsen_with_same_id']} of the {nl['nielsen_donors']} Nielsen donors, the 5,000 nearest runs of the query include a sample of Li under the same identifier. The pipeline counts such a sample as another person, so it can have higher similarity to the query than the donor's other sample. The {nl['nielsen_outranked']} Nielsen donors whose other sample is at rank 2 or more are only in part the {nl['nielsen_outranked']} donors with a sample in Li. Of these {nl['nielsen_outranked']} donors, {nl['nielsen_best_other_cohort']} have a run of Li as the nearest run of other people, {nl['nielsen_best_same_id']} of them under the donor's own identifier. We then counted a sample of the other cohort that carries the donor's identifier as a sample of the donor and not as another person. With this counting, {nl['nielsen_rank1_after']} of the {nl['nielsen_donors']} Nielsen donors are recovered instead of {nl['nielsen_rank1_before']}, and {nl['li_rank1_after']} of the {nl['li_donors']} Li donors instead of {nl['li_rank1_before']}. Nielsen also shares {nl['shared_nielsen_lechatelier']} identifiers with LeChatelier 2013, a cohort of the same consortium that is not in the analysis and whose runs the pipeline cannot assign to a person. These numbers are therefore lower bounds. We left both cohorts out rather than merge identifiers across cohorts that we cannot all resolve.

These recovered links nevertheless indicate that a person's samples can be recovered across two studies that differ in instrument, read length and depth (below).

@@LINKS@@

\subsubsection*{{Vincent}}

Many samples of this cohort hold mostly human DNA. We measured human DNA for the {vx['pool']['pool_donors']} donors of the {len(by)} cohorts kept and of Vincent (the {vx['others']['n']} donors of the cohorts kept, before the two rules below, and the {vx['vincent']['n']} donors of Vincent). The human share of a sample is the percentage of the hashes of its sketch that match a human genome, the sketch being recomputed from the public unitigs of the run (Methods). The human share of a donor is the higher of its two samples. It has a median of {vx['vincent']['share_median']:.0f}\% over the {vx['vincent']['n']} donors of Vincent (up to {vx['vincent']['share_max']:.0f}\%), compared with {vx['others']['share_median']:.2f}\% over the {vx['others']['n']} donors of the {len(by)} cohorts kept. It is 5\% or more for {vx['vincent']['n_over5']} of the {vx['vincent']['n']} donors of Vincent, compared with {vx['others']['n_over5']} of the {vx['others']['n']} (Supplementary Fig. 3A). Over the {vx['pool']['pool_donors']} donors the share recovered falls as the human share rises. It is {pct(vx['strata'][0]['linked'], vx['strata'][0]['donors'])}\% below 0.5\% human hashes ({vx['strata'][0]['donors']} donors), {pct(vx['strata'][1]['linked'], vx['strata'][1]['donors'])}\% from 0.5 to 5\% ({vx['strata'][1]['donors']}), {pct(vx['strata'][2]['linked'], vx['strata'][2]['donors'])}\% from 5 to 50\% ({vx['strata'][2]['donors']}) and {pct(vx['strata'][3]['linked'], vx['strata'][3]['donors'])}\% at 50\% or more ({vx['strata'][3]['donors']}; Supplementary Fig. 3B). We then removed the hashes that match a human genome and searched each query among the {vx['pool']['pool_runs']:,} runs (the two runs of each of the {vx['pool']['pool_donors']} donors). For the {vx['pool']['poor_removed']['n']} donors under 5\%, the share whose other sample is the nearest of these runs does not change, {vx['pool']['poor_all']['percent']:.0f}\% before the removal and after it (Supplementary Fig. 3C). The removal restores part of the deficit of the {vx['pool']['rich_all']['n']} donors with 5\% or more. Their other sample is the nearest of the {vx['pool']['pool_runs']:,} runs for {vx['pool']['rich_all']['top1']} of them before the removal and for {vx['pool']['rich_removed']['top1']} after. For the {vx['vincent']['n']} donors of Vincent these counts are {vx['vincent']['pool_all']} before and {vx['vincent']['pool_removed']} after, and {vx['vincent']['linked_stock']} of the {vx['vincent']['n']} are recovered in the HGM set. For {vx['vincent']['discordant']} donors of Vincent one sample holds under 20\% human hashes and the other over 60\%, and {vx['vincent']['discordant_linked']} of them are recovered. The cohort was left out for its human content. Among the cohorts kept, {vx['others']['n_over5_linked']} of the {vx['others']['n_over5']} donors with 5\% or more human hashes are recovered.
{stat_par}
\subsection*{{Rules of contamination and re-deposit}}

Two rules remove {exs['contamination'] + exs['re_deposit']} of the {exs['candidates']} donors of the cohorts kept. The first removes the samples studied that take part, as source or target, in a CroCoDeEL event between two people at a rate above 5\% (below). It removes {exs['samples_removed']} samples, and {exs['contamination']} donors leave the analysis, {exs['contamination_query']} whose query sample is removed and {exs['contamination_no_other']} that have no other sample left. Of these {exs['contamination']} donors, {exs['recovered_contamination']} would have been recovered. The second removes the donors whose nearest run of other people has a similarity to the query above the highest similarity of the other sample of any donor that the first rule leaves. They are the {NUMW[len(seven_redep)].lower()} donors of Costea described below that the first rule leaves (the first rule removes the {'seventh' if len(seven_conta) == 1 else 'others'}). Their nearest run of other people carries the reads of their own sample under another subject code, and none of them was recovered. In HMP IBD, the public table of the consortium gives to its participants runs of two other projects of the consortium, PRJNA306874 and PRJNA389280. Such a run counts as the participant's own, neither its other sample nor a run of other people ({hl['runs']} runs of {hl['participants']} participants, {hl['donors']} of them donors of the analysis). This changes the rank of the other sample of no donor. It changes the rank of the first sample alone for {hle['rank_v1_narrow']} donors and that of the samples at least 30 days from the query for {hle['rank_far_narrow']}. It also removes the deposits of the donor's own samples from the search for an absent person (Fig. 1H).

@@RULES@@
"""

    kp, bz = f2["kept"], f2["by_cohort"]
    top_c = max([c for c in bz if not bz[c]["left_out"]], key=lambda c: bz[c]["percent_samples_in_other"])
    pk = {c: bz[c]["percent_samples_in_other"] for c in bz if not bz[c]["left_out"]}
    top_pct = max(round(v, 1) for v in pk.values())
    top_txt = " and ".join(f"{NAMES[c]} ({bz[c]['samples_in_other']} of {bz[c]['samples']})" for c in sorted([c for c in pk if round(pk[c], 1) == top_pct], key=lambda c: -pk[c]))
    ex_ = f2["exposed"]
    kept_txt = (f"The last, a donor of {' and '.join(NAMES[c] for c in ex_['kept_other_cohorts'])}, keeps its most similar other sample, which takes part in an event at {100 * ex_['kept_other_max_rate']:.1f}\\%, and stays in the analysis." if ex_['kept_other'] == 1 and ex_['kept_query'] == 0 else
                f"The other {ex_['kept_other'] + ex_['kept_query']} stay in the analysis, in events at rates up to {100 * ex_['kept_other_max_rate']:.1f}\\%.")
    fs = f2["studied"]; fu, fo, bs = fs["kept"], fs["left_out"], fs["by_cohort"]        # Supplementary Note 1 is limited to the events in which a sample studied takes part
    bu = {c: v for c, v in bs.items() if not v["left_out"]}
    pu = {c: bu[c]["percent_studied_in_other"] for c in bu}
    top_pct_u = max(round(v, 1) for v in pu.values())
    top_txt_u = " and ".join(f"{NAMES[c]} ({bu[c]['studied_in_other']} of {bu[c]['samples_studied']})" for c in sorted([c for c in pu if round(pu[c], 1) == top_pct_u], key=lambda c: NAMES[c]))
    zero_u = sorted([c for c in bu if bu[c]["studied_in_other"] == 0], key=lambda c: NAMES[c])
    low_u = min(round(v, 1) for v in pu.values())
    low_txt_u = ("none of the samples studied of " + " and ".join(NAMES[c] for c in zero_u)) if zero_u else (f"{low_u:.1f}\\% of those of " + " and ".join(NAMES[c] for c in sorted(bu, key=lambda c: NAMES[c]) if round(pu[c], 1) == low_u))
    range_u = ("none" if zero_u else f"{low_u:.1f}\\%") + f" to {top_pct_u:.1f}\\%"
    peh_s, nie_s, li_s, vin_s = (bs[c] for c in ("PehrssonE_2016", "NielsenHB_2014", "LiJ_2014", "VincentC_2016"))
    assert all(xc[k]["events_studied"] == xc[k]["events"] for k in ("all", "same", "other")), "a sample studied takes part in every event between Nielsen and Li, so the numbers of these events are those of the samples studied"
    xws = xc["within_studied"]
    ofn = lambda k, n: f"all {n}" if k == n else f"{k} of the {n}"                  # 71 of the 71 reads as all 71
    assert ex_["samples_used"] == fu["studied_in_other"] and ex_["query_linked_stock"] == 0 and ex_["query_linked_stock"] + ex_["other_linked_stock"] == ex_["linked_stock"] and ex_["query_donors"] > 1 and ex_["other_donors"] > 1 and ex_["query_donors"] + ex_["other_donors"] == ex_["n"] and ex_["excluded_query"] + ex_["excluded_other"] == sc["excluded"]["contamination"] and ex_["candidates"] == sc["excluded"]["candidates"], "Supplementary Note 1 counts the samples and the donors of the two groups and says that none of the donors with the query sample in an event with another person is recovered"
    note2 = rf"""
\subsubsection*{{Cross-sample contamination}}

We ran CroCoDeEL version 1.3.1 on the species table of each of {len(bs)} cohorts, the {len(bu)} kept and the {len(bs) - len(bu)} left out of the analysis (above), with all the samples of the cohort, the default settings of the tool and no abundance filter. We report only the events in which a sample studied takes part, as source or target. The samples studied are the query sample, the first sample and the most similar other sample of each donor (Supplementary Fig. 1B). They are {fu['samples_studied']:,} samples for the {fu['donors']} donors of the cohorts kept, before the two rules above, and {fo['samples_studied']} for the {fo['donors']} donors of the cohorts left out. The tool treats cross-sample contamination as a dilution, in which a source sample adds species to a target sample. These species keep in the target the relative abundances that they have in the source, multiplied by one constant, the rate. On log scales of source and target abundances, these species lie on a line of slope one, the contamination line. A random forest scores each pair that has a contamination line, and pairs with a probability of at least 0.5 are reported as events. Every pair of samples is tested in both directions, so it can give two events. An event is between two people when its two samples belong to different subjects of the curated metadata, and between two samples of one person otherwise. Indeed, CroCoDeEL can flag as contaminated two longitudinal samples of one individual, especially over short intervals, among the causes of false positives\cite{{goulet2026}} (README of the tool).

Events between two people are few among the samples studied of the cohorts kept (Supplementary Fig. 4 and Supplementary Table 8). These {fu['samples_studied']:,} samples take part in {fu['events']} events, {fu['other']} of them ({pct(fu['other'], fu['events'])}\%) between two people. Of the {fu['samples_studied']:,} samples, {fu['studied_in_other']} ({100 * fu['studied_in_other'] / fu['samples_studied']:.1f}\%) take part in at least one of the {fu['other']} events between two people, from {low_txt_u} to {top_pct_u:.1f}\% of those of {top_txt_u}. The rate of these {fu['other']} events has a median of {rp(fu['median_rate_other'])}\% (range {rp(fu['min_rate_other'])}\% to {rp(fu['max_rate_other'])}\%). In Pehrsson, {ofn(peh_s['other'], peh_s['events'])} events in which a sample studied takes part are between two people, and {peh_s['studied_in_other']} of its {peh_s['samples_studied']} samples studied ({peh_s['percent_studied_in_other']:.0f}\%) take part in at least one of these {peh_s['other']} events (above). In Nielsen, Li and Vincent, also left out, {nie_s['studied_in_other']} of {nie_s['samples_studied']}, {li_s['studied_in_other']} of {li_s['samples_studied']} and {vin_s['studied_in_other']} of {vin_s['samples_studied']} samples studied take part in at least one event between two people ({nie_s['percent_studied_in_other']:.1f}\%, {li_s['percent_studied_in_other']:.1f}\% and {vin_s['percent_studied_in_other']:.1f}\%), compared with {range_u} in the cohorts kept.

Nielsen and Li carry the same subject identifiers for {lk['shared_subject_ids']} patients (above), so we also ran CroCoDeEL on their two species tables merged ({xc['samples']['merged']} samples; Supplementary Fig. 4C) to find events between the two cohorts. The tool reports {xc['all']['events']} events between a Nielsen and a Li sample, and a sample studied takes part in each of them. Their median rate is {rp(xc['all']['median_rate'])}\% (range {100 * xc['all']['min_rate']:.1f}\% to {rp(xc['all']['max_rate'])}\%). These {xc['all']['events']} events join {xc['all']['pairs']} pairs of samples. Of the {xc['samples']['nielsen_studied']} Nielsen samples studied, {xc['all']['studied_nielsen']} take part in at least one of them, as do {xc['all']['studied_li']} of the {xc['samples']['li_studied']} Li samples studied. {NUMW[xc['same']['events']]} of these events join two samples that carry the same subject identifier (rates {rp(xc['same']['min_rate'])}\% to {rp(xc['same']['max_rate'])}\%) and {xc['other']['events']} join samples with different identifiers (median rate {rp(xc['other']['median_rate'])}\%). The {xc['all']['studied_nielsen']} Nielsen samples studied that take part in these {xc['all']['events']} events belong to {xc['all']['donors_nielsen_studied']} donors, too few to account for the {nl['nielsen_donors'] - nl['nielsen_linked_all']} Nielsen donors that are not recovered.

The {ex_['samples_used']} samples studied of the cohorts kept that take part in an event between two people belong to {ex_['n']} donors. The query sample is one of them for {ex_['query_donors']} donors, in events at rates of {rp(ex_['query_min_rate'])}\% to {rp(ex_['query_max_rate'])}\%, and none of these donors is recovered. Another sample studied is one of them for {ex_['other_donors']} donors, in events at rates of {rp(ex_['other_min_rate'])}\% to {rp(ex_['other_max_rate'])}\%. The rule of contamination removes the samples in such an event at a rate above 5\%, so {ex_['excluded_query'] + ex_['excluded_other']} of these donors leave the analysis. {kept_txt}

Contamination between samples of a cohort, as CroCoDeEL detects it, does not account for the recoveries of the {N} donors of the analysis. The query sample took part in no event for {A2['none/stock'][1]} of the {N} donors, and {p_('none/stock')}\% of these {A2['none/stock'][1]} donors are recovered in the HGM set, compared with {pct1(n_stock, N)}\% of all {N} donors. The 30-day subset (Methods) holds the {len(rows_far)} donors with another sample at least 30 days from the query, and the most similar of those samples takes the place of their other sample. The query sample took part in no event for {A2['none/30d'][1]} of these {len(rows_far)} donors, and {p_('none/30d')}\% of these {A2['none/30d'][1]} donors are recovered, compared with {pct1(sum(r['rank_far_narrow'] == 1 for r in rows_far), len(rows_far))}\% of all {len(rows_far)}.
"""
    fig2 = rf"""\begin{{figure}}[t]\centering
\includegraphics[width=\textwidth]{{supplementary_figure4}}
\caption{{\textbf{{CroCoDeEL events between two people, by cohort.}} Events in which a sample studied takes part, as source or target, in the {len(by)} cohorts kept (grey) and the {len(bz) - len(by)} cohorts left out of the analysis (red). The samples studied are the query sample, the first sample and the most similar other sample of each donor, in the cohorts kept those of the {sc['excluded']['candidates']} donors before the two rules of Supplementary Note 1. CroCoDeEL was run on the species table of each cohort, with all its samples. Each event that the tool reports has a source sample and a target sample in the same species table (Supplementary Note 1), and it is between two people when its two samples belong to different subjects of the curated metadata. \textbf{{A}}, Share of the samples studied of each cohort that take part in at least one event between two people (samples concerned over samples studied), hatched according to whether the sample is a source (right-leaning), a target (left-leaning), or both (cross-hatched). \textbf{{B}}, Rate of each of these events, the share of the target that comes from the source, in percent (a dot per event, the bar is the median), with the number of events in which a sample studied takes part and the number of them between two people. Dashed line, 5\%, the rate above which the rule of contamination of Supplementary Note 1 removes the samples studied of an event between two people in the cohorts kept. Cohorts are ordered by the share in A within each group. \textbf{{C}}, Nielsen and Li carry the same subject identifiers for {lk['shared_subject_ids']} patients, so CroCoDeEL was also run on their two species tables merged, and the panel shows the events between a Nielsen sample and a Li sample in which a sample studied takes part, all together and by whether the two samples carry the same subject identifier (a dot per event, the bar is the median), with the number of events and the number of samples studied of each cohort that take part in one.}}
\end{{figure}}
"""
tab67 = ""
if FIG2 and peh:
    bcz = f2["by_cohort"]
    kept6 = sorted([c for c in bcz if not bcz[c]["left_out"]], key=lambda c: (-bcz[c]["percent_samples_in_other"], c))
    out6 = sorted([c for c in bcz if bcz[c]["left_out"]], key=lambda c: (-bcz[c]["percent_samples_in_other"], c))
    fmtrate = lambda med_, lo_, hi_: f"{rp(med_)} ({rp(lo_)}--{rp(hi_)})" if med_ is not None else "--"
    bs6, fs6 = f2["studied"]["by_cohort"], f2["studied"]                              # the events in which a sample studied takes part
    kept6 = sorted([c for c in bs6 if not bs6[c]["left_out"]], key=lambda c: (-bs6[c]["percent_studied_in_other"], c))
    out6 = sorted([c for c in bs6 if bs6[c]["left_out"]], key=lambda c: (-bs6[c]["percent_studied_in_other"], c))
    row6 = lambda lab, b: f"{lab} & {b['samples']:,} ({b['samples_studied']:,}) & {b['events']} & {b['other']} ({pct(b['other'], b['events']) if b['events'] else 0}) & {b['studied_in_other']} ({100 * b['studied_in_other'] / b['samples_studied']:.1f}) & {fmtrate(b['median_rate_other'], b['min_rate_other'], b['max_rate_other'])} \\\\"
    T6 = "\\multicolumn{6}{l}{\\textit{Cohorts kept}} \\\\\n" + "\n".join(row6(NAMES[c], bs6[c]) for c in kept6) + "\n\\midrule\n" + row6("All cohorts kept", fs6["kept"])
    T6 += "\n\\midrule\n\\multicolumn{6}{l}{\\textit{Cohorts left out}} \\\\\n" + "\n".join(row6(NAMES[c], bs6[c]) for c in out6) + "\n\\midrule\n" + row6("All cohorts left out", fs6["left_out"])
    L7 = {"141118_ECP_Lane1": "141118, lane 1", "141118_ECP_Lane2": "141118, lane 2", "150102_ECP_Lane1": "150102, lane 1", "150102_ECP_Lane2": "150102, lane 2", "150102_ECP_Lane3": "150102, lane 3",
          "150102_ECP_Lane4": "150102, lane 4", "140805_ECP": "140805 (no lane name)"}
    pl = peh["per_lane"]
    T7 = "\n".join(f"{L7[c]} & {pl[c]['samples']} & {pl[c]['pairs12']} & {pl[c]['near12']} & {pl[c]['near_all']} \\\\" for c in L7)
    T7 += f"\n\\midrule\nAll rows above ({sum(v['samples'] for v in pl.values())} of {peh['samples']} samples) & {sum(v['samples'] for v in pl.values())} & {sum(v['pairs12'] for v in pl.values())} & {sum(v['near12'] for v in pl.values())} & {sum(v['near_all'] for v in pl.values())} \\\\"
    tab67 = rf"""
\begin{{table}}[!ht]\centering\small
\caption{{\textbf{{Events between two people, by cohort.}} CroCoDeEL events in which a sample studied takes part, as source or target, in each of the {len(bs6)} cohorts (the {len(kept6)} kept and the {len(out6)} left out); Supplementary Fig. 4 draws them. \textit{{Samples}}, samples of the cohort in the curated tables and, in parentheses, the samples studied, namely the query sample, the first sample and the most similar other sample of each donor, in the cohorts kept those of the {sc['excluded']['candidates']} donors before the two rules of Supplementary Note 1. \textit{{Events}}, events in which a sample studied takes part. \textit{{Between two people}}, number of these events whose two samples belong to different subjects of the curated metadata, and share of the events (\%). \textit{{In such an event}}, samples studied that take part in at least one event between two people, and share of the samples studied (\%). \textit{{Rate}}, median and range of the rate of the events between two people, the share of the target that comes from the source, in percent.}}
\setlength{{\tabcolsep}}{{5pt}}
\begin{{tabular}}{{lrrrrr}}
\toprule
Cohort & Samples (studied) & Events & Between two people & In such an event & Rate (\%) \\
\midrule
{T6}
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[!ht]\centering\small
\caption{{\textbf{{Pairs of samples of different participants with near-identical profiles in Pehrsson, by sequencing lane.}} \textit{{Lane}}, the date and the lane named in the SRA alias of the unpaired run of each sample. The position of a sample in a lane is the rank of the accession of that run among the samples of the lane. \textit{{Samples}}, samples of the cohort in the lane. \textit{{Pairs 12 apart}}, pairs of samples of different participants 12 positions apart; \textit{{below 0.10}}, those with a Bray-Curtis dissimilarity below 0.10. \textit{{Any distance}}, pairs of samples of different participants with a Bray-Curtis dissimilarity below 0.10 inside the lane, whatever the number of positions between them.}}
\begin{{tabular}}{{lrrrr}}
\toprule
Lane & Samples & Pairs 12 apart & Below 0.10 & Any distance \\
\midrule
{T7}
\bottomrule
\end{{tabular}}
\end{{table}}
"""

note4 = tab8 = fig3 = fig6 = ""
if EXT:
    ex = json.load(open(EXT)); ho = ex["holdout"]
    def sci(v):
        if v is None or not (v < float("inf")):
            return "none"
        m, e = f"{v:.1e}".split("e")                                                   # rounds the mantissa before the exponent is read
        return f"${m}\\times10^{{{int(e)}}}$"
    f1 = lambda v: f"{v:.1f}"

    assert FIG3 and QUAD, "Supplementary Note 3 and Supplementary Tables 15 to 17 need the values of supplementary_figure5.py and of figure1_brief_communication.py (seventh and eighth arguments)"
    f3 = json.load(open(FIG3)); f3m = f3["models"]; f3g = f3["grid"]; f3meta = f3["meta"]; WINS = f3meta["windows"]
    hw = ex["holdout_windows"]                                                       # the test on thinned HGM sets of the fifteen fits (three models x five windows)
    sg = lambda v: "0.0" if abs(v) < 0.05 else f"{v:+.1f}".replace("-", "$-$")
    rmse = lambda c_, f: math.hypot(hw[f]["error"][c_], hw[f]["sd"][c_])            # root mean square error over the thinnings, in percentage points
    DRAWN = {k_: ho["0.01"]["windows"][k_] for k_ in sorted(ho["0.01"]["windows"], key=lambda k_: -f3m[k_]["first_8e9"])}   # the fits drawn in Fig. 1B, in the order of its legend, the highest value first
    assert all(w_ == 1000 for w_ in DRAWN.values()), "the three fits drawn in Fig. 1B are those for N >= 1,000, where the measured share and its upper bound differ by at most 0.2 points"
    NUM = {2: "two", 3: "three", 4: "four", 5: "five", 15: "fifteen"}                                     # the counts of models, of windows and of fits, in words
    geN = lambda w_: f"$N \\ge {w_:,}$".replace(",", "{,}")                           # a window, in the text
    NAME9 = {"pycmb": "Pitman--Yor correctness", "line": r"Straight line in $\log N$", "miss": r"Power law of the miss rate $1 - a N^{b}$"}
    assert all(f3m[k_]["window"] == w_ for k_, w_ in DRAWN.items()), "the windows of the fits drawn in Fig. 1B"
    assert all(abs(hw[f]["error"][f"{k_}|{w_}"] - ho[f]["error"][k_]) < 1e-9 for f in ho for k_, w_ in DRAWN.items()), "the test of the fits drawn is that of the models of Fig. 1B"
    NAME4 = {"pycmb": "the Pitman--Yor model", "line": "the straight line", "miss": "the power law of the miss rate"}
    sgr = lambda k_: (lambda a_, b_: sg(a_) if sg(a_) == sg(b_) else f"{sg(a_)} to {sg(b_)}")(min(ho[f]["error"][k_] for f in ho), max(ho[f]["error"][k_] for f in ho))   # the range of the mean errors of a fit drawn
    grp76 = {}
    for k_ in DRAWN:                                                                 # the fits drawn grouped by the range of their mean errors, in the order of the legend of Fig. 1B
        grp76.setdefault(sgr(k_), []).append(NAME4[k_])
    rng_pts = lambda a, b: f"{f1(a)} percentage points" if f1(a) == f1(b) else f"{f1(a)} to {f1(b)} percentage points"                        # a range of errors, or one value when the two bounds round to the same
    sp = ex["small_pools"]; sps = sp["share"]; sp8 = sp["at_8e9"]                                                          # the measured share and its upper bound for small pools
    gap = lambda n_: sps[str(n_)]["upper"] - sps[str(n_)]["measured"]
    assert gap(1000) < 0.25 and gap(10000) < 0.1 and gap(100000) < 0.01, "Note 3 says that the measured curve is exact to within 0.2 percentage points from 1,000 runs on"
    sp8w = sp["at_8e9_windows"]; LO0 = WINS[0]                                                                             # the same fits for each window: the window from 2 runs shows what a fit that starts in small pools owes to the convention
    assert all(abs(sp8w[c_][f"{k_}|{DRAWN[k_]}"] - sp8[c_][k_]) < 1e-9 for c_ in sp8 for k_ in DRAWN), "the fits of the five windows include those drawn"
    assert sp8w["upper"][f"miss|{LO0}"] < 0 < min(sp8w["upper"][f"{k_}|{LO0}"] for k_ in ("pycmb", "line")), "Note 3 says that the power law fitted from 2 runs to the upper bound gives a miss rate above 100%"
    g = lambda v: "0" if v < 0.05 else (f"{v:.1f}" if v < 1 else f"{v:.0f}")         # a share in percent: whole numbers, one decimal below 1
    hdr11 = " & ".join("$\\times$" + f"{v['factor']:.0f}" for v in ho.values())                    # distance beyond the thinned HGM set
    fs11 = ", ".join(f"{float(f_):g}" for f_ in ho); mp11 = ", ".join(f"{v['Mprime']:,}" for v in ho.values())
    v20 = [v_["first_8e9"] for v_ in f3g.values()]
    e_drawn = [ho[f]["error"][k_] for f in ho for k_ in DRAWN]; e_all = [hw[f]["error"][c_] for f in ho for c_ in hw[f]["error"]]
    # the test set against a prediction without any fit: the share measured in the thinned HGM set itself. Its error is the decline of the recovery rate from the thinned HGM set to the HGM set, exact from the ranks:
    # a donor with k runs more similar than the other sample is first in the thinned HGM set with probability (1 - f)^k
    kk_all = np.array([-1 if r["rank_any_narrow"] is None else r["rank_any_narrow"] - 1 for r in rows]); kf_ = kk_all[kk_all >= 0].astype(float)
    rate_ = 100 * float((kk_all == 0).sum()) / N
    assert all(abs(rate_ - ho[f]["measured"]) < 1e-4 for f in ho), "the ranks and the thinning test give the same measured recovery rate"
    keep_ = {f: (1 - float(f)) ** kf_ for f in ho}
    decl = {f: 100 * float(keep_[f].sum()) / N - rate_ for f in ho}                  # in percentage points
    decl_sd = {f: 100 * math.sqrt(float((keep_[f] * (1 - keep_[f])).sum())) / N for f in ho}
    base4 = {f: math.hypot(decl[f], decl_sd[f]) for f in ho}                         # root mean square error of the prediction without a fit
    better_all = [c_ for c_ in hw["0.01"]["error"] if all(math.hypot(hw[f]["error"][c_], hw[f]["sd"][c_]) < base4[f] for f in ho)]
    jn = lambda L: " and ".join([", ".join(L[:-1]), L[-1]] if len(L) > 1 else L)
    beat = {f: [k_ for k_ in DRAWN if math.hypot(ho[f]["error"][k_], ho[f]["sd"][k_]) < base4[f]] for f in ho}      # the fits drawn that do better than no fit, at each distance
    d3 = [decl[f] for f in ho]
    assert len(d3) == 3
    r100_drawn = [math.hypot(ho["0.01"]["error"][k_], ho["0.01"]["sd"][k_]) for k_ in DRAWN]
    assert not beat["0.1"] and beat["0.03"] and len(beat["0.01"]) == len(DRAWN), "the three drawn fits gain their advantage over the constant-rate baseline as the thinning factor grows"
    starts = [f"{max(2, round(ho[x_]['Mprime'] * 1000 / 291881)):,}" for x_ in ho]                                      # where the window from 1,000 runs starts in each thinned HGM set
    assert all(int(s_.replace(",", "")) <= 100 for s_ in starts) and gap(10) > gap(100) > gap(1000), "Note 3 says that the two bounds differ more where the window from 1,000 runs starts in the thinned HGM sets"
    assert set(DRAWN.values()) == {1000} and LO0 == 2, "Note 3 describes the fits drawn as fitted for N >= 1,000, and the fits from 2 runs as those that start in the smallest pools"
    dw_ = max(abs(sp8w["measured"][f"{k_}|1000"] - sp8w["upper"][f"{k_}|1000"]) for k_ in DRAWN)          # how much a fit drawn moves between the measured share and its upper bound, in points
    rel = lambda k_: [100 * ho[f]["error"][k_] / decl[f] for f in ho]                # error of a fit drawn in Fig. 1B, in percent of the decline
    pc3 = lambda v: "{:.0f}\\%, {:.0f}\\% and {:.0f}\\%".format(*v)

    # ---- Supplementary Table 17: the recovery rate and the three fits drawn in Fig. 1B with one cohort left out, and for all the 497 donors before the two rules of Supplementary Note 1
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                  # python -I does not put the script's folder on the path
    import pyc_mb as PY                                                             # Pitman-Yor correctness model of Rocher et al. 2025
    from scipy.special import gammaln
    M_HGM = 291881; NS14 = np.unique(np.round(np.logspace(np.log10(2), np.log10(M_HGM), 70)).astype(int))      # the grid of Fig. 1B

    def refit(kk):
        """Recovery rate and the three fits drawn in Fig. 1B at the population of the world, for the donors whose numbers of runs more similar than the other sample are kk (-1, other sample not returned)."""
        kv = kk[kk >= 0].astype(float); x_ = NS14.astype(float); R_ = []
        for n_ in NS14:
            m_ = n_ - 1
            with np.errstate(all="ignore"):
                lg = gammaln(M_HGM - kv + 1) - gammaln(M_HGM - kv - m_ + 1) - (gammaln(M_HGM + 1) - gammaln(M_HGM - m_ + 1))
            R_.append(float(np.where(M_HGM - kv >= m_, np.exp(lg), 0.0).sum() / len(kk)))
        R_ = np.array(R_); o_ = dict(n=len(kk), rec=int((kk == 0).sum()))
        for k_, lo_ in DRAWN.items():
            s_ = NS14 >= lo_
            if k_ == "line":
                v_ = np.polyval(np.polyfit(np.log2(x_[s_]), R_[s_], 1), np.log2(8e9))
            elif k_ == "pycmb":
                v_ = np.atleast_1d(PY.kappa(np.asarray(8e9, float), PY.fit(x_[s_], R_[s_])))[0]
            else:                                                                      # the miss rate 1 - R = a N^b, least squares on the logarithms
                v_ = 1 - np.exp(np.polyval(np.polyfit(np.log(x_[s_]), np.log(np.clip(1 - R_[s_], 1e-12, None)), 1), np.log(8e9)))
            o_[k_] = 100 * float(np.clip(v_, 0, 1))
        return o_

    all14 = refit(kk_all)
    assert all(abs(all14[k_] - f3m[k_]["first_8e9"]) < 1e-6 for k_ in DRAWN), "the fits made here on all the donors are those of Supplementary Fig. 5"
    kv_all = kk_all[kk_all >= 0].astype(float); x_all = NS14.astype(float)
    with np.errstate(all="ignore"):
        R_all = np.array([float(np.where(M_HGM - kv_all >= (n_ - 1), np.exp(gammaln(M_HGM - kv_all + 1) - gammaln(M_HGM - kv_all - (n_ - 1) + 1) - (gammaln(M_HGM + 1) - gammaln(M_HGM - (n_ - 1) + 1))), 0.0).sum() / len(kk_all)) for n_ in NS14])
    s_py = NS14 >= DRAWN["pycmb"]
    py_fit = PY.fit(x_all[s_py], R_all[s_py])
    assert abs(100 * float(np.atleast_1d(PY.kappa(np.asarray(8e9, float), py_fit))[0]) - f3m["pycmb"]["first_8e9"]) < 1e-6 and not py_fit.bounds_hit
    py_h_nats, py_h_bits, py_gamma, py_d = py_fit.h, py_fit.h / math.log(2), py_fit.gamma, py_fit.d
    st_arr = np.array([r["study"] for r in rows])
    loo = {st: refit(kk_all[st_arr != st]) for st in order}
    rate14 = {st: 100 * loo[st]["rec"] / loo[st]["n"] for st in loo}
    lo14, hi14 = min(rate14, key=rate14.get), max(rate14, key=rate14.get)              # the recovery rate is lowest and highest without these cohorts
    big14 = max(by, key=lambda s: len(by[s]))
    fail14 = sum(1 for r in by[hi14] if r["rank_any_narrow"] != 1)                  # donors not recovered of the cohort whose absence gives the highest rate
    lo14_txt = f"{NAMES[lo14]}" + (f", the cohort with the most donors ({len(by[lo14])} of {N})" if lo14 == big14 else f" ({len(by[lo14])} of {N} donors)")
    rng14 = lambda k_: f"{g(min(v_[k_] for v_ in loo.values()))}\\% to {g(max(v_[k_] for v_ in loo.values()))}\\%"
    # ---- the extreme-value model of the similarity tail (evt_model.py), the fourth model of Fig. 1B: model subsection + unified evaluation snippets, Supplementary Table 18 and a column of Supplementary Table 17
    evt_model_sec = evt_sens_intro = evt_sens_txt = evt_thin_txt = evt_loo_txt = evt_diag_txt = evt_lim_txt = evt_tab = ""; NEV = 0; evt_cell = lambda st: ""
    if EVT:
        assert EVS, "the checks of evt_similarities.py come with the output of evt_model.py"
        ev = json.load(open(EVT)); evs = json.load(open(EVS)); em = ev["models"]["free"]; et = ev["thinning"]; dr = str(ev["drawn"]); ed = em[dr]; NEV = 1
        assert ev["donors"] == N == evs["returned"]["donors"] + evs["not_returned"]["donors"] and ev["recovered"] == n_stock == evs["recovered"] and ev["hgm_runs"] == M_HGM == evs["hgm_runs"], "the donors and the HGM set of the extreme-value model are those of the analysis"
        assert list(et) == list(ho) and ev["thinnings"] == list(ho.values())[0]["thinnings"] and ev["min_excesses"] == 10, "the model is tested at the same distances and on as many thinnings as the three others"
        assert evs["not_returned"]["donors"] == sp["not_found"] and min(evs["not_returned"]["k"]) > 0 and evs["service"]["max_abs_difference"] < 1e-5 and evs["top"] == ev["similarities_stored"]
        assert set(ed["leave_one_out"]) == set(order) and list(ev["ranges"]) == ["free", "bounded", "unbounded"] and ev["ranges"]["free"] == [-0.5, 0.5] and ev["ranges"]["bounded"] == [-0.5, 0.0] and ev["ranges"]["unbounded"] == [0.0, 1.0]
        sig2 = lambda v: f"{v:.{max(0, 1 - int(math.floor(math.log10(v))))}f}"         # two significant digits
        sh = lambda r_: sig2(ev["tail_share"][r_])                                       # the share of the runs of other people that a tail holds, in percent
        neg = lambda v: f"{v:.2f}".replace("-", "$-$")
        tails = [str(r_) for r_ in ev["tails"]]; evloo = ed["leave_one_out"]; F10, F33, F100 = list(et)
        tested = {f: [r_ for r_ in tails if r_ in et[f]] for f in et}                    # the tails fitted at each distance
        assert tested[F10] == tails and tested[F33] == tested[F100] == tails[-len(tested[F33]):] and len(tested[F33]) < len(tails), "all the tails are tested at 10 times, the largest ones at 33 and 100 times"
        e2 = lambda v: f"{v:.2f}"                                                       # an error as it is printed
        best_r = {f: min(tested[f], key=lambda r_: et[f][r_]["rmse"]) for f in (F10, F33, F100)}
        assert all(best_r[f] == dr for f in (F10, F33, F100)), "Fig. 1B draws the tail with the lowest thinning RMSE at all three distances"
        assert et[F100][tested[F100][0]]["excesses"] == ev["min_excesses"] and et[F100][dr]["excesses"] >= ev["min_excesses"], "at 100 times the smallest tested tail keeps min_excesses similarities"
        assert all(et[f][dr]["rmse"] < et[f]["no_fit"]["rmse"] for f in (F10, F33, F100)), "the tail drawn improves on the constant-rate baseline at all three thinning factors"
        assert all(abs(et[f][r_]["given_to_not_recovered"] + et[f][r_]["taken_from_recovered"] - et[f][r_]["error"]) < 1e-9 and et[f][r_]["given_to_not_recovered"] >= 0 >= et[f][r_]["taken_from_recovered"] for f in et for r_ in tested[f]), "the two parts of an error"
        vt = [em[r_]["at_world"] for r_ in tails]
        kd = sorted(b_ - a_ for a_, b_ in evs["returned"]["k_lists_and_here"])
        assert len(kd) == 2 and kd[0] > 0 and evs["returned"]["same_k"] + len(kd) == evs["returned"]["donors"]
        evt_cell = lambda st: f" & {(ed['at_world'] if st is None else evloo[st]):.1f}"
        evb = lambda x_, on_: f"\\textbf{{{x_}}}" if on_ else x_
        e2 = lambda v: f"{v:.2f}"; sg2 = lambda v: f"{v:+.2f}".replace("-", "$-$")
        T18 = []
        for r_, v_ in em.items():
            on_ = r_ == dr; ok_ = lambda f: r_ in et[f]
            cells = [f"{int(r_):,}", sh(r_), neg(v_["shape_mean"]), f"{v_['at_1e7']:.1f}", f"{v_['at_world']:.1f}", f"{v_['interval_world'][0]:.0f}--{v_['interval_world'][1]:.0f}"]
            cells += [sg2(et[f][r_]["error"]) if ok_(f) else "--" for f in et] + [e2(et[f][r_]["rmse"]) if ok_(f) else "--" for f in et]
            T18.append(" & ".join(evb(x_, on_) for x_ in cells) + " \\\\")
        T18.append("\\midrule\n\\multicolumn{3}{l}{Constant-rate baseline} & -- & -- & -- & " + " & ".join(sg2(et[f]["no_fit"]["error"]) for f in et) + " & " + " & ".join(e2(et[f]["no_fit"]["rmse"]) for f in et) + " \\\\")
        dr_tex = f"{ev['drawn']:,}".replace(",", "{,}")
        max_tail_tex = f"{int(tails[-1]):,}".replace(",", "{,}")
        tails_list_txt = ", ".join(f"{int(t_):,}" for t_ in tails[:-1]) + f" or {int(tails[-1]):,}"
        evt_model_sec = rf"""\subsubsection*{{Extreme-value model of the similarity tail}}

The {NUM[len(DRAWN)]} curve-fitted models above reduce the information available in the data in two ways: (i) they collapse each donor's similarities into a binary outcome (recovered or not recovered) at each pool size; and (ii) they fit only the resulting aggregate curve, discarding how far each recovered donor's other sample lies above the runs of other people. To make the most of the available information, we model the distribution of similarities directly, building on biometric identification models that express rank-one recovery in a gallery of size $N$ as the probability that a person's other sample exceeds the maximum of $N - 1$ similarities to other people\cite{{grother2004,bolle2005}}, and on the large-gallery extreme-value analysis of Baveja \textit{{et al.}}\cite{{baveja2010}}. Whereas Baveja \textit{{et al.}} apply a Type I (Gumbel) limit to a pooled (or quality-stratified) similarity distribution and expand the resulting integral to obtain the first-order straight line in $\log N$ fitted above, we relax both assumptions: we condition on each donor's observed other-sample similarity and fit a donor-specific upper-tail model with an unconstrained shape parameter.

Formally, for each donor $i \in \{{1, \dots, {N}\}}$, let $v_i$ denote the observed similarity between donor $i$'s query and its other sample, and let $X_{{i,1}}, \dots, X_{{i,n_i}} \stackrel{{\mathrm{{i.i.d.}}}}{{\sim}} F_i$ denote the similarities to the $n_i$ runs of other people in the HGM set, with observed order statistics $x_{{i,(1)}} \ge x_{{i,(2)}} \ge \dots \ge x_{{i,(n_i)}}$. Extending the HGM set of size $M$ by $N - M$ new independent runs of other people with similarities $X'_{{i,1}}, \dots, X'_{{i,N-M}} \stackrel{{\mathrm{{i.i.d.}}}}{{\sim}} F_i$, let $E_i(N)$ be the event that donor $i$ is recovered in this enlarged pool of size $N \ge M$. If there is already a sample from another person with higher similarity than $v_i$ in the HGM set ($k_i > 0$), then $\mathbb{{P}}(E_i(N)) = 0$. If donor $i$ is recovered ($k_i = 0$, so $v_i > x_{{i,(1)}}$), recovery holds if and only if none of the added runs exceeds $v_i$, i.e. $E_i(N) = \{{\max_{{1 \le j \le N - M}} X'_{{i,j}} \le v_i\}}$, which conditional on the single-run exceedance probability $\lambda_i = \mathbb{{P}}(X'_{{i,1}} > v_i)$ has probability
\[
    \mathbb{{P}}\!\left(E_i(N) \mid \lambda_i\right) = (1 - \lambda_i)^{{N - M}}.
\]
Because only the highest similarities among runs of other people can overtake $v_i$, the bulk of $F_i$ is uninformative for extrapolation and only its upper tail matters, leading directly to a peaks-over-threshold extreme-value model.

Following standard peaks-over-threshold extreme-value practice\cite{{coles2001}}, we set each donor's threshold to an upper empirical quantile $u_i = x_{{i,(r+1)}}$ (retaining the $r$ highest similarities, or top $r/n_i$ fraction of donor $i$'s similarities to other people, and assessing sensitivity across a range of $r$ below), and model the resulting excesses $y_{{i,j}} = x_{{i,(j)}} - u_i > 0$ ($j = 1, \dots, r$) separately from the threshold exceedance probability $\zeta_i = \mathbb{{P}}(X'_{{i,1}} > u_i)$. By the Pickands--Balkema--de Haan theorem\cite{{coles2001}}, excesses $Y_{{i,j}} = (X - u_i) \mid (X > u_i)$ over a high threshold follow a generalised Pareto distribution (GPD) with scale $\sigma_i > 0$ and shape $\xi_i \in \mathbb{{R}}$, while conditional on $r$ exceedances out of $n_i$ runs of other people, assigning $\zeta_i$ the standard Jeffreys' conjugate $\mathrm{{Beta}}(1/2, 1/2)$ prior yields the exact posterior $\zeta_i \mid r, n_i \sim \mathrm{{Beta}}(r + 1/2,\, n_i - r + 1/2)$ independently of the excesses. Under the standard scale-invariant reference prior $p(\sigma_i \mid \xi_i) \propto \sigma_i^{{-1}}\mathbf{{1}}_{{\{{\sigma_i > \sigma_{{\min,i}}(\xi_i)\}}}}$, where $\sigma_{{\min,i}}(\xi_i) = \max(0,\, -\xi_i \max_j y_{{i,j}})$ is the lower support bound required for positivity of the GPD density, the conditional posterior for $\sigma_i$ is proper for every $\xi_i \in [-0.5, 0.5]$ and $r \ge 1$. Because $\xi_i$ is poorly identified from a single donor's tail alone, we pool shape information across donors via a hierarchical prior\cite{{cooley2007}}:
\begin{{align*}}
    Y_{{i,j}} \mid \sigma_i, \xi_i &\stackrel{{\mathrm{{i.i.d.}}}}{{\sim}} \mathrm{{GPD}}(\sigma_i, \xi_i), \qquad j = 1, \dots, r,\; i = 1, \dots, {N}, \\
    p(\sigma_i \mid \xi_i) &\propto \sigma_i^{{-1}}\mathbf{{1}}_{{\{{\sigma_i > \max(0,\, -\xi_i \max_j y_{{i,j}})\}}}}, \\
    \xi_i \mid a, b &\stackrel{{\mathrm{{i.i.d.}}}}{{\sim}} \mathrm{{Beta}}_{{[-0.5,\,0.5]}}(a, b),
\end{{align*}}
where $\mathrm{{Beta}}_{{[-0.5,\,0.5]}}(a, b)$ denotes a $\mathrm{{Beta}}(a, b)$ prior rescaled to $[-0.5, 0.5]$ following Martins and Stedinger\cite{{martins2000}} (see also Coles and Dixon\cite{{coles1999}}). Constraining $\xi_i \in [-0.5, 0.5]$ regularises shape estimates in finite tails while permitting both bounded ($\xi_i < 0$) and unbounded ($\xi_i \ge 0$) upper tails. Conditional on $(\zeta_i, \sigma_i, \xi_i)$, the exceedance probability at $v_i > u_i$ is
\[
    \lambda_i(\zeta_i, \sigma_i, \xi_i) = \zeta_i\left(1 + \xi_i\,\frac{{v_i - u_i}}{{\sigma_i}}\right)_{{+}}^{{-1/\xi_i}}, \qquad (z)_{{+}} = \max(z, 0),
\]
interpreted as $\zeta_i \exp(-(v_i - u_i)/\sigma_i)$ when $\xi_i = 0$; when $\xi_i < 0$, the fitted tail has a finite upper endpoint at $u_i - \sigma_i / \xi_i$, beyond which $\lambda_i = 0$. Averaging the joint posterior predictive recovery probability over all {N} donors yields
\[
    R(N) = \frac{{1}}{{{N}}}\sum_{{i\,:\,k_i = 0}} \mathrm{{E}}\!\left[\left(1 - \lambda_i(\zeta_i, \sigma_i, \xi_i)\right)^{{N - M}} \;\middle|\; r, n_i, \mathbf{{y}}_i\right], \qquad N \ge M,
\]
where $\mathbf{{y}}_i = (y_{{i,1}}, \dots, y_{{i,r}})$ and $R(M)$ equals the measured HGM recovery rate. Intuitively, donors whose other sample lies well beyond their fitted upper tail retain recovery probability close to one even at $N = 8\times10^9$, whereas donors near the tail edge gradually lose rank-one status as $N$ grows.

We fit this hierarchical model using empirical Bayes (maximum marginal likelihood) across all {N} donors. Because donors are conditionally independent given $(a, b)$ and the threshold exceedance probability $\zeta_i$ factorises cleanly from $(\sigma_i, \xi_i)$ and depends only on $(r, n_i)$, both parameter estimation and prediction reduce to low-dimensional deterministic integrals evaluated with standard numerical quadrature.

First, for each donor $i$ and candidate shape $\xi_i \in [-0.5, 0.5]$, we marginalise out the nuisance scale parameter $\sigma_i$ under its reference prior $p(\sigma_i \mid \xi_i) \propto \sigma_i^{{-1}}\mathbf{{1}}_{{\{{\sigma_i > \sigma_{{\min,i}}(\xi_i)\}}}}$:
\begin{{align*}}
    L_i(\xi_i) &= \int_{{\sigma_{{\min,i}}(\xi_i)}}^{{\infty}} p(\mathbf{{y}}_i \mid \sigma_i, \xi_i)\, p(\sigma_i \mid \xi_i)\, d\sigma_i \\
               &= \int_{{ -\infty }}^{{ \infty }} p\!\left(\mathbf{{y}}_i \mid \sigma_{{\min,i}}(\xi_i) + e^\phi,\, \xi_i\right)\, p\!\left(\sigma_{{\min,i}}(\xi_i) + e^\phi \mid \xi_i\right)\, e^\phi\, d\phi,
\end{{align*}}
where $\phi = \log(\sigma_i - \sigma_{{\min,i}}(\xi_i))$ maps the support onto $(-\infty, \infty)$. Shifting the log-integrand by its mode $\hat{{\phi}}_{{i,\xi}}$ for floating-point stability, we evaluate this one-dimensional integral directly over $(-\infty, \infty)$ by adaptive Gauss--Kronrod quadrature, obtaining $\ell_i(\xi_i) = \log L_i(\xi_i)$ up to an additive constant common across shapes.

Second, integrating over the shape prior $p(\xi \mid a, b) = \mathrm{{Beta}}_{{[-0.5,\,0.5]}}(a, b)$ yields the log marginal likelihood of the hyperparameters:
\[
    \ell(a, b) = \sum_{{i=1}}^{{{N}}} \log \int_{{-0.5}}^{{0.5}} L_i(\xi)\, p(\xi \mid a, b)\, d\xi.
\]
We discretise $[-0.5, 0.5]$ on a uniform grid of 80 equal-width cells $C_g$ with midpoints $\xi_g$ (cell width $\Delta\xi = 0.0125$, roughly four to nine times finer than the fitted shape prior standard deviation across all five tails; doubling to 160 cells alters extrapolated recovery rates by less than $10^{{-3}}$ percentage points in synthetic checks), weighting $L_i(\xi_g)$ by the exact $\mathrm{{Beta}}_{{[-0.5,\,0.5]}}(a, b)$ cell probability mass $\mathbb{{P}}(\xi \in C_g \mid a, b)$; this decouples the inner scale integrals from hyperparameter optimisation while remaining numerically stable even for sharply concentrated priors. Maximising $\ell(a, b)$ over $(a, b)$ gives the unpenalised empirical Bayes estimates $(\hat{{a}}, \hat{{b}})$ and each donor's posterior weights $p(\xi_g \mid \mathbf{{y}}_i, \hat{{a}}, \hat{{b}}) \propto L_i(\xi_g)\,\mathbb{{P}}(\xi \in C_g \mid \hat{{a}}, \hat{{b}})$.

Finally, for each recovered donor ($k_i = 0$) and pool size $N \ge M$, the posterior predictive recovery probability factorises into three nested one-dimensional integrals:
\begin{{align*}}
    &\mathrm{{E}}\!\left[\left(1 - \lambda_i(\zeta_i, \sigma_i, \xi_i)\right)^{{N - M}} \;\middle|\; r, n_i, \mathbf{{y}}_i, \hat{{a}}, \hat{{b}}\right] \\
    &\qquad = \sum_{{g=1}}^{{80}} p(\xi_g \mid \mathbf{{y}}_i, \hat{{a}}, \hat{{b}}) \int_{{\sigma_{{\min,i}}(\xi_g)}}^{{\infty}} p(\sigma_i \mid \mathbf{{y}}_i, \xi_g) \left[\int_{{0}}^{{1}} \left(1 - F_{{\mathrm{{Beta}}}}^{{-1}}(u)\, q_i(\sigma_i, \xi_g)\right)^{{N - M}}\! du\right] d\sigma_i,
\end{{align*}}
where $q_i(\sigma_i, \xi_g) = (1 + \xi_g (v_i - u_i)/\sigma_i)_{{+}}^{{-1/\xi_g}}$ and $p(\sigma_i \mid \mathbf{{y}}_i, \xi_g) = p(\mathbf{{y}}_i \mid \sigma_i, \xi_g)\, p(\sigma_i \mid \xi_g) / L_i(\xi_g)$, and the inner integral rewrites $\zeta_i \sim \mathrm{{Beta}}(r + 1/2,\, n_i - r + 1/2)$ via its exact quantile function $\zeta_i = F_{{\mathrm{{Beta}}}}^{{-1}}(u)$ on $u \in (0, 1)$. For every shape cell, the inner integral over $u \in (0, 1)$ is evaluated by 25-point Gauss--Legendre quadrature (avoiding nested adaptive quadrature callbacks while differing from 100-point quadrature by less than $10^{{-6}}$ for the main $r = {dr_tex}$ model and less than $5\times 10^{{-5}}$ across all thinned tails), nested inside direct adaptive Gauss--Kronrod quadrature over $\phi = \log(\sigma_i - \sigma_{{\min,i}}(\xi_g)) \in (-\infty, \infty)$.

Because the search service returns at most 5,000 runs per query and stops at similarity 0.05 (leaving some queries with far fewer than {ev['drawn']:,} HGM runs of other people, and preventing uniform tail sizes up to $r = {max_tail_tex}$ or uncensored thinned-pool validation), we computed the similarity of each query to every run in the HGM set directly from the hypervector database searched by the Faure \textit{{et al.}} service (which contains {evs['database']['accessions']:,} runs, including {M_HGM - evs['hgm_runs_absent']:,} of the {M_HGM:,} HGM runs, and agrees with the returned query--HGM similarities to within $10^{{-5}}$ across all {evs['service']['pairs']:,} pairs). This yields the exact count $k_i$ for every donor, matching the returned lists for {evs['returned']['same_k']} of the {evs['returned']['donors']} donors with a returned partner (exceeding them by {kd[0]} and {kd[1]} runs for the remaining two), ranging from {min(evs['not_returned']['k']):,} to {max(evs['not_returned']['k']):,} across $\mathcal{{C}}$ and leaving the same {n_stock} donors recovered. To select the threshold size $r$ on empirical grounds, we choose the tail that minimises thinned-pool extrapolation RMSE across all three validation distances (detailed below), namely the $r = {dr_tex}$ highest similarities per donor ({sh(dr)}\% of the runs of other people) for the main curve in Fig. 1B. This choice reflects a principled bias--variance compromise: smaller tails retain fewer excesses under heavy thinning and suffer elevated estimation variance, whereas larger tails extend further into the bulk of the similarity distribution and introduce upward bias. At $r = {dr_tex}$, the fitted hierarchical model yields a mean shape of {neg(ed['shape_mean'])} with negative posterior mean shapes for {100 * ed['shape_negative']:.0f}\% of the recovered donors; sensitivity across all {NUM[len(tails)]} tail sizes ($r = {int(tails[0]):,}$ to ${max_tail_tex}$) is examined below.

"""
        evt_sens_intro = (
            f" and refitted the extreme-value model across {NUM[len(tails)]} tail sizes $r$ = {tails_list_txt} similarities per donor "
            f"({sh(tails[0])}\\% to {sh(tails[-1])}\\% of the runs of other people; bottom row)"
        )
        evt_sens_txt = (
            f" Across all {NUM[len(tails)]} extreme-value tail sizes $r$, world-population predictions range similarly from "
            f"{g(min(vt))}\\% to {g(max(vt))}\\%, with larger $r$ flattening slightly earlier because they reach further into the bulk of the similarity distribution."
        )
        evt_censor_txt = (
            f" The extreme-value model avoids this censoring ambiguity altogether by using the exact hypervector counts $k_i$ across all {N} donors."
        )
        evt_thin_txt = (
            f" Third, using donor-level similarity tails yields the highest validation accuracy: fitted with exact counts $k_i$ and the same tail fraction in each thinned set "
            f"(excluding tails with fewer than ten excesses, so that all {NUM[len(tested[F10])]} tails are tested at 10 times and the {NUM[len(tested[F33])]} largest at 33 and 100 times), "
            f"the $r = {dr_tex}$ extreme-value model achieves the lowest RMSE of all tested tails and improves on both the constant-rate baseline and the curve-fitted models "
            f"at all three thinning distances, with RMSE {jn([e2(et[f][dr]['rmse']) for f in (F10, F33, F100)])} percentage points "
            f"(mean errors {jn([sg2(et[f][dr]['error']) for f in (F10, F33, F100)])}, retaining {et[F100][dr]['excesses']} excesses per donor at 100-fold thinning) "
            f"versus baseline RMSE {jn([e2(et[f]['no_fit']['rmse']) for f in (F10, F33, F100)])} percentage points, directly motivating the choice of $r = {dr_tex}$ for Fig. 1B. "
            f"Decomposing the 33-fold thinning error of the $r = {dr_tex}$ model into probability erroneously retained for unrecovered donors and probability "
            f"withdrawn from recovered donors yields {sg2(et[F33][dr]['given_to_not_recovered'])} and {sg2(et[F33][dr]['taken_from_recovered'])} percentage points, respectively."
        )
        evt_loo_txt = f" and between {g(min(evloo.values()))}\\% and {g(max(evloo.values()))}\\% for the extreme-value model ($r = {dr_tex}$)"
        evt_diag_txt = f"Full numerical diagnostics across all {NUM[len(tails)]} tails of the extreme-value model appear in Supplementary Table 18. "
        evt_lim_clauses = (
            " (including, for the extreme-value model, the tendency of negative GPD shape parameters to induce finite upper endpoints that can arise from finite-threshold samples of unbounded distributions such as Gaussian tails)",
            "; for the extreme-value model they include the re-estimation of the distribution of tail shapes on every resample, but not any dependence among the non-match similarities to a given query, which its likelihood treats as conditionally independent"
        )
        evt_tab = rf"""
\begin{{table}}[!ht]\centering\small\setlength{{\tabcolsep}}{{3.5pt}}
\caption{{\textbf{{Extreme-value model of the similarity tail for five tails and test on thinned HGM sets.}} Generalised Pareto shape $\xi_i \in [-0.5, 0.5]$. \textit{{Tail}}, the number of highest similarities of the runs of other people used for each donor, and their share of these runs; bold, the tail drawn in Fig. 1B. \textit{{Mean shape}}, the mean of the distribution of shapes. \textit{{Recovered}}, share of the {N} donors whose other sample is first according to the model in a pool of $10^7$ runs and of $8\times10^9$ runs, the size of the population of the world, with its 95\% interval on {ev['resamples']:,} resamples of cohorts, then donors. Test, each run of other people is kept with probability $f$ ({fs11}) and the model is fitted on the thinned data alone with the same share of the runs of other people in the tail. \textit{{Mean error}} and \textit{{RMSE}}, mean and root mean square, over {ev['thinnings']} thinnings, of the recovery rate that the model predicts at 291,881 runs minus the measured share ({f1(list(ho.values())[0]['measured'])}\%), in percentage points, at 10, 33 and 100 times the thinned HGM set. --, a tail of fewer than ten similarities in the thinned HGM set, which is not fitted. Last row, constant-rate baseline, the share measured in the thinned HGM set itself taken as the prediction.}}
\begin{{tabular}}{{rrrrrcrrrrrr}}
\toprule
\multicolumn{{2}}{{c}}{{Tail}} & & \multicolumn{{3}}{{c}}{{Recovered (\%)}} & \multicolumn{{3}}{{c}}{{\makecell[cb]{{Mean error\\(percentage points)}}}} & \multicolumn{{3}}{{c}}{{\makecell[cb]{{RMSE\\(percentage points)}}}} \\
\cmidrule(lr){{1-2}}\cmidrule(lr){{4-6}}\cmidrule(lr){{7-9}}\cmidrule(lr){{10-12}}
Runs & \% & \makecell[cb]{{Mean\\shape}} & At $10^7$ & At $8\times10^9$ & Interval & {hdr11} & {hdr11} \\
\midrule
""" + "\n".join(T18) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    better_txt = (
        f"Across all three validation distances, the thinning experiment yields {'three' if NEV else 'two'} main conclusions. "
        f"First, all {'four' if NEV else NUM[len(DRAWN)]} models predict held-out full-HGM recovery accurately over up to two orders of magnitude without systematically overstating risk "
        f"(for the {NUM[len(DRAWN)]} curve-fitted fits drawn in Fig. 1B, mean errors span {sg(min(e_drawn))} to {sg(max(e_drawn))} percentage points, and {sg(min(e_all))} to {sg(max(e_all))} across all {NUM[len(v20)]} fits). "
        f"Second, modelling the decline clearly outperforms a naive constant-rate baseline predicting no further drop beyond the thinned pool "
        f"(which overestimates full-HGM recovery by {f1(d3[0])}, {f1(d3[1])} and {f1(d3[2])} percentage points at 10, 33 and 100 times): "
        f"although the true decline over a 10-fold step is small enough that the {NUM[len(DRAWN)]} $N \\ge 1{{,}}000$ curve-fitted models do not yet beat a flat baseline in root mean square error (RMSE), "
        f"{jn([NAME4[k_] for k_ in beat['0.03']])} outperform it at 33-fold growth, and all {NUM[len(DRAWN)]} reduce RMSE substantially at 100-fold growth "
        f"({min(r100_drawn):.1f} to {max(r100_drawn):.1f} versus {base4['0.01']:.1f} percentage points)."
    )
    row14c = lambda name, o_, e_="": f"{name} & {o_['n']} & {o_['rec']} ({pct1(o_['rec'], o_['n'])}) & " + " & ".join(f"{o_[k_]:.1f}" for k_ in DRAWN) + e_ + " \\\\"
    T14c = [row14c("None", all14, evt_cell(None))] + [row14c(st.replace("_", chr(92) + "_"), loo[st], evt_cell(st)) for st in order]
    six_cap = ""
    SHORT14 = {"pycmb": "the Pitman--Yor correctness model", "line": r"the straight line in $\log N$", "miss": "the power law of the miss rate"}
    HEAD14 = {"pycmb": "Pitman--Yor", "line": "Straight line", "miss": "Power law"}
    ex497_txt = ""
    if rows_ex:                                                                   # the rate before the two rules of Supplementary Note 1 stays visible
        kk497 = np.concatenate([kk_all, np.array([-1 if e["rank_any_narrow"] is None else e["rank_any_narrow"] - 1 for e in rows_ex])])
        all497 = refit(kk497); ev497 = " & --" * NEV; ev497_prose = ""
        if NEV and EVX:                                                           # the extreme-value model made again on these donors (evt_removed_donors.py)
            exr = json.load(open(EVX))
            assert exr["donors"] == all497["n"] == N + len(rows_ex) and exr["recovered"] == all497["rec"] and exr["removed"] == len(rows_ex) and exr["donors_of_the_analysis"] == N and exr["hgm_runs"] == M_HGM, "the donors of evt_removed_donors.py are those of the last row"
            assert exr["tail"] == ev["drawn"] and exr["range"] == ev["ranges"]["free"] and abs(exr["analysis"]["at_world"] - ed["at_world"]) < 1e-9 and abs(exr["analysis"]["a"] - ed["a"]) < 1e-9, "evt_removed_donors.py found the model of evt_model.py again before it fitted these donors"
            ev497 = f" & {exr['at_world']:.1f}"
            ev497_prose = f", while the extreme-value model predicts {exr['at_world']:.1f}\\%"
        T14c.append("\\midrule\n" + row14c("None, with the donors", all497, ev497) + "\nremoved by the two rules" + " &" * (2 + len(DRAWN) + NEV) + " \\\\")
        six_cap = f" Last row, all the {len(kk497)} donors of the cohorts kept, with the {len(rows_ex)} that the rules of contamination and re-deposit remove (Supplementary Note 1)."
        ex497_txt = (f" Even when all {len(kk497)} donors are retained prior to applying the contamination and re-deposit rules, "
                     f"the measured recovery rate stands at {pct1(all497['rec'], all497['n'])}\\% and the {NUM[len(DRAWN)]} curve-fitted extrapolations still range from "
                     f"{g(min(all497[k_] for k_ in DRAWN))}\\% to {g(max(all497[k_] for k_ in DRAWN))}\\% at $8\\times10^9$ runs{ev497_prose}.")
    LONG4 = {"pycmb": r"the Pitman--Yor correctness model of Rocher \textit{et al.}\cite{rocher2025}", "line": r"the straight line in $\log N$", "miss": r"the power law of the miss rate $1 - a N^{b}$\cite{grother2019}"}
    EVROW = (f"Last row, the extreme-value model of the similarity tail, which is not fitted to the measured curve and starts from the measured share at 291,881 runs; "
             f"columns, the tail size $r$, using the {', '.join(f'{int(t_):,}' for t_ in tails[:-1])} or {int(tails[-1]):,} highest similarities of each donor (Supplementary Table 18). ") if NEV else ""
    fig3 = rf"""\begin{{figure}}[t]\centering
\includegraphics[width=\textwidth,height=0.745\textheight,keepaspectratio]{{supplementary_figure5}}
\caption{{\textbf{{Extrapolations of the recovery curve of Fig. 1B for five fit windows{' and five tails' if NEV else ''}.}} Share of the {f3meta['donors']} donors whose other sample is first in a pool of $N$ runs, measured in pools drawn from the 291,881 runs of the HGM set (thick line, darker over the pool sizes used by the fit), given by the {'model' if NEV else 'fit'} up to $10^{{12}}$ runs (dashed line). {'First three rows' if NEV else 'Rows'}, the three models of Fig. 1B{' fitted to the measured curve,' if NEV else ''} in the order of its legend{'; columns' if NEV else '. Columns'}, the fit window, from $N$ = {', '.join(f'{w_:,}' for w_ in WINS[:-1])} or {WINS[-1]:,} to 291,881. {EVROW}Shaded, the fit {'and the tail ' if NEV else ''}drawn in Fig. 1B. Cross and value, the {'model' if NEV else 'fit'} at the population of the world, $8\times10^9$, with in parentheses{', for the first three rows,' if NEV else ''} its root mean square error on thinned HGM sets at 100 times, in percentage points (Supplementary Table 16); at 100 times the windows from 10 and 100 runs give the same thinned fit. Dotted lines, the HGM set and the population of the world.}}
\end{{figure}}
"""
    qd = json.load(open(QUAD)); qf = qd["first"]; WK = str(int(8e9)); K7 = str(int(1e7))
    assert abs(qf[WK]["est"] - f3m["line"]["first_8e9"]) < 1e-6      # Fig. 1B and Supplementary Fig. 5 draw the same fit
    bt = json.load(open(BOOT)) if BOOT else None
    assert bt is not None, "Supplementary Note 3 and Supplementary Fig. 6 need extrapolation_bootstrap.json (nineteenth argument)"
    bm = bt["models"]; bt_meta = bt["meta"]
    assert all(abs(bm[k_]["at_N"][WK]["est"] - f3m[k_]["first_8e9"]) < 1e-6 for k_ in DRAWN)
    biv = lambda k_: f"{bm[k_]['at_N'][WK]['lo']:.1f}--{bm[k_]['at_N'][WK]['hi']:.1f}"
    bnum = lambda v: format(v, ",").replace(",", "{,}")                               # a number of resamples, written for math mode
    bstr = lambda k_: f"{biv(k_)}\\% for {NAME4[k_]} ($B = {bnum(bm[k_]['B'])}$)"
    decades = f"{math.log10(8e9 / 291881):.1f}"
    min_lo_bt = min(bm[k_]["at_N"][WK]["lo"] for k_ in DRAWN)
    evt_censor_txt = "" if not NEV else evt_censor_txt
    evt_lim_clauses = ("", "") if not NEV else evt_lim_clauses
    bt_strs = [bstr(k_) for k_ in DRAWN] + ([f"{ed['interval_world'][0]:.1f}--{ed['interval_world'][1]:.1f}\\% for the extreme-value model ($B = {bnum(ev['resamples'])}$)"] if NEV else [])
    bt_txt = (
        f"To quantify sampling uncertainty across cohorts and donors, we employ a two-stage cohort-then-donor resampling scheme (Methods), which could also be called a bootstrap method. "
        f"In each resample, we draw the {len(by)} cohorts with replacement and, for each drawn cohort of size $n_c$, draw $n_c$ donors with replacement from that cohort; "
        f"we then take the 2.5th and 97.5th percentiles across resamples as the 95\\% interval (for the three curve-fitted models we recompute the measured curve and refit each model on every resample, whereas for the extreme-value model we estimate the distribution of tail shapes, $\\mathrm{{Beta}}(a, b)$, again from the donors drawn, each weighted by the number of times it is drawn, which updates the posterior of the shape of every donor, and average the resulting donor probabilities). "
        f"Resampling cohorts first accounts for between-cohort heterogeneity that an interval treating the {N} donors as independent would ignore. "
        f"For every model, we select the number of resamples $B$ adaptively by starting at $B = {bnum(bm['pycmb']['doubling_trace'][0]['B'])}$ and doubling $B$ until the Monte Carlo standard error of both interval bounds (estimated from {bt_meta['sub_bootstrap_draws']:,} draws with replacement among the $B$ resamples) falls below {100 * bt_meta['se_over_width_threshold']:.0f}\\% of the interval width. "
        f"This {100 * bt_meta['se_over_width_threshold']:.0f}\\% threshold provides a practical compromise: looser stopping rules leave noisy tail quantiles, whereas far tighter rules multiply refitting cost without materially altering the bounds. "
        f"Supplementary Fig. 6 illustrates both the resulting uncertainty bands ({'Panels A--C for the curve-fitted models and Panel D for the extreme-value model' if EVR else 'Panels A--C'}, showing point extrapolations as dashed lines, 95\\% bands as shaded regions and world-population intervals as vertical error bars) and the convergence diagnostics across doubling steps ({'Panels E--F' if EVR else 'Panels D--E'}, showing the relative Monte Carlo standard error falling below the 0.05 threshold and the bounds stabilising). "
        f"At $8\\times10^9$ runs, the resulting 95\\% intervals are {jn(bt_strs)}. "
        f"The power law of the miss rate displays a noticeably wider lower tail and requires a larger $B$ because uncertainty in its log--log slope compounds multiplicatively across {decades} orders of magnitude, yet even its 95\\% lower bound ({min_lo_bt:.1f}\\%) remains well above the 9--33\\% re-identification risk thresholds."
    )
    fig6 = rf"""\begin{{figure}}[t]\centering
\includegraphics[width=\textwidth,height=0.745\textheight,keepaspectratio]{{supplementary_figure6}}
\caption{{\textbf{{Cohort-then-donor resampling intervals for the {NUM[len(DRAWN)]} models{' fitted to the measured curve' if NEV else ''} and adaptive selection of $B$.}} \textbf{{A--C}}, Recovery rate against pool size $N$ measured up to 291,881 HGM runs (orange curve and shaded 95\% interval) and extrapolated by the Pitman--Yor correctness model (\textbf{{A}}), a straight line in $\log N$ (\textbf{{B}}) and a power law of the miss rate (\textbf{{C}}), each fitted for {geN(qd['fit_from'])}. Dashed lines, point extrapolations; shaded bands, 95\% intervals across $B$ cohort-then-donor resamples; vertical error bars, 95\% intervals at $8\times10^9$ runs. Blue area and lines, re-identification risk thresholds (9\%, 20\% and 33\%). \textbf{{D}}, Convergence of the maximum Monte Carlo standard error (from {bt_meta['sub_bootstrap_draws']:,} draws with replacement among the $B$ resamples) of the 2.5th and 97.5th percentiles at $8\times10^9$ runs relative to interval width as $B$ doubles from {bm['pycmb']['doubling_trace'][0]['B']:,}; ringed points mark the stopping $B$ below the {bt_meta['se_over_width_threshold']:.2f} threshold. \textbf{{E}}, Lower (2.5th percentile, solid) and upper (97.5th percentile, dashed) interval bounds at $8\times10^9$ runs with $\pm 1$ Monte Carlo standard error across doubling steps.}}
\end{{figure}}
"""
    if EVR:                                                                      # the extreme-value model in Supplementary Fig. 6: panel D and its points in the two lower panels
        er = json.load(open(EVR)); erm = er["models"]["evt"]; erw = erm["at_N"][str(int(8e9))]
        assert NEV and er["meta"]["tail"] == ev["drawn"] and er["meta"]["resamples"] == ev["resamples"] == erm["B"] and er["meta"]["resamples_drawn"] == ev["resamples_drawn"] and er["meta"]["sub_draws"] == bt_meta["sub_bootstrap_draws"], "the band is that of the tail drawn in Fig. 1B"
        assert abs(erw["est"] - ed["at_world"]) < 1e-6 and abs(erw["lo"] - ed["interval_world"][0]) < 1e-6 and abs(erw["hi"] - ed["interval_world"][1]) < 1e-6, "the interval of Supplementary Fig. 6D is that of Supplementary Table 18"
        assert er["meta"]["se_over_width_threshold"] == bt_meta["se_over_width_threshold"] and erm["doubling_trace"][0]["B"] == bm["pycmb"]["doubling_trace"][0]["B"] and erm["converged"] and erm["doubling_trace"][-1]["passed"] and erm["doubling_trace"][-1]["B"] == erm["B"], "the same rule as for the three other models, and the number of resamples that it selects"
        for old_, new_ in ((" for the three models fitted to the measured curve and adaptive selection of $B$.}", " for the four models of Fig. 1B and adaptive selection of $B$.}"),
                           (r"\textbf{A--C}, Recovery rate", r"\textbf{A--D}, Recovery rate"),
                           (r"(\textbf{C}), each fitted for " + geN(qd["fit_from"]) + ". Dashed lines", r"(\textbf{C}), each fitted for " + geN(qd["fit_from"]) + r", and by the extreme-value model of the similarity tail (\textbf{D}). Dashed lines"),
                           (r"vertical error bars, 95\% intervals at $8\times10^9$ runs. Blue area",
                            r"vertical error bars, 95\% intervals at $8\times10^9$ runs. In \textbf{A--C} each model is fitted again on each resample. In \textbf{D} the distribution of tail shapes is estimated again on each resample from the donors drawn, and the resample averages the resulting probabilities of these donors. Blue area"),
                           (r"\textbf{D}, Convergence", r"\textbf{E}, Convergence"), (r"\textbf{E}, Lower (2.5th percentile", r"\textbf{F}, Lower (2.5th percentile")):
            assert fig6.count(old_) == 1, (old_, fig6.count(old_)); fig6 = fig6.replace(old_, new_)
    iv = lambda d: f"{d['lo']:.0f}--{d['hi']:.0f}"                                       # 95% interval of the fit over cohorts
    fmt_m = lambda v: f"{int(round(v)):,}".replace(",", "{,}")
    med_s = int(round(sp["list_runs"]["median"]))
    p_med_1000 = math.exp(gammaln(M_HGM - med_s + 1) - gammaln(M_HGM - med_s - 999 + 1) - (gammaln(M_HGM + 1) - gammaln(M_HGM - 999 + 1)))
    s_short, s_full = sp["list_runs"]["short"], sp["list_runs"]["full"]              # runs of other people in the lists that stop at the similarity of 0.05, and in the lists of 5,000 runs
    assert p_med_1000 < 1e-7 and s_short and s_full and len(s_short) + len(s_full) == sp["not_found"] and max(s_short + s_full) < 5000
    assert min(s_short + s_full) == sp["list_runs"]["min"] and max(s_short + s_full) == sp["list_runs"]["max"]
    in_words = lambda n: {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}.get(n, str(n))
    s_short_tex = "$s_i = " + ("$ and $".join(["$, $".join(fmt_m(v) for v in s_short[:-1]), fmt_m(s_short[-1])]) if len(s_short) > 1 else fmt_m(s_short[0])) + "$"
    K6 = str(int(1e6)); K11 = str(int(1e11))
    pred_w_vals = [f1(bm[k_]["at_N"][WK]["est"]) + r"\%" for k_ in DRAWN] + ([f1(ed["at_world"]) + r"\%"] if NEV else [])
    pred_7_vals = [f1(bm[k_]["at_N"][K7]["est"]) + r"\%" for k_ in DRAWN] + ([f1(ed["at_1e7"]) + r"\%"] if NEV else [])
    pred_names = [NAME4[k_] for k_ in DRAWN] + (["the extreme-value model"] if NEV else [])
    evt_traj_txt = (
        f", while the extreme-value model ({f1(erm['at_N'][K6]['est'])}\\% to {f1(erm['at_N'][K11]['est'])}\\%) closely matches the Pitman--Yor and log-linear trajectories at world population ({f1(ed['at_world'])}\\%) and flattens smoothly as donors whose other sample sits well above their fitted upper tail retain recovery probabilities close to one"
        if (NEV and EVR) else (
            f", while the extreme-value model predicts {f1(ed['at_world'])}\\% at world population and flattens smoothly as donors whose other sample sits well above their fitted upper tail retain recovery probabilities close to one"
            if NEV else ""
        )
    )
    evt_censor_bridge = ""
    if NEV:
        k_nr = np.asarray(evs["not_returned"]["k"], dtype=float)
        def _add_exact(n_):
            m_ = n_ - 1
            with np.errstate(all="ignore"):
                lg = gammaln(M_HGM - k_nr + 1) - gammaln(M_HGM - k_nr - m_ + 1) - (gammaln(M_HGM + 1) - gammaln(M_HGM - m_ + 1))
            return 100.0 * float(np.where(M_HGM - k_nr >= m_, np.exp(lg), 0.0).sum() / N)
        d_ex_2, d_ex_1k, d_ex_10k = _add_exact(2), _add_exact(1000), _add_exact(10000)
        r_ex_2, r_ex_1k = sps["2"]["measured"] + d_ex_2, sps["1000"]["measured"] + d_ex_1k
        assert sps["2"]["measured"] <= r_ex_2 <= sps["2"]["upper"] and sps["1000"]["measured"] <= r_ex_1k <= sps["1000"]["upper"]
        evt_censor_bridge = (
            f"\n\nIn principle, this right-censoring could be removed: for the donor-level extreme-value model introduced below, "
            f"we recompute query similarities directly from the offline hypervector database underlying the Faure \\textit{{et al.}} service, "
            f"which yields exact counts $k_i$ for all {sp['not_found']} unreturned donors in $\\mathcal{{C}}$. "
            f"However, we deliberately retain the right-censored lower bound $R_{{\\mathrm{{low}}}}(N)$ for the measured curve in Fig. 1B and the {NUM[len(DRAWN)]} curve-fitted models. "
            f"Doing so keeps those results faithful to what the deployed search service actually returns when queried directly across all {M_HGM:,} HGM runs, "
            f"whereas the offline hypervector snapshot omits {evs['hgm_runs_absent']} HGM runs and is needed only to extract deep similarity tails beyond the service's 5,000-run limit. "
            f"As an independent check, substituting the exact offline counts $k_i$ for $i \\in \\mathcal{{C}}$ places the uncensored curve within $[R_{{\\mathrm{{low}}}}(N), R_{{\\mathrm{{up}}}}(N)]$ across all pool sizes: "
            f"it gives {r_ex_2:.1f}\\% at $N = 2$ and {r_ex_1k:.1f}\\% at $N = 1{{,}}000$, differing from $R_{{\\mathrm{{low}}}}(N)$ by only {d_ex_1k:.2f} percentage points at $N = 1{{,}}000$ and {d_ex_10k:.3f} percentage points at $N = 10{{,}}000$."
        )
    note4 = rf"""
\section*{{Supplementary Note 3. Extrapolation beyond the HGM set}}

In Fig. 1B, we measure the recovery rate in pools of up to 291,881 runs, the size of the HGM set, and extrapolate it to the population of the world using two complementary approaches: {NUM[len(DRAWN)]} aggregate curve-fitted laws from the privacy and biometrics literature and{' a fourth,' if NEV else ''} donor-level extreme-value model of the similarity upper tail. This Note details how this was done and is organised as follows. First, we derive exact bounds on the measured recovery curve to account for right-censoring in the search output. Next, we specify the curve-fitted scaling laws{' and the donor-level extreme-value tail model' if NEV else ''}. We then compare their population-scale predictions across pool sizes, before evaluating thinned-pool validation accuracy, sensitivity to fit windows, tail size and censoring, resampling uncertainty and cohort robustness across all {'four' if NEV else NUM[len(DRAWN)]} models together. Finally, we summarise the scope and shared methodological limitations of these extrapolations.

\subsection*{{Right-censoring bounds on the measured recovery curve}}

The measured curve in Fig. 1B is a lower bound on the true recovery rate $R(N)$. For each donor $i \in \{{1,\dots,{N}\}}$, denote by $k_i$ the number of runs of other people in the HGM set ($M = {fmt_m(M_HGM)}$) that are more similar to the query than donor $i$'s other sample (Methods). In a random pool of size $N$, donor $i$ is recovered with probability
\[
    p(k_i, N) = \frac{{\binom{{M - k_i}}{{N - 1}}}}{{\binom{{M}}{{N - 1}}}} \quad (\text{{or }} 0 \text{{ if }} M - k_i < N - 1), \qquad R(N) = \frac{{1}}{{{N}}}\sum_{{i=1}}^{{{N}}} p(k_i, N).
\]
Because the search returns at most 5,000 runs across all {INDEX_RUNS} indexed runs and stops at a similarity of 0.05, $k_i$ is right-censored for a set $\mathcal{{C}}$ of {sp['not_found']} donors whose other sample is not returned by the query. For each $i \in \mathcal{{C}}$, let $s_i$ denote the number of runs of other people in the HGM set present in the returned list; since every such run is more similar to the query than donor $i$'s other sample, we have $s_i \le k_i \le M$. Here $s_i < 5{{,}}000$ for all $i \in \mathcal{{C}}$, both because the search stops early at similarity 0.05 for {in_words(len(s_short))} donor{'s' if len(s_short) > 1 else ''} ({s_short_tex}) and because non-HGM runs and the donor's own runs occupy part of the 5,000 returned entries for the other {len(s_full)} donors ($s_i = {fmt_m(min(s_full))}$ to ${fmt_m(max(s_full))}$).

Since $p(k, N)$ decreases monotonically with $k$, replacing $k_i$ for $i \in \mathcal{{C}}$ by the two endpoints $M$ and $s_i$ brackets the true recovery rate:
\[
    R_{{\mathrm{{low}}}}(N) \le R(N) \le R_{{\mathrm{{up}}}}(N).
\]
Setting $k_i = M$ yields $p(M, N) = 0$ for all $N \ge 2$ and defines the measured lower bound $R_{{\mathrm{{low}}}}(N)$ plotted in Fig. 1B; conversely, setting $k_i = s_i$ assumes that no unreturned HGM run exceeds donor $i$'s other sample and defines the upper bound $R_{{\mathrm{{up}}}}(N)$. Because donors with known $k_i$ contribute identically to both bounds, their gap reduces to
\[
    \Delta R(N) = R_{{\mathrm{{up}}}}(N) - R_{{\mathrm{{low}}}}(N) = \frac{{1}}{{{N}}}\sum_{{i \in \mathcal{{C}}}} p(s_i, N).
\]
As $N$ grows, each term $p(s_i, N)$ decays rapidly towards zero (for the median $s_i = {fmt_m(med_s)}$, $p(s_i, 1{{,}}000) < 10^{{-7}}$), implying that the two bounds converge quickly. Specifically, we find $\Delta R(2) = {f1(gap(2))}$ percentage points ($R_{{\mathrm{{low}}}}(2) = {sps['2']['measured']:.1f}\%$, $R_{{\mathrm{{up}}}}(2) = {sps['2']['upper']:.1f}\%$), $\Delta R(100) = {f1(gap(100))}$ percentage point ($R_{{\mathrm{{low}}}}(100) = {sps['100']['measured']:.1f}\%$, $R_{{\mathrm{{up}}}}(100) = {sps['100']['upper']:.1f}\%$), and $\Delta R(1{{,}}000) = {f1(gap(1000))}$ percentage points ($R_{{\mathrm{{low}}}}(1{{,}}000) = {sps['1000']['measured']:.1f}\%$, $R_{{\mathrm{{up}}}}(1{{,}}000) = {sps['1000']['upper']:.1f}\%$), with $\Delta R(N) \le {f1(gap(10000))}$ percentage points for $N \ge 10{{,}}000$.{evt_censor_bridge}

\subsection*{{Extrapolation models}}

\subsubsection*{{Curve-fitted scaling laws}}

To extrapolate beyond $M = {fmt_m(M_HGM)}$ runs, we evaluate the measured curve on a grid of {f3meta['grid_sizes']} distinct pool sizes between 2 and 291,881 ({f3meta['grid_values']} log-spaced values rounded to whole numbers; Methods) and fit {NUM[len(DRAWN)]} complementary laws from the privacy and biometrics literature: (i) the Pitman--Yor correctness model (PYC-MB) of Rocher \textit{{et al.}}\cite{{rocher2025}}, a two-parameter Bayesian scaling law developed specifically in the privacy literature to forecast population-scale re-identification risk from finite galleries; (ii) a straight line in $\log N$, the first-order asymptotic form of rank-one identification against gallery size in biometrics\cite{{baveja2010}}; and (iii) a power law of the miss rate, $1 - a N^{{b}}$, where the miss rate is the share of unrecovered donors, as reported by the National Institute of Standards and Technology for face identification\cite{{grother2019}}. In Fig. 1B we fit all {NUM[len(DRAWN)]} models to $R_{{\mathrm{{low}}}}(N)$ over the {qd['fit_points']} grid sizes satisfying {geN(qd['fit_from'])}. This fit window represents a reasoned compromise: starting in smaller pools exposes the fit to right-censoring distortion, whereas starting later (e.g., at $N \ge 10{{,}}000$) restricts the fit to fewer grid points and inflates extrapolation variance. At {geN(qd['fit_from'])}, the censoring gap is already at most {f1(gap(1000))} percentage points and the large-pool scaling regime applies.

Among curve-fitted laws, the Pitman--Yor correctness model of Rocher \textit{{et al.}}\cite{{rocher2025}} directly targets the inference problem considered here: forecasting the fraction of individuals correctly singled out in a large population from a smaller empirical benchmark. Rather than fitting an ad-hoc curve, it derives a closed-form expected correctness $\kappa(N) = \mathrm{{E}}[K_N]/N$ from a generative Pitman--Yor process $\mathrm{{PY}}(d, \alpha)$ partitioning a gallery of size $N$ into anonymity sets, reparametrised by two interpretable quantities, the entropy $h$ (effective distinguishing information) and the tail complexity $\gamma$ (equivalently the discount $d$), and fitted by log-weighted least squares. Rocher \textit{{et al.}} validated this two-parameter law across 476 exact, sparse and robust identification curves. In this taxonomy, k-mer sketch matching of gut metagenomes acts as a robust matching technique. Fitted to our measured curve for {geN(qd['fit_from'])}, PYC-MB yields a unique interior optimum with $h = {py_h_nats:.0f}$~nats ($\approx {py_h_bits:.0f}$~bits), $\gamma = {py_gamma:.2f}$ and discount $d = {py_d:.3f}$. Mechanistically, because $\kappa(N)$ decays asymptotically as $N^{{d-1}}$ under a Pitman--Yor process, a discount close to unity ($d - 1 \approx {py_d - 1:.3f}$) directly accounts for the gentle decline of recovery across orders of magnitude. Nevertheless, because the underlying theory is derived for exact discrete matching and extended empirically to similarity-based matching (and still relies solely on the aggregate curve $R(N)$), we cross-check it against the two biometric laws and against the donor-level extreme-value model below.

Over {geN(qd['fit_from'])}, the straight line in $\log N$ stays within {f1(qd['max_abs_residual_first'])} percentage points of the measured curve and declines by {f1(-100 * qd['coef_first'][0])} percentage points per doubling of the pool. By contrast, the power law of the miss rate ($1 - a N^{{b}}$) models the unrecovered share as growing multiplicatively with $N$, producing a steeper decline at extreme pool sizes.

{evt_model_sec}\subsection*{{Population-scale predictions across all {'four' if NEV else NUM[len(DRAWN)]} models}}

Extrapolated to {decades} orders of magnitude beyond the HGM set to the world population ($8\times10^9$ runs), all {'four' if NEV else NUM[len(DRAWN)]} models remain in close agreement, predicting recovery rates of {jn(pred_w_vals)} for {jn(pred_names)}, respectively ({jn(pred_7_vals)} at $10^7$ runs; Supplementary Table 15). Across pool sizes from $10^6$ to $10^{{11}}$ runs, the Pitman--Yor correctness model ({f1(bm['pycmb']['at_N'][K6]['est'])}\% to {f1(bm['pycmb']['at_N'][K11]['est'])}\%) and the straight line in $\log N$ ({f1(bm['line']['at_N'][K6]['est'])}\% to {f1(bm['line']['at_N'][K11]['est'])}\%) track one another closely, whereas the power law of the miss rate declines more steeply ({f1(bm['miss']['at_N'][K6]['est'])}\% to {f1(bm['miss']['at_N'][K11]['est'])}\%){evt_traj_txt}.

\subsection*{{Validation, sensitivity, uncertainty and robustness across all {'four' if NEV else NUM[len(DRAWN)]} models}}

To evaluate extrapolation accuracy against a known ground truth, we thinned the HGM set by factors of 10, 33 and 100 and used fits on the thinned sets alone to predict the recovery rate measured on the full HGM set (Supplementary Table 16). Specifically, we kept each run of another person independently with probability $f$ ({fs11}), yielding thinned sets of size $M' \approx f M$; for the curve-fitted models we scaled the lower bound of each fit window $[N_{{\min}}, M]$ proportionally to $\max(2, f N_{{\min}})$ so that the window covers the same relative range of pool sizes (starting at {jn(starts)} runs for the 1,000-run window at $f$ = {fs11}) and extrapolated 10-, 33- or 100-fold to 291,881 runs. {better_txt}{evt_thin_txt}

To assess sensitivity to fit windows{', tail size $r$' if NEV else ''} and right-censoring, we refitted all {NUM[len(DRAWN)]} curve-fitted models across {NUM[len(WINS)]} fit windows starting from $N$ = {', '.join(f'{w_:,}' for w_ in WINS[:-1])} or {WINS[-1]:,} up to 291,881 runs (Supplementary Fig. 5, first three rows; thick lines denote measured rates, darker over fitted windows, dashed lines show extrapolations to $10^{{12}}$ runs, and crosses mark $8\times10^9$){evt_sens_intro}. Across all {NUM[len(v20)]} curve-fitted specifications, world-population predictions cluster tightly between {g(min(v20))}\% and {g(max(v20))}\%, with the Pitman--Yor and log-linear trajectories tracking one another closely.{evt_sens_txt} We also repeated each curve-fitted specification on the upper censoring bound $R_{{\mathrm{{up}}}}(N)$ ($k_i = s_i$ for $i \in \mathcal{{C}}$): when fitted from {LO0} runs to $R_{{\mathrm{{up}}}}(N)$ rather than $R_{{\mathrm{{low}}}}(N)$, the Pitman--Yor model and the straight line shift to {g(sp8w['upper'][f'pycmb|{LO0}'])}\% and {g(sp8w['upper'][f'line|{LO0}'])}\% at $8\times10^9$ runs (versus {g(sp8w['measured'][f'pycmb|{LO0}'])}\% and {g(sp8w['measured'][f'line|{LO0}'])}\%) and the power law exceeds a 100\% miss rate, whereas fits starting at 1,000 runs yield {jn([f1(sp8w['upper'][k_ + '|1000']) + chr(92) + '%' for k_ in DRAWN])} on $R_{{\mathrm{{up}}}}(N)$ versus {jn([f1(sp8w['measured'][k_ + '|1000']) + chr(92) + '%' for k_ in DRAWN])} on $R_{{\mathrm{{low}}}}(N)$ (shifting by at most {f1(dw_)} percentage points).{evt_censor_txt}

{bt_txt}

Beyond thinning accuracy, tuning sensitivity and sampling uncertainty, we verified across all {'four' if NEV else NUM[len(DRAWN)]} models that our extrapolations are driven neither by any single cohort nor by the exclusion of the {len(rows_ex)} donors removed under the contamination and re-deposit rules, by recomputing the measured rate and the Fig. 1B fits after leaving each cohort out in turn as well as with all {N + len(rows_ex)} donors retained (Supplementary Table 17). Across leave-one-out subsets, the measured recovery rate varies narrowly between {min(rate14.values()):.1f}\% and {max(rate14.values()):.1f}\%, while world-population extrapolations remain between {g(min(v_[k_] for v_ in loo.values() for k_ in DRAWN))}\% and {g(max(v_[k_] for v_ in loo.values() for k_ in DRAWN))}\% across the {NUM[len(DRAWN)]} curve-fitted models{evt_loo_txt}; the measured rate is lowest without {lo14_txt}, and highest without {NAMES[hi14]}, which accounts for {fail14} of the {N - n_stock} unrecovered donors.{ex497_txt}

{evt_diag_txt}Taken together, these validation, sensitivity, resampling and cohort-removal checks confirm that the high extrapolated recovery rates are consistent across model families and robust to sampling variation and data selection, remaining well above the 9--33\% re-identification risk thresholds in every scenario tested.

\subsection*{{Scope and limitations of the extrapolations}}

Four methodological limitations apply across these extrapolations and warrant caution when interpreting predictions at population scale. First, every model assumes that newly added runs drawn into larger pools follow the same similarity distribution as the HGM set; outliers far above the empirical distribution, such as close relatives or untracked re-deposits of the donor, are not represented. Second, extrapolating to the world population spans roughly {round(8e9 / M_HGM, -3):,.0f} times the HGM set, whereas thinned-pool validation covers at most 100-fold growth (a regime in which the smallest tails of the extreme-value model retain too few excesses to be fitted). Third, behaviour far beyond $8\times10^9$ runs remains sensitive to parametric tail assumptions{evt_lim_clauses[0]}, so these curves should not be extrapolated arbitrarily beyond human population scale. Fourth, the reported 95\% intervals capture cohort and donor sampling uncertainty only{evt_lim_clauses[1]}.

"""
    NLAB = [(1e6, "$10^6$"), (1e7, "$10^7$"), (1e8, "$10^8$"), (1e9, "$10^9$"), (8e9, "$8\\times10^9$"), (1e10, "$10^{10}$"), (1e11, "$10^{11}$")]
    pyb, lnb, msb = bm["pycmb"]["at_N"], bm["line"]["at_N"], bm["miss"]["at_N"]
    assert bm["pycmb"]["B"] == qd["resamples"] and all(str(int(n_)) in d_ for d_ in (pyb, lnb, msb, qf) for n_, _ in NLAB)
    assert all(abs(lnb[str(int(n_))][f_] - qf[str(int(n_))][f_]) < 1e-9 for n_, _ in NLAB for f_ in ("est", "lo", "hi"))
    evb15 = (erm["at_N"] if EVR else {K7: dict(est=ed["at_1e7"], lo=ed["interval_1e7"][0], hi=ed["interval_1e7"][1]), WK: dict(est=ed["at_world"], lo=ed["interval_world"][0], hi=ed["interval_world"][1])}) if NEV else {}
    evt_miss15 = NEV and any(str(int(n_)) not in evb15 for n_, _ in NLAB)
    c15 = lambda d_, k_: f"{d_[k_]['est']:.1f} & {iv(d_[k_])}" if k_ in d_ else "-- & --"
    T8 = "\n".join(" & ".join([lab] + [c15(d_, str(int(n_))) for d_ in ([pyb, lnb, msb] + ([evb15] if NEV else []))]) + " \\\\" for n_, lab in NLAB)
    n9_ = lambda v: "above $10^{40}$" if not (v is not None and v < float("inf")) else sci(v)
    bf = lambda s_, on: f"\\textbf{{{s_}}}" if on else s_
    T9 = []
    for i_, (k_, w_) in enumerate(DRAWN.items()):
        for j_, lo_ in enumerate(WINS):
            c_ = f"{k_}|{lo_}"; on_ = lo_ == w_; v_ = f3g[c_]
            up_ = v_["first_8e9_upper"]; assert abs(up_ - sp8w["upper"][c_]) < 1e-6, c_
            cells = [NAME9[k_] if j_ == 0 else "", f"{lo_:,}"] + [f"{v_['first_8e9']:.1f}", f"{up_:.1f}" if up_ >= 0 else "$<$0", n9_(v_["N9"])] + [sg(hw[f]["error"][c_]) for f in ho] + [f"{rmse(c_, '0.01'):.1f}"]
            T9.append(" & ".join(cells[:1] + [bf(x_, on_) for x_ in cells[1:]]) + " \\\\" + ("\n\\addlinespace" if j_ == len(WINS) - 1 and i_ < len(DRAWN) - 1 else ""))
    if NEV:                                                                      # the extreme-value model, one row per tail, with the errors of its own test (Supplementary Table 18)
        T9[-1] += "\n\\addlinespace"
        for j_, r_ in enumerate(tails):
            e_ = lambda f_: sg(et[f_][r_]["error"]) if r_ in et[f_] else "--"
            cells = ["Extreme-value model" if j_ == 0 else "", f"{int(r_):,}", f"{em[r_]['at_world']:.1f}", "--", "--", e_(F10), e_(F33), e_(F100), f"{et[F100][r_]['rmse']:.1f}" if r_ in et[F100] else "--"]
            T9.append(" & ".join(cells[:1] + [bf(x_, r_ == dr) for x_ in cells[1:]]) + " \\\\")
    EV16 = (r"For the extreme-value model the value at $8\times10^9$ is given under $R_{\mathrm{low}}$, as its curve starts from the measured recovery rate, and the errors are those of its own test, "
            r"in which the count $k_i$ of every donor is used and the tail keeps the same share of the runs of other people (Supplementary Table 18). --, not computed for this model, or a tail of fewer than ten similarities in the thinned HGM set, which is not fitted. ") if NEV else ""
    EV16B = " This baseline is that of the test of the three models fitted to the measured curve, and Supplementary Table 18 gives that of the test of the extreme-value model." if NEV else ""
    WHDR = r"\makecell[lb]{Window\\or tail}" if NEV else "Window"                # the second column: the fit window of a model fitted to the measured curve, or the tail of the extreme-value model
    T9.append("\\midrule\n\\multicolumn{2}{l}{Constant-rate baseline (thinned HGM set)} & -- & -- & -- & " + " & ".join(sg(decl[f]) for f in ho) + f" & {math.hypot(decl['0.01'], decl_sd['0.01']):.1f} \\\\")
    ev_cap15 = (rf", and the extreme-value model of the similarity tail fitted to the {int(dr):,} highest similarities of each donor" if NEV else "")
    ev_bt15 = (rf", and $B = {bnum(ev['resamples'])}$ for the extreme-value model (whose distribution of tail shapes is estimated again on each resample)" if NEV else "")
    ev_head15 = ' & \\multicolumn{2}{c}{Extreme value}' if NEV else ''             # kept out of the f-string below: Python 3.11, the version of environment.txt, refuses a backslash in an f-string expression
    ev_rule15 = '\\cmidrule(lr){8-9}' if NEV else ''
    ev_cols15 = ' & \\% & Interval' if NEV else ''
    tab8 = rf"""
\begin{{table}}[!ht]\centering\small\setlength{{\tabcolsep}}{{4pt}}
\caption{{\textbf{{Extrapolation of the recovery curve of Fig. 1B across all {'four' if NEV else 'three'} models.}} Share of the {ex['donors']} donors whose other sample is first in a pool of $N$ runs, in the order of the legend of Fig. 1B: the Pitman--Yor correctness model (PYC-MB) of Rocher \textit{{et al.}}\cite{{rocher2025}}, a straight line in $\log N$ and a power law of the miss rate\cite{{grother2019}}, each fitted to the measured curve at the {qd['fit_points']} pool sizes from {qd['fit_from']:,} to 291,881 (measured share at $N$ = 291,881, {f1(qd['measured_first_at_M'])}\%){ev_cap15}. \textit{{Interval}}, 95\% interval on cohort-then-donor resamples with adaptive $B$ (Supplementary Fig. 6): $B = {bnum(bm['pycmb']['B'])}$ for Pitman--Yor, $B = {bnum(bm['line']['B'])}$ for the straight line and $B = {bnum(bm['miss']['B'])}$ for the power law{ev_bt15}.{' --, not yet computed at this pool size.' if evt_miss15 else ''}}}
\begin{{tabular}}{{l{'rc' * (3 + NEV)}}}
\toprule
 & \multicolumn{{2}}{{c}}{{Pitman--Yor}} & \multicolumn{{2}}{{c}}{{Straight line}} & \multicolumn{{2}}{{c}}{{Power law}}{ev_head15} \\
\cmidrule(lr){{2-3}}\cmidrule(lr){{4-5}}\cmidrule(lr){{6-7}}{ev_rule15}
$N$ & \% & Interval & \% & Interval & \% & Interval{ev_cols15} \\
\midrule
{T8}
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[!ht]\centering\small\setlength{{\tabcolsep}}{{3.8pt}}
\caption{{\textbf{{Predictions of the {'four' if NEV else 'three'} models of Fig. 1B, for five fit windows{' or five tails,' if NEV else ''} and test on thinned HGM sets.}} The {NUM[len(v20)]} fits{' and the ' + NUM[len(tails)] + ' tails' if NEV else ''} of Supplementary Fig. 5, in the order of the legend of Fig. 1B; bold, the fit{' and the tail' if NEV else ''} drawn in Fig. 1B. \textit{{Window{' or tail' if NEV else ''}}}, {'for the three models fitted to the measured curve, ' if NEV else ''}the fit uses the pool sizes of the measured curve from $N \ge$ the value given to 291,881{'; for the extreme-value model, the number of highest similarities of each donor that make the tail' if NEV else ''}. \textit{{At $8\times10^9$}}, extrapolated recovery rate at the population of the world ($8\times10^9$ runs) when fitted to the measured lower bound ($R_{{\mathrm{{low}}}}$) or to the upper bound ($R_{{\mathrm{{up}}}}$; Supplementary Note 3; $<$0, a miss rate above 100\%). $N_9$, pool size at which the recovery curve of the fit falls to 9\%, the re-identification risk threshold of Fig. 1B (above $10^{{40}}$, the curve is still above 9\% at $10^{{40}}$ runs, the largest pool computed). Test, each run of another person in the HGM set is kept with probability $f$ ({fs11}, HGM sets of {mp11} runs) and each fit is made on the thinned set alone by scaling the lower bound of its window $[N_{{\min}}, M]$ proportionally to $\max(2, f N_{{\min}})$. At 100 times the windows from 2, 10 and 100 runs therefore all start at 2 runs. In this test the window from 2 runs is fitted on more pool sizes than the other windows, so its errors differ slightly from those of the windows from 10 and 100 runs. \textit{{Mean error}}, the recovery rate that the fit predicts at 291,881 runs, that is 10, 33 and 100 times the thinned HGM set for the three values of $f$, minus the measured share ({f1(list(ho.values())[0]['measured'])}\%), in percentage points, as a mean over {list(ho.values())[0]['thinnings']} thinnings; positive, the fit overestimates. \textit{{RMSE}}, root mean square error over the thinnings of the prediction at 100 times, in percentage points, given in parentheses in {'the first three rows of ' if NEV else ''}Supplementary Fig. 5. {EV16}Last row, constant-rate baseline assuming no further decline beyond the thinned HGM set (predicting the share measured in the thinned HGM set itself, as an expected value over the thinnings).{EV16B}}}
\begin{{tabular}}{{llrrrrrrr}}
\toprule
 & & \multicolumn{{2}}{{c}}{{At $8\times10^9$ (\%)}} & & \multicolumn{{3}}{{c}}{{\makecell[cb]{{Mean error\\(percentage points)}}}} & RMSE \\
\cmidrule(lr){{3-4}}\cmidrule(lr){{6-8}}
Model & {WHDR} & $R_{{\mathrm{{low}}}}$ & $R_{{\mathrm{{up}}}}$ & $N_9$ & {hdr11} & $\times$100 \\
\midrule
""" + "\n".join(T9) + r"""
\bottomrule
\end{tabular}
\end{table}
""" + rf"""
\begin{{table}}[!ht]\centering\small\setlength{{\tabcolsep}}{{{3.4 if NEV else 4.5}pt}}
\caption{{\textbf{{Recovery rate and extrapolations with one cohort left out.}} Cohorts are named by their curatedMetagenomicData identifiers, in the order of Supplementary Table 11. \textit{{Donors}}, donors kept. \textit{{Recovered}}, donors whose other sample is at rank 1 in the HGM set. \textit{{At $8\times10^9$}}, the {NUM[len(DRAWN)]} fits drawn in Fig. 1B (bold in Supplementary Table 16), namely {jn([SHORT14[k_] for k_ in DRAWN])}, made again on the measured curve of the donors kept{', and the extreme-value model of the similarity tail (bold in Supplementary Table 18), made again on the donors kept' if NEV else ''}, as shares in \%. \textit{{None}}, all the {N} donors.{six_cap}{' --, not computed.' if NEV and rows_ex and not EVX else ''}}}
\begin{{tabular}}{{lrr{'r' * (len(DRAWN) + NEV)}}}
\toprule
 & & & \multicolumn{{{len(DRAWN) + NEV}}}{{c}}{{At $8\times10^9$ (\%)}} \\
\cmidrule(lr){{4-{3 + len(DRAWN) + NEV}}}
Cohort left out & Donors & Recovered (\%) & {' & '.join([HEAD14[k_] for k_ in DRAWN] + ['Extreme value'] * NEV)} \\
\midrule
""" + "\n".join(T14c) + r"""
\bottomrule
\end{tabular}
\end{table}
""" + evt_tab


SA_MORE = sc["sample_arrangements"]["more"] if sc else []                           # cohort_selection.py: the donors with more stool samples than Supplementary Fig. 1 lists one by one
assert all(m_["studies"] == ["HMP_2019_ibdmdb"] for m_ in SA_MORE), "the donors with more than four stool samples are all of HMP IBD"
SA_LO, SA_HI = (min(m_["min_samples"] for m_ in SA_MORE), max(m_["max_samples"] for m_ in SA_MORE)) if SA_MORE else (0, 0)
figsel = rf"""\begin{{figure}}[t]\centering
\includegraphics[width=\textwidth,height=0.79\textheight,keepaspectratio]{{supplementary_figure1}}
\caption{{\textbf{{Selection of the cohorts, the donors and their samples.}} \textbf{{A}}, From curatedMetagenomicData to the donors of the analysis, each step with the cohorts or studies and the samples or donors that it keeps (left, as in Supplementary Table 1) and what it removes (right, red; Supplementary Note 1 and Supplementary Table 2). \textbf{{B}}, The samples of each donor as circles along time, large for the samples studied, namely its query sample, its first sample and its most similar other sample (Supplementary Note 1), and small for its other samples. A line is one arrangement of these samples with its number of donors, left for the donors recovered and right for the others. The last line of each side stands for the donors of HMP IBD, who have {SA_LO} to {SA_HI} stool samples. A donor is a subject identifier of the curated metadata. A curated run of a sample is a run that the curated metadata list for it; HMP IBD, whose curated metadata list none, takes the runs of its project in the European Nucleotide Archive (ENA), and its first and last samples are the earliest and the latest with a run in the index.}}
\end{{figure}}
"""
fig4 = rf"""\begin{{figure}}[t]\centering
\includegraphics[width=\textwidth,height=0.71\textheight,keepaspectratio]{{supplementary_figure3}}
\caption{{\textbf{{Human DNA in the samples of each cohort and share of donors recovered.}} \textbf{{A}}, Share of the hashes of the sketch of each sample that match a human genome (T2T-CHM13v2.0 or GRCh38.p14), for the {vx['pool']['pool_donors']} donors of {len(by) + 1} cohorts, namely the {len(by)} cohorts kept and Vincent, which was left out (red). A dot is a sample and a thin line joins the two samples of a donor; the bar is the median, over the donors of the cohort, of the more human of their two samples; the dotted line marks 5\%; shares under 0.001\% are drawn at 0.001\%. \textbf{{B}}, Share of the {vx['pool']['pool_donors']} donors whose other sample is the nearest of the 291,881 runs of the HGM set, by the human share of the more human of their two samples (95\% Wilson intervals; the number of donors is under each bar; red, the two classes at 5\% or more). \textbf{{C}}, Share of donors whose other sample is the nearest of the {vx['pool']['pool_runs']:,} runs of the {vx['pool']['pool_donors']} donors (the two samples of each donor), with all hashes, without the hashes that match a human genome and with these hashes only, for the {vx['pool']['poor_all']['n']} donors whose more human sample holds under 5\% human hashes and the {vx['pool']['rich_all']['n']} donors with 5\% or more. The sketches were recomputed from the public unitigs of the runs (Methods).{stat_cap}}}
\end{{figure}}
"""

# ---- Supplementary Note 1, Supplementary Fig. 2 and Supplementary Tables 6 and 7: the links between Nielsen and Li
note5 = fig5 = tab1516 = ""
if lk and FIG2:
    F, R, SP, SQ, CU, MD = lk["forward"], lk["reverse"], lk["species"], lk["sequencing"], lk["curation"], lk["metadata_cmd"]
    xc5 = f2["nielsen_li_cross"]; nshared = lk["shared_subject_ids"]; bcs, aits, spm = SP["Bray-Curtis"], SP["Aitchison"], SP["Spearman"]; wn = SP["within_nielsen_bc"]
    assert SQ["N"]["first_created"] == ["2012-12"] and SQ["L"]["first_created"] == ["2014-01", "2014-05"] and lk["libraries_shared"] == 0, "the text of Note 1 gives these dates and says that no library is shared"
    nl_before, nl_after = nl["nielsen_rank1_before"], nl["nielsen_rank1_after"]
    OB = lk["outranked_by"]; QR = lk["queried_runs"]; SE = lk["sensitivity"]; ac5, mg5, kp5 = SE["with_two_cohorts_as_counted"], SE["with_two_cohorts_merged"], SE["kept_first"]
    p1 = lambda a, b: f"{100 * a / b:.1f}"
    note5 = rf"""
The {nshared} patients all have inflammatory bowel disease and come from the Spanish MetaHIT collection. For {lk['time_points']['two_nielsen_one_li']} of the {nshared} patients Nielsen holds two samples (time points 0 and 2, by their names) and Li one (time point 1). The other {lk['time_points']['others']} patients have a single Nielsen sample and one Li sample, and one of them has a second Li sample. The dates of the Li samples are not recorded. The samples of the two cohorts were deposited in two studies (ENA projects PRJEB1220 for Nielsen and PRJEB5224 for Li), which share no run accession, no sample accession and no library name. The supplementary table of Li lists the Nielsen samples of these individuals under the same individual code. Of the {F['donors']} Nielsen donors, {F['with_li_sample']} are among these patients.

When a Nielsen donor queries the index, the Li sample with the same identifier is the nearest run of the HGM set for {F['rank1']} of the {F['with_li_sample']} donors, apart from the donor's own Nielsen samples. It is among the 5,000 nearest runs that the Faure \textit{{et al.}}\cite{{faure2026}} service returns for {F['in_list']} of the {F['with_li_sample']} (Supplementary Fig. 2A). In the other direction, when the {R['queries']} Li samples of the {R['patients']} patients query the index, {R['rank1']} of them find a Nielsen sample of their identifier first. For {R['mutual_rank1_among_donor_patients']} of the {F['with_li_sample']} donors the sample of the other cohort is first in both directions. When the {F['with_li_sample']} Nielsen donors query the index, the Li sample is more similar to the query than the donor's own other sample for {F['li_closer_than_within']} of them. The median similarity to the query is {F['median_p']:.2f} for the Li sample, {F['median_within_p']:.2f} for the donor's own other sample and {F['median_best_impostor']:.2f} for the nearest run of other people when the Li sample of the same identifier is not counted as one (Supplementary Fig. 2B). The highest similarity between a query and a sample of the same identifier in the other cohort, in either direction, is {lk['max_similarity_to_same_identifier_sample']:.2f}, {costea_compare}far from the 0.99 at which the pipeline removes a run as a quasi-exact re-deposit. With the Li sample counted as another person, the donor's other sample is first for only {F['within_rank1']} of these {F['with_li_sample']} donors.

The queried runs are paired runs of Genome Analyzer IIx instruments for the {QR['N']['n']} Nielsen donors and of HiSeq 2000 instruments for the {QR['L']['n']} Li donors (median {QR['N']['gbp_median']:.1f} and {QR['L']['gbp_median']:.1f} Gbp). Over the whole cohorts, the paired runs of Nielsen have a median of {SQ['N']['mate_median']:.0f} bp per mate, compared with {SQ['L']['mate_median']:.0f} bp for Li (Supplementary Table 6). The link thus holds across a change of instrument, read length and depth, against runs of other people that include those of the other patients of Nielsen.

\textbf{{Resemblance between samples with the same identifier.}} We looked for resemblance measured without the sketch (Supplementary Table 7). Species profiles find a Nielsen sample of the same identifier first among the {SQ['N']['samples']} Nielsen samples for {bcs['li_to_nielsen_rank1']} of the {bcs['li_to_nielsen_n']} Li samples with the Bray-Curtis dissimilarity, {aits['li_to_nielsen_rank1']} with the Aitchison distance and {spm['li_to_nielsen_rank1']} with the Spearman correlation. They separate the pairs of a Nielsen and a Li sample of the same identifier from all the other pairs of a Nielsen and a Li sample with an area under the curve of {bcs['auc']:.2f}, {aits['auc']:.3f} and {spm['auc']:.2f} (Supplementary Table 7; Supplementary Fig. 2C for the Bray-Curtis dissimilarity). The Bray-Curtis dissimilarity of a Nielsen and a Li sample of the same identifier (median {bcs['median_matched']:.2f}) is close to that of two samples of one Nielsen patient ({wn['median_same_patient']:.2f}) and below that of the samples of two Nielsen patients ({wn['median_different_patients']:.2f}). The sex, age and body mass index recorded for the two cohorts agree (sex for {MD['sex']} of {MD['patients']} patients, age within one year for {MD['age_le1']}, body mass index within 0.5 kg per m$^2$ for {MD['bmi_le05']}), but the supplementary table of Li gives no clinical variable and the values curated for Li equal those of the published table of Nielsen.

We keep Nielsen and Li out of the main analysis for two reasons. (1) The links cannot be separated from the identifiers that define them (previous paragraph). (2) We can resolve the overlap of identifiers only in part (above), and the recovery rate of Nielsen depends on how the samples under a shared identifier are counted. Counting a Li sample under the donor's identifier as the donor's own, and no longer as another person's, raises the recovery rate of Nielsen from {pct(nl_before, nl['nielsen_donors'])}\% to {pct(nl_after, nl['nielsen_donors'])}\% ({nl_before} and {nl_after} of {nl['nielsen_donors']} donors; Supplementary Fig. 2D,E).
"""
    fig5 = rf"""\begin{{figure}}[t]\centering
\includegraphics[width=\textwidth]{{supplementary_figure2}}
\caption{{\textbf{{Links between the samples of Nielsen and Li.}} \textbf{{A}}, Rank, among the 291,881 runs of the HGM set, of the most similar sample of the other cohort with the same subject identifier, when a donor of Nielsen queries the index (the {F['with_li_sample']} of the {F['donors']} Nielsen donors whose identifier has a sample in Li) and when a Li sample queries it (the {R['queries']} Li samples whose identifier has a sample in Nielsen); the other runs under that identifier, in either cohort, are not counted as runs of other people. \textbf{{B}}, Similarity to the query, for the {F['with_li_sample']} Nielsen donors of A, of the nearest run of other people among the runs of the HGM set (a sample of Li under the donor's identifier is not counted as another person here), of the donor's other sample in Nielsen and of the sample of the same identifier in Li (a line joins these last two samples of a donor when both are drawn; bar, median). For {sum(v == 0 for v in F['within_p'])} donors the other sample is not among the 5,000 nearest runs, and it is drawn at 0 and enters the median as 0. The sample of Li is among the 5,000 nearest runs for {F['in_list']} of the {F['with_li_sample']} donors, and only these {F['in_list']} are drawn and enter the median. \textbf{{C}}, Bray-Curtis dissimilarity between species profiles, for a Nielsen and a Li sample with the same subject identifier (a dot per pair), two samples of one Nielsen patient (a dot per pair), and all the other pairs of a Nielsen and a Li sample and of two Nielsen samples (box, interquartile range; whiskers, 5th to 95th percentile). The bar is the median of each of the four groups. \textbf{{D}}, Share of the {nl['nielsen_donors']} Nielsen donors and of the {nl['li_donors']} Li donors whose other sample is first among the 291,881 runs, when a sample of the other cohort under the donor's identifier counts as another person (as counted) and when it counts as the donor's own (identifiers merged); the band is the range of this share across the {len(by)} cohorts kept. \textbf{{E}}, Outcome for the {OB['first'] + OB['li_same_identifier'] + OB['li_other_identifier'] + OB['nielsen_other_person'] + OB['other_run'] + OB['not_found']} Nielsen donors when Li samples count as other people: the other sample is first, another sample has higher similarity (classed by the nearest run of other people: a sample of Li with the donor's identifier, a sample of Li with another identifier, another patient of Nielsen or a run outside the two cohorts), or it is not among the 5,000 nearest runs.}}
\end{{figure}}
"""
    ev = lambda m: f"{m['li_to_nielsen_rank1']} of {m['li_to_nielsen_n']}; AUC {m['auc']:.3f}"
    T15 = "\n".join([
        f"Sketch index, a Nielsen donor queries ({F['with_li_sample']} donors) & Li sample of the same identifier first for {F['rank1']} & Same labels and same index in both directions \\\\",
        f"Sketch index, a Li sample queries ({R['queries']} samples) & Nielsen sample of the same identifier first for {R['rank1']}; both directions for {R['mutual_rank1_among_donor_patients']} of {F['with_li_sample']} donors & Same labels and same index; deeper runs queried \\\\",
        f"Species profiles, Bray-Curtis & Nielsen sample first among {SQ['N']['samples']} for {ev(bcs)} & Alike composition, not identity \\\\",
        f"Species profiles, Aitchison & First for {ev(aits)} & As above; depends on a pseudocount \\\\",
        f"Species profiles, Spearman & First for {ev(spm)} & As above \\\\",
        f"Bray-Curtis dissimilarity, median & {bcs['median_matched']:.2f} (same identifier, Nielsen and Li), {wn['median_same_patient']:.2f} (two samples of one Nielsen patient), {wn['median_different_patients']:.2f} (two Nielsen patients), {bcs['median_other_cross']:.2f} (other pairs of Nielsen and Li) & Same consortium and disease \\\\",
        f"CroCoDeEL on the species tables of the two cohorts together & {xc5['all']['events']} events between a Nielsen and a Li sample, {xc5['same']['events']} with the same identifier and {xc5['other']['events']} with different identifiers & Alike profiles; resemblance not limited to samples of one identifier \\\\",
        f"Curated sex, age, body mass index & Sex {MD['sex']} of {MD['patients']}, age within one year {MD['age_le1']}, body mass index within 0.5 kg per m$^2$ {MD['bmi_le05']} & The values of Li equal those of the published table of Nielsen \\\\",
        f"Run accessions, sample accessions or library names in common & None; highest similarity to a sample of the same identifier {lk['max_similarity_to_same_identifier_sample']:.2f} & Reads not compared, so a re-deposit is not excluded \\\\",
        "Host genotype, strain tracking, family & Not done; no family information in the metadata & Identity not confirmed \\\\",
    ])
    fm = lambda v: f"{v[1]:.0f}--{v[2]:.0f}"
    T16 = "\n".join([
        "ENA project & PRJEB1220 & PRJEB5224 \\\\",
        f"Runs; samples & {SQ['N']['runs']:,}; {SQ['N']['samples']} & {SQ['L']['runs']:,}; {SQ['L']['samples']} \\\\",
        f"Instrument of the paired runs & Genome Analyzer IIx ({SQ['N']['models_paired']['Illumina Genome Analyzer IIx']} of {SQ['N']['paired_runs']}) & HiSeq 2000 ({SQ['L']['models_paired']['Illumina HiSeq 2000']} of {SQ['L']['paired_runs']}) \\\\",
        f"Instrument of all the runs & Genome Analyzer IIx {SQ['N']['models_all']['Illumina Genome Analyzer IIx']:,}, Genome Analyzer II {SQ['N']['models_all']['Illumina Genome Analyzer II']}, HiSeq 2000 {SQ['N']['models_all']['Illumina HiSeq 2000']} & HiSeq 2000 {SQ['L']['models_all']['Illumina HiSeq 2000']} \\\\",
        f"Bases per mate of the paired runs, median (5th to 95th percentile) & {SQ['N']['mate_median']:.0f} ({SQ['N']['mate_p5_p95'][0]:.0f}--{SQ['N']['mate_p5_p95'][1]:.0f}) & {SQ['L']['mate_median']:.0f} ({SQ['L']['mate_p5_p95'][0]:.0f}--{SQ['L']['mate_p5_p95'][1]:.0f}) \\\\",
        f"Gbp per sample, all runs, median (interquartile range) & {SQ['N']['gbp_median']:.1f} ({SQ['N']['gbp_iqr'][0]:.1f}--{SQ['N']['gbp_iqr'][1]:.1f}) & {SQ['L']['gbp_median']:.1f} ({SQ['L']['gbp_iqr'][0]:.1f}--{SQ['L']['gbp_iqr'][1]:.1f}) \\\\",
        f"Queried runs of the donors & {QR['N']['n']} paired Genome Analyzer IIx runs, median {QR['N']['gbp_median']:.1f} Gbp & {QR['L']['n']} paired HiSeq 2000 runs, median {QR['L']['gbp_median']:.1f} Gbp \\\\",
        f"Library names & {SQ['N']['libraries']} & {SQ['L']['libraries']} ({lk['libraries_shared']} in common) \\\\",
        "Record created in the ENA & December 2012 & January and May 2014 \\\\",
        "Centre label & BGI & Beijing Genome Institute; BGI \\\\",
        "DNA extraction & Not recorded (none in the curated metadata) & Not recorded (MoBio kit in the curated metadata) \\\\",
    ])
    tab1516 = rf"""
\begin{{table}}[!ht]\centering\footnotesize\setlength{{\tabcolsep}}{{4pt}}
\caption{{\textbf{{Resemblance between the samples of Nielsen and Li with the same identifier.}} A pair is a Nielsen sample and a Li sample; the pairs of the same identifier are those of the {nshared} subject identifiers present in both cohorts. First means rank 1 among the 291,881 runs of the HGM set in the two lines on the sketch index, where the runs of the other samples of the same identifier are not counted as runs of other people, and among the 396 Nielsen samples in the three lines on species profiles. AUC, area under the curve of the dissimilarity for telling the pairs of the same identifier from the other pairs of Nielsen and Li.}}
\begin{{tabular}}{{>{{\raggedright\arraybackslash}}p{{4.3cm}}>{{\raggedright\arraybackslash}}p{{6.4cm}}>{{\raggedright\arraybackslash}}p{{4.4cm}}}}
\toprule
Line of evidence & Result & Limit \\
\midrule
{T15}
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[!ht]\centering\small
\caption{{\textbf{{Sequencing of the samples of Nielsen and Li.}} Run-level records of the ENA for the two projects, all stool samples of the two cohorts. \textit{{Bases per mate}}, the mean read length of a paired run, for the reads as deposited after the processing of the submitters.}}
\begin{{tabular}}{{>{{\raggedright\arraybackslash}}p{{6.0cm}}>{{\raggedright\arraybackslash}}p{{4.5cm}}>{{\raggedright\arraybackslash}}p{{4.5cm}}}}
\toprule
 & Nielsen 2014 & Li 2014 \\
\midrule
{T16}
\bottomrule
\end{{tabular}}
\end{{table}}
"""

# ---- Supplementary Note 1 and Supplementary Table 9: the reads of seven Costea donors deposited under another subject code
note6 = tab17 = ""
if cs7:
    D7 = cs7["donors"]; C7 = cs7["cohort"]; SE7 = cs7["sensitivity"]; OS7 = cs7["open_set"]; NP7 = cs7["no_partner"]; PR7 = cs7["projects"]
    by_ps7 = {r["subject"]: r for r in rows + rows_ex}; assert len(by_ps7) == len(rows) + len(rows_ex), "a subject identifier names two donors"                               # the seven are among the donors removed by the two rules (Supplementary Note 1)
    for d in D7:                                                                # the local aggregate must agree with the table of the donors
        r = by_ps7[d["subject"]]
        assert r["study"] == "CosteaPI_2017" and r["rank_any_narrow"] == d["rank_stock"] and r["rank_any_narrow"] != 1 and abs((r["partner_any"] or 0) - d["partner"]) < 1e-9 and abs(r["max_imp"] - d["closest_other"]) < 1e-9 and r["days"] == d["days"]
    cand7 = rows + rows_ex; cos7 = [r for r in cand7 if r["study"] == "CosteaPI_2017"]
    assert SE7["all_as_counted"] == [sum(r["rank_any_narrow"] == 1 for r in cand7), len(cand7)] and SE7["costea_as_counted"][1] == len(cos7) and SE7["costea_as_counted"][0] == sum(r["rank_any_narrow"] == 1 for r in cos7) and C7["failures_german"] == len(D7) == C7["german"]
    assert all(d["same_year"] and d["same_site"] and d["third_deposit_same_day"] and d["age_gap_years"] <= 1 and d["curated_days_in_third_deposit"][0] == d["curated_days_in_third_deposit"][1] and d["median_length_matched"] == 94 and d["reads_contained_first_tile"] == d["reads_found_first_tile"] for d in D7)
    with_p = [d for d in D7 if d["partner"] > 0]; relab = [d for d in with_p if d["of_a_relabelled_sample"]]
    assert all(d["stock_runs_above_partner"] == d["of_the_subject_code"] + d["of_a_new_deposit"] and d["rank_without_the_code"] == 1 for d in with_p) and len(with_p) == SE7["regained"] == 6 and len(relab) == 1 and D7 and len(D7) - len(with_p) == 1
    assert [d["date_gap_days"] > 0 for d in D7].count(True) == 1 and [d["sex"] for d in D7].count("same") == 6 and [d["sex"] for d in D7].count("not recorded") == 1
    found = [d["reads_found_first_tile"] / d["reads_queried_first_tile"] for d in D7]; ratio = [d["reads_matched_first_tile"] / d["reads_queried_first_tile"] for d in D7]
    lenq = [d["median_length_queried"] for d in D7]; nq = [d["reads_queried_first_tile"] for d in D7]; series = [d["curated_days_in_third_deposit"][1] for d in D7]; days7 = sorted({d["sampling_day"] for d in D7})
    ownr = [x for d in D7 if d["other_runs_of_queried_sample"] for x in d["other_runs_of_queried_sample"]]; matr = [x for d in D7 for x in d["other_runs_of_matched_sample"]]
    TF = [d["trim_filter"] for d in D7]                                           # costea_trim_filter.py: trimmed and filtered? (first tile of read 1 of each run)
    assert all(t["relation"]["identical"] == t["relation"]["matched_inside_queried"] == t["relation"]["neither"] == t["relation"]["ambiguous_position"] == 0 and t["relation"]["queried_inside_matched"] == t["common"] for t in TF)
    assert all(t["cut_5prime"]["share_of_common_reads"] == 1 for t in TF) and min(t["length_matched"]["min"] for t in TF) == max(t["length_matched"]["max"] for t in TF)
    rg = lambda xs, nd=0: f"{min(xs):.{nd}f}" if f"{min(xs):.{nd}f}" == f"{max(xs):.{nd}f}" else f"{min(xs):.{nd}f} to {max(xs):.{nd}f}"        # a range, or the single value when it is the same for all
    n_common = sum(t["common"] for t in TF); same_q = sum(t["same_quality_on_kept_bases"]["same"] for t in TF); checked_q = sum(t["same_quality_on_kept_bases"]["checked"] for t in TF)
    len_m = TF[0]["length_matched"]["min"]; len_q = (min(t["length_queried"]["min"] for t in TF), max(t["length_queried"]["max"] for t in TF))
    above7 = [d["stock_runs_above_partner"] for d in with_p]; ceil2 = lambda x: f"{math.ceil(x * 100) / 100:.2f}"; pr = lambda p: PR7[p]["first_public"]
    assert all(PR7[p]["centre"] == "EMBL" for p in PR7) and pr("PRJEB41102").startswith("2020-11") and pr("PRJEB50977").startswith("2023-10") and days7 == [30, 918] and NP7["other_code_days"] == [0, 7] and len(NP7["other_code_scores"]) == 2
    NUMW7 = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six", 7: "Seven"}
    p7 = {d["subject"] for d in D7}
    n_redep7 = sum(1 for e in rows_ex if e["rule"] == "re-deposit" and e["subject"] in p7); n_conta7 = sum(1 for e in rows_ex if e["rule"] == "contamination" and e["subject"] in p7)
    assert n_redep7 + n_conta7 == len(D7) == 7 and n_conta7 == 1, "Supplementary Note 1 says that the rule of re-deposit removes six of the seven and the rule of contamination the seventh"
    note6 = rf"""
\subsubsection*{{Reads of seven Costea donors deposited under another subject code}}

For {len(D7)} of the {SE7['all_as_counted'][1]} donors of the cohorts kept the nearest run of other people has a similarity of {min(d['closest_other'] for d in D7):.2f} to {max(d['closest_other'] for d in D7):.2f}. Among the {SE7['all_as_counted'][1] - len(D7)} other donors its highest similarity is {cs7['largest_closest_other_of_the_490']:.2f}. The seven similarities are also higher than the highest similarity of the other sample of any donor, which is {cs7['largest_partner_of_all']:.2f}. The seven are the {C7['german']} donors of Costea with a German subject identifier. None of them is recovered, and they make up {C7['failures_german']} of the {C7['failures']} donors of the cohort that are not recovered. For each of the seven we looked up the nearest run of other people, which we call the matched run, in the public records of the ENA and of BioSamples. We then compared its reads with those of the queried run. Supplementary Table 9 gives the results donor by donor.

For all seven, the matched run belongs to one project of the ENA, PRJEB41102, first public in November 2020, on strain persistence in the adults and infants of families. The queried runs belong to PRJEB17632, the project of Costea, and the same centre (EMBL) deposited both. PRJEB41102 records the matched sample under a subject code of its own. The curated metadata contain no link between the two projects, so the pipeline could have counted the matched run as belonging to another person.

We compared read 1 of the two runs on the first tile of the flow cell, from the public FASTQ files ({min(nq):,} to {max(nq):,} queried reads per pair). In all seven pairs the two runs share the instrument, lane and tile. In each pair, {100 * min(found):.1f}\% or more of the queried reads have a read at the same position on the tile (same x and y coordinates) in the matched run ({min(t['only_in_q'] for t in TF):,} to {max(t['only_in_q'] for t in TF):,} queried reads per pair have none). Over the seven pairs this makes {n_common:,} reads found in both runs. For every one of them the queried read is a strictly shorter part of the matched read (none is identical, none is longer). For {same_q:,} of the {checked_q:,} this part has the same quality string in the two runs. The matched reads all have {len_m} bases, the queried reads {len_q[0]} to {len_q[1]}. In every read found in both runs, the queried run lacks the first bases of the matched read (a median of {rg([t['cut_5prime']['median_bases'] for t in TF])} bases depending on the pair). In {rg([100 * t['cut_3prime']['share_of_common_reads'] for t in TF])}\% of these reads the queried run also lacks a tail (a median of {rg([t['cut_3prime']['median_bases'] for t in TF])} bases depending on the pair). On the tile the matched run holds {min(ratio):.1f} to {max(ratio):.1f} times as many reads as the queried run. The queried run is thus a trimmed and filtered version of the reads of the matched run. Each queried run also has an exact re-deposit in a third project of the same centre (PRJEB24041), a run with the same number of bases that has a similarity of 0.99 or more, and the pipeline removed it (Methods). These comparisons use the scripts \texttt{{costea\_reads\_compare.py}} and \texttt{{costea\_trim\_filter.py}} of the reproduction package.

The subject code of the matched run also names samples in a fourth project of the same centre (PRJEB50977, first public in October 2023), whose sample names give the day of the sample. For all seven, this project holds a sample of the code on the day given in the Costea name of the query sample (day {days7[0]} or {days7[1]}). It also holds a sample of the code on all the other sampling days of the donor's curated series ({min(series)} to {max(series)} sampling days per donor in all). The PRJEB50977 sample of the queried day has the collection date of the query sample, or a date {max(d['date_gap_days'] for d in D7)} days away for 1 donor. The record of the matched sample gives the collection year of the query sample, the same sex as the curated metadata (not recorded for 1 donor) and an age within one year of the curated age plus the time from day 0 of the series to the query sample.

Six of the seven donors have their other sample among the 5,000 nearest runs. For each of them we looked at every run of the HGM set that the pipeline counts as another person's and that is more similar to the query than the other sample ({min(above7)} to {max(above7)} runs per donor). Each of these runs either carries the subject code of the matched run or is a run of PRJEB24041, with a similarity below 0.99, whose sample name gives the same person as the Costea name of the query sample. For one of the six, {relab[0]['of_a_relabelled_sample']} of these runs belong to a sample of PRJEB41102 labelled with another subject code. The record of this sample (sex, age, family role) is that of the subject code of the matched run, and its sample number is consecutive with those of that code. We therefore count these {relab[0]['of_a_relabelled_sample']} runs under the subject code of the matched run. With all the runs above the other sample counted as the donor's own, the other sample would rank first for each of the six donors. The number of donors recovered would then rise from {SE7['costea_as_counted'][0]} of {SE7['costea_as_counted'][1]} ({pct(*SE7['costea_as_counted'])}\%) to {SE7['costea_with_regained'][0]} of {SE7['costea_with_regained'][1]} ({pct(*SE7['costea_with_regained'])}\%) in Costea and, over all the donors of the cohorts kept, from {SE7['all_as_counted'][0]} of {SE7['all_as_counted'][1]} ({pct1(*SE7['all_as_counted'])}\%) to {SE7['all_with_regained'][0]} of {SE7['all_with_regained'][1]} ({pct1(*SE7['all_with_regained'])}\%). The seventh donor would still not be recovered. In the open-set test of the {SE7['all_as_counted'][1]} donors of the cohorts kept, the search for an absent person was accepted for all seven.

The seven are removed from the analysis (above), {NUMW7[n_redep7].lower()} by the rule of re-deposit and the one whose query sample takes part in a contamination event by the rule of contamination. We cannot list the runs of a person that sit in other projects under other identifiers, so other donors may have such runs at lower similarities, which the pipeline counts as other people. This could possibly explain why the Costea cohort has a lower recovery rate than the other cohorts (see below).
"""
    T7 = []
    for d in D7:
        other = "Not in 5,000" if d["partner"] == 0 else f"{d['partner']:.2f}, rank {d['rank_stock']}"
        above = "--" if d["partner"] == 0 else f"{d['stock_runs_above_partner']} ({d['of_the_subject_code'] + d['of_a_new_deposit']})"
        without = "--" if d["partner"] == 0 else str(d["rank_without_the_code"])
        sex = "Same" if d["sex"] == "same" else "n.r."
        assert "queried_run" not in d and "matched_run" not in d, "the accessions of the seven pairs stay out of the supplement"
        T7.append(f"{d['subject'].replace('_', chr(92) + '_')} & {d['days']:.0f} & {other} & {d['closest_other']:.2f} & {d['reads_queried_first_tile']:,} & {100 * d['reads_found_first_tile'] / d['reads_queried_first_tile']:.1f} & {d['sampling_day']} & {d['date_gap_days']} & {sex} & {above} & {without} \\\\")
    T7 = "\n".join(T7)
    tab17 = rf"""
\begin{{table}}[!ht]\centering\footnotesize\setlength{{\tabcolsep}}{{2.5pt}}
\caption{{\textbf{{Queried and matched runs of the seven Costea donors.}} \textit{{Donor}}, subject identifier of curatedMetagenomicData; the reproduction package gives the runs of each donor. \textit{{Days}}, between the first and the last sample of the donor. \textit{{Other sample}}, its similarity to the query and its rank among the 291,881 runs of the HGM set, before the removal of the seven donors (Supplementary Note 1). \textit{{Match}}, similarity of the nearest run of other people, a run of PRJEB41102. \textit{{Reads}}, read 1 of the queried run on the first tile of the flow cell. \textit{{Found}}, share of these reads (\%) with a read at the same position (instrument, lane, tile, x, y) in the matched run. \textit{{Day}}, day of the query sample in the series of its donor, as given by its Costea name; PRJEB50977 holds a sample whose name gives that day and the subject code of the matched run. \textit{{Gap}}, days between the collection dates of the query sample and of that sample of PRJEB50977. \textit{{Day}} is counted from day 0 of the series, so it exceeds \textit{{Days}} for the one donor whose first curated sample is later. \textit{{Sex}}, record of the matched sample against the curated metadata (n.r., not recorded). \textit{{Above}}, number of runs of the HGM set that the pipeline counts as another person's and that are more similar than the donor's other sample and, in parentheses, how many of them carry the subject code of the matched run or are runs of PRJEB24041 whose sample name gives the same person as the Costea name of the query sample (the runs of one sample labelled with another subject code are counted under the subject code of the matched run; Supplementary Note 1). \textit{{Rank}}, rank of the other sample without those runs. --, the other sample lies beyond the 5,000 nearest runs.}}
\begin{{tabular}}{{lrlrrrrrlrr}}
\toprule
Donor & Days & Other sample & Match & Reads & Found & Day & Gap & Sex & Above & Rank \\
\midrule
{T7}
\bottomrule
\end{{tabular}}
\end{{table}}
"""

# ---- Supplementary Note 4: the thresholds of re-identification risk (the 9% line of Fig. 1B); no number of the analysis is used here
note7 = rf"""
\section*{{Supplementary Note 4. Thresholds of re-identification risk}}

The main text compares the recovery rate with the thresholds of re-identification risk advised for the release of clinical data. This Note explains how this risk is measured, how our recovery rate relates to this measure, and where the values of the thresholds come from.

The thresholds come from the methods used to anonymise clinical data before their release, in the Institute of Medicine report\cite{{iom2015}} (appendix B by El Emam and Malin) and the guidance of the European Medicines Agency\cite{{ema}}. These methods give each record of a data set a probability of being re-identified. Records that share the same values of the indirectly identifying variables cannot be told apart and form a group. The re-identification probability of a record is 1 divided by the number of individuals in its group. For example, a record in a group of 11 has a probability of 1 in 11, or approximately 0.09, while a record alone in its group has a probability of 1. For the whole data set, the maximum risk is generally the highest probability among its records but the average of these probabilities (the average risk) is also considered. Both texts mentioned above retain the maximum risk for a public release. According to the report of the Institute of Medicine, the reason is that no contract binds the recipient of public data and that an adversary can target the record that is easiest to re-identify.

In our analysis the microbiome plays the role of the indirectly identifying variables. A recovered donor is analogous to a record alone in its group, since no record from another individual in the set is more similar to the query than the donor's other sample. This sample can thus be unambiguously singled out. The recovery links the two samples of this donor with a probability of 1, and a donor that is not recovered counts for 0. The recovery rate is the mean of these values over the donors and is therefore an average risk. This holds for the rate measured with the full HGM set, {n_stock} of {N} donors or {100 * n_stock / N:.1f}\%, and for the rates extrapolated to the world population (Supplementary Note 3). As mentioned, the maximum re-identification risk is the highest individual probability of re-identification, so recovering even a single donor is theoretically sufficient to bring it to 1. Our analysis using the full HGM set yields a maximum risk of 1 and an average risk of {n_stock / N:.2f}, both exceeding the 0.09 threshold (see below). Although thresholds for public data release are generally based on maximum risk, in this article, we report the average across donors because a maximum risk of 1 could arise from a single donor, whereas the average demonstrates that this risk is a widespread property of adult gut microbiomes.

A threshold on the maximum risk is the inverse of the smallest group that a data custodian accepts, called the minimum cell size. Data custodians, in particular the government agencies that report statistics, have long applied such minimum sizes to what they release, and the appendix of the report of the Institute of Medicine reviews them. The most common minimum is 5, a probability of 1 in 5 or 0.2, but it was not recommended. Some public health offices use 3 (0.33). Larger minima are in use, such as 11 in the United States (0.09) and 20 in Canada (0.05), and the largest that this review found is 25 (0.04). The review states that no threshold is agreed upon, that the lowest probabilities are used for public releases, and that a minimum of 20 is appropriate for them.

The minimum of 11 is the rule of the Centers for Medicare \& Medicaid Services, which allow no cell smaller than 11 to be displayed from their data\cite{{cms2020}}. The European Medicines Agency advises a threshold of 0.09 on the maximum risk for the clinical reports that it publishes. It bases this value on the recommendations made in the report of the Institute of Medicine and on the precedents for the public release of health data, and an applicant may justify another value. Health Canada encourages the same threshold, which it equates to a smallest group of 11 patients who share the same indirectly identifying values\cite{{healthcanada2019}}.

Thresholds used in practice range from 33\% to 4\%, and the traditional 20\% is used more often for trusted recipients than for public release\cite{{iom2015,elemam2011b,elemam2013}}. The fading area of Fig. 1B spans 33\% to 9\%.

One of the main principles of many privacy laws is that identification does not have to be proved and that it is anonymisation that has to be proved (principle of privacy by default, recognised in various data protection frameworks, including those of the European Union, the United Kingdom, Switzerland and Kenya). Recovering a donor establishes a link between two metagenomic samples from the same individual, even when they are deposited under different identifiers. This ability, known as \textit{{linkability}}, is one of the three risks that anonymisation must address\cite{{wp29}} to consider data as anonymous. The European Data Protection Board is working on an updated version of the criteria. These new guidelines\cite{{edpb2026}} also rule the absence of linkage between records as a criterion for assessing anonymity. This criterion is particularly important for publicly released data, which can be linked by anyone with the necessary means (in our case effortlessly thanks to the Faure \textit{{et al.}} web service). Our results demonstrate that such linkage is already achievable using publicly accessible tools, directly challenging the assumption that publicly released gut metagenomic sequences can be considered anonymous in the European Union and in countries sharing the same kind of rule.
"""

# ---- Supplementary Tables 1 and 2: how the cohorts were chosen (Supplementary Note 1)
tab1314 = ""
if sc:
    NAMES14 = dict(NAMES, HMP_2012="HMP 2012", HMP_2019_t2d="HMP 2019 T2D", IjazUZ_2017="Ijaz")
    cat_, lg, lb, ac, at = sc["catalogue"], sc["longitudinal"], sc["labels"], sc["adult_cohorts"], sc["adult_totals"]
    ex_b = sc["biosample_only_in_index"]
    n_left = sum(c["adult_dna"] for c in ac.values() if c["adult_dna"] and not c["analysed"]); k_left = sum(1 for c in ac.values() if c["adult_dna"] and not c["analysed"])
    n_idx = sum(1 for c in ac.values() if c["indexed"]); n_dna = sum(1 for c in ac.values() if c["adult_dna"]); n_ana = sum(1 for c in ac.values() if c["analysed"])
    T13 = "\n".join([
        f"curatedMetagenomicData & {cat_['studies']} studies & {cat_['samples']:,} samples \\\\",
        f"\\quad stool samples & {cat_['stool_studies']} studies & {cat_['stool']:,} samples \\\\",
        f"\\quad donors with two or more stool samples & {sc['two_stool']['studies']} studies & {sc['two_stool']['donors']:,} donors \\\\",
        "\\midrule",
        f"Longitudinal cohorts ({lg['min_donors']} or more such donors) & {lg['cohorts']} cohorts & {lg['donors']:,} donors \\\\",
        f"\\quad labelled infants & {lb['infants']['cohorts']} & {lb['infants']['donors']} \\\\",
        f"\\quad labelled mixed & {lb['mixed']['cohorts']} & {lb['mixed']['donors']} \\\\",
        f"\\quad not labelled (fewer than 70\\% of donors with an age) & {lb['not recorded']['cohorts']} & {lb['not recorded']['donors']} \\\\",
        f"\\quad labelled adults & {lb['adults']['cohorts']} & {lb['adults']['donors']} \\\\",
        "\\midrule",
        f"Adult cohorts, both samples in the index (before Dec. 2023) & {n_idx} cohorts & {at['indexed']} donors \\\\",
        f"Adults (18 years or more, or Mehta) with a DNA query & {n_dna} cohorts & {at['adult_dna']} donors \\\\",
        f"Without the {k_left} cohorts left out & {n_ana} cohorts & {sc['excluded']['candidates']} donors \\\\",
        f"Without the {at['contamination'] + at['re_deposit']} donors removed by two rules & {n_ana} cohorts & {at['analysed']} donors \\\\",
    ])
    zd = lambda v: str(v) if v else "--"
    rm14 = lambda c_: f"{c_['contamination']} + {c_['re_deposit']}" if c_.get("contamination") or c_.get("re_deposit") else "--"     # removed by the rule of contamination, then by that of re-deposit
    dm_ = sc["demographics"]                                                       # cohort_selection.py: sex and country of the donors of the analysis, from the curated metadata
    CNAME = dict(USA="United States", KAZ="Kazakhstan", DNK="Denmark", SWE="Sweden", CAN="Canada", LUX="Luxembourg")
    assert set(dm_["total"]["countries"]) <= set(CNAME), "a country without a name in Supplementary Table 2"
    assert dm_["total"]["donors"] == at["analysed"] and dm_["total"]["sex_not_recorded"] == 0 and dm_["total"]["women"] + dm_["total"]["men"] == at["analysed"], "the caption of Supplementary Table 2 says every donor of the analysis has a recorded sex"

    def demo14(c):
        if not ac[c]["analysed"]:
            assert c not in dm_["by_cohort"]
            return "-- & \\multicolumn{1}{c}{--}"                                   # no donor of the analysis; the dash of the country is centred, apart from that of the column before
        d_ = dm_["by_cohort"][c]
        assert d_["donors"] == ac[c]["analysed"]
        ctry = ", ".join(CNAME[k_] if len(d_["countries"]) == 1 else f"{CNAME[k_]} ({n_})" for k_, n_ in d_["countries"].items())
        return f"{d_['women']}/{d_['men']} & {ctry}"
    row14 = lambda c: f"{NAMES14[c]} & {ac[c]['donors']} & {ac[c]['indexed']} & {zd(ac[c]['under18'])} & {zd(ac[c]['no_age'])} & {zd(ac[c]['rna'])} & {ac[c]['adult_dna']} & {rm14(ac[c])} & {ac[c]['analysed']} & {demo14(c)} \\\\"
    kept14 = sorted([c for c in ac if ac[c]["analysed"]], key=lambda c: (-ac[c]["donors"], c))
    left14 = sorted([c for c in ac if ac[c]["adult_dna"] and not ac[c]["analysed"]], key=lambda c: (-ac[c]["donors"], c))
    none14 = sorted([c for c in ac if not ac[c]["adult_dna"]], key=lambda c: (-ac[c]["donors"], c))
    T14 = "\\multicolumn{11}{l}{\\textit{Cohorts kept}} \\\\\n" + "\n".join(row14(c) for c in kept14)
    T14 += "\n\\midrule\n\\multicolumn{11}{l}{\\textit{Cohorts left out (Supplementary Note 1)}} \\\\\n" + "\n".join(row14(c) for c in left14)
    T14 += "\n\\midrule\n\\multicolumn{11}{l}{\\textit{Cohorts without a donor to analyse}} \\\\\n" + "\n".join(row14(c) for c in none14)
    T14 += f"\n\\midrule\nAll {len(ac)} adult cohorts & {at['donors']} & {at['indexed']} & {at['under18']} & {at['no_age']} & {at['rna']} & {at['adult_dna']} & {at['contamination']} + {at['re_deposit']} & {at['analysed']} & {dm_['total']['women']}/{dm_['total']['men']} & {len(dm_['total']['countries'])} countries \\\\"
    tab1314 = rf"""
\begin{{table}}[!ht]\centering\small
\caption{{\textbf{{Selection of the cohorts and donors of the analysis from curatedMetagenomicData.}} Studies and samples are those of the metadata tables of curatedMetagenomicData as of 1 October 2026. A longitudinal cohort is a study with at least {lg['min_donors']} donors who have two or more stool samples. The pair of a donor is its first and its last stool sample, and only the last one is queried. The labels (infants, mixed, adults, not labelled) belong to the cohorts and not to their donors.}}
\begin{{tabular}}{{lrr}}
\toprule
Step & Studies or cohorts & Samples or donors \\
\midrule
{T13}
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[!ht]\centering\footnotesize\setlength{{\tabcolsep}}{{3pt}}
\caption{{\textbf{{Selection of the donors of the analysis in each of the {len(ac)} adult cohorts.}} \textit{{Donors}}, donors with two or more stool samples in the curated metadata. \textit{{In the index}}, each of the two samples of the donor has a run in the 4,791,208 runs of the index (Logan v1.1 freeze of the Sequence Read Archive, last updated in December 2023, with a new update in progress), among the runs that the curated metadata list for the sample (for HMP IBD, among the runs of its ENA project; Supplementary Note 1). \textit{{Under 18}}, the lowest recorded age of the donor is under 18. \textit{{No age}}, no age is recorded; Mehta records none, but its donors are not counted here (Supplementary Note 1). \textit{{RNA}}, the queried run is named as an RNA library. \textit{{Adults, DNA}} (donors in the index that are in none of these three columns). \textit{{Removed}}, donors of the cohorts kept removed by the rule of contamination and by that of re-deposit, in this order (Supplementary Note 1). \textit{{In the analysis}}, the donors of the {n_ana} cohorts kept that the two rules leave. \textit{{Women/men}}, the women and the men among the donors of the analysis, by the sex that the curated metadata record. \textit{{Country}}, the country that the curated metadata record for these donors.}}
\begin{{tabular}}{{lrrrrrrrrrl}}
\toprule
Cohort & Donors & \makecell[rb]{{In the\\index}} & \makecell[rb]{{Under\\18}} & \makecell[rb]{{No\\age}} & RNA & \makecell[rb]{{Adults,\\DNA}} & Removed & \makecell[rb]{{In the\\analysis}} & \makecell[rb]{{Women/\\men}} & Country \\
\midrule
{T14}
\bottomrule
\end{{tabular}}
\end{{table}}
"""

# ---- Supplementary Tables 3 and 4: the HGM set of human gut metagenomes, and what the numbers become with other sets of runs
fmt_i = lambda v: f"{v:,}"
lab_ = st_["neighbours_added_by_label"]; left_ = st_["neighbours_left_out"]; gl_ = st_["gut_labels_added_by_label"]
rows17 = lambda L, names: "\n".join(f"\\quad {names[k_]} & {fmt_i(L[k_])} \\\\" for k_ in names if k_ in L)
N17 = {"gut or feces metagenome": "gut or feces metagenome", "human metagenome": "human metagenome", "metagenome, human host declared": "metagenome, human host declared", "metagenome, no host declared": "metagenome, no host declared", "Homo sapiens": "\\textit{Homo sapiens}", "human gut or feces metagenome, library other than WGS": "human gut or feces metagenome, library other than WGS", "human gut or feces metagenome, WGS": "human gut or feces metagenome, WGS"}
assert sum(lab_.values()) == st_["neighbours_added"] and sum(gl_.values()) == st_["gut_labels_added"] and st_["labelled"] + st_["curated_added"] + st_["neighbours_added"] + st_["gut_labels_added"] == st_["stock"]
row18 = lambda name, key: f"{name} & {fmt_i(se_[key]['runs'])} & {se_[key]['first']} & {se_[key]['threshold']:.4f} & {100 * se_[key]['absent'] / N:.1f} & {se_[key]['other_at_or_above_first']} ({se_[key]['other_at_or_above_first_pct']:.1f}) & {se_[key]['line_8e9']:.1f} \\\\"
sc18 = se_["scale_of_the_pools"]
tabgut = rf"""
\begin{{table}}[!ht]\centering\small
\caption{{\textbf{{Composition of the HGM set, the runs of the index taken as human gut metagenomes.}} \textit{{Labelled}}, human gut or faecal shotgun metagenome in the Logan metadata. These runs include none of the indexed runs of Costea, Hansen and Raymond, which are deposited under \textit{{Homo sapiens}}, feces metagenome and gut metagenome respectively. \textit{{Curated}}, the other runs that curatedMetagenomicData lists for its stool samples or that belong to the BioSamples of the stool samples of the 7 cohorts, kept when the ENA gives them a metagenomic library of strategy WGS, WGA or OTHER. \textit{{Neighbours}}, the other runs that are among the 5,000 nearest runs of at least one of the {st_['queries']} stool queries (those of the {sc['adult_totals']['indexed']} donors in the index and of {sum(sc['queried_outside_the_donors'].values())} more; Supplementary Note 1) and that have a metagenomic library of strategy WGS, WGA or OTHER, one of the organism labels of the rows below, no declared non-human host and no other part of the body as declared source (ENA, November 2026). \textit{{Other gut or feces labels}}, the other runs of the index whose organism label at the ENA is gut metagenome, feces metagenome, human gut metagenome or human feces metagenome, with the same conditions on the library, the host and the source as the neighbours. These runs are in no list of neighbours, so they change no rank, only the size of the pools.}}
\begin{{tabular}}{{lr}}
\toprule
Runs & Number \\
\midrule
Labelled & {fmt_i(st_['labelled'])} \\
Curated & {fmt_i(st_['curated_added'])} \\
Neighbours & {fmt_i(st_['neighbours_added'])} \\
{rows17(lab_, N17)}
Other gut or feces labels & {fmt_i(st_['gut_labels_added'])} \\
{rows17(gl_, N17)}
\midrule
All & {fmt_i(st_['stock'])} \\
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[!ht]\centering\footnotesize\setlength{{\tabcolsep}}{{4pt}}
\caption{{\textbf{{Sensitivity of the results to the set of runs counted as other people.}} Each row counts only the runs of the set as other people. \textit{{Recovered}}, donors for whom no run of the set has higher similarity to the query than their other sample. \textit{{Threshold}}, the similarity at which as many searches for an absent person are accepted as donors present are rejected (equal error rate, Fig. 1H). \textit{{Error}}, the share of the {N} searches for an absent person that are accepted at the threshold, equal to the share of the {N} donors present that are rejected, in \%. \textit{{Accepted, recovered}} (donors whose other sample reaches the threshold and who are recovered, \% of {N}). \textit{{Straight line}}, the value at $8\times10^9$, in \%, of the straight line of Fig. 1B fitted again with the ranks among the runs of the set. The pools are drawn from $M$ = {fmt_i(se_['stock']['runs'])} runs, the size of the HGM set, whatever the size of the set, except in the last row ($M$ = 172,075, the labelled runs). Four projects each hold people of one of the 7 cohorts. The runs of PRJNA306874 and PRJNA389280 that the public table of the HMP consortium gives to a donor count as that donor's own (Supplementary Note 1), and the other runs of the four projects count as other people in the first row. The three rows that remove the four projects keep the runs of PRJEB41102 (Supplementary Note 1), and the rule of re-deposit of Supplementary Note 1 removes the donors whose nearest run of other people is one of these runs. PRJEB24041 deposits again the reads of part of Costea (338 of its 633 runs carry the name of a Costea sample), PRJEB28237 holds cultures of the stool of the donors of Raymond (same subject codes), PRJNA306874 (``iHMP ibdmdb'') has sample titles that carry the curated sample identifiers of HMP IBD, and PRJNA389280 is a second deposit of the metagenomes of HMP IBD, whose sample aliases the table of the consortium names.}}
\begin{{tabular}}{{lrrrrrr}}
\toprule
Runs & Number & Recovered & Threshold & Error & \makecell[rb]{{Accepted,\\recovered}} & Straight line \\
\midrule
{row18('The HGM set (the paper)', 'stock')}
{row18('Without PRJEB24041 and PRJEB28237', 'stock_without_two_projects')}
{row18('Also without PRJNA306874', 'stock_without_three_projects')}
{row18('Also without PRJNA389280', 'stock_without_four_projects')}
{row18('All the runs of the index', 'index')}
{row18('Runs labelled human gut or faecal', 'labelled')}
\bottomrule
\end{{tabular}}
\end{{table}}
"""
# the tables other than the five of Supplementary Note 2 in the order of their first citation in the Notes, so that their numbers follow one another in the document (place_floats then puts each next to that citation)
TABLE_ORDER = ["Selection of the cohorts and donors of the analysis", "Selection of the donors of the analysis in each", "Composition of the HGM set", "Sensitivity of the results",
               "Pairs of samples of different participants", "Sequencing of the samples of Nielsen and Li", "Resemblance between the samples of Nielsen and Li", "Events between two people, by cohort.",
               "Queried and matched runs", "Extrapolation of the recovery curve", "Predictions of the ", "Recovery rate and extrapolations with one cohort left out", "Extreme-value model of the similarity tail"]
TABLES_OF_NOTE_1 = TABLE_ORDER.index("Extrapolation of the recovery curve")       # the tables cited first in Supplementary Note 1 come before the five tables of Supplementary Note 2, the others after
_tabs = tab67 + tab8 + tab1314 + tab1516 + tab17 + tabgut
_blocks = re.findall(r"\\begin\{table\}.*?\\end\{table\}\n", _tabs, re.S)
assert not re.sub(r"\\begin\{table\}.*?\\end\{table\}\n", "", _tabs, flags=re.S).strip(), "only tables in these blocks"
_rank = lambda b: [i for i, k in enumerate(TABLE_ORDER) if re.search(r"\\caption\{\\textbf\{" + re.escape(k), b)]
assert all(len(_rank(b)) == 1 for b in _blocks), "each table has one place in TABLE_ORDER"
tables_before = "\n".join(b for b in sorted(_blocks, key=lambda b: _rank(b)[0]) if _rank(b)[0] < TABLES_OF_NOTE_1)
tables_after = "\n".join(b for b in sorted(_blocks, key=lambda b: _rank(b)[0]) if _rank(b)[0] >= TABLES_OF_NOTE_1)

# Supplementary Note 1, on the cohorts, takes at the two marks left in its text the three texts that document its steps: the links between Nielsen and Li under the heading
# of these two cohorts, then cross-sample contamination and the reads of the seven Costea donors under the heading of the two rules
assert note3.count("@@LINKS@@\n") == note3.count("@@RULES@@\n") == (1 if note3 else 0), "the Note on the cohorts holds one mark for each text that goes into it"
note_cohorts = note3.replace("@@LINKS@@\n", note5.strip("\n") + "\n").replace("@@RULES@@\n", note2.strip("\n") + "\n\n" + note6.strip("\n") + "\n")
assert "@@" not in note_cohorts

doc = rf"""\documentclass[11pt,a4paper]{{article}}
\usepackage[utf8]{{inputenc}}
\usepackage[T1]{{fontenc}}
\usepackage{{helvet}}
\renewcommand{{\familydefault}}{{\sfdefault}}
\usepackage[margin=2.2cm]{{geometry}}
\usepackage{{booktabs}}
\usepackage{{graphicx}}
\usepackage{{caption}}
\usepackage{{placeins}}
\usepackage{{makecell}}
\usepackage{{amsmath}}
\usepackage{{amssymb}}
\usepackage{{xurl}}
\usepackage[super,sort&compress,comma]{{natbib}}
\providecommand{{\nolinkurl}}[1]{{\url{{#1}}}}
\providecommand{{\urlprefix}}{{URL }}
\providecommand{{\doi}}[1]{{https://doi.org/#1}}
\providecommand{{\bibcommenthead}}{{}}
\captionsetup{{labelfont=bf,labelsep=space,justification=raggedright,singlelinecheck=false,font=small}}
\renewcommand{{\tablename}}{{Supplementary Table}}
\renewcommand{{\figurename}}{{Supplementary Figure}}
\setlength{{\parindent}}{{0pt}}\setlength{{\parskip}}{{5.5pt}}\emergencystretch=2em
\renewcommand{{\topfraction}}{{0.9}}\renewcommand{{\bottomfraction}}{{0.8}}\renewcommand{{\textfraction}}{{0.08}}\renewcommand{{\floatpagefraction}}{{0.75}}\setcounter{{topnumber}}{{4}}\setcounter{{bottomnumber}}{{3}}\setcounter{{totalnumber}}{{6}}
\makeatletter\setlength{{\@fptop}}{{0pt}}\setlength{{\@fpsep}}{{12pt}}\setlength{{\@fpbot}}{{0pt plus 1fil}}\makeatother
\begin{{document}}
{{\large\bfseries Supplementary Information}}

{{\bfseries {title}}}

{AUTHORS}
{note_cohorts}{note}{note4}{note7}
{figsel}{fig5}{fig4}{fig2}{fig3}{fig6}\clearpage
{tables_before}
\begin{{table}}[!ht]\centering\footnotesize\setlength{{\tabcolsep}}{{2.2pt}}
\caption{{\textbf{{Donors by rank of the other sample.}} \textit{{Rank 2 or more}}, the other sample is among the 5,000 nearest runs returned, but at least one sample of another person has higher similarity to the query. \textit{{Not in 5,000}}, it is not among the 5,000 nearest runs returned by the Faure \textit{{et al.}} service. Similarities are medians of the similarities returned by the Faure \textit{{et al.}} service, with the interquartile range (IQR) in parentheses; the nearest run of other people is the most similar run of the whole index that is not assigned to the donor. We assigned runs to subjects only in the 12 adult cohorts that have donors in the index (Supplementary Table 2). Number of donors whose nearest run of other people is assigned to another subject of the donor's cohort (\textit{{Same cohort}}), to a subject of another of these 12 cohorts (\textit{{Other cohort}}) or to no subject of these 12 cohorts (\textit{{Outside the 12}}).}}
\begin{{tabular}}{{lrrrrrrr}}
\toprule
 & & & \multicolumn{{2}}{{c}}{{Median similarity (IQR)}} & \multicolumn{{3}}{{c}}{{Nearest run of other people from}} \\
\cmidrule(lr){{4-5}}\cmidrule(lr){{6-8}}
 & Donors & \% of {N} & Other sample & Other people & Same cohort & Other cohort & Outside the 12 \\
\midrule
{T1}
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[!ht]\centering\footnotesize\setlength{{\tabcolsep}}{{4pt}}
\caption{{\textbf{{Cohorts, ordered by the share of donors recovered.}} Cohorts are named by their curatedMetagenomicData identifiers. \textit{{Recovered}}, donors whose other sample is at rank 1 in the HGM set, the measure of Fig. 1D (also at rank 1 among all the runs of the index). \textit{{Not in 5,000}}, donors whose other sample is not among the 5,000 nearest runs returned, counted as not recovered. \textit{{Interval}}, median number of days between the query and the first sample, with the interquartile range (IQR) in parentheses. \textit{{Query run}}, median number of bases of the queried run in the ENA, in Gbp, with the IQR in parentheses. \textit{{Runs}}, median number of runs that the curated metadata list for the query sample (-- for HMP\_2019\_ibdmdb, whose curated metadata list none, and for all cohorts together; Supplementary Note 1).}}
\begin{{tabular}}{{lrrrrrr}}
\toprule
Cohort & Donors & Recovered (\%) & Not in 5,000 & Interval (days) & Query run (Gbp) & Runs \\
\midrule
{T2}
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[!ht]\centering\small
\caption{{\textbf{{Interval between the query and the first sample.}} The query is the last sample of the donor. The first sample is the earliest one (for HMP IBD, the earliest with a run in the index; Supplementary Note 1). \textit{{Interval}}, number of days between the first and the last sample. Under \textit{{Most similar other sample}}, the other sample of a donor is the most similar of all its other samples, which is not always the first sample (as in Supplementary Table 11); \textit{{Recovered}} counts the donors for whom it is at rank 1 in the HGM set and \textit{{Not in 5,000}} those for whom it is not among the 5,000 nearest runs returned. Under \textit{{First sample alone}}, the first sample takes the place of the other sample, and \textit{{Recovered}} counts the donors for whom it is at rank 1 in the HGM set, the measure of Fig. 1E.}}
\begin{{tabular}}{{lrrrr}}
\toprule
 & & \multicolumn{{2}}{{c}}{{Most similar other sample}} & \multicolumn{{1}}{{c}}{{First sample alone}} \\
\cmidrule(lr){{3-4}}\cmidrule(lr){{5-5}}
Interval (days) & Donors & Recovered (\%) & Not in 5,000 & Recovered (\%) \\
\midrule
{T3}
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[!ht]\centering\small
\caption{{\textbf{{Bases of the queried run.}} \textit{{Gbp}}, from the ENA, before any filtering. \textit{{Recovered}}, as in Supplementary Table 11 (rank 1), the measure of Fig. 1F, where the two classes under 0.5 Gbp are drawn as one.}}
\begin{{tabular}}{{lrr}}
\toprule
Gbp & Donors & Recovered (\%) \\
\midrule
{T4}
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[!ht]\centering\small
\caption{{\textbf{{Bases of the queried run within each cohort.}} Cohorts are named by their curatedMetagenomicData identifiers. Each cohort is split at the median number of bases of its queried runs, into the queries at or below this median and the queries above it. RaymondF\_2016 is not split, because all its queried runs hold the same number of bases. \textit{{At or below}}, \textit{{Above}}, donors recovered over donors of the group, with the share in percent.}}
\begin{{tabular}}{{lrrr}}
\toprule
Cohort & Median (Gbp) & At or below (\%)& Above (\%) \\
\midrule
{T5}
\bottomrule
\end{{tabular}}
\end{{table}}

{tables_after}\bibliographystyle{{sn-nature}}
\bibliography{{references}}
\end{{document}}
"""
CITED = re.compile(r"Supplementary (Tables?|Figs?\.|Figures?) ((?:\d+[A-Z]?(?:--\d+[A-Z]?)?)(?:(?:, and |, | and | to )\d+[A-Z]?)*)")


def cited_numbers(s):
    """The numbers of a citation such as '4 and 5', '8, 9 and 11', '1 to 4' or '3D'."""
    out = []
    for part in re.split(r", and |, | and ", s):
        n = [int(x) for x in re.findall(r"\d+", part)]
        out += list(range(n[0], n[1] + 1)) if len(n) == 2 and (" to " in part or "--" in part) else n[:1]
    return out




def place_floats(doc):
    """Puts each figure and each table next to the paragraph of the Notes that cites it first, and keeps its number.
    The floats are written after the Notes, in the order of their numbers. Each one is taken out and put back after the paragraph (one line of the file) that holds its first citation in the Notes; a float
    that the Notes do not cite follows the float of the number before it. The floats are then no longer in the order of their numbers, so a \setcounter before each one keeps its number, and a \FloatBarrier
    (placeins) before each Note keeps the floats of a Note inside that Note."""
    pat = re.compile(r"\\begin\{(figure|table)\}.*?\\end\{\1\}\n", re.S)
    blocks, count = [], collections.Counter()
    for m in pat.finditer(doc):
        count[m.group(1)] += 1
        blocks.append((m.group(1), count[m.group(1)], m.group(0)))
    lines = pat.sub("", doc).split("\n")
    end = next(i for i, l in enumerate(lines) if "\\bibliographystyle" in l)
    first = {}
    for i, l in enumerate(lines[:end]):
        for m in CITED.finditer(l):
            kind = "table" if m.group(1).startswith("Table") else "figure"
            for n in cited_numbers(m.group(2)):
                first.setdefault((kind, n), (i, m.start()))
    after = collections.defaultdict(list)
    for kind, n, block in blocks:
        where = first.get((kind, n))
        k = n
        while where is None and k > 1:                                           # not cited in the Notes: it follows the float of the number before it
            k -= 1
            where = first.get((kind, k))
        assert where is not None, f"Supplementary {kind} {n} is cited nowhere in the Notes and has no float before it"
        first.setdefault((kind, n), (where[0], where[1] + 0.5))
        block = block.replace("\\begin{figure}[t]", "\\begin{figure}[!htbp]").replace("\\begin{table}[!ht]", "\\begin{table}[!htbp]")
        after[where[0]].append((kind, n, first[(kind, n)][1], f"\\setcounter{{{kind}}}{{{n - 1}}}\n{block.rstrip()}"))
    out = []
    for i, l in enumerate(lines):
        if i == end or (l.startswith("\\section*{Supplementary Note") and i > 0):
            out.append("\\FloatBarrier")
        out.append(l)
        here = after.get(i, [])
        lead = {k: min(pos for kk, _, pos, _ in here if kk == k) for k in {kk for kk, _, _, _ in here}}       # in one paragraph: the kind cited first comes first, then the floats of a kind in the order of their numbers
        for _, _, _, b in sorted(here, key=lambda f: (lead[f[0]], f[1])):
            out += ["", b]
    assert sum(len(v) for v in after.values()) == len(blocks)
    for kind in ("table", "figure"):                                             # the numbers follow one another in the document
        seq = [int(m.group(1)) + 1 for m in re.finditer(r"\\setcounter\{" + kind + r"\}\{(\d+)\}", "\n".join(out))]
        assert seq == sorted(seq), f"the Supplementary {kind}s are not in the order of their numbers: {seq}"
    return "\n".join(out)


def add_contents(doc):
    """Adds Contents, Supplementary Figures and Supplementary Tables lists, TOC entries for every heading, and short caption titles."""
    doc = doc.replace(
        "\\providecommand{\\bibcommenthead}{}\n",
        "\\providecommand{\\bibcommenthead}{}\n\\renewcommand{\\bibsection}{\\section*{\\refname}\\addcontentsline{toc}{section}{\\refname}}\n",
        1,
    )
    doc = doc.replace(
        "\\FloatBarrier\n\\section*{Supplementary Note 1. The cohorts}",
        "{\\small\\setcounter{tocdepth}{3}\n\\renewcommand{\\contentsname}{Contents}\\tableofcontents\n\\renewcommand{\\listfigurename}{Supplementary Figures}\\listoffigures\n\\renewcommand{\\listtablename}{Supplementary Tables}\\listoftables}\n\\clearpage\n\n\\FloatBarrier\n\\section*{Supplementary Note 1. The cohorts}",
        1,
    )
    doc = re.sub(
        r"^\\((?:sub){0,2}section)\*\{([^\n]+)\}$",
        lambda m: f"\\{m.group(1)}*{{{m.group(2)}}}\\addcontentsline{{toc}}{{{m.group(1)}}}{{{m.group(2)}}}",
        doc,
        flags=re.M,
    )
    doc = re.sub(
        r"\\caption\{\\textbf\{([^}]+?)\.\}",
        lambda m: f"\\caption[{m.group(1)}]{{{chr(92)}textbf{{{m.group(1)}.}}",
        doc,
    )
    return doc


# The text types as literals a few numbers that depend on the data (the size of the HGM set, the adult cohorts with donors in the index, the cohorts kept). They are checked here, so that a change of the data cannot leave them stale.
assert st_ is None or st_["stock"] == 291881, "the text types the size of the HGM set as 291,881"
assert sc is None or sum(1 for c in sc["adult_cohorts"].values() if c["indexed"]) == 12, "the text types 12 adult cohorts with donors in the index"
assert len(by) == 7, "the text types 7 cohorts kept"
open(os.path.join(OVL, "supplementary.tex"), "w", encoding="utf-8").write(add_contents(place_floats(doc)))
print("wrote supplementary.tex;", len(ok), "recovered,", len(out), "outranked,", len(ab), "not in the 5,000;", len(small), "small queries")
