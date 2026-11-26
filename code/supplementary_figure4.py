#!/usr/bin/env python3
"""Supplementary Figure 4: CroCoDeEL events between two people in which a sample studied takes part, cohort by cohort, for the 7 cohorts kept and the 4 cohorts left out.

    python3 -I supplementary_figure4.py ROWS.json OUTPUT_PREFIX [NIELSEN_LI.conta.tsv [CROCODEEL_DIR]]     # writes OUTPUT_PREFIX.png (600 dpi), .pdf and _stats.json

ROWS.json is data/donors.json: the donors kept, with the cohorts kept and the donors of the statistics of Supplementary Note 1.
CroCoDeEL results (<cohort>.conta.tsv, version 1.3.1, default settings) and the sample metadata (<cohort>.metadata.tsv: sample, subject, days from the first collection) are read from CROCODEEL_DIR (default data/crocodeel/);
sample and subject names are those of curatedMetagenomicData (make_inputs.py).
An event joins two samples of one cohort, a source and a target, with a rate, the share of the target that comes from the source (CroCoDeEL). It is "between two people" when the two samples belong to two different subjects of the
curated metadata, "same donor" when they belong to one subject. The samples studied are the query sample, the first sample and the most similar other sample of each donor, for the donors
of the 7 cohorts kept (ROWS.json and donors_excluded.json, read from the folder of ROWS.json: the 497 donors before the two rules) and of the 4 cohorts left out (donors_left_out.json, read from the folder of ROWS.json); an event counts when one of them is its source
or its target, and CroCoDeEL itself was run on all the samples of each cohort. A, share of the samples studied of each cohort that take part in at least one event between two people, hatched by the role of the sample in its events (source, target or both); B, rate of each of these
events. C (with the third argument): Nielsen and Li carry the same subject identifiers for 44 patients, so CroCoDeEL was also run on their two tables merged; events between a Nielsen and a Li sample in which a sample studied takes part,
by subject identifier.
The _stats.json file holds under "studied" the statistics limited to the samples studied (Supplementary Note 1, Supplementary Fig. 4 and Supplementary Table 8), and under "by_cohort" those of all the samples
of the 11 cohorts, which Supplementary Note 1 quotes, with the status of the query sample and the donors whose samples take part in an event between two people (Supplementary Note 1).
Palette and type follow Figure 1; grey = cohort kept, red = cohort left out.
"""
import csv, json, math, os, statistics as S, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.patches
from matplotlib.transforms import blended_transform_factory

ROWS, PREFIX = sys.argv[1], sys.argv[2]
NL_CONTA = sys.argv[3] if len(sys.argv) > 3 else None                  # CroCoDeEL run on the species tables of Nielsen and Li merged (crocodeel_nielsen_li.py): panel C
CRO = sys.argv[4] if len(sys.argv) > 4 else "data/crocodeel"
INK, INK2, INK3, RULE, GUIDE = "#1B2227", "#3F4A52", "#6B777F", "#C9D0D4", "#E3E7EA"
C_KEPT, C_OUT = "#6B777F", "#C8553D"
import os
from matplotlib import font_manager
FONTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fonts")                  # Liberation Sans shipped with the package (licences in fonts/), so that the drawing does not depend on the machine
font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if f.name != "Liberation Sans"]
for fn in sorted(os.listdir(FONTS)):
    if fn.lower().endswith(".ttf"):
        font_manager.fontManager.addfont(os.path.join(FONTS, fn))
plt.rcParams.update({"font.family": ["Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42})
NAMES = {"MehtaRS_2018": "Mehta", "CosteaPI_2017": "Costea", "HMP_2019_ibdmdb": "HMP IBD", "HansenLBS_2018": "Hansen", "Bengtsson-PalmeJ_2015": "Bengtsson-Palme", "RaymondF_2016": "Raymond",
         "Heitz-BuschartA_2016": "Heintz-Buschart", "PehrssonE_2016": "Pehrsson", "NielsenHB_2014": "Nielsen", "LiJ_2014": "Li", "VincentC_2016": "Vincent"}
LEFT_OUT = ["PehrssonE_2016", "NielsenHB_2014", "LiJ_2014", "VincentC_2016"]

rows = [r for r in json.load(open(ROWS)) if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]
STUDIES = sorted({r["study"] for r in rows})                      # the cohorts kept
assert set(STUDIES) | set(LEFT_OUT) == set(NAMES) and not set(STUDIES) & set(LEFT_OUT), set(STUDIES) ^ (set(NAMES) - set(LEFT_OUT))
ROWS_LO = os.path.join(os.path.dirname(os.path.abspath(ROWS)), "donors_left_out.json")     # the 183 donors of the 4 cohorts left out, in the format of ROWS.json, for the samples studied of these cohorts
rows_lo = [r for r in json.load(open(ROWS_LO)) if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]
assert {r["study"] for r in rows_lo} == set(LEFT_OUT), "donors_left_out.json must hold the donors of the 4 cohorts left out"
ROWS_EX = os.path.join(os.path.dirname(os.path.abspath(ROWS)), "donors_excluded.json")      # the donors of the cohorts kept, before the two rules, that the rules of contamination and re-deposit remove (Supplementary Note 1)
rows_ex = json.load(open(ROWS_EX)) if os.path.exists(ROWS_EX) else []
assert all(r["study"] in STUDIES for r in rows_ex)
cand = rows + rows_ex                                                             # the donors of the cohorts kept, before the two rules, (497): the samples studied and the evidence of the rule of contamination
ev, subj, nsamples = {}, {}, {}
for st in list(NAMES):
    md = {m["sample_id"]: m for m in csv.DictReader(open(f"{CRO}/{st}.metadata.tsv"), delimiter="\t")}
    nsamples[st] = len(md)
    subj[st] = {k: m["subject_id"] for k, m in md.items()}
    ev[st] = [(d["source"], d["target"], float(d["rate"])) for d in csv.DictReader([l for l in open(f"{CRO}/{st}.conta.tsv") if not l.startswith("#")], delimiter="\t")]


def known(st, k):
    v = subj[st].get(k)
    return v if v not in (None, "", "NA") else None


def kind(st, s, t):
    a, b = known(st, s), known(st, t)
    return "unknown" if a is None or b is None else ("same" if a == b else "other")


def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


# ------------------------------------------------------------------ the numbers of Supplementary Note 1 (the 7 cohorts kept)
def status(r):                                              # the query sample, against every event of its cohort
    st, Q = r["study"], r["v2s"]; same = other = False
    for s, t, _ in ev[st]:
        if Q in (s, t):
            sx = subj[st].get(t if s == Q else s)
            if sx == r["subject"]:
                same = True
            elif sx:
                other = True
    return "other" if other else "same" if same else "none"


for r in rows:
    r["status"] = status(r)
far = [r for r in rows if r["far_possible"]]
CATS = ["none", "same", "other"]
stats_A = {}
for key in CATS:
    for lab, grp, ok in (("all", rows, lambda r: r["rank_any"] == 1), ("stock", rows, lambda r: r["rank_any_narrow"] == 1), ("30d", far, lambda r: r["rank_far_narrow"] == 1)):
        g = [r for r in grp if r["status"] == key]; stats_A[f"{key}/{lab}"] = [sum(ok(r) for r in g), len(g)]
# donors whose first sample, query sample or most similar other sample takes part in an event with another person: the most that contamination between two people could change
oth_samples = {st: {x for s, t, _ in ev[st] if kind(st, s, t) == "other" for x in (s, t)} for st in STUDIES}
exposed = [r for r in cand if any(x in oth_samples[r["study"]] for x in (r["v1s"], r["v2s"], r["best_any"]) if x)]        # among the 497 donors of the cohorts kept
exp_ids = {id(r) for r in exposed}; rest_ = [r for r in cand if id(r) not in exp_ids]
stats_exposed = dict(n=len(exposed), linked=sum(r["rank_any"] == 1 for r in exposed), rest_n=len(rest_), rest_linked=sum(r["rank_any"] == 1 for r in rest_),
                     linked_stock=sum(r["rank_any_narrow"] == 1 for r in exposed), rest_linked_stock=sum(r["rank_any_narrow"] == 1 for r in rest_))
# the same donors, sample by sample: the samples that take part in an event between two people, and which of them the analysis uses (the query sample, the first sample or the most similar other sample of a donor)
used_samples = {(r["study"], x) for r in cand for x in (r["v1s"], r["v2s"], r["best_any"]) if x and x in oth_samples[r["study"]]}
q_don = [r for r in exposed if r["v2s"] in oth_samples[r["study"]]]                       # the query sample is one of them
o_don = [r for r in exposed if r["v2s"] not in oth_samples[r["study"]]]                   # only the first sample or the most similar other sample is one of them
assert len(q_don) + len(o_don) == len(exposed) and len(q_don) == sum(status(r) == "other" for r in cand), "the donors with the query sample in an event with another person must be those of the status 'other'"


def rates_in(st, samples):                                                                  # the rates of the events between two people in which one of these samples takes part
    return [rt for s, t, rt in ev[st] if kind(st, s, t) == "other" and (s in samples or t in samples)]


q_rates = [rt for r in q_don for rt in rates_in(r["study"], {r["v2s"]})]
o_rates = [rt for r in o_don for rt in rates_in(r["study"], {x for x in (r["v1s"], r["best_any"]) if x in oth_samples[r["study"]]})]
stats_exposed.update(samples_in_other=sum(len(oth_samples[st]) for st in STUDIES), samples_used=len(used_samples),
                     query_donors=len(q_don), query_linked_stock=sum(r["rank_any_narrow"] == 1 for r in q_don), query_min_rate=min(q_rates), query_max_rate=max(q_rates),
                     other_donors=len(o_don), other_linked_stock=sum(r["rank_any_narrow"] == 1 for r in o_don), other_min_rate=min(o_rates), other_max_rate=max(o_rates))
ex_keys = {(r["study"], r["subject"]) for r in rows_ex}
out_ = lambda r: (r["study"], r["subject"]) in ex_keys
o_ex, o_kept = [r for r in o_don if out_(r)], [r for r in o_don if not out_(r)]
o_ex_rates = [rt for r in o_ex for rt in rates_in(r["study"], {x for x in (r["v1s"], r["best_any"]) if x in oth_samples[r["study"]]})]
o_kept_rates = [rt for r in o_kept for rt in rates_in(r["study"], {x for x in (r["v1s"], r["best_any"]) if x in oth_samples[r["study"]]})]
stats_exposed.update(excluded_query=sum(out_(r) for r in q_don), excluded_query_linked_stock=sum(r["rank_any_narrow"] == 1 for r in q_don if out_(r)),
                     excluded_other=len(o_ex), excluded_other_linked_stock=sum(r["rank_any_narrow"] == 1 for r in o_ex), excluded_other_cohorts=sorted({r["study"] for r in o_ex}),
                     excluded_other_min_rate=min(o_ex_rates) if o_ex_rates else None, excluded_other_max_rate=max(o_ex_rates) if o_ex_rates else None,
                     kept_other=len(o_kept), kept_other_linked_stock=sum(r["rank_any_narrow"] == 1 for r in o_kept), kept_other_cohorts=sorted({r["study"] for r in o_kept}),
                     kept_other_max_rate=max(o_kept_rates) if o_kept_rates else None, kept_query=sum(not out_(r) for r in q_don), candidates=len(cand))
allev =[(s, t, rt, kind(st, s, t)) for st in STUDIES for s, t, rt in ev[st]]
same_r = [e[2] for e in allev if e[3] == "same"]; other_r = [e[2] for e in allev if e[3] != "same"]      # an event that is not between two samples of one donor is counted as between two people (no unknown subject in these cohorts)
assert all(e[3] != "unknown" for e in allev)
allr = same_r + other_r

# ------------------------------------------------------------------ events between two people, cohort by cohort (the 11 cohorts)
studied = {st: {x for r in cand + rows_lo if r["study"] == st for x in (r["v1s"], r["v2s"], r["best_any"]) if x} for st in NAMES}     # the samples studied: the query sample, the first sample and the most similar other sample of each donor (in the cohorts kept, those of the 497 donors before the two rules)
assert all(studied[st] and studied[st] <= set(subj[st]) for st in NAMES), "the samples studied must be samples of the curated tables of their cohort"
by_cohort = {}
for st in NAMES:
    kinds = [(kind(st, s, t), rt, s, t) for s, t, rt in ev[st]]
    oth = [(rt, s, t) for k, rt, s, t in kinds if k == "other"]
    in_other = {x for _, s, t in oth for x in (s, t)}
    r_ = [x[0] for x in oth]
    by_cohort[st] = dict(samples=nsamples[st], events=len(kinds), same=sum(1 for k, *_ in kinds if k == "same"), other=len(oth), unknown=sum(1 for k, *_ in kinds if k == "unknown"),
                         median_rate_other=S.median(r_) if r_ else None, min_rate_other=min(r_) if r_ else None, max_rate_other=max(r_) if r_ else None,
                         q25_rate_other=float(np.percentile(r_, 25)) if r_ else None, q75_rate_other=float(np.percentile(r_, 75)) if r_ else None,
                         samples_in_other=len(in_other), percent_samples_in_other=100 * len(in_other) / nsamples[st], left_out=st in LEFT_OUT, rates_other=r_)
# ------------------------------------------------------------------ the same events, limited to those in which a sample studied takes part (Supplementary Note 1, Supplementary Fig. 4 and Supplementary Table 8)
by_st = {}
for st in NAMES:
    U = studied[st]
    kinds = [(kind(st, s, t), rt, s, t) for s, t, rt in ev[st] if s in U or t in U]                  # the sample studied is the source or the target
    oth = [(rt, s, t) for k, rt, s, t in kinds if k == "other"]
    r_ = [x[0] for x in oth]
    in_other = {x for _, s, t in oth for x in (s, t) if x in U}
    by_st[st] = dict(samples=nsamples[st], samples_studied=len(U), events=len(kinds), same=sum(1 for k, *_ in kinds if k == "same"), other=len(oth), unknown=sum(1 for k, *_ in kinds if k == "unknown"),
                     median_rate_other=S.median(r_) if r_ else None, min_rate_other=min(r_) if r_ else None, max_rate_other=max(r_) if r_ else None,
                     studied_in_other=len(in_other), percent_studied_in_other=100 * len(in_other) / len(U), left_out=st in LEFT_OUT, rates_other=r_)


def total(cs):
    r_ = [x for c in cs for x in by_st[c]["rates_other"]]
    return dict(samples=sum(nsamples[c] for c in cs), samples_studied=sum(len(studied[c]) for c in cs), events=sum(by_st[c]["events"] for c in cs), same=sum(by_st[c]["same"] for c in cs),
                other=len(r_), studied_in_other=sum(by_st[c]["studied_in_other"] for c in cs), donors=sum(1 for r in cand + rows_lo if r["study"] in cs),
                median_rate_other=S.median(r_) if r_ else None, min_rate_other=min(r_) if r_ else None, max_rate_other=max(r_) if r_ else None)


studied_kept, studied_left_out = total(STUDIES), total(LEFT_OUT)
ev_u = [(st, s, t, rt, kind(st, s, t)) for st in STUDIES for s, t, rt in ev[st] if s in studied[st] or t in studied[st]]     # the events of the cohorts kept in which a sample studied takes part
assert all(e[4] != "unknown" for e in ev_u) and len(ev_u) == studied_kept["events"]
assert studied_kept["studied_in_other"] == stats_exposed["samples_used"] and studied_kept["min_rate_other"] == min(stats_exposed["query_min_rate"], stats_exposed["other_min_rate"]) \
    and studied_kept["max_rate_other"] == max(stats_exposed["query_max_rate"], stats_exposed["other_max_rate"]), "in the cohorts kept, the samples studied in an event between two people, and the range of the rates of these events, must be those of the donors that Supplementary Note 1 counts"
# ------------------------------------------------------------------ Nielsen and Li analysed together
cross = None
if NL_CONTA:
    sj = {f"{tag}|{k}": v for tag, st in (("N", "NielsenHB_2014"), ("L", "LiJ_2014")) for k, v in subj[st].items()}
    nl_ev = [(d["source"], d["target"], float(d["rate"])) for d in csv.DictReader([l for l in open(NL_CONTA) if not l.startswith("#")], delimiter="\t")]
    assert all(s_ in sj and t_ in sj for s_, t_, _ in nl_ev)
    xe = [(s_, t_, rt) for s_, t_, rt in nl_ev if s_[0] != t_[0]]                               # a source in one cohort, a target in the other
    cls = lambda s_, t_: "same" if sj[s_] == sj[t_] else "other"
    groups = {"all": xe, "same": [e for e in xe if cls(e[0], e[1]) == "same"], "other": [e for e in xe if cls(e[0], e[1]) == "other"]}

    def summ(L):
        r_ = [e[2] for e in L]
        return dict(events=len(L), pairs=len({frozenset(e[:2]) for e in L}), samples_nielsen=len({x for e in L for x in e[:2] if x[0] == "N"}), samples_li=len({x for e in L for x in e[:2] if x[0] == "L"}),
                    median_rate=S.median(r_) if r_ else None, min_rate=min(r_) if r_ else None, max_rate=max(r_) if r_ else None, rates=r_, source_nielsen=sum(e[0][0] == "N" for e in L))
    cross = {k: summ(v) for k, v in groups.items()}
    old = {st: {(d["source"], d["target"]) for d in csv.DictReader([l for l in open(f"{CRO}/{st}.conta.tsv") if not l.startswith("#")], delimiter="\t")} for st in ("NielsenHB_2014", "LiJ_2014")}
    new = {st: {(s_[2:], t_[2:]) for s_, t_, _ in nl_ev if s_[0] == t_[0] == tag} for tag, st in (("N", "NielsenHB_2014"), ("L", "LiJ_2014"))}
    cross["within"] = {st: dict(separate_run=len(old[st]), merged_run=len(new[st]), common=len(old[st] & new[st])) for st in old}
    cross["samples"] = dict(nielsen=nsamples["NielsenHB_2014"], li=nsamples["LiJ_2014"], merged=len(sj), nielsen_studied=len(studied["NielsenHB_2014"]), li_studied=len(studied["LiJ_2014"]))
    is_st = lambda x: x[2:] in studied["NielsenHB_2014" if x[0] == "N" else "LiJ_2014"]          # a sample of the merged run (N| or L| and the sample) that is a sample studied
    for k, v in groups.items():
        L_ = [e for e in v if is_st(e[0]) or is_st(e[1])]
        cross[k].update(events_studied=len(L_), rates_studied=[e[2] for e in L_], median_rate_studied=S.median([e[2] for e in L_]) if L_ else None,
                        studied_nielsen=len({x for e in L_ for x in e[:2] if x[0] == "N" and is_st(x)}), studied_li=len({x for e in L_ for x in e[:2] if x[0] == "L" and is_st(x)}),
                        donors_nielsen_studied=len({sj[x] for e in L_ for x in e[:2] if x[0] == "N" and is_st(x)}))
    st_old = {st: {e for e in old[st] if e[0] in studied[st] or e[1] in studied[st]} for st in old}
    st_new = {st: {e for e in new[st] if e[0] in studied[st] or e[1] in studied[st]} for st in new}
    cross["within_studied"] = {st: dict(separate_run=len(st_old[st]), merged_run=len(st_new[st]), common=len(st_old[st] & st_new[st])) for st in st_old}
json.dump(dict(A=stats_A, exposed=stats_exposed, studied=dict(kept=studied_kept, left_out=studied_left_out, by_cohort={k: {a: b for a, b in v.items() if a != 'rates_other'} for k, v in by_st.items()}), nielsen_li_cross={k: ({a_: b_ for a_, b_ in v.items() if a_ not in ('rates', 'rates_studied')} if isinstance(v, dict) and 'rates' in v else v) for k, v in cross.items()} if cross else None, by_cohort={k: {a: b for a, b in v.items() if a != "rates_other"} for k, v in by_cohort.items()},
               n_same=len(same_r), n_other=len(other_r), median_rate_other=S.median(other_r), n_events=len(allr),
               kept=dict(events=len(allr), same=len(same_r), other=len(other_r), median_rate_other=S.median(other_r), q25_rate_other=float(np.percentile(other_r, 25)), q75_rate_other=float(np.percentile(other_r, 75)),
                         min_rate_other=min(other_r), max_rate_other=max(other_r),
                         samples=sum(nsamples[s] for s in STUDIES), samples_in_other=sum(by_cohort[s]["samples_in_other"] for s in STUDIES)),
               left_out=dict(events=sum(by_cohort[s]["events"] for s in LEFT_OUT), same=sum(by_cohort[s]["same"] for s in LEFT_OUT), other=sum(by_cohort[s]["other"] for s in LEFT_OUT),
                             median_rate_other=S.median([x for c in LEFT_OUT for x in by_cohort[c]["rates_other"]]), min_rate_other=min(x for c in LEFT_OUT for x in by_cohort[c]["rates_other"]), max_rate_other=max(x for c in LEFT_OUT for x in by_cohort[c]["rates_other"]),
                             samples=sum(nsamples[s] for s in LEFT_OUT), samples_in_other=sum(by_cohort[s]["samples_in_other"] for s in LEFT_OUT))), open(PREFIX + "_stats.json", "w"))

# ------------------------------------------------------------------ the figure: the samples studied, in the cohorts kept (grey) and left out (red)
order_k = sorted(STUDIES, key=lambda c: (-by_st[c]["percent_studied_in_other"], c))
order_o = sorted(LEFT_OUT, key=lambda c: (-by_st[c]["percent_studied_in_other"], c))
GAP = 1.5
ypos = {}
n_k, n_o = len(order_k), len(order_o)
for i, c in enumerate(order_k):
    ypos[c] = n_o + GAP + (n_k - 1 - i)
for i, c in enumerate(order_o):
    ypos[c] = n_o - 1 - i
YMAX = n_o + GAP + n_k - 1
W, H_AB = 6.5, 4.7
H_C = 1.75 if cross else 0.0                                              # height of panel C in inches
H = H_AB + H_C
fig = plt.figure(figsize=(W, H), dpi=300)
LEFT, BOT, HH = 0.205, (H_C + 0.105 * H_AB) / H, 0.775 * H_AB / H
axa = fig.add_axes((LEFT, BOT, 0.215, HH)); axb = fig.add_axes((0.495, BOT, 0.255, HH))
for ax in (axa, axb):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK3)
    ax.tick_params(length=2.5, width=0.5, color=INK3, labelcolor=INK2, labelsize=6.5)
    ax.set_ylim(-0.7, YMAX + 0.9); ax.grid(axis="x", color=GUIDE, linewidth=0.45, zorder=0); ax.set_axisbelow(True)
    ax.axhline(n_o + GAP / 2 - 0.5, color=RULE, lw=0.6, zorder=1)
axa.set_xlim(0, 100); axa.set_xticks([0, 25, 50, 75, 100])
axb.set_xlim(0, 100); axb.set_xticks([0, 25, 50, 75, 100]); axb.set_yticks([])
axa.set_yticks([ypos[c] for c in order_k + order_o]); axa.set_yticklabels([NAMES[c] for c in order_k + order_o], fontsize=6.8, color=INK)
axa.set_xlabel("Samples studied in an event between two people (%)", fontsize=6.8, color=INK, labelpad=3)
axb.set_xlabel("Rate (% of the target from the source)", fontsize=6.8, color=INK, labelpad=3)
ROLE = {}                                                                               # the samples studied in an event between two people, by their role in these events: source only, both, target only
for c in order_k + order_o:
    U = studied[c]; oth = [(s, t) for s, t, rt in ev[c] if (s in U or t in U) and kind(c, s, t) == "other"]
    S_, T_ = {s for s, t in oth if s in U}, {t for s, t in oth if t in U}
    assert len(S_ | T_) == by_st[c]["studied_in_other"], "a sample studied in an event between two people is its source, its target or both"
    ROLE[c] = (len(S_ - T_), len(S_ & T_), len(T_ - S_))
HATCH = (("As source", "//////"), ("As both", "xxxxxx"), ("As target", "\\\\\\\\\\\\"))   # two textures, crossed for a sample that is the source of one event and the target of another
plt.rcParams["hatch.linewidth"] = 0.55
rng = np.random.default_rng(5)
trB = blended_transform_factory(axb.transAxes, axb.transData)
for c in order_k + order_o:
    y = ypos[c]; col = C_OUT if c in LEFT_OUT else C_KEPT; b = by_st[c]
    left_ = 0.0
    for (lab_, h_), n_ in zip(HATCH, ROLE[c]):                                          # one segment per role, white hatching on the colour of the cohort
        if n_:
            w_ = 100 * n_ / b["samples_studied"]
            axa.barh(y, w_, left=left_, height=0.62, color=col, edgecolor="white", linewidth=0, hatch=h_, zorder=2); left_ += w_
    assert abs(left_ - b["percent_studied_in_other"]) < 1e-9, "the three roles make the bar"
    axa.text(b["percent_studied_in_other"] + 1.5, y, f"{b['studied_in_other']}/{b['samples_studied']}", ha="left", va="center", fontsize=6.0, color=INK2, zorder=4)
    r_ = b["rates_other"]
    if r_:
        axb.scatter(100 * np.asarray(r_), y + rng.uniform(-0.26, 0.26, len(r_)), s=6.0, color=col, alpha=0.60, linewidths=0, zorder=3)
        axb.plot([100 * b["median_rate_other"]] * 2, [y - 0.36, y + 0.36], color=INK, lw=1.4, solid_capstyle="butt", zorder=5)
    axb.text(1.20, y, f"{b['events']}", transform=trB, ha="right", va="center", fontsize=6.2, color=INK2)
    axb.text(1.78, y, f"{b['other']} ({100 * b['other'] / b['events']:.0f}%)" if b["events"] else "0", transform=trB, ha="right", va="center", fontsize=6.2, color=INK2)
for x_, s_ in ((1.20, "Events"), (1.78, "Between\ntwo people")):
    axb.text(x_, YMAX + 0.55, s_, transform=trB, ha="right", va="bottom", fontsize=5.8, color=INK3, clip_on=False, linespacing=1.1)
axa.text(-0.90, n_o + GAP + n_k - 0.2 + 0.55, "Cohorts kept", transform=blended_transform_factory(axa.transAxes, axa.transData), fontsize=7.0, fontweight="bold", color=INK, ha="left", va="bottom", clip_on=False)
axa.text(-0.90, n_o - 0.25, "Cohorts left out", transform=blended_transform_factory(axa.transAxes, axa.transData), fontsize=7.0, fontweight="bold", color=C_OUT, ha="left", va="bottom", clip_on=False)
for t_, x0 in (("A", 0.006), ("B", 0.452)):
    fig.text(x0, (H_C + 0.985 * H_AB) / H, t_, fontsize=10, fontweight="bold", color=INK, ha="left", va="top")
fig.text(0.495, (H_C + 0.930 * H_AB) / H, "Each dot is an event between two people; bar, median", fontsize=6.2, color=INK3, ha="left", va="bottom")
fig.legend(handles=[matplotlib.patches.Patch(facecolor=INK3, edgecolor="white", linewidth=0, hatch=h_, label=lab_) for lab_, h_ in HATCH], loc="lower left",
           bbox_to_anchor=(LEFT, (H_C + 0.930 * H_AB) / H), ncol=3, frameon=False, fontsize=6.2, labelcolor=INK3, handlelength=1.5, handleheight=1.0,
           columnspacing=0.9, handletextpad=0.4, borderaxespad=0, borderpad=0)   # the role of the sample studied in panel A
CONTA_RULE = 5                                                                          # the rule of contamination (Supplementary Note 1): in the cohorts kept, a sample studied in an event between two people at a rate above 5% is removed
axb.plot([CONTA_RULE, CONTA_RULE], [n_o + GAP / 2 - 0.5, YMAX + 0.62], color=INK, lw=0.8, ls=(0, (3, 2)), zorder=4)
axb.text(CONTA_RULE + 1.5, YMAX + 0.62, f"{CONTA_RULE}%, removal threshold", fontsize=5.8, color=INK, ha="left", va="center", zorder=6)
if cross:
    rows_c = [("all", "All events between the two cohorts"), ("same", "Same identifier"), ("other", "Different identifiers")]
    axc = fig.add_axes((0.495, 0.49 / H, 0.255, 0.80 / H))
    for s_ in ("top", "right"):
        axc.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"):
        axc.spines[s_].set_color(INK3)
    axc.tick_params(length=2.5, width=0.5, color=INK3, labelcolor=INK2, labelsize=6.5)
    axc.set_xlim(0, 100); axc.set_xticks([0, 25, 50, 75, 100]); axc.set_yticks([])
    axc.set_ylim(-0.6, 2.6); axc.grid(axis="x", color=GUIDE, linewidth=0.45, zorder=0); axc.set_axisbelow(True)
    axc.set_xlabel("Rate (% of the target from the source)", fontsize=6.8, color=INK, labelpad=3)
    trC = blended_transform_factory(axc.transAxes, axc.transData)
    rngc = np.random.default_rng(7)
    for i, (key, lab) in enumerate(rows_c):
        y = 2 - i; g = cross[key]
        if g["rates_studied"]:
            axc.scatter(100 * np.asarray(g["rates_studied"]), y + rngc.uniform(-0.22, 0.22, len(g["rates_studied"])), s=6.0, color=C_OUT, alpha=0.60, linewidths=0, zorder=3)
            axc.plot([100 * g["median_rate_studied"]] * 2, [y - 0.34, y + 0.34], color=INK, lw=1.4, solid_capstyle="butt", zorder=5)
        axc.text(-0.03, y, lab, transform=trC, ha="right", va="center", fontsize=6.8, color=INK, fontweight="bold" if key == "all" else "normal")
        axc.text(1.12, y, f"{g['events_studied']}", transform=trC, ha="right", va="center", fontsize=6.2, color=INK2)
        axc.text(1.50, y, f"{g['studied_nielsen']}/{cross['samples']['nielsen_studied']}", transform=trC, ha="right", va="center", fontsize=6.2, color=INK2)
        axc.text(1.82, y, f"{g['studied_li']}/{cross['samples']['li_studied']}", transform=trC, ha="right", va="center", fontsize=6.2, color=INK2)
    axc.axhline(1.5, color=RULE, lw=0.6, zorder=1)
    for x_, s_ in ((1.12, "Events"), (1.50, "Nielsen\nstudied"), (1.82, "Li\nstudied")):
        axc.text(x_, 2.75, s_, transform=trC, ha="right", va="bottom", fontsize=5.8, color=INK3, clip_on=False, linespacing=1.1)
    fig.text(0.006, (0.49 + 0.80 + 0.30) / H, "C", fontsize=10, fontweight="bold", color=INK, ha="left", va="bottom")
    fig.text(0.045, (0.49 + 0.80 + 0.30) / H, "Nielsen and Li analysed together: events between a Nielsen sample and a Li sample", fontsize=6.8, color=INK, ha="left", va="bottom")
fig.savefig(PREFIX + ".png", dpi=600)
fig.savefig(PREFIX + ".pdf")
print("wrote", PREFIX + ".png and .pdf")
print("samples studied")
print(f"{'cohort':16} {'samples':>7} {'studied':>7} {'events':>6} {'same':>5} {'other':>5} {'studied in other':>16} {'%':>5} {'median rate':>11} {'min-max':>10}")
for c in order_k + order_o:
    b = by_st[c]
    med_s = "--" if b["median_rate_other"] is None else "%.2f" % b["median_rate_other"]
    rng_s = "--" if b["min_rate_other"] is None else "%.2f-%.2f" % (b["min_rate_other"], b["max_rate_other"])
    print(f"{NAMES[c]:16} {b['samples']:7d} {b['samples_studied']:7d} {b['events']:6d} {b['same']:5d} {b['other']:5d} {b['studied_in_other']:16d} {b['percent_studied_in_other']:5.1f} {med_s:>11} {rng_s:>10}")
print("kept", studied_kept)
print("left out", studied_left_out)
if cross:
    print("Nielsen and Li merged", {k: {a: b for a, b in v.items() if a not in ("rates", "rates_studied")} for k, v in cross.items() if k in ("all", "same", "other", "within_studied", "samples")})
print(f"A: {stats_A}\nkept: {len(allr)} events, {len(same_r)} same, {len(other_r)} other, median rate of the events between two people {S.median(other_r):.2f}")
