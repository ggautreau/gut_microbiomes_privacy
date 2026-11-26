#!/usr/bin/env python3
"""Supplementary Figure 1: the selection of the cohorts, the donors and their samples, from curatedMetagenomicData to the donors of the analysis (Supplementary Note 1 and
Supplementary Tables 1 and 2). The file is named supplementary_figure1; build_supplementary.py places the figure first.

    python3 -I supplementary_figure1.py COHORT_SELECTION.json DONORS.json DONORS_EXCLUDED.json OUTPUT_PREFIX     # writes OUTPUT_PREFIX.png (600 dpi) and .pdf

COHORT_SELECTION.json  output of cohort_selection.py (aggregates/cohort_selection.json): the steps of Supplementary Table 1, the donors of each adult cohort (Supplementary Table 2)
                       and the arrangements of the samples of the donors along time (sample_arrangements), drawn in panel B
DONORS.json            the donors of the analysis (data/donors.json); DONORS_EXCLUDED.json the donors of the cohorts kept that the rules of contamination and re-deposit remove (data/donors_excluded.json)
The samples studied of a donor are its query sample, its first sample and its most similar other sample (v2s, v1s and best_any of the rows), each counted once, two or three per donor.
Every number is read from these files. The steps are those of Supplementary Table 1; the last two are named by what they keep. Palette and type follow Figure 1; red, what a step removes.
"""
import json, sys, textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

SEL, DON, EXC, PREFIX = sys.argv[1:5]
INK, INK2, INK3, RULE, FILL, C_OUT = "#1B2227", "#3F4A52", "#6B777F", "#C9D0D4", "#F2F4F5", "#C8553D"
C_SK, C_SK_INK, G_F1 = "#D55E00", "#8A3D00", "#5C646C"                           # the orange, its dark shade for words and the grey of the samples of Fig. 1A
import os
from matplotlib import font_manager
FONTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fonts")                  # Liberation Sans shipped with the package (licences in fonts/), so that the drawing does not depend on the machine
font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if f.name != "Liberation Sans"]
for fn in sorted(os.listdir(FONTS)):
    if fn.lower().endswith(".ttf"):
        font_manager.fontManager.addfont(os.path.join(FONTS, fn))
plt.rcParams.update({"font.family": ["Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42})
NAMES = {"PehrssonE_2016": "Pehrsson", "NielsenHB_2014": "Nielsen", "LiJ_2014": "Li", "VincentC_2016": "Vincent", "CosteaPI_2017": "Costea"}

sc = json.load(open(SEL)); donors = json.load(open(DON)); exc = json.load(open(EXC))
cat, two, lg, lb, ac, at, ex = sc["catalogue"], sc["two_stool"], sc["longitudinal"], sc["labels"], sc["adult_cohorts"], sc["adult_totals"], sc["excluded"]
n_idx = sum(1 for c in ac.values() if c["indexed"]); n_dna = sum(1 for c in ac.values() if c["adult_dna"]); n_kept = sum(1 for c in ac.values() if c["analysed"])
none_idx = sum(1 for c in ac.values() if c["donors"] and not c["indexed"])
lost_dna = [st for st, c in ac.items() if c["indexed"] and not c["adult_dna"]]      # cohorts whose donors in the index are all removed at the step of the adults
assert all(ac[st]["under18"] == ac[st]["indexed"] for st in lost_dna), "a cohort that the step of the adults removes loses its donors because they are under 18"
bio = sc["biosample_only_in_index"]                                              # donors whose indexed run of a sample is in its BioSample but not among its curated accessions
assert list(bio) == ["MehtaRS_2018"], "the donors with an indexed run only in the BioSample are all of Mehta"
left = {st: c["adult_dna"] for st, c in ac.items() if c["adult_dna"] and not c["analysed"]}
studied_set = lambda rows: {(r["study"], x) for r in rows for x in (r["v1s"], r["v2s"], r["best_any"]) if x}
S_ana, S_exc = studied_set(donors), studied_set(exc)
s_cand, s_ana, s_exc = len(S_ana | S_exc), len(S_ana), len(S_exc)
assert len(donors) == at["analysed"] and len(donors) + len(exc) == ex["candidates"] == at["adult_dna"] - sum(left.values()), "the donors, the donors removed and the steps of cohort_selection.py must add up"
assert not S_ana & S_exc and s_cand - s_exc == s_ana, "the donors removed take away samples studied of their own only"
assert {e["study"] for e in exc if e["rule"] == "re-deposit"} == {"CosteaPI_2017"}, "the rule of re-deposit removes donors of Costea only"
f = lambda n: f"{n:,}"
n_before = ex["candidates"]                                                      # the donors of the cohorts kept, before the two rules (key of cohort_selection.py)
count = lambda n, one, many: f"{n} {one if n == 1 else many}"
no_stool = cat["studies"] - cat["stool_studies"]                                 # studies without a stool sample
no_two = cat["stool_studies"] - two["studies"]                                   # studies whose donors all have a single stool sample

# each step: its label (that of Supplementary Table 1, except the last two, named by what they keep), what it keeps, and what it removes (right, in red; paragraphs separated by a blank line)
STEPS = [
    ("curatedMetagenomicData", f"{cat['studies']} studies, {f(cat['samples'])} samples", None),
    ("Stool samples", f"{cat['stool_studies']} studies, {f(cat['stool'])} samples",
     f"{f(cat['samples'] - cat['stool'])} samples of other body sites, among them all those of {count(no_stool, 'study', 'studies')}"),
    ("Donors with two or more stool samples", f"{two['studies']} studies, {f(two['donors'])} donors",
     f"Donors with a single stool sample, among them\nall those of {count(no_two, 'study', 'studies')}"),
    (f"Longitudinal cohorts ({lg['min_donors']} or more such donors)", f"{lg['cohorts']} cohorts, {f(lg['donors'])} donors",
     f"{lg['below_threshold_studies']} studies with fewer such donors ({lg['below_threshold_donors']} donors)"),
    ("Cohorts labelled adults", f"{lb['adults']['cohorts']} cohorts, {at['donors']} donors",
     f"{lb['infants']['cohorts']} cohorts labelled infants ({lb['infants']['donors']} donors), {lb['mixed']['cohorts']} mixed ({lb['mixed']['donors']}) and {lb['not recorded']['cohorts']} not labelled for lack of ages ({lb['not recorded']['donors']})"),
    ("Adult cohorts, both samples in the index", f"{n_idx} cohorts, {at['indexed']} donors",
     f"{at['donors'] - at['indexed']} donors without a curated run in the index for each of their two samples, "
     f"among them all those of {count(none_idx, 'cohort', 'cohorts')} and {bio['MehtaRS_2018']} of Mehta whose indexed run is only in the BioSample of the sample"),
    ("Adults (18 years or more, or Mehta) with a DNA query", f"{n_dna} cohorts, {at['adult_dna']} donors",
     f"{at['no_age']} donors without a recorded age, {at['rna']} whose query is an RNA library and {at['under18']} under 18"
     + (f", among them all those left in {count(len(lost_dna), 'cohort', 'cohorts')}" if lost_dna else "")),
    ("Donors of the cohorts kept", f"{n_kept} cohorts, {n_before} donors, {f(s_cand)} samples studied",
     f"{len(left)} cohorts left out, {sum(left.values())} donors (" + ", ".join(f"{NAMES[st]} {left[st]}" for st in sorted(left, key=lambda s: -left[s])) + ")"),
    ("Donors of the analysis", f"{n_kept} cohorts, {at['analysed']} donors, {f(s_ana)} samples studied",
     f"{ex['samples_removed']} samples studied that take part in a CroCoDeEL event between two people at a rate above 5%, leaving {ex['contamination']} donors without a query or another sample\n\n"
     f"{ex['re_deposit']} donors of Costea whose nearest run of other people carries the reads of their own sample under another subject code\n\n"
     f"In all, {ex['contamination'] + ex['re_deposit']} donors and their {s_exc} samples studied"),
]
# the samples of each donor, drawn as circles along time with the signs of Fig. 1A: the samples studied large, the others small; a donor has two or three samples studied
per_donor = [len({x for x in (r["v1s"], r["v2s"], r["best_any"]) if x}) for r in donors]
n3, n2 = per_donor.count(3), per_donor.count(2)
n2_none = sum(1 for r in donors if not r["best_any"]); n2_first = sum(1 for r in donors if r["best_any"] and r["best_any"] == r["v1s"])
assert n2 + n3 == len(donors) and n2 == n2_none + n2_first and sum(per_donor) == s_ana, "two or three samples studied per donor, none shared between donors"
# one line per arrangement of the samples of a donor along time, apart for the donors recovered and for the others (sample_arrangements of cohort_selection.py).
# Roles: query; first (the first sample); first_and_most_similar (the first sample is also the most similar other sample); first_none_returned (the first sample of a donor
# with no other sample among the runs returned); most_similar (the most similar other sample); other (another sample)
sa = sc["sample_arrangements"]
BLOCKS = {True: [], False: []}                                                   # recovered or not: (roles along time, donors, note)
for a_ in sa["listed"]:
    BLOCKS[a_["recovered"]].append((tuple(a_["roles"]), a_["donors"], None))
for a_ in sa["more"]:                                                           # the donors with more samples than are listed: a line of circles cut by an ellipsis, and their numbers of samples
    BLOCKS[a_["recovered"]].append((("first", "other", "other", "...", "other", "most_similar", "other", "query") if a_["most_similar_returned"] else ("first_none_returned", "other", "other", "...", "other", "other", "query"),
                                    a_["donors"], f"{a_['min_samples']} samples" if a_["min_samples"] == a_["max_samples"] else f"{a_['min_samples']} to {a_['max_samples']} samples"))
N_REC, N_NOT = sum(b[1] for b in BLOCKS[True]), sum(b[1] for b in BLOCKS[False])
with_role = lambda role: sum(b[1] for k_ in BLOCKS for b in BLOCKS[k_] if role in b[0])
assert N_REC + N_NOT == len(donors) and N_REC == sum(1 for r in donors if r["rank_any_narrow"] == 1), "the arrangements of cohort_selection.py must hold every donor, the donors recovered apart"
assert with_role("first_none_returned") == n2_none and 0 <= n2_first - with_role("first_and_most_similar") <= sum(a_["donors"] for a_ in sa["more"] if a_["most_similar_returned"]), "the arrangements must agree with the donors on the other samples returned"
N_LINES = max(len(BLOCKS[True]), len(BLOCKS[False]))
DEFS = [("query", "Query sample", "the last stool sample, queried through its largest curated run in the index"),
        ("first", "First sample", "the earliest stool sample, marked 1, ranked in place of the most similar other sample in Fig. 1E"),
        ("most_similar", "Most similar\n(nearest)\nother sample", "of the other stool samples of the donor, the one most similar to the query, marked 1 when it is also the first sample"),
        ("first_none_returned", "Not returned", f"the first sample of the {n2_none} donors with no other sample among the runs returned"),
        ("other", "Other sample", "any other stool sample of the donor, neither the first nor the most similar")]

# layout in inches, from the top: the steps in a column on the left, what each step removes on the right, beside the arrow that leads to the step, and under the chart,
# across the whole width, the samples of a donor (circles on the left, the words on the right)
W = 6.5
XL, WL = 0.06, 3.0                                                              # the column of the steps
XR, WR = 3.36, 3.08                                                             # the column on the right
LH, LS, LO = 0.122, 0.112, 0.100                                                # height of a line at 7 pt, at 6.6 pt in the samples box, at 6.6 pt in a red box
BOX = 2 * LH + 0.08                                                             # a step: its name and what it keeps
wrap = lambda s_, n: "\n".join(textwrap.fill(q, n, break_on_hyphens=False) for q in s_.split("\n"))     # never at the hyphen of re-deposit
outs = [None if o is None else wrap(o, 68) for _, _, o in STEPS]
paras = [None if o is None else o.split("\n\n") for o in outs]
hout = [None if p is None else sum(q.count("\n") + 1 for q in p) * LO + (len(p) - 1) * 0.5 * LO + 0.075 for p in paras]     # the lines as drawn at 6.6 pt, a blank line counting half
TERM = 0.88                                                                     # room for the terms of the samples box
defs = [(k_, t_, textwrap.fill(d_, 35, break_on_hyphens=False)) for k_, t_, d_ in DEFS]
XD = XL + 3.56                                                                  # the words, right of the circles
R_B, R_Q, R_S = 0.046, 0.042, 0.017                                              # radius of a disc, of a ring, of another sample
ys, y = [], 0.20                                                               # under the letter of panel A
for i in range(len(STEPS)):
    if i:
        y += 0.20 if hout[i] is None else hout[i] + 0.08
    ys.append(y)
    y += BOX
last = len(STEPS) - 1
y_box = ys[-1] + BOX + 0.24                                                     # the samples box, under the chart
XB0, WB0 = XL, XR + WR - XL
PITCH = 0.168                                                                   # between two lines of circles
y_h = y_box + 0.08 + LH + 0.13                                                  # the headings of the two blocks
y_r = [y_h + 0.20 + PITCH * i_ for i_ in range(N_LINES)]                         # the centres of the lines of circles
y_bar = y_r[-1] + 0.105                                                         # the bar under the counts
y_tot = y_bar + 0.095                                                           # the total
y_t = y_tot + 0.17                                                              # the arrow of time
h_pict = y_t + 0.10 - y_box
lines = lambda t_, d_: max(t_.count("\n"), d_.count("\n")) + 1                   # lines of a definition, its term included
h_defs = 0.08 + LH + 0.05 + sum(lines(t_, d_) * LS + 0.04 for _, t_, d_ in defs) + 0.02
hbox = max(h_pict, h_defs) + 0.04
H = y_box + hbox + 0.06
fig = plt.figure(figsize=(W, H), dpi=300)
ax = fig.add_axes((0, 0, 1, 1)); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
frames = []                                                                     # (text, its box): every text must stay inside its box
ax.text(XL - 0.01, 0.03, "A", fontsize=10, fontweight="bold", color=INK, ha="left", va="top", zorder=3)                   # the two panels: A the steps, B the samples of each donor
ax.text(XL - 0.01, y_box - 0.19, "B", fontsize=10, fontweight="bold", color=INK, ha="left", va="top", zorder=3)
mid = XL + WL / 2
for i, (title, keep, _) in enumerate(STEPS):
    end = i == last
    y = ys[i]
    ax.add_patch(Rectangle((XL, y), WL, BOX, facecolor=FILL if end else "white", edgecolor=INK if end else INK3, lw=0.9 if end else 0.6, zorder=2))
    t1 = ax.text(XL + 0.08, y + 0.045, title, fontsize=7.0, fontweight="bold", color=INK, ha="left", va="top", zorder=3)
    t2 = ax.text(XL + 0.08, y + 0.045 + LH, keep, fontsize=7.0, color=INK2, ha="left", va="top", zorder=3)
    frames += [(t1, (XL, y, WL, BOX)), (t2, (XL, y, WL, BOX))]
    if i:
        ax.add_patch(FancyArrowPatch((mid, ys[i - 1] + BOX), (mid, y), arrowstyle="-|>", mutation_scale=6, color=INK3, lw=0.6, zorder=1, shrinkA=0, shrinkB=0))
    if outs[i]:
        yc = (ys[i - 1] + BOX + y) / 2                                           # the middle of the arrow that leads to the step
        yr = yc - hout[i] / 2
        ax.add_patch(Rectangle((XR, yr), WR, hout[i], facecolor="white", edgecolor=C_OUT, lw=0.6, zorder=2))
        yy = yr + 0.04
        for par in paras[i]:                                                     # a paragraph per rule
            t3 = ax.text(XR + 0.08, yy, par, fontsize=6.6, color=INK2, ha="left", va="top", linespacing=1.22, zorder=3)
            frames.append((t3, (XR, yr, WR, hout[i])))
            yy += (par.count("\n") + 1.5) * LO
        ax.add_patch(FancyArrowPatch((mid, yc), (XR, yc), arrowstyle="-|>", mutation_scale=6, color=C_OUT, lw=0.6, zorder=1, shrinkA=0, shrinkB=0))
ax.add_patch(Rectangle((XB0, y_box), WB0, hbox, facecolor=FILL, edgecolor=INK3, lw=0.6, zorder=2))
FR_P, FR_D = (XB0, y_box, XD - 0.08 - XB0, hbox), (XD - 0.08, y_box, XB0 + WB0 - XD + 0.08, hbox)        # the frames of the circles and of the words
t5 = ax.text(XB0 + 0.08, y_box + 0.08, "The samples of each donor", fontsize=7.0, fontweight="bold", color=INK, ha="left", va="top", zorder=3)
frames.append((t5, FR_P))


def sign(kind, x_, y_, col=None):
    """The signs of Fig. 1A: the query is a white disc with an orange ring; a disc is orange for a donor recovered and grey otherwise."""
    if kind == "query":
        ax.add_patch(Circle((x_, y_), R_Q, facecolor="white", edgecolor=C_SK, lw=1.0, zorder=3))
    elif kind in ("first", "most_similar", "first_and_most_similar"):
        ax.add_patch(Circle((x_, y_), R_B, facecolor=col or C_SK, edgecolor="none", zorder=3))
    elif kind == "first_none_returned":                                          # no other sample among the runs returned: a grey ring
        ax.add_patch(Circle((x_, y_), R_Q, facecolor=FILL, edgecolor=G_F1, lw=1.0, zorder=3))
    else:                                                                        # another sample: a small black dot, neither orange nor grey
        assert kind == "other", kind
        ax.add_patch(Circle((x_, y_), R_S, facecolor=INK, edgecolor="none", zorder=3))


def one(x_, y_, frame):
    """The subscript of the first sample."""
    frames.append((ax.text(x_ + R_B * 0.80, y_ + R_B * 0.95, "1", fontsize=4.8, fontweight="bold", color=INK, ha="left", va="center", zorder=4), frame))


for rec_, xo_ in ((True, XB0 + 0.08), (False, XB0 + 1.74)):
    col_, ink_ = (C_SK, C_SK_INK) if rec_ else (G_F1, INK2)
    sign("most_similar", xo_ + 0.10, y_h, col_)                                  # the heading carries the colour, as at the bottom of Fig. 1A
    frames.append((ax.text(xo_ + 0.20, y_h, "Recovered" if rec_ else "Not recovered", fontsize=6.8, fontweight="bold", color=ink_, ha="left", va="center", zorder=3), FR_P))
    xa, xb = xo_ + 0.10, xo_ + 1.02                                              # the line of time of a donor
    for (roles_, n_, note_), yc_ in zip(BLOCKS[rec_], y_r):
        ax.plot([xa - 0.06, xb + 0.04], [yc_, yc_], color=RULE, lw=0.6, zorder=2.5, solid_capstyle="round")
        for j_, r_ in enumerate(roles_):
            x_ = xa + (xb - xa) * j_ / (len(roles_) - 1)
            if r_ == "...":                                                      # more samples than are drawn
                frames.append((ax.text(x_, yc_ - 0.012, "\u2026", fontsize=6.6, color=INK3, ha="center", va="center", zorder=3), FR_P))
                continue
            sign(r_, x_, yc_, col_)
            if r_.startswith("first"):
                one(x_, yc_, FR_P)
        frames.append((ax.text(xb + 0.33, yc_, f"{n_}", fontsize=6.6, fontweight="bold", color=INK, ha="right", va="center", zorder=3), FR_P))
        if note_:
            frames.append((ax.text((xa + xb) / 2, yc_ + 0.075, note_, fontsize=5.6, color=INK3, ha="center", va="top", zorder=3), FR_P))
    ax.plot([xb + 0.10, xb + 0.35], [y_bar, y_bar], color=INK, lw=0.6, solid_capstyle="butt", zorder=3)        # the total of the block, under a bar
    frames.append((ax.text(xb + 0.33, y_tot, f"{N_REC if rec_ else N_NOT}", fontsize=6.6, fontweight="bold", color=INK, ha="right", va="center", zorder=3), FR_P))
    frames.append((ax.text(xb + 0.37, y_tot, "donors", fontsize=6.6, color=INK2, ha="left", va="center", zorder=3), FR_P))
    ax.add_patch(FancyArrowPatch((xa - 0.06, y_t), (xb + 0.04, y_t), arrowstyle="-|>", mutation_scale=5, color=INK3, lw=0.5, zorder=3, shrinkA=0, shrinkB=0))
    frames.append((ax.text(xb + 0.10, y_t, "time", fontsize=6.2, color=INK3, ha="left", va="center", zorder=3), FR_P))
yy = y_box + 0.08 + LH + 0.05
for k_, t_, d_ in defs:
    sign(k_, XD + 0.03, yy + 0.05)
    if k_.startswith("first"):
        one(XD + 0.03, yy + 0.05, FR_D)
    t6 = ax.text(XD + 0.19, yy, t_, fontsize=6.6, fontweight="bold", color=INK2, ha="left", va="top", linespacing=1.22, zorder=3)
    t4 = ax.text(XD + 0.19 + TERM, yy, d_, fontsize=6.6, color=INK2, ha="left", va="top", linespacing=1.22, zorder=3)
    frames += [(t6, FR_D), (t4, FR_D)]
    yy += lines(t_, d_) * LS + 0.04
rend = fig.canvas.get_renderer(); inv = ax.transData.inverted()
for t_, (x0, y0, w0, h0) in frames:                                           # in inches of the figure, y downwards
    (a0, b0), (a1, b1) = inv.transform(t_.get_window_extent(rend).get_points())
    assert x0 + 0.03 <= a0 and a1 <= x0 + w0 - 0.03 and y0 + 0.01 <= min(b0, b1) and max(b0, b1) <= y0 + h0 - 0.01, (t_.get_text()[:50], a1 - x0 - w0, max(b0, b1) - y0 - h0)
fig.savefig(PREFIX + ".png", dpi=600)
fig.savefig(PREFIX + ".pdf")
print("wrote", PREFIX + ".png and .pdf", f"({W} x {H:.2f} in)")
print({"donors_before_the_two_rules": ex["candidates"], "samples_studied_before_the_two_rules": s_cand, "donors": at["analysed"], "samples_studied_donors": s_ana, "samples_studied_removed": s_exc, "donors_with_3": n3, "donors_with_2": n2, "other_is_first": n2_first, "other_not_returned": n2_none})
