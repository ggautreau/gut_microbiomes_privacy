#!/usr/bin/env python3
"""How the numbers move with the runs that count as other people (Methods, "HGM set"; caption of Supplementary Table 4). Aggregates only.

    python3 -I stock_sensitivity.py DATA_DIR LABELLED.txt.gz ENA.tsv ROWS_V4.json OUT.json [--hmp LINKS.json]

DATA_DIR      the DATA folder of rows_with_partner.py (results.jsonl, logan_index_ids.txt.gz, logan_gut_wgs_narrow_ids.txt.gz = the HGM set, biosample_runs.json, curated_sample_biosample.json, runs_hmp2.json, ena_PRJNA389280_runs.tsv, matrices/, and mehta_sibling_runs.json when present)
LABELLED      the 172,075 runs labelled human gut or faecal shotgun metagenome in the Logan metadata
ENA.tsv       public metadata of the other runs considered for the HGM set (build_stock.py)
ROWS_V4.json  the 483 donors of the analysis, with the visit codes of HMP IBD (rows_v4.json, exclude_donors.py)

For each set of runs: the donors first among them (no run of another person above the other sample, as in rows_with_partner.py), the threshold of the open-set test (the equal error rate: the similarity at which as many searches of an absent person
are accepted as donors present are rejected, as in figure1_brief_communication.py), the donors whose other sample reaches it and is the top hit, and the straight line of Fig. 1B at the population of the world. Sets: the HGM set (the paper); the same without the projects
that hold the same people as a cohort, whose runs count as other people (the Mehta runs that the curated BioSamples do not name are identified by their alias, mehta_siblings.py): PRJEB24041 deposits the reads of part of Costea again (its sample names are those of Costea with a final underscore),
PRJEB28237 holds cultures of the stool of the patients of Raymond (same patient codes), PRJNA306874 ("iHMP ibdmdb"; its sample titles carry the curated sample identifiers of HMP IBD), PRJNA389280 (a second deposit of the metagenomes of HMP IBD, paired with its curated samples by read counts);
every run of the index; the labelled runs alone. Size of the pools: M = the HGM set, and 172,075 for the labelled runs; scale_of_the_pools gives the line of the paper for other values of M.
"""
import csv, gzip, json, math, os, sys
import numpy as np
from scipy.special import gammaln

DATA, LABELLED, ENA, ROWS, OUT = sys.argv[1:6]
LINKS = json.load(open(sys.argv[sys.argv.index('--hmp') + 1])) if '--hmp' in sys.argv else {}   # option --hmp LINKS.json (hmp_sibling_links.py): the runs that public records give to the donor are its own
DUP = 0.99
ids = lambda p: set(gzip.open(p, "rt").read().split())
labelled, stock, index = ids(LABELLED), ids(f"{DATA}/logan_gut_wgs_narrow_ids.txt.gz"), ids(f"{DATA}/logan_index_ids.txt.gz")
M_GUT = len(stock)
csv.field_size_limit(10**9)
project = {p: set() for p in ("PRJEB24041", "PRJEB28237", "PRJNA306874")}
hmp_second = {r["run_accession"] for r in csv.DictReader(open(f"{DATA}/ena_PRJNA389280_runs.tsv"), delimiter="\t", quoting=csv.QUOTE_NONE)}      # PRJNA389280: its runs are labelled, so ENA.tsv does not list them (mehta_siblings.py fetch)
for r in csv.DictReader(open(ENA), delimiter="\t", quoting=csv.QUOTE_NONE):
    if r["study_accession"] in project:
        project[r["study_accession"]].add(r["run_accession"])
rows = json.load(open(ROWS)); studies = {r["study"] for r in rows}
stool = {(st, m["sample_id"]): m["subject_id"] for st in studies for m in csv.DictReader(open(f"{DATA}/matrices/{st}.metadata.tsv"), delimiter="\t") if m["body_site"] == "stool" and m["subject_id"]}
bs = {b: tuple(k.split("|", 1)) for k, b in json.load(open(f"{DATA}/curated_sample_biosample.json")).items() if tuple(k.split("|", 1)) in stool}
owner = {run: (bs[v["biosample"]][0], stool[bs[v["biosample"]]]) for run, v in json.load(open(f"{DATA}/biosample_runs.json")).items() if v["biosample"] in bs}
owner.update({run: (v[0], v[1]) for run, v in json.load(open(f"{DATA}/runs_hmp2.json")).items() if v[0] == "HMP_2019_ibdmdb"})
# runs of the Mehta project that the curated BioSamples do not name (DNA and RNA libraries of one sample sit in two BioSamples), identified by their ENA alias (mehta_siblings.py)
if os.path.exists(f"{DATA}/mehta_sibling_runs.json"):
    owner.update({run: (v[0], v[1]) for run, v in json.load(open(f"{DATA}/mehta_sibling_runs.json")).items() if run not in owner})
answers = {}
for line in open(f"{DATA}/results.jsonl"):
    q = json.loads(line)
    answers.setdefault((q["study"], q["subject"]), q)
other = []                                                                            # for each donor, the runs of other people in its list
for r in rows:
    key = (r["study"], r["subject"])
    other.append({a: s for a, s in answers[key]["results"].items() if s < DUP and owner.get(a) != key and tuple(LINKS.get(a, ())) != key})
partner = np.array([r["partner_any"] or 0.0 for r in rows]); study = np.array([r["study"] for r in rows]); D = len(rows)
XW = math.log2(8e9)


def line_at_world(k, M):
    Ns = np.unique(np.round(np.logspace(np.log10(2), np.log10(M), 70)).astype(int)); Ns = Ns[Ns >= 1000]; f = k >= 0; y = []
    for N in Ns:
        ok = f & (M - k >= N - 1)
        y.append(np.exp(gammaln(M - k[ok] + 1) - gammaln(M - k[ok] - N + 2) - gammaln(M + 1) + gammaln(M - N + 2)).sum() / D)
    return 100 * float(np.clip(np.polyval(np.polyfit(np.log2(Ns.astype(float)), y, 1), XW), 0, 1))


out = {}
if LINKS:
    keys = {(r["study"], r["subject"]) for r in rows}
    out["hmp_links"] = dict(runs=len(LINKS), participants=len({v[1] for v in LINKS.values()}), donors=len({tuple(v) for v in LINKS.values()} & keys),
                            runs_in_the_lists=sum(1 for r in rows for a, s_ in answers[(r["study"], r["subject"])]["results"].items() if s_ < DUP and tuple(LINKS.get(a, ())) == (r["study"], r["subject"])))
two = project["PRJEB24041"] | project["PRJEB28237"]
three = two | project["PRJNA306874"]
for name, runs, M in (("stock", stock, M_GUT), ("stock_without_two_projects", stock - two, M_GUT), ("stock_without_three_projects", stock - three, M_GUT), ("stock_without_four_projects", stock - three - hmp_second, M_GUT),
                      ("index", index, M_GUT), ("labelled", labelled, len(labelled))):
    best = np.array([max((s for a, s in o.items() if a in runs), default=0.0) for o in other])
    k = np.array([sum(s > p for a, s in o.items() if a in runs) if p > 0 else -1 for o, p in zip(other, partner)])
    cand = np.unique(np.concatenate([partner[partner > 0], best]))
    far = np.array([np.mean(best >= c) for c in cand]); frr = np.array([np.mean(partner < c) for c in cand])           # absent persons accepted; donors present rejected (a donor whose other sample is beyond the list counts as rejected)
    t = float(cand[int(np.argmin(np.abs(far - frr)))])                                                                 # equal error rate, as in figure1_brief_communication.py
    out[name] = dict(runs=len(runs), first=int((k == 0).sum()), threshold=round(t, 5), absent=int((best >= t).sum()), other_at_or_above=int((partner >= t).sum()), other_at_or_above_first=int(((partner >= t) & (partner >= best)).sum()),
                     other_at_or_above_first_pct=100 * float(((partner >= t) & (partner >= best)).mean()), runs_above=int(k[k > 0].sum()), line_8e9=line_at_world(k.astype(float), M))
    print(name, out[name])
k_stock = np.array([sum(s > p for a, s in o.items() if a in stock) if p > 0 else -1 for o, p in zip(other, partner)], float)
out["scale_of_the_pools"] = {str(M): line_at_world(k_stock, M) for M in (172075, 220000, 255120, M_GUT, 350000)}          # the same ranks; only the number of human gut metagenomes changes
print(out["scale_of_the_pools"])
json.dump(out, open(OUT, "w"), indent=1)
