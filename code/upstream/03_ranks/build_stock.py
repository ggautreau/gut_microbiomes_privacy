#!/usr/bin/env python3
"""The HGM set: the runs of the index taken as human gut metagenomes (Methods, "HGM set"; Supplementary Table 3). Every rank, every similarity and every pool of the paper is computed on them.

    python3 -I build_stock.py fetch  DATA_DIR LABELLED.txt.gz CMD_META.tsv ENA.tsv                                  # asks the ENA for the runs it lacks, accession by accession (resumes)
    python3 -I build_stock.py select DATA_DIR LABELLED.txt.gz CMD_META.tsv ENA.tsv OUT_STOCK.txt.gz OUT_SUMMARY.json

DATA_DIR      the DATA folder of rows_with_partner.py: results.jsonl (the 5,000 nearest runs of each query), logan_index_ids.txt.gz, biosample_runs.json, curated_sample_biosample.json, runs_hmp2.json, matrices/
LABELLED      the runs labelled human gut or faecal shotgun metagenome in the Logan metadata (172,075 runs; the query that made the list was not kept)
CMD_META      the metadata tables of curatedMetagenomicData in one file (body_site, NCBI_accession)
ENA.tsv       public metadata of the other runs considered for the HGM set: run_accession, scientific_name, tax_id, library_strategy, library_source, host, host_tax_id, isolation_source, study_accession, sample_title

Labels miss human gut metagenomes: Costea, Hansen and Raymond have no run among the labelled runs, because their runs were deposited under Homo sapiens, feces metagenome and gut metagenome. A run enters the HGM set when it is
  1. labelled (LABELLED), or
  2. a run of a curated stool sample (the accessions of curatedMetagenomicData and the runs of the BioSamples of the 7 cohorts of the paper), whatever its label, if the ENA gives a metagenomic library of strategy WGS, WGA or OTHER, or
  3. among the 5,000 nearest runs of at least one of the stool queries, with a metagenomic library of strategy WGS, WGA or OTHER, an organism label among LABELS, no declared non-human host
     and no other part of the body as declared source, or
  4. labelled gut or feces at the ENA (GUT_TAXA, read by taxon), close to a query or not, with the same conditions on the library, the host and the source. Such a run is in no list of neighbours, so it changes no rank, only the size of the pools; a run with no declared host can come from an animal.
The condition on the neighbours is what keeps the generic labels: on the labels alone the same rule would admit about 350,000 runs, among them sputum, skin or soil metagenomes.
Requests hold public run accessions only, with no link to a donor.
"""
import collections, csv, gzip, json, os, re, sys, time, urllib.parse, urllib.request

MODE, DATA, LABELLED, CMD, ENA = sys.argv[1:6]
FIELDS = "run_accession,scientific_name,tax_id,library_strategy,library_source,host,host_tax_id,isolation_source,study_accession,sample_title"
GUT_TAXA = {"human gut metagenome": 408170, "human feces metagenome": 2705415, "feces metagenome": 1861841, "gut metagenome": 749906}
COHORTS = ["Bengtsson-PalmeJ_2015", "CosteaPI_2017", "HMP_2019_ibdmdb", "HansenLBS_2018", "Heitz-BuschartA_2016", "MehtaRS_2018", "RaymondF_2016"]
LABELS = {"human gut metagenome", "human feces metagenome", "gut metagenome", "feces metagenome", "fecal metagenome", "intestinal metagenome", "homo sapiens", "human metagenome", "metagenome"}
STRATEGY = {"WGS", "WGA", "OTHER"}
HOST_HUMAN = {"homo sapiens", "homo", "human", "human gut", "h. sapiens"}
HOST_OTHER = {"mus musculus", "mouse", "porcine", "sus scrofa", "piglet", "pig", "canis lupus familiaris", "dog", "rat", "rattus norvegicus", "chicken", "gallus gallus", "bos taurus", "cow", "macaca mulatta"}
OTHER_SITE = re.compile(r"sputum|conjunctiv|bronchoalv|\bbal\b|balf|\bblood\b|plasma|serum|breast ?milk|\bmilk\b|urine|clean-catch|\bskin\b|forehead|\bnares\b|nasal|nasoph|vagin|cervi|\blung\b|\boral\b|mouth|dental|plaque|saliva|tongue|throat|tonsil|torso|\bface\b|osteolog|\btalus\b|\bbone\b|\bteeth\b|\btooth\b|calculus", re.I)
GUT_WORD = re.compile(r"stool|feces|faeces|fecal|faecal|facel|gut\b|intestin|colon|rect|caec|cecum|ileum|ileal|bowel", re.I)
NO_SOURCE = re.compile(r"^(|missing|not applicable|not collected|not provided|na|n/a|unknown|non|none|metagenome)$", re.I)
ids = lambda p: set(gzip.open(p, "rt").read().split())

labelled, index = ids(LABELLED), ids(f"{DATA}/logan_index_ids.txt.gz")
assert labelled <= index
queries = {}
for fn in ("results.jsonl", "results_pass2.jsonl"):
    if os.path.exists(f"{DATA}/{fn}"):
        for l in open(f"{DATA}/{fn}"):
            try:
                q = json.loads(l)
            except ValueError:
                continue
            if "results" in q and (q["study"], q["subject"]) not in queries:
                queries[(q["study"], q["subject"])] = q
neighbours = set().union(*(set(q["results"]) for q in queries.values()))
assert neighbours <= index

csv.field_size_limit(10**9)
curated = set()
for r in csv.DictReader(open(CMD), delimiter="\t"):
    if r["body_site"] == "stool":
        curated |= {a.strip() for a in (r.get("NCBI_accession") or "").replace(",", ";").split(";") if a.strip() in index}
curated_listed = set(curated); n_curated_listed = len(curated_listed)
stool = {(st, r["sample_id"]) for st in COHORTS for r in csv.DictReader(open(f"{DATA}/matrices/{st}.metadata.tsv"), delimiter="\t") if r["body_site"] == "stool" and r["subject_id"]}
biosample = {bs for k, bs in json.load(open(f"{DATA}/curated_sample_biosample.json")).items() if tuple(k.split("|", 1)) in stool}
cohort_runs = {run for run, v in json.load(open(f"{DATA}/biosample_runs.json")).items() if v["biosample"] in biosample}
cohort_runs |= {run for run, v in json.load(open(f"{DATA}/runs_hmp2.json")).items() if v[0] == "HMP_2019_ibdmdb"}
cohort_runs &= index
curated |= cohort_runs

if MODE == "fetch":
    done = {l.split("\t", 1)[0] for l in open(ENA)} if os.path.exists(ENA) else set()
    todo = sorted((neighbours | curated) - labelled - done)
    print(len(queries), "queries;", len(neighbours), "neighbours;", len(curated), "curated runs;", len(todo), "to ask", flush=True)
    with open(ENA, "a") as fo:
        if not done:
            fo.write(FIELDS.replace(",", "\t") + "\n")
        for i in range(0, len(todo), 400):
            body = urllib.parse.urlencode(dict(result="read_run", query=" OR ".join(f'run_accession="{a}"' for a in todo[i:i + 400]), fields=FIELDS, format="tsv", limit=0)).encode()
            for t in range(6):
                try:
                    with urllib.request.urlopen(urllib.request.Request("https://www.ebi.ac.uk/ena/portal/api/search", body), timeout=180) as f:
                        lines = f.read().decode().split("\n")
                    break
                except OSError as e:
                    print("retry", i, repr(e)[:80], flush=True); time.sleep(5 * (t + 1))
            else:
                raise SystemExit(f"no answer for the batch at {i}")
            for l in lines[1:]:
                if l.strip():
                    fo.write(l.rstrip("\n") + "\n")
            fo.flush(); time.sleep(0.1)
        done = {l.split("\t", 1)[0] for l in open(ENA)}
        for name, tax in GUT_TAXA.items():
            query = f'tax_eq({tax}) AND library_source="METAGENOMIC" AND (library_strategy="WGS" OR library_strategy="WGA" OR library_strategy="OTHER")'
            body = urllib.parse.urlencode(dict(result="read_run", query=query, fields=FIELDS, format="tsv", limit=0)).encode()
            with urllib.request.urlopen(urllib.request.Request("https://www.ebi.ac.uk/ena/portal/api/search", body), timeout=500) as f:
                lines = f.read().decode().split("\n")
            new = [l for l in lines[1:] if l.strip() and l.split("\t", 1)[0] in index and l.split("\t", 1)[0] not in done and l.split("\t", 1)[0] not in labelled]
            fo.write("".join(l.rstrip("\n") + "\n" for l in new)); fo.flush(); done |= {l.split("\t", 1)[0] for l in new}
            print(name, len(new), "new runs of the index", flush=True)
    sys.exit(0)

OUT, SUMMARY = sys.argv[6:8]
ena = {r["run_accession"]: r for r in csv.DictReader(open(ENA), delimiter="\t", quoting=csv.QUOTE_NONE)}


def library_ok(r):
    return r["library_source"] == "METAGENOMIC" and r["library_strategy"] in STRATEGY


def neighbour_reason(r):
    """None when a neighbour of another label enters the HGM set, else the reason it does not."""
    host = r["host"].strip().lower()
    if not library_ok(r):
        return "another library"
    if r["host_tax_id"] not in ("", "9606") or host in HOST_OTHER:
        return "a non-human host"
    if r["scientific_name"].lower() not in LABELS:
        return "another organism label"
    if OTHER_SITE.search(r["isolation_source"]) and not GUT_WORD.search(r["isolation_source"]):
        return "another part of the body as source"
    return None


cand = (neighbours | curated) - labelled
unknown = cand - set(ena)
add_cur = {a for a in curated - labelled if a in ena and library_ok(ena[a])}
assert not cohort_runs - labelled - add_cur, "every indexed run of the 7 cohorts is in the HGM set"
reasons = {a: neighbour_reason(ena[a]) for a in neighbours - labelled - add_cur if a in ena}
add_nei = {a for a, why in reasons.items() if why is None}
add_gut = {a for a, r in ena.items() if a in index and a not in labelled and a not in add_cur and a not in add_nei and r["scientific_name"].lower() in GUT_TAXA and neighbour_reason(r) is None}
stock = labelled | add_cur | add_nei | add_gut
assert stock <= index and len(stock) == len(labelled) + len(add_cur) + len(add_nei) + len(add_gut)
with open(OUT, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as fo:      # the same bytes whatever the name and the date
    fo.write(("\n".join(sorted(stock)) + "\n").encode())


def label_of(a):
    r = ena[a]; name = r["scientific_name"].lower(); host = r["host"].strip().lower()
    if name == "metagenome":
        return "metagenome, human host declared" if r["host_tax_id"] == "9606" or host in HOST_HUMAN else "metagenome, no host declared"
    if name in ("human gut metagenome", "human feces metagenome"):
        return "human gut or feces metagenome, library other than WGS" if r["library_strategy"] != "WGS" else "human gut or feces metagenome, WGS"
    return {"gut metagenome": "gut or feces metagenome", "feces metagenome": "gut or feces metagenome", "fecal metagenome": "gut or feces metagenome", "intestinal metagenome": "gut or feces metagenome",
            "homo sapiens": "Homo sapiens", "human metagenome": "human metagenome"}.get(name, "other labels")


declared = [a for a in add_nei if not NO_SOURCE.match(ena[a]["isolation_source"].strip())]
summary = dict(queries=len(queries), neighbours=len(neighbours), neighbours_outside_labelled=len(neighbours - labelled), not_known_to_the_ena=len(unknown), stock=len(stock), labelled=len(labelled), curated_added=len(add_cur), neighbours_added=len(add_nei), gut_labels_added=len(add_gut),
               neighbours_left_out=dict(collections.Counter(w for w in reasons.values() if w)), curated_added_by_label=dict(collections.Counter(label_of(a) for a in add_cur)), neighbours_added_by_label=dict(collections.Counter(label_of(a) for a in add_nei)), gut_labels_added_by_label=dict(collections.Counter(label_of(a) for a in add_gut)), gut_labels_added_no_host_gut_or_feces=sum(label_of(a) == "gut or feces metagenome" and ena[a]["host_tax_id"] != "9606" and ena[a]["host"].strip().lower() not in HOST_HUMAN for a in add_gut),
               gut_labels_added_host_none=sum(label_of(a) == "gut or feces metagenome" and ena[a]["host_tax_id"] != "9606" and ena[a]["host"].strip().lower() not in HOST_HUMAN and ena[a]["host"].strip().lower() in ("", "missing") for a in add_gut),
               gut_labels_added_host_text=sum(label_of(a) == "gut or feces metagenome" and ena[a]["host_tax_id"] != "9606" and ena[a]["host"].strip().lower() not in HOST_HUMAN and ena[a]["host"].strip().lower() not in ("", "missing") for a in add_gut),
               added_projects=len({ena[a]["study_accession"] for a in add_cur | add_nei | add_gut}), neighbours_source_declared=len(declared), neighbours_source_gut=sum(bool(GUT_WORD.search(ena[a]["isolation_source"])) for a in declared),
               curated_stool_listed=n_curated_listed, curated_stool_listed_labelled=len(curated_listed & labelled), cohort_runs_indexed=len(cohort_runs), cohort_runs_labelled=len(cohort_runs & labelled))
json.dump(summary, open(SUMMARY, "w"), indent=1)
print(json.dumps(summary, indent=1))
