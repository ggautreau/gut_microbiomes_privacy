#!/usr/bin/env python3
"""Re-runs the steps of code/upstream/ from the files of local/, in a scratch copy of the reproduction package, and compares what they write with the files of the package (local/, data/, aggregates/, reference/).

    python3 -I rerun_upstream.py [--fast] [--extra DIR] [--table OUT.md] [SCRATCH_DIR]          # from the root of the package; about 5 h 20 min on 32 cores, 16 min with --fast

SCRATCH_DIR   the copy is made there (default: a temporary folder, removed at the end)
--fast        leaves out extrapolation_tables.py, extrapolation_bootstrap.py and evt_model.py (the slow extrapolation steps)
--extra DIR   a folder with the table that is not in the package, hmp2_metadata_2018-08-20.csv (https://ibdmdb.org); the step that reads it is then run too
--table       writes the table of the steps (provenance/VERIFIED_LOCALLY.md holds the one of the paper)
The steps run with the network cut when the system allows it (unshare -rn); the first line printed says so. The steps that query the service or the ENA, run CroCoDeEL, sketch runs or read the hypervector database of the service are not repeated
(code/upstream/README.md lists them).
"""
import csv, filecmp, hashlib, json, math, os, shutil, subprocess, sys, tempfile, time
import numpy as np

args = sys.argv[1:]
FAST = "--fast" in args
EXTRA = args[args.index("--extra") + 1] if "--extra" in args else None
OUTMD = args[args.index("--table") + 1] if "--table" in args else None
rest = [a for i, a in enumerate(args) if a not in ("--fast", "--extra", "--table") and (i == 0 or args[i - 1] not in ("--extra", "--table"))]
ARCH = os.getcwd()
assert os.path.isdir(f"{ARCH}/local/adultes") and os.path.isdir(f"{ARCH}/code/upstream"), "run from the root of the package"
TEST = os.path.abspath(rest[0]) if rest else tempfile.mkdtemp(prefix="rerun_upstream_")
KEEP = bool(rest)
EXTRA = os.path.abspath(EXTRA) if EXTRA else None
OUTMD = os.path.abspath(OUTMD) if OUTMD else None

os.environ["OPENBLAS_NUM_THREADS"] = os.environ["OMP_NUM_THREADS"] = "1"                 # extrapolation_tables.py and evt_model.py run in several processes
NET = ["unshare", "-rn"] if shutil.which("unshare") and subprocess.run(["unshare", "-rn", "true"], capture_output=True).returncode == 0 else []
print("network cut for every step (unshare -rn)" if NET else "the network could not be cut on this system (no unshare -rn); no step needs it", flush=True)
shutil.rmtree(TEST, ignore_errors=True)
shutil.copytree(ARCH, TEST, ignore=shutil.ignore_patterns("output", ".git", "__pycache__", "w", "adultes_data_v2"))
os.makedirs(f"{TEST}/w")
p = subprocess.run(NET + ["bash", "code/upstream/prepare_local.sh"], cwd=TEST, capture_output=True, text=True)
assert p.returncode == 0, "prepare_local.sh failed: " + p.stderr[-600:]
manifest = {r["path_in_local_tree"]: r for r in csv.DictReader(open(f"{TEST}/provenance/MANIFEST_upstream.tsv"), delimiter="\t")}
results = []


def run(a, step):
    t = time.time()
    p = subprocess.run(NET + ["python3", "-I"] + a, cwd=TEST, capture_output=True, text=True)
    if p.returncode:
        print(step, "FAILED:", p.stderr[-900:], flush=True)
    return p.returncode == 0, time.time() - t


ABS_TOL, REL_TOL = 1e-3, 1e-9            # as check_reproduction.py: 3 decimals, or 1e-9 (relative) for the very large numbers
NEAR = [False]                            # set when a JSON file of the step is equal only to 3 decimals (another Python or numpy moves the optimum of a fit); record() reads it and clears it


def close(a, b):
    """As diff_json of check_reproduction.py: 2 when equal, 1 when the only differences are floats equal to 3 decimals (or to REL_TOL relative), 0 otherwise."""
    if isinstance(a, dict) and isinstance(b, dict):
        return 0 if set(a) != set(b) else min((close(a[k], b[k]) for k in a), default=2)
    if isinstance(a, list) and isinstance(b, list):
        return 0 if len(a) != len(b) else min((close(x, y) for x, y in zip(a, b)), default=2)
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
        if a == b:
            return 2
        return 1 if (isinstance(a, float) or isinstance(b, float)) and math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL) else 0
    return 2 if a == b else 0


def same_json(a, b):
    c = close(json.load(open(a)), json.load(open(b)))
    NEAR[0] = NEAR[0] or c == 1
    return c > 0


def same_npz(a_, b_):
    a = np.load(a_, allow_pickle=True); b = np.load(b_, allow_pickle=True)
    return sorted(a.files) == sorted(b.files) and all((len(a[k]) == len(b[k]) and all(np.array_equal(x, y) for x, y in zip(a[k], b[k]))) if k == "imp_nar" else np.array_equal(a[k], b[k]) for k in a.files)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def record(step, cmd, ok, same, secs, note=""):
    near, NEAR[0] = NEAR[0], False
    results.append((step, cmd, ("equal to 3 decimals" if near else "identical") if ok and same else ("RUN FAILED" if not ok else "DIFFERENT"), f"{secs:.0f} s", note))
    print(f"{results[-1][2]:10s} {step} ({secs:.0f} s) {note}", flush=True)


def skipped(step, cmd, why):
    results.append((step, cmd, "not run", "", why))
    print(f"{'not run':10s} {step}: {why}", flush=True)


U = "code/upstream"; T = TEST; LA = f"{T}/local/adultes"; HC = "local/humancheck/all558"
# 0. donor pairs from the curated metadata and the ENA table of HMP IBD (the steps before the queries that need no network)
shutil.copy(f"{T}/{U}/02_pairs_and_queries/adult_cohorts.txt", f"{T}/adult_cohorts.txt")
ok, t = run([f"{U}/02_pairs_and_queries/build_pairs.py"], "pairs")
record("build_pairs.py", "local/crocodeel_valles/matrices -> pairs.json (first and last stool sample of each donor of the 14 adult cohorts)", ok, ok and same_json(f"{T}/pairs.json", f"{LA}/pairs.json"), t, "compared with local/adultes/pairs.json")
shutil.copy(f"{LA}/ibdmdb_runs_ena.tsv", f"{T}/ibdmdb_runs_ena.tsv")
shutil.copy(f"{LA}/pairs_hmp2_v1_backup.json", f"{T}/pairs_hmp2.json")                              # build_hmp2_v2.py updates the pairs that build_hmp2.py (first version, ENA) wrote
ok, t = run([f"{U}/02_pairs_and_queries/build_hmp2_v2.py"], "hmp2")
record("build_hmp2_v2.py", "ibdmdb_runs_ena.tsv (ENA) and the pairs of build_hmp2.py -> pairs_hmp2.json, runs_hmp2.json, runs_hmp2_dates.json", ok,
       ok and same_json(f"{T}/pairs_hmp2.json", f"{LA}/pairs_hmp2_v2_backup.json") and all(same_json(f"{T}/{f}", f"{LA}/{f}") for f in ("runs_hmp2.json", "runs_hmp2_dates.json")), t,
       "compared with the files of local/adultes/ (the pairs before the rebuild of 8 November 2026)")
os.makedirs(f"{T}/w/v3data"); os.makedirs(f"{T}/w/v3out")
for f in os.listdir(f"{T}/local/adultes_data_v2"):
    if f != "pairs_hmp2.json": os.symlink(f"{T}/local/adultes_data_v2/{f}", f"{T}/w/v3data/{f}")
shutil.copy(f"{T}/pairs_hmp2.json", f"{T}/w/v3data/pairs_hmp2.json")                               # the pairs that build_hmp2_v2.py just wrote
ok, t = run([f"{U}/02_pairs_and_queries/build_hmp2_v3.py", "local/rows_v3.json", "w/v3data", "local/adultes/ibdmdb_runs_ena.tsv", "w/v3out"], "hmp3")
record("build_hmp2_v3.py", "the pairs of build_hmp2_v2.py, the curated metadata and ibdmdb_runs_ena.tsv -> pairs_hmp2.json (first and last sample of the 48 adult donors of HMP IBD by the curated link), requery_todo.json (the 14 donors whose query changes)", ok,
       ok and same_json(f"{T}/w/v3out/pairs_hmp2.json", f"{LA}/pairs_hmp2.json") and same_json(f"{T}/w/v3out/requery_todo.json", f"{T}/local/work/requery_todo_hmp_ibd_v3.json"), t, "compared with the files of local/")
ok, t = run([f"{U}/02_pairs_and_queries/merge_results.py", "local/adultes", "w/merged"], "merge")
record("merge_results.py", "first answers + the 198 + 14 queries with the largest run -> adultes_data_v2/results.jsonl", ok, ok and sha(f"{T}/w/merged/results.jsonl") == manifest["local/adultes_data_v2/results.jsonl"]["sha256"], t, "the sha256 of the file of the paper (MANIFEST_upstream.tsv)")
ok, t = run([f"{U}/01_cohort_tables/cohort_age.py", "local/crocodeel_valles/matrices", "w/cohort_age.csv"], "age")
record("cohort_age.py", "local/crocodeel_valles/matrices -> aggregates/cohort_age.csv (age label of the 30 longitudinal cohorts)", ok, ok and filecmp.cmp(f"{T}/w/cohort_age.csv", f"{T}/aggregates/cohort_age.csv", shallow=False), t, "byte for byte")
# 1. the HGM set, and the ranks from the answers of the service
ok, t = run([f"{U}/03_ranks/build_stock.py", "select", "local/adultes_data_v2", "local/adultes/logan_gut_wgs_labelled_172075_ids.txt.gz", "local/cmd_meta/cmd_meta.tsv", "local/adultes/ena_stock_candidates.tsv", "w/human_gut.txt.gz", "w/stock_summary.json"], "stock")
record("build_stock.py select", "-> reference/human_gut_stock_291881_runs.txt.gz, aggregates/stock_summary.json", ok, ok and filecmp.cmp(f"{T}/w/human_gut.txt.gz", f"{T}/reference/human_gut_stock_291881_runs.txt.gz", shallow=False) and same_json(f"{T}/w/stock_summary.json", f"{T}/aggregates/stock_summary.json"), t, "from the table of the ENA of local/adultes/")
ok, t = run([f"{U}/03_ranks/mehta_siblings.py", "map", "local/adultes_data_v2", "local/adultes/ena_PRJNA354235_runs.tsv", "w/mehta_sibling_runs.json"], "mehta")
record("mehta_siblings.py map", "-> mehta_sibling_runs.json (runs of the Mehta project that the curated BioSamples do not name, by ENA alias)", ok, ok and filecmp.cmp(f"{T}/w/mehta_sibling_runs.json", f"{LA}/mehta_sibling_runs.json", shallow=False), t, "byte for byte")
ok, t = run([f"{U}/03_ranks/rows_with_partner.py", "local/adultes_data_v2", "w/rows_v2_full_without_links.json"], "rows0")
record("rows_with_partner.py", "local/adultes_data_v2 -> rows_v2_full_without_links.json (without --hmp)", ok, ok and same_json(f"{T}/w/rows_v2_full_without_links.json", f"{T}/local/rows_v2_full_without_links.json"), t, "831 donors")
HMP2 = "hmp2_metadata_2018-08-20.csv"
if EXTRA and os.path.exists(f"{EXTRA}/{HMP2}"):
    shutil.copy(f"{EXTRA}/{HMP2}", f"{T}/local/{HMP2}")
    ok, t = run([f"{U}/03_ranks/hmp_sibling_links.py", "w/rows_v2_full_without_links.json", f"local/{HMP2}", "local/adultes/ena_PRJNA389280_runs.tsv", "local/adultes/ena_stock_candidates.tsv", "w/hmp_sibling_links.json"], "links")
    record("hmp_sibling_links.py", "the public table of the HMP consortium -> hmp_sibling_links.json (runs of PRJNA306874 and PRJNA389280 given to a donor)", ok, ok and same_json(f"{T}/w/hmp_sibling_links.json", f"{T}/local/hmp_sibling_links.json"), t, f"with {HMP2}, which is not in the package")
else:
    skipped("hmp_sibling_links.py", "the public table of the HMP consortium -> hmp_sibling_links.json", f"reads {HMP2} (https://ibdmdb.org), which is not in the package; the steps below read local/hmp_sibling_links.json")
    shutil.copy(f"{T}/local/hmp_sibling_links.json", f"{T}/w/hmp_sibling_links.json")
ok, t = run([f"{U}/03_ranks/rows_with_partner.py", "local/adultes_data_v2", "w/rows_v2_full.json", "--hmp", "w/hmp_sibling_links.json"], "rows")
record("rows_with_partner.py --hmp", "-> rows_v2_full.json", ok, ok and same_json(f"{T}/w/rows_v2_full.json", f"{T}/local/rows_v2_full.json"), t, "831 donors")
ok, t = run([f"{U}/03_ranks/make_rows_v4.py", "w/rows_v2_full.json", "w/rows_v4_497.json"], "v4")
record("make_rows_v4.py", "-> rows_v4_497.json (the 497 donors of the 7 cohorts kept, before the two rules)", ok, ok and same_json(f"{T}/w/rows_v4_497.json", f"{T}/local/rows_v4_497.json"), t)
ok, t = run([f"{U}/03_ranks/make_rows_v4.py", "w/rows_v2_full_without_links.json", "w/rows_v3.json", "--keep-vincent"], "v3")
record("make_rows_v4.py --keep-vincent", "-> rows_v3.json (558 donors, 8 cohorts, for the human-DNA check, without --hmp)", ok, ok and same_json(f"{T}/w/rows_v3.json", f"{T}/local/rows_v3.json"), t)
ok, t = run([f"{U}/03_ranks/donors_left_out.py", "w/rows_v2_full_without_links.json", "w/donors_left_out.json"], "leftout")
record("donors_left_out.py", "-> donors_left_out.json (the 183 donors of the 4 cohorts left out)", ok, ok and same_json(f"{T}/w/donors_left_out.json", f"{T}/local/donors_left_out.json"), t)
ok, t = run([f"{U}/05_human_dna/all558_join.py", f"{HC}/final_table.csv", "w/rows_v3.json", "w/join_final"], "join")
record("all558_join.py", "-> join_final_table.csv, join_final_summary.json (the 558 donors: share of human DNA and rank)", ok, ok and filecmp.cmp(f"{T}/w/join_final_table.csv", f"{T}/{HC}/join_final_table.csv", shallow=False) and same_json(f"{T}/w/join_final_summary.json", f"{T}/{HC}/join_final_summary.json"), t, "from the table of the human-DNA check")
ok, t = run([f"{U}/05_human_dna/all558_pairs.py", "local", "w/pairs_service.json"], "pairs558")
record("all558_pairs.py", "-> pairs_service.json (the 558 pairs of the human-DNA check: query run, run of the other sample)", ok, ok and same_json(f"{T}/w/pairs_service.json", f"{T}/{HC}/pairs_service.json"), t, "also identical to the list of runs sketched")
ok, t = run([f"{U}/03_ranks/extract_lists.py", "local/adultes_data_v2", "w/rows_v4_497.json", "w/lists_v4_497.npz", "--hmp", "w/hmp_sibling_links.json"], "lists")
record("extract_lists.py --hmp", "-> lists_v4_497.npz (similarities of the runs of other people of each of these 497 donors)", ok, ok and same_npz(f"{T}/w/lists_v4_497.npz", f"{T}/local/lists_v4_497.npz"), t)
ok, t = run([f"{U}/03_ranks/hmp_curated_ids.py", "w/rows_v4_497.json", "local/adultes", "local/crocodeel_valles/matrices", "w/rows_v4_497_curated.json"], "curated")
record("hmp_curated_ids.py", "-> rows_v4_497_curated.json (the curated sample identifiers of HMP IBD)", ok, ok and same_json(f"{T}/w/rows_v4_497_curated.json", f"{T}/local/rows_v4_497_curated.json"), t)
ok, t = run([f"{U}/03_ranks/exclude_donors.py", "w/rows_v4_497.json", "w/rows_v4_497_curated.json", "w/lists_v4_497.npz", "local/crocodeel_valles", "w"], "exclude")
record("exclude_donors.py", "-> rows_v4.json, rows_v4_curated.json, lists_v4.npz (the 483 donors) and donors_excluded.json (the 14 donors removed)", ok, ok and all(same_json(f"{T}/w/{f}", f"{T}/local/{f}") for f in ("rows_v4.json", "rows_v4_curated.json", "donors_excluded.json")) and same_npz(f"{T}/w/lists_v4.npz", f"{T}/local/lists_v4.npz"), t)
ok, t = run([f"{U}/03_ranks/make_open_set_scores.py", "w/lists_v4_497.npz", "w/open_set_scores_497.csv"], "open497")
record("make_open_set_scores.py (the 497 donors)", "-> open_set_scores_497.csv (for costea_same_people.py)", ok, ok and filecmp.cmp(f"{T}/w/open_set_scores_497.csv", f"{T}/local/open_set_scores_497.csv", shallow=False), t, "byte for byte")
ok, t = run([f"{U}/03_ranks/make_open_set_scores.py", "w/lists_v4.npz", "w/open_set_scores.csv"], "open")
record("make_open_set_scores.py", "-> data/open_set_scores.csv (panel H of Figure 1)", ok, ok and filecmp.cmp(f"{T}/w/open_set_scores.csv", f"{T}/data/open_set_scores.csv", shallow=False), t, "byte for byte")
ok, t = run([f"{U}/02_pairs_and_queries/query_sizes.py", "local/adultes_data_v2/results.jsonl", "w/query_gbp.json", "local/adultes/ena2_all.tsv", "local/adultes/ena_query_runs_bases.tsv"], "gbp")
g1 = json.load(open(f"{T}/w/query_gbp.json")) if ok else {}; g2 = json.load(open(f"{T}/local/query_gbp.json"))
record("query_sizes.py", "-> query_gbp.json (from the two tables of the ENA of local/adultes/)", ok, ok and all(g1.get(k) == g2[k] for k in g2) and set(g2) <= set(g1), t, f"the {len(g2)} donors of local/query_gbp.json; no request to the ENA")
ok, t = run([f"{U}/03_ranks/make_adult_ranks.py", "w/rows_v4.json", "local/query_gbp.json", "w/adult_ranks.csv"], "ranks")
record("make_adult_ranks.py", "-> data/adult_ranks.csv (the table of Figure 1)", ok, ok and filecmp.cmp(f"{T}/w/adult_ranks.csv", f"{T}/data/adult_ranks.csv", shallow=False), t, "byte for byte")
if FAST:
    skipped("extrapolation_tables.py", "-> aggregates/extrapolation_summary.json", "left out by --fast")
    skipped("extrapolation_bootstrap.py", "-> aggregates/extrapolation_bootstrap.json", "left out by --fast")
else:
    ok, t = run([f"{U}/03_ranks/extrapolation_tables.py", "w/lists_v4.npz", "w/extrapolation_summary.json"], "extrap")
    record("extrapolation_tables.py", "-> aggregates/extrapolation_summary.json", ok, ok and same_json(f"{T}/w/extrapolation_summary.json", f"{T}/aggregates/extrapolation_summary.json"), t, "400 thinnings at each fraction")
    ok, t = run([f"{U}/03_ranks/extrapolation_bootstrap.py", "w/adult_ranks.csv", "w/extrapolation_bootstrap.json"], "boot")
    record("extrapolation_bootstrap.py", "-> aggregates/extrapolation_bootstrap.json", ok, ok and same_json(f"{T}/w/extrapolation_bootstrap.json", f"{T}/aggregates/extrapolation_bootstrap.json"), t, "adaptive doubling bootstrap")
if FAST:
    skipped("evt_model.py", "-> aggregates/evt_summary.json", "left out by --fast")
else:
    ok, t = run([f"{U}/03_ranks/evt_model.py", "local/evt_similarities.npz", "w/evt_summary.json", "400"], "evt")
    record("evt_model.py", "-> aggregates/evt_summary.json (the extreme-value model of the similarity tail, Supplementary Table 18)", ok, ok and same_json(f"{T}/w/evt_summary.json", f"{T}/aggregates/evt_summary.json"), t, "400 thinnings at each fraction")
ok, t = run([f"{U}/03_ranks/evt_resampling.py", "local/evt_similarities.npz", "aggregates/evt_summary.json", "w/evt_resampling.json"], "evt band")      # one process, a few seconds: also run with --fast
record("evt_resampling.py", "-> aggregates/evt_resampling.json (the resampling band of the extreme-value model, Supplementary Fig. 6)", ok, ok and same_json(f"{T}/w/evt_resampling.json", f"{T}/aggregates/evt_resampling.json"), t, "finds the curve and the intervals of evt_model.py again before it writes")
ok, t = run([f"{U}/03_ranks/evt_removed_donors.py", "local/evt_similarities_497.npz", "local/rows_v4_497.json", "local/rows_v4.json", "local/evt_similarities.npz", "aggregates/evt_summary.json", "w/evt_removed_donors.json"], "evt 497")      # one process, a few seconds: also run with --fast
record("evt_removed_donors.py", "-> aggregates/evt_removed_donors.json (the extreme-value model with the donors that the two rules remove, last row of Supplementary Table 17)", ok, ok and same_json(f"{T}/w/evt_removed_donors.json", f"{T}/aggregates/evt_removed_donors.json"), t, "finds the model of evt_model.py again on the 483 donors before it writes")
ok, t = run([f"{U}/03_ranks/stock_sensitivity.py", "local/adultes_data_v2", "local/adultes/logan_gut_wgs_labelled_172075_ids.txt.gz", "local/adultes/ena_stock_candidates.tsv", "w/rows_v4.json", "w/stock_sensitivity.json", "--hmp", "w/hmp_sibling_links.json"], "sensitivity")
record("stock_sensitivity.py", "-> aggregates/stock_sensitivity.json", ok, ok and same_json(f"{T}/w/stock_sensitivity.json", f"{T}/aggregates/stock_sensitivity.json"), t)
# 2. aggregates of the supplement
ok, t = run([f"{U}/04_aggregates/cohort_selection.py", "w/rows_v2_full.json", "w/rows_v4.json", "w/cohort_selection.json", "w/donors_excluded.json", "w/rows_v2_full_without_links.json"], "sel")
record("cohort_selection.py", "-> aggregates/cohort_selection.json", ok, ok and same_json(f"{T}/w/cohort_selection.json", f"{T}/aggregates/cohort_selection.json"), t)
ok, t = run([f"{U}/04_aggregates/pehrsson_lanes.py", "w/pehrsson_summary.json", "w/rows_v2_full.json", "local/adultes_data_v2"], "peh")
record("pehrsson_lanes.py", "-> aggregates/pehrsson_summary.json", ok, ok and same_json(f"{T}/w/pehrsson_summary.json", f"{T}/aggregates/pehrsson_summary.json"), t)
ok, t = run([f"{U}/04_aggregates/nielsen_li_overlap.py", "local/adultes_data_v2", "w/rows_v2_full.json", "w/nielsen_li_overlap.json"], "nlo")
record("nielsen_li_overlap.py", "-> aggregates/nielsen_li_overlap.json", ok, ok and same_json(f"{T}/w/nielsen_li_overlap.json", f"{T}/aggregates/nielsen_li_overlap.json"), t)
ok, t = run([f"{U}/04_aggregates/nielsen_li_links.py", "w/nielsen_li_links.json", "w/rows_v4.json", "w/nielsen_li_overlap.json"], "nll")
record("nielsen_li_links.py", "-> aggregates/nielsen_li_links.json", ok, ok and same_json(f"{T}/w/nielsen_li_links.json", f"{T}/aggregates/nielsen_li_links.json"), t)
ok, t = run([f"{U}/04_aggregates/raymond_antibiotic.py", "w/raymond_antibiotic.json", "w/rows_v4.json", "local/cmd_meta/cmd_meta.tsv"], "ray")
record("raymond_antibiotic.py", "-> aggregates/raymond_antibiotic.json (panel G of Figure 1)", ok, ok and same_json(f"{T}/w/raymond_antibiotic.json", f"{T}/aggregates/raymond_antibiotic.json"), t)
ok, t = run([f"{U}/04_aggregates/costea_same_people.py", "w/costea_same_people.json", "w/rows_v4_497.json", "w/open_set_scores_497.csv"], "cst")
record("costea_same_people.py", "-> aggregates/costea_same_people.json", ok, ok and same_json(f"{T}/w/costea_same_people.json", f"{T}/aggregates/costea_same_people.json"), t)
ok, t = run([f"{U}/05_human_dna/human_dna_exclusion.py", f"{HC}/join_final_table.csv", f"{HC}/pool_final_ranks.csv", f"{HC}/join_final_summary.json", f"{HC}/pool_final_summary.json", "w/human_dna_exclusion.json"], "hum")
record("human_dna_exclusion.py", "-> aggregates/human_dna_exclusion.json", ok, ok and same_json(f"{T}/w/human_dna_exclusion.json", f"{T}/aggregates/human_dna_exclusion.json"), t, "from the tables of the human-DNA check")
# 3. inputs of reproduce.sh
os.makedirs(f"{T}/w/inputs")
ok, t = run([f"{U}/06_inputs/make_inputs.py", "w/rows_v4_497_curated.json", "local/query_gbp.json", "local/crocodeel_valles", "local/adultes", "local/crocodeel_nielsen_li/nielsen_li.conta.tsv",
             f"{HC}/join_final_table.csv", f"{HC}/pairs_service.json", f"{HC}/run_projects.json", f"{HC}/stat_per_run.json", "w/inputs", "w/donors_left_out.json", "w/donors_excluded.json"], "inputs")
files = ["donors.json", "donors_left_out.json", "donors_excluded.json", "query_gbp.json", "accessions.tsv", "human_dna/donors_human_share.csv"] + [f"crocodeel/{f}" for f in sorted(os.listdir(f"{T}/data/crocodeel"))]
eq = ok and all(filecmp.cmp(f"{T}/w/inputs/{f}", f"{T}/data/{f}", shallow=False) for f in files)
record("make_inputs.py", f"-> data/ ({len(files)} files)", ok, eq, t, "byte for byte")

if OUTMD:
    with open(OUTMD, "w") as fh:
        fh.write("| Script | Writes | Result against the file of the package | Time | Note |\n|---|---|---|---|---|\n")
        for step, cmd, res, secs, note in results:
            fh.write(f"| `{step}` | {cmd} | **{res}** | {secs} | {note} |\n")
ran = [r for r in results if r[2] != "not run"]; bad = [r for r in ran if r[2] not in ("identical", "equal to 3 decimals")]; near = sum(r[2] == "equal to 3 decimals" for r in ran)
print("\n%d of %d steps identical%s; %d not run" % (len(ran) - len(bad), len(ran), f" ({near} of them equal to 3 decimals)" if near else "", len(results) - len(ran)))
if not KEEP:
    shutil.rmtree(TEST, ignore_errors=True)
sys.exit(1 if bad else 0)
