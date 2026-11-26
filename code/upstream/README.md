# `code/upstream/`: from public data to the tables of the reproduction package

These scripts turn public data (curatedMetagenomicData, the ENA, the public similarity service) into the tables that `reproduce.sh` reads. They read the files of `local/`: the tables of curatedMetagenomicData and of the ENA,
the lists of the 5,000 nearest runs of every query, and the tables that each step writes for the next. Two tables are not in the package (see `../../README.md`, "What is not here"). The scripts are
the ones that produced the paper's tables; the only edits are the replacement of absolute paths by paths under `./local/`, the correction of some docstrings, and one line of `extract_all.R`,
which writes the column `infant_age` of the 15 studies with family records as their tables have it.

Run every command from the root of the package, after `bash code/upstream/prepare_local.sh`, which unpacks the six files that the package holds compressed and builds `local/adultes_data_v2/`.
`python3 -I rerun_upstream.py` re-runs, in a scratch copy and without network, the steps marked verified below, and compares what they write with the files of the package.
`provenance/MANIFEST_upstream.tsv` gives the size and sha256 of each file of `local/`, and `provenance/VERIFIED_LOCALLY.md` the steps re-run for the paper and what they gave.

## The `./local/` tree

```
local/
  sm.rda                                     sampleMetadata of curatedMetagenomicData (R: save(sampleMetadata, file = "local/sm.rda"))
  hmp2_metadata_2018-08-20.csv               public metadata table of the HMP consortium; not in the package (https://ibdmdb.org)
  hmp_sibling_links.json                     output of hmp_sibling_links.py: the runs that this table gives to a donor of HMP IBD
  crocodeel_valles/
    matrices/<study>.metadata.tsv            93 studies: sample, subject, body site, days from first collection, reads, run accessions, age ...
    matrices/<study>.species_relab.tsv       species-level relative abundances of curatedMetagenomicData, for the 11 cohorts of the 680 donors (extract_all.R writes the 93)
    results/<study>.conta.tsv                CroCoDeEL 1.3.1 events, default settings, for the same 11 cohorts
    reproduce/cohort_age.csv                 age label of the 30 longitudinal cohorts (aggregates; also in aggregates/)
    pehrsson/pehrsson_run_alias.json         run -> SRA run alias, PRJNA300541
  cmd_meta/cmd_meta.tsv                      the 93 metadata tables of curatedMetagenomicData in one file
  adultes/                                   pairs.json, pairs_hmp2.json and the pairs before the rebuilds (pairs_hmp2_v1_backup.json, pairs_hmp2_v2_backup.json), runs_hmp2.json, runs_hmp2_dates.json, biosample_runs.json,
                                             curated_sample_biosample.json, the lists of runs logan_index_ids.txt.gz, logan_gut_wgs_labelled_172075_ids.txt.gz, logan_gut_wgs_broad_ids.txt.gz and logan_gut_wgs_narrow_ids.txt.gz
                                             (the HGM set; the labelled list and the HGM set are also in reference/), the tables of the ENA (ena_stock_candidates.tsv, ibdmdb_runs_ena.tsv, ena_PRJNA354235_runs.tsv,
                                             ena_PRJNA389280_runs.tsv, ena_query_runs_bases.tsv, ena2_all.tsv), mehta_sibling_runs.json, and the first answers of the service (results.jsonl, results_pass2.jsonl)
  work/                                      results_largest.jsonl (the answers of the queries with the largest run), requery_todo.json, requery_todo_hmp_ibd_v3.json
  adultes_data_v2/                           the DATA folder of rows_with_partner.py, built by prepare_local.sh: the files of adultes/, matrices/ (a copy of crocodeel_valles/matrices), and
                                             results.jsonl (final queries: the 5,000 nearest runs of each donor's query run, merged by merge_results.py)
  rows_v2_full_without_links.json            output of rows_with_partner.py
  rows_v2_full.json                          output of rows_with_partner.py --hmp
  rows_v3.json  rows_v4_497.json             output of make_rows_v4.py --keep-vincent and of make_rows_v4.py
  rows_v4_497_curated.json                   output of hmp_curated_ids.py
  lists_v4_497.npz  open_set_scores_497.csv  output of extract_lists.py and of make_open_set_scores.py for the 497 donors
  rows_v4.json  rows_v4_curated.json  lists_v4.npz  donors_excluded.json      output of exclude_donors.py
  evt_similarities.npz  evt_similarities_497.npz  output of evt_similarities.py for the 483 donors of the analysis and for the 497 donors of the cohorts kept
  donors_left_out.json  query_gbp.json       output of donors_left_out.py and of query_sizes.py
  crocodeel_nielsen_li/nielsen_li.conta.tsv  CroCoDeEL on the species tables of Nielsen and Li merged
  nielsen_li_validation/                     ena_PRJEB1220.tsv, ena_PRJEB5224.tsv (ENA), shared_donors_runs.tsv, reverse/prep.json, reverse/answers.jsonl
  costea_same_people/                        ena_PRJEB*.tsv, ena_studies.tsv (ENA), biosamples.json (BioSamples), pairs.tsv, reads_compare_first_tile.json, trim_filter_first_tile.json, relabel.json
  humancheck/all558/                         outputs of the human-DNA check (final_table.csv, join_final_table.csv, pool_final_ranks.csv, *_summary.json), all558_runs.json, pairs_service.json (all558_pairs.py),
                                             run_projects.json (run_projects.py) and stat_per_run.json (the counts of the NCBI taxonomic analysis, STAT, of the 1,116 runs of the check; R. Faure, 8 November 2026)
```

## Steps

"Verified" means: re-run for this package from the files of `local/` (`rerun_upstream.py`), and the file written is identical to the one of the package.
"Not re-run" means: the step needs the network, a long computation or a tool that is not in the package; the script is the one that was used.

### 1. Cohort tables (`01_cohort_tables/`)

| Script | Command | Writes | Status |
|---|---|---|---|
| `extract_all.R` | `Rscript code/upstream/01_cohort_tables/extract_all.R` | `local/crocodeel_valles/matrices/` (species table and metadata of each study) | not re-run for the package (R, ExperimentHub download); run on 15 November 2026 from the download cache, into an empty folder, it rewrote the 186 tables of 1 October 2026 byte for byte |
| `run_crocodeel.sh` | `bash code/upstream/01_cohort_tables/run_crocodeel.sh` | `local/crocodeel_valles/results/<study>.conta.tsv` (CroCoDeEL 1.3.1, `search_conta`, default settings) | not re-run (hours; results of 1 October 2026) |
| `cohort_age.py` | `python3 -I code/upstream/01_cohort_tables/cohort_age.py local/crocodeel_valles/matrices OUT.csv` | the age label of the 30 longitudinal cohorts (`aggregates/cohort_age.csv`) | verified, byte for byte |

### 2. Donor pairs, the index and the queries (`02_pairs_and_queries/`)

The scripts of this step read and write in the current folder. `build_pairs.py` runs from the root of the package, after `cp code/upstream/02_pairs_and_queries/adult_cohorts.txt .` (`local/crocodeel_valles/matrices` is relative to the root); `build_hmp2_v2.py` needs `ibdmdb_runs_ena.tsv` and `pairs_hmp2.json` in the current folder. `build_hmp2_v3.py` takes its paths as arguments (usage in its docstring). `requery_todo.json` (the 198 donors of the first `requery_largest.py`) and `all558_runs.json` (the pairs of the human-DNA check and the size of each unitig file, asked to the Logan bucket with HEAD requests) were written by commands that were not kept; the 14 donors of the second `requery_largest.py` are in the `requery_todo.json` that `build_hmp2_v3.py` writes, kept here as `local/work/requery_todo_hmp_ibd_v3.json`.

| Script | Writes | Status |
|---|---|---|
| `build_pairs.py` (reads `adult_cohorts.txt`, the 14 cohorts labelled as adults) | `pairs.json`: first and last stool sample of each donor, with run accessions | verified |
| `build_hmp2.py` (ENA filereport of PRJNA398089 and PRJNA281230) | first version of `pairs_hmp2.json` and `runs_hmp2.json`: first and last stool metagenome of each participant of HMP IBD and HMP 2019 T2D, from the deposits (the curated metadata hold no run accession) | not re-run (network) |
| `build_hmp2_v2.py` (reads `ibdmdb_runs_ena.tsv`, the filereport of PRJNA398089: `run_accession, sample_accession, sample_alias, collection_date, library_strategy, library_source, read_count, first_public`) | updates `pairs_hmp2.json` and writes `runs_hmp2.json`, `runs_hmp2_dates.json`: a sample is a participant-visit, whatever the library suffix of its alias | verified (from the pairs of `build_hmp2.py`; the pairs it wrote are `pairs_hmp2_v2_backup.json`, the version before the rebuild of 8 November 2026) |
| `build_hmp2_v3.py` (reads the pairs of `build_hmp2_v2.py`, `runs_hmp2.json`, `ibdmdb_runs_ena.tsv`, the curated metadata of HMP IBD and a rows file that lists the adult donors, such as `rows_v3.json`) | `pairs_hmp2.json` with the first and last sample of the 48 adult donors of HMP IBD set by the curated link (a visit counts as a curated sample when a curated sample of the donor has twice the read count of one of its ENA runs, the link of `hmp_curated_ids.py`; first and last = earliest and latest such visit with a run in the index; query = deepest run of the last visit), and `requery_todo.json`, the 14 donors whose queried run changes | verified (the file written is identical to the one used) |
| `resolve_biosamples.py` (ENA) | `biosample_runs.json`, `curated_sample_biosample.json`: every run of every curated BioSample | not re-run (network) |
| `fetch_index_ids.py` | `logan_index_ids.txt.gz`: the 4,791,208 runs the index listed on 9 October 2026 | not re-run (network; the index changes) |
| `run_adults.py`, `run_extra.py` | `results.jsonl` and `results_pass2.jsonl`: the 5,000 nearest runs of each donor's last sample, first and second pass (9 October 2026) | not re-run (service) |
| `requery_largest.py` | `results_largest.jsonl`: the 198 donors whose queried run was not the largest indexed run of the sample are queried again with that run (21 October 2026); the 14 donors of HMP IBD whose query `build_hmp2_v3.py` changes are queried the same way (8 November 2026) and their answers are appended to the same file | not re-run (service) |
| `merge_results.py` | `adultes_data_v2/results.jsonl`: the first answers, replaced by the answers of the largest run where there is one | verified, byte for byte |
| `query_sizes.py` | `query_gbp.json`: bases of the run queried for each donor (ENA `base_count`) | verified (from tables of the ENA already downloaded; the script asks the ENA for the runs it lacks) |

The service is `https://similarity.logan-search.org/search`, request `{"sample_id": <run accession>, "sync": true, "topk": 5000}`: a request holds one public run accession. The ranks and the pools are computed on the HGM set, the 291,881 runs of the index taken as human gut metagenomes (`reference/human_gut_stock_291881_runs.txt.gz`, written by `03_ranks/build_stock.py`). It starts from the 172,075 runs labelled human
gut or faecal shotgun in the Logan metadata (`reference/logan_gut_wgs_172075_runs.txt.gz`; the query that made that list was not kept).

### 3. Ranks and extrapolation (`03_ranks/`)

```
python3 -I code/upstream/03_ranks/build_stock.py fetch local/adultes_data_v2 local/adultes/logan_gut_wgs_labelled_172075_ids.txt.gz local/cmd_meta/cmd_meta.tsv local/adultes/ena_stock_candidates.tsv   # public metadata of the runs considered for the HGM set (ENA)
python3 -I code/upstream/03_ranks/build_stock.py select local/adultes_data_v2 local/adultes/logan_gut_wgs_labelled_172075_ids.txt.gz local/cmd_meta/cmd_meta.tsv local/adultes/ena_stock_candidates.tsv local/adultes/logan_gut_wgs_narrow_ids.txt.gz aggregates/stock_summary.json   # the HGM set (= reference/human_gut_stock_291881_runs.txt.gz), Supplementary Table 3
python3 -I code/upstream/03_ranks/mehta_siblings.py fetch local/adultes/ena_PRJNA354235_runs.tsv                    # public run table of the Mehta project (ENA; the request holds the project accession only)
python3 -I code/upstream/03_ranks/mehta_siblings.py fetch local/adultes/ena_PRJNA389280_runs.tsv PRJNA389280        # and of a second deposit of HMP IBD (Supplementary Table 4)
python3 -I code/upstream/03_ranks/mehta_siblings.py map local/adultes_data_v2 local/adultes/ena_PRJNA354235_runs.tsv local/adultes/mehta_sibling_runs.json   # the runs of the Mehta project that the curated BioSamples do not name (the other library of a sample), by ENA alias; verified, byte for byte
python3 -I code/upstream/03_ranks/hmp_sibling_links.py local/rows_v2_full_without_links.json local/hmp2_metadata_2018-08-20.csv local/adultes/ena_PRJNA389280_runs.tsv local/adultes/ena_stock_candidates.tsv local/hmp_sibling_links.json   # the runs of PRJNA306874 and PRJNA389280 that the public table of the consortium (ibdmdb.org) gives to a donor of HMP IBD
python3 -I code/upstream/03_ranks/rows_with_partner.py local/adultes_data_v2 local/rows_v2_full.json --hmp local/hmp_sibling_links.json   # rank of the other sample of every donor (831 donors); the runs given back by the table are neither the donor's other sample nor another person (without --hmp: local/rows_v2_full_without_links.json, the same rows with these runs counted as runs of other people)
python3 -I code/upstream/03_ranks/make_rows_v4.py local/rows_v2_full.json local/rows_v4_497.json              # the 497 donors of the 7 cohorts kept, before the two rules (local/rows_v3.json, the 558 donors of 8 cohorts of the human-DNA check, comes from the same script with --keep-vincent on local/rows_v2_full_without_links.json)
python3 -I code/upstream/03_ranks/hmp_curated_ids.py local/rows_v4_497.json local/adultes local/crocodeel_valles/matrices local/rows_v4_497_curated.json   # the curated sample identifiers of HMP IBD, which join the CroCoDeEL tables
python3 -I code/upstream/03_ranks/donors_left_out.py local/rows_v2_full_without_links.json local/donors_left_out.json   # the 183 donors of the 4 cohorts left out, for the samples studied of Supplementary Note 1 (written to data/ by make_inputs.py)
python3 -I code/upstream/03_ranks/extract_lists.py local/adultes_data_v2 local/rows_v4_497.json local/lists_v4_497.npz --hmp local/hmp_sibling_links.json # similarities of the runs of other people of each of these donors in the HGM set
python3 -I code/upstream/03_ranks/exclude_donors.py local/rows_v4_497.json local/rows_v4_497_curated.json local/lists_v4_497.npz local/crocodeel_valles local   # the 483 donors of the analysis: removes the donors hit by contamination between two people above 5% (8) and by a re-deposit of Costea (6); writes rows_v4.json, rows_v4_curated.json, lists_v4.npz and donors_excluded.json
python3 -I code/upstream/03_ranks/make_adult_ranks.py local/rows_v4.json local/query_gbp.json data/adult_ranks.csv     # the table of Figure 1
python3 -I code/upstream/03_ranks/extrapolation_tables.py local/lists_v4.npz aggregates/extrapolation_summary.json   # Supplementary Note 3 (thinning test, small pools), 400 thinnings at each fraction
python3 -I code/upstream/03_ranks/extrapolation_bootstrap.py data/adult_ranks.csv aggregates/extrapolation_bootstrap.json   # resampling intervals of the three models fitted to the measured curve, the number of resamples doubled until the bounds are stable (Supplementary Fig. 6)
python3 -I code/upstream/03_ranks/evt_similarities.py HYPERVECTOR_DB local/adultes_data_v2 local/rows_v4.json local/lists_v4.npz local/hmp_sibling_links.json local/evt_similarities.npz aggregates/evt_similarities.json   # the similarity of each query to every run of the HGM set, from the hypervector database of the service (HYPERVECTOR_DB: the folder big_db_reordered of the file named in the script, 18 Gb unpacked, not in the package); keeps the 25,000 highest of the runs of other people of each donor
python3 -I code/upstream/03_ranks/evt_model.py local/evt_similarities.npz aggregates/evt_summary.json          # the extreme-value model of the similarity tail, the fourth model of Fig. 1B (Supplementary Note 3, Supplementary Tables 17 and 18), 400 thinnings at each fraction; its intervals use the first 400 of its 100,000 resamples
python3 -I code/upstream/03_ranks/evt_resampling.py local/evt_similarities.npz aggregates/evt_summary.json aggregates/evt_resampling.json   # the number of resamples of this model, by the stopping rule of extrapolation_bootstrap.py, and its resampling band for the tail drawn in Fig. 1B (Supplementary Fig. 6); it stops if the rule does not select the number that evt_model.py uses, and it finds the curve and the intervals of evt_model.py again before it writes
python3 -I code/upstream/03_ranks/evt_similarities.py HYPERVECTOR_DB local/adultes_data_v2 local/rows_v4_497.json local/lists_v4_497.npz local/hmp_sibling_links.json local/evt_similarities_497.npz aggregates/evt_similarities_497.json   # the same similarities for the 497 donors of the cohorts kept, before the two rules
python3 -I code/upstream/03_ranks/evt_removed_donors.py local/evt_similarities_497.npz local/rows_v4_497.json local/rows_v4.json local/evt_similarities.npz aggregates/evt_summary.json aggregates/evt_removed_donors.json   # the extreme-value model made again with the 14 donors that the two rules remove (last row of Supplementary Table 17); it checks that the 483 donors of the analysis have the same similarities in the two files and finds the model of evt_model.py again before it writes
python3 -I code/upstream/03_ranks/make_open_set_scores.py local/lists_v4.npz data/open_set_scores.csv              # panel H of Figure 1: the other sample against the nearest run of other people
python3 -I code/upstream/03_ranks/make_open_set_scores.py local/lists_v4_497.npz local/open_set_scores_497.csv    # the same two similarities for the 497 donors of the cohorts kept, read by costea_same_people.py
python3 -I code/upstream/03_ranks/stock_sensitivity.py local/adultes_data_v2 local/adultes/logan_gut_wgs_labelled_172075_ids.txt.gz local/adultes/ena_stock_candidates.tsv local/rows_v4.json aggregates/stock_sensitivity.json --hmp local/hmp_sibling_links.json   # Supplementary Table 4
```
The three `fetch` commands ask the ENA and were not re-run. The others are verified, byte for byte or as the same data (the summary of `build_stock.py select`, `hmp_sibling_links.py`, `rows_with_partner.py`, `make_rows_v4.py`, `hmp_curated_ids.py`, `donors_left_out.py`, `extract_lists.py`, `exclude_donors.py`, `extrapolation_tables.py`, `extrapolation_bootstrap.py`, `evt_model.py`, `evt_resampling.py`, `evt_removed_donors.py`, `stock_sensitivity.py`), against the files of the paper. `evt_similarities.py` was not re-run, since the hypervector database is not in the package.

### 4. Aggregates of the supplement (`04_aggregates/`)

```
python3 -I code/upstream/04_aggregates/cohort_selection.py local/rows_v2_full.json local/rows_v4.json aggregates/cohort_selection.json local/donors_excluded.json local/rows_v2_full_without_links.json   # Note 1, Tables 1 and 2
python3 -I code/upstream/04_aggregates/pehrsson_lanes.py aggregates/pehrsson_summary.json local/rows_v2_full.json local/adultes_data_v2                                  # Note 1 and Table 5
python3 -I code/upstream/04_aggregates/raymond_antibiotic.py aggregates/raymond_antibiotic.json local/rows_v4.json local/cmd_meta/cmd_meta.tsv                 # panel G of Figure 1 (Raymond, cephalosporins)
python3 -I code/upstream/04_aggregates/nielsen_li_overlap.py local/adultes_data_v2 local/rows_v2_full.json aggregates/nielsen_li_overlap.json      # Note 1
python3 -I code/upstream/04_aggregates/nielsen_li_links.py aggregates/nielsen_li_links.json local/rows_v4.json aggregates/nielsen_li_overlap.json # Note 1, Fig. 2, Tables 6 and 7
python3 -I code/upstream/04_aggregates/costea_same_people.py aggregates/costea_same_people.json local/rows_v4_497.json local/open_set_scores_497.csv   # Note 1, Table 9 (the 497 donors of the cohorts kept, before the two rules)
python3 -I code/upstream/04_aggregates/costea_reads_compare.py local/costea_same_people local/costea_same_people/reads_compare_first_tile.json local/costea_same_people/pairs.tsv   # first tile of the seven pairs of runs (public FASTQ files)
python3 -I code/upstream/04_aggregates/costea_trim_filter.py local/costea_same_people local/costea_same_people/trim_filter_first_tile.json local/costea_same_people/pairs.tsv   # the same tile: is the queried run a trimmed and filtered version of the matched run? (Note 1; aggregates only)
python3 -I code/upstream/04_aggregates/crocodeel_nielsen_li.py local/crocodeel_nielsen_li                                                          # CroCoDeEL on the merged tables of Nielsen and Li (Supplementary Fig. 4C)
python3 -I code/upstream/04_aggregates/reverse_queries_prep.py     # the 45 Li samples of the patients that carry the same identifier in Nielsen: local/nielsen_li_validation/reverse/prep.json
python3 -I code/upstream/04_aggregates/reverse_queries_query.py    # their 5,000 nearest runs (25 October 2026): local/nielsen_li_validation/reverse/answers.jsonl
```
The first six are verified. `costea_reads_compare.py` and `costea_trim_filter.py` need the network, `crocodeel_nielsen_li.py` needs CroCoDeEL 1.3.1 (set `CROCODEEL` to its executable) and the reverse queries need the service: not re-run. `reverse_queries_prep.py` reads
`local/nielsen_li_validation/ena_PRJEB5224.tsv` and `shared_donors_runs.tsv`.

### 5. Human DNA in the samples (`05_human_dna/`)

Sketches (FracMinHash, k = 31, one hash in 1,000, sourmash 4.9.4) of the 1,116 runs of the 558 donors are recomputed from the public unitigs of Logan (`https://s3.amazonaws.com/logan-pub/u/<run>/<run>.unitigs.fa.zst`) and
compared with the sketches of T2T-CHM13v2.0 and GRCh38.p14.

```
python3 -I code/upstream/05_human_dna/all558_pairs.py local pairs_service.json                 # the 558 pairs (query run, other sample) with the similarity of the service
python3 -I code/upstream/05_human_dna/run_projects.py pairs_service.json run_projects.json       # the BioProject of each of their runs (ENA; a request holds run accessions only)
python3 -I code/upstream/05_human_dna/all558_sketch.py run WORKDIR local/humancheck/all558/all558_runs.json 16               # sketch of each run, with the k-mer behind each hash (streamed, nothing else kept)
bash code/upstream/05_human_dna/sketch_grch38.sh WORKDIR                                      # sketch of GRCh38.p14 (the T2T sketch, 2,511,062 hashes, was made with sketch_stream.py from GCF_009914755.1; the command line was not kept)
bash code/upstream/05_human_dna/build_t2t_index.sh WORKDIR                                    # bowtie2 index of T2T for the near-human k-mers
python3 -I code/upstream/05_human_dna/all558_variants.py near WORKDIR SKDIR pairs_service.json IDX THREADS OUT_PREFIX
python3 -I code/upstream/05_human_dna/all558_analysis.py WORKDIR SKDIR pairs_service.json OUT_PREFIX [NEAR.npz]      # share of human hashes, similarity with and without them
python3 -I code/upstream/05_human_dna/all558_pool.py WORKDIR SKDIR pairs_service.json OUT_PREFIX [NEAR.npz]          # rank of the other sample in the pool of the 1,116 runs (its summary is also aggregates/human_dna_pool_summary.json)
python3 -I code/upstream/05_human_dna/all558_join.py ANALYSIS_table.csv local/rows_v3.json OUT_PREFIX        # join with the rank in the HGM set
python3 -I code/upstream/05_human_dna/human_dna_exclusion.py JOIN_TABLE POOL_RANKS JOIN_SUMMARY POOL_SUMMARY aggregates/human_dna_exclusion.json
```
`all558_pairs.py`, `all558_join.py` and the last step (aggregates of the exclusion of Vincent) were re-run for the package and are identical to the files of the paper; the sketches took 2 h 50 min on 16 cores. On 8 November 2026, after `build_hmp2_v3.py` changed 14 pairs of HMP IBD and after `all558_pairs.py` was given the sibling runs of Mehta (`mehta_sibling_runs.json`, as in `rows_with_partner.py`), which changed the other sample of 3 Mehta donors, the 25 runs that the new pairs add were sketched (22 of them with sourmash 4.4.3 and again with 4.9.4, with identical hashes and 31-mers, the 3 others with 4.9.4), and `all558_variants.py near`, `all558_analysis.py`, `all558_pool.py`, `all558_join.py` and `human_dna_exclusion.py` were run again on the 1,116 runs of the 558 donors (the 1,091 other runs keep their sketches of 24 October 2026). `stat_per_run.json` holds, for each of the 1,116 runs, the two counts of the taxonomic analysis that the NCBI runs on every run of the SRA (STAT; human hits and species-level hits; R. Faure, 8 November 2026) and the hashes and human hashes of the sketch of the run; the human share of STAT is the ratio of the two counts.

### 6. Inputs of `reproduce.sh` (`06_inputs/`)

`make_inputs.py` writes `data/donors.json`, `data/donors_left_out.json`, `data/donors_excluded.json`, `data/query_gbp.json`, `data/crocodeel/`, `data/human_dna/` and `data/accessions.tsv` (the run accessions and BioProjects of the two runs of each donor) from the files of `local/` (verified, byte for byte). Usage in its docstring.
