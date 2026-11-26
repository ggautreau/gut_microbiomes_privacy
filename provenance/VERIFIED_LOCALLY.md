# Steps re-run for this reproduction package

On 25 and 26 November 2026, on the machine where the paper was prepared, `rerun_upstream.py` ran the scripts of `code/upstream/` that need neither the network nor a long computation (all but `reverse_queries_prep.py`, below), in a scratch copy of this package,
from the files of `local/` (their sizes and sha256 are in `MANIFEST_upstream.tsv`) and with the network cut. What each one writes was compared with the file of the package (`local/`, `data/`, `aggregates/`, `reference/`).
"Identical" means the same bytes (JSON files: the same data); "equal to 3 decimals", JSON files whose floats are equal to 3 decimals (or to 1e-9 relative) but not in their last digits (at most 3e-7 relative here), as `check_reproduction.py` also accepts them: another Python or numpy moves the optimum of a fit. Every step below gave identical files or files equal to 3 decimals. One of them reads a table that is not in the package and was run with it (`--extra`); without this table
`rerun_upstream.py` leaves that step out.

| Script | Writes | Result against the file of the package | Time | Note |
|---|---|---|---|---|
| `build_pairs.py` | local/crocodeel_valles/matrices -> pairs.json (first and last stool sample of each donor of the 14 adult cohorts) | **identical** | 0 s | compared with local/adultes/pairs.json |
| `build_hmp2_v2.py` | ibdmdb_runs_ena.tsv (ENA) and the pairs of build_hmp2.py -> pairs_hmp2.json, runs_hmp2.json, runs_hmp2_dates.json | **identical** | 0 s | compared with the files of local/adultes/ (the pairs before the rebuild of 8 November 2026) |
| `build_hmp2_v3.py` | the pairs of build_hmp2_v2.py, the curated metadata and ibdmdb_runs_ena.tsv -> pairs_hmp2.json (first and last sample of the 48 adult donors of HMP IBD by the curated link), requery_todo.json (the 14 donors whose query changes) | **identical** | 1 s | compared with the files of local/ |
| `merge_results.py` | first answers + the 198 + 14 queries with the largest run -> adultes_data_v2/results.jsonl | **identical** | 2 s | the sha256 of the file of the paper (MANIFEST_upstream.tsv) |
| `cohort_age.py` | local/crocodeel_valles/matrices -> aggregates/cohort_age.csv (age label of the 30 longitudinal cohorts) | **identical** | 0 s | byte for byte |
| `build_stock.py select` | -> reference/human_gut_stock_291881_runs.txt.gz, aggregates/stock_summary.json | **identical** | 4 s | from the table of the ENA of local/adultes/ |
| `mehta_siblings.py map` | -> mehta_sibling_runs.json (runs of the Mehta project that the curated BioSamples do not name, by ENA alias) | **identical** | 0 s | byte for byte |
| `rows_with_partner.py` | local/adultes_data_v2 -> rows_v2_full_without_links.json (without --hmp) | **identical** | 3 s | 831 donors |
| `hmp_sibling_links.py` | the public table of the HMP consortium -> hmp_sibling_links.json (runs of PRJNA306874 and PRJNA389280 given to a donor) | **identical** | 0 s | with hmp2_metadata_2018-08-20.csv, which is not in the package |
| `rows_with_partner.py --hmp` | -> rows_v2_full.json | **identical** | 4 s | 831 donors |
| `make_rows_v4.py` | -> rows_v4_497.json (the 497 donors of the 7 cohorts kept, before the two rules) | **identical** | 0 s |  |
| `make_rows_v4.py --keep-vincent` | -> rows_v3.json (558 donors, 8 cohorts, for the human-DNA check, without --hmp) | **identical** | 0 s |  |
| `donors_left_out.py` | -> donors_left_out.json (the 183 donors of the 4 cohorts left out) | **identical** | 0 s |  |
| `all558_join.py` | -> join_final_table.csv, join_final_summary.json (the 558 donors: share of human DNA and rank) | **identical** | 0 s | from the table of the human-DNA check |
| `all558_pairs.py` | -> pairs_service.json (the 558 pairs of the human-DNA check: query run, run of the other sample) | **identical** | 1 s | also identical to the list of runs sketched |
| `extract_lists.py --hmp` | -> lists_v4_497.npz (similarities of the runs of other people of each of these 497 donors) | **identical** | 4 s |  |
| `hmp_curated_ids.py` | -> rows_v4_497_curated.json (the curated sample identifiers of HMP IBD) | **identical** | 0 s |  |
| `exclude_donors.py` | -> rows_v4.json, rows_v4_curated.json, lists_v4.npz (the 483 donors) and donors_excluded.json (the 14 donors removed) | **identical** | 0 s |  |
| `make_open_set_scores.py (the 497 donors)` | -> open_set_scores_497.csv (for costea_same_people.py) | **identical** | 0 s | byte for byte |
| `make_open_set_scores.py` | -> data/open_set_scores.csv (panel H of Figure 1) | **identical** | 0 s | byte for byte |
| `query_sizes.py` | -> query_gbp.json (from the two tables of the ENA of local/adultes/) | **identical** | 1 s | the 839 donors of local/query_gbp.json; no request to the ENA |
| `make_adult_ranks.py` | -> data/adult_ranks.csv (the table of Figure 1) | **identical** | 0 s | byte for byte |
| `extrapolation_tables.py` | -> aggregates/extrapolation_summary.json | **identical** | 548 s | 400 thinnings at each fraction |
| `extrapolation_bootstrap.py` | -> aggregates/extrapolation_bootstrap.json | **equal to 3 decimals** | 68 s | adaptive doubling bootstrap |
| `evt_model.py` | -> aggregates/evt_summary.json (the extreme-value model of the similarity tail, Supplementary Table 18) | **equal to 3 decimals** | 17663 s | 400 thinnings at each fraction |
| `evt_resampling.py` | -> aggregates/evt_resampling.json (the resampling band of the extreme-value model, Supplementary Fig. 6) | **equal to 3 decimals** | 440 s | finds the curve and the intervals of evt_model.py again before it writes |
| `evt_removed_donors.py` | -> aggregates/evt_removed_donors.json (the extreme-value model with the donors that the two rules remove, last row of Supplementary Table 17) | **identical** | 489 s | finds the model of evt_model.py again on the 483 donors before it writes |
| `stock_sensitivity.py` | -> aggregates/stock_sensitivity.json | **identical** | 6 s |  |
| `cohort_selection.py` | -> aggregates/cohort_selection.json | **identical** | 2 s |  |
| `pehrsson_lanes.py` | -> aggregates/pehrsson_summary.json | **identical** | 1 s |  |
| `nielsen_li_overlap.py` | -> aggregates/nielsen_li_overlap.json | **identical** | 1 s |  |
| `nielsen_li_links.py` | -> aggregates/nielsen_li_links.json | **identical** | 2 s |  |
| `raymond_antibiotic.py` | -> aggregates/raymond_antibiotic.json (panel G of Figure 1) | **identical** | 0 s |  |
| `costea_same_people.py` | -> aggregates/costea_same_people.json | **identical** | 1 s |  |
| `human_dna_exclusion.py` | -> aggregates/human_dna_exclusion.json | **identical** | 0 s | from the tables of the human-DNA check |
| `make_inputs.py` | -> data/ (29 files) | **identical** | 0 s | byte for byte |

## Steps not re-run

All but `reverse_queries_prep.py` need the network (ENA, the similarity service, the Logan unitigs), R and curatedMetagenomicData, CroCoDeEL, the hypervector database of the service, or hours of computation. The scripts are the ones that produced the files of the paper.

| Script | Why |
|---|---|
| `extract_all.R` | R with ExperimentHub, whose download cache is not in the package; run on 15 November 2026 from that cache, into an empty folder, it rewrote the 186 tables of 1 October 2026 byte for byte |
| `run_crocodeel.sh`, `crocodeel_nielsen_li.py` | CroCoDeEL 1.3.1 (about 2.5 hours for the 93 tables on 1 October 2026; under a minute for the merged Nielsen and Li table) |
| `build_hmp2.py`, `resolve_biosamples.py`, `fetch_index_ids.py`, `run_projects.py` | ENA and service listing (the index changes over time); the BioProject of each run of the human-DNA check, asked to the ENA on 8 November 2026 (a request holds run accessions only) |
| `costea_reads_compare.py`, `costea_trim_filter.py` | public FASTQ files of the ENA, streamed on 31 October and on 8 November 2026 (first tile of the seven pairs of runs; only the counts are kept; the second script tests whether the queried run is a trimmed and filtered version of the matched run) |
| `run_adults.py`, `run_extra.py`, `requery_largest.py`, `reverse_queries_query.py` | queries of the public similarity service of 9, 21 and 25 October and 8 November 2026 (one request at a time, 1 s apart) |
| `reverse_queries_prep.py` | no network; from files of `local/` it prepares the reverse queries of the 45 Li samples and writes `local/nielsen_li_validation/reverse/prep.json`; it is not among the steps of `rerun_upstream.py` |
| `build_stock.py fetch`, `mehta_siblings.py fetch` | public metadata and run tables of the ENA (1 November 2026; the records of the 68 runs that the queries of 8 November 2026 add to the table of `build_stock.py fetch` were asked on that day) |
| `evt_similarities.py` | the hypervector database of the service (a file of 10 GB, 18 GB unpacked, downloaded on 23 November 2026), which is not in the package; from it the script computed the similarity of each query to every run of the HGM set and wrote `local/evt_similarities.npz` and `aggregates/evt_similarities.json` (on the 2,284,827 pairs that the service had returned, the two similarities differ by less than 0.00001); run again on `local/rows_v4_497.json` and `local/lists_v4_497.npz` on 24 November 2026, it wrote `local/evt_similarities_497.npz` and `aggregates/evt_similarities_497.json` for the 497 donors of the cohorts kept (2,352,301 pairs, the same largest difference), and `evt_removed_donors.py` checks that the 483 donors of the analysis have the same similarities in the two files |
| `all558_sketch.py`, `sketch_grch38.sh`, `build_t2t_index.sh`, `all558_variants.py`, `all558_analysis.py`, `all558_pool.py` | the human-DNA check: sketches of 1,116 runs recomputed from the Logan unitigs (2 h 50 min on 16 cores; 25 of them on 8 November 2026, after the rebuild of the HMP IBD pairs and the correction of the other sample of 3 Mehta donors; 22 with sourmash 4.4.3 and again with 4.9.4, identical, the 3 others with 4.9.4), bowtie2 on T2T. The steps from `all558_variants.py` to `all558_join.py` were run again on 8 November 2026 on the 1,116 runs, and the files they wrote are the ones read above; only the aggregation of their results (`human_dna_exclusion.py`) and the join are in the table |

Every script of `code/upstream/` compiles (`python3 -m py_compile`) and every shell script parses (`bash -n`).
