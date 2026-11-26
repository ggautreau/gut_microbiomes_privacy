#!/usr/bin/env bash
# Regenerates every figure, every number of the text and every supplementary table of the Brief Communication from the tables of this reproduction package,
# then compares each file with the one used in the manuscript (expected/) and recomputes the headline numbers from the definitions. About 40 seconds on a machine with several cores, about a minute on one core (with one thread, figure1_quadratic.json can differ from the expected file in the last digit, which is reported as NEAR).
set -euo pipefail
cd "$(dirname "$0")"
# the files of the package are the ones that were checked in (set SKIP_SHA=1 to run on inputs you changed)
if [ -z "${SKIP_SHA:-}" ]; then
    if sha256sum -c SHA256SUMS --quiet; then echo "SHA256SUMS: every file is as shipped"; else echo "ERROR: a file differs from the one that was shipped (SKIP_SHA=1 runs anyway)"; exit 1; fi
fi
OUT=output
rm -rf "$OUT"; mkdir -p "$OUT/manuscript"

# Supplementary Fig. 5, first: the three models of Fig. 1B fitted to the measured curve, on five windows (errors of the thinning test from aggregates/extrapolation_summary.json), and the extreme-value model for its five tails (aggregates/evt_summary.json); Figure 1 checks the fits that it draws against its values
python3 -I code/supplementary_figure5.py data/adult_ranks.csv aggregates/extrapolation_summary.json "$OUT/supplementary_figure5" aggregates/evt_summary.json

# Supplementary Fig. 6: cohort-then-donor resampling intervals for the four models and adaptive selection of B (the extreme-value model from aggregates/evt_resampling.json)
python3 -I code/supplementary_figure6.py data/adult_ranks.csv aggregates/extrapolation_bootstrap.json "$OUT/supplementary_figure6" aggregates/evt_resampling.json

# Figure 1: data/adult_ranks.csv (483 donors), aggregates/raymond_antibiotic.json (panel G), aggregates/extrapolation_summary.json (the windows of the fits drawn), the values of the three models fitted to the measured curve and aggregates/evt_summary.json (the extreme-value model of the similarity tail) -> figure1.png/.pdf and the straight lines of panel B
python3 -I code/figure1_brief_communication.py . "$OUT/figure1" "$OUT/supplementary_figure5_values.json"

# Supplementary Fig. 1: the choice of the cohorts, the donors and their samples, from aggregates/cohort_selection.json and the donors
python3 -I code/supplementary_figure1.py aggregates/cohort_selection.json data/donors.json data/donors_excluded.json "$OUT/supplementary_figure1"

# Supplementary Fig. 2: links between the samples of Nielsen and Li
python3 -I code/supplementary_figure2.py aggregates/nielsen_li_links.json aggregates/nielsen_li_overlap.json data/donors.json "$OUT/supplementary_figure2"

# Supplementary Fig. 3: human DNA in the samples and linking (558 donors, Vincent included), with the human share of the NCBI taxonomic analysis of the same runs (panels D and E)
python3 -I code/supplementary_figure3.py data/human_dna/donors_human_share.csv aggregates/human_dna_pool_summary.json "$OUT/supplementary_figure3"

# Supplementary Fig. 4: CroCoDeEL events between two people, cohort by cohort, and between Nielsen and Li (panel C)
python3 -I code/supplementary_figure4.py data/donors.json "$OUT/supplementary_figure4" data/crocodeel/nielsen_li_merged.conta.tsv data/crocodeel

# Every number of the main text and the Methods that depends on the donors
python3 -I code/reproduce_numbers.py data/donors.json data/query_gbp.json aggregates/extrapolation_summary.json "$OUT/supplementary_figure4_stats.json" "$OUT/figure1_quadratic.json" "$OUT/numbers_in_text.json"

# Supplementary Notes 1 to 4 and Tables 1 to 18 (supplementary.tex)
cp data/manuscript_title.txt "$OUT/manuscript/title.txt"
python3 -I code/build_supplementary.py data/donors.json data/query_gbp.json "$OUT/manuscript" "$OUT/supplementary_figure4_stats.json" aggregates/pehrsson_summary.json \
    aggregates/extrapolation_summary.json "$OUT/supplementary_figure5_values.json" "$OUT/figure1_quadratic.json" aggregates/nielsen_li_overlap.json \
    aggregates/human_dna_exclusion.json aggregates/cohort_selection.json aggregates/nielsen_li_links.json aggregates/costea_same_people.json aggregates/stock_summary.json aggregates/stock_sensitivity.json \
    "$OUT/supplementary_figure3_values.json" aggregates/evt_summary.json aggregates/evt_similarities.json aggregates/extrapolation_bootstrap.json aggregates/evt_resampling.json aggregates/evt_removed_donors.json

python3 -I check_reproduction.py
python3 -I check_independent.py
