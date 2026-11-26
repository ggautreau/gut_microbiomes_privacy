#!/usr/bin/env bash
# Prepares local/ for the scripts of code/upstream/. Run once from the root of the package:  bash code/upstream/prepare_local.sh
#   1. unpacks the six files that the package holds compressed (the scripts read them uncompressed);
#   2. builds local/adultes_data_v2/, the folder that the scripts of 03_ranks/ and 04_aggregates/ read: the files of local/adultes/, the tables of
#      local/crocodeel_valles/matrices/, and the answers of the service merged by merge_results.py (the first answers, replaced by the answers of the largest run where there is one).
# No network. The files written here are listed in .gitignore.
set -euo pipefail
[ -d local/adultes ] && [ -d code/upstream ] || { echo "run from the root of the package" >&2; exit 1; }
for f in local/adultes/results.jsonl local/adultes/results_pass2.jsonl local/work/results_largest.jsonl local/adultes/ena_stock_candidates.tsv local/adultes/biosample_runs.json \
         local/nielsen_li_validation/reverse/answers.jsonl; do
    [ -f "$f" ] || gzip -dc "$f.gz" > "$f"
done
D=local/adultes_data_v2
mkdir -p $D
for f in biosample_runs.json curated_sample_biosample.json ena_PRJNA389280_runs.tsv logan_gut_wgs_broad_ids.txt.gz logan_gut_wgs_narrow_ids.txt.gz logan_index_ids.txt.gz mehta_sibling_runs.json \
         pairs.json pairs_hmp2.json runs_hmp2.json runs_hmp2_dates.json; do
    cp -p local/adultes/$f $D/$f
done
rm -rf $D/matrices && cp -r local/crocodeel_valles/matrices $D/matrices
python3 -I code/upstream/02_pairs_and_queries/merge_results.py local/adultes $D
echo "local/ is ready"
