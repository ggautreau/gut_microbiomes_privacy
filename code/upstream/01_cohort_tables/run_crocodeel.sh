#!/bin/bash
E=${CROCODEEL_ENV:-$HOME/miniconda3/envs/CroCoDeEL_env}                  # a conda environment with CroCoDeEL 1.3.1 (pip install crocodeel==1.3.1)
cd local/crocodeel_valles
mkdir -p results
for f in matrices/*.species_relab.tsv; do
  s=$(basename "$f" .species_relab.tsv)
  [ -f "results/${s}.conta.tsv" ] && { echo "== $s déjà fait"; continue; }
  n=$(head -1 "$f" | awk -F'\t' '{print NF-1}')
  echo "== $s ($n échantillons)"
  "$E/bin/crocodeel" search_conta -s "$f" -c "results/${s}.conta.tsv" 2>&1 \
    | grep -E "contamination events detected|samples contaminated|Search completed"
done
echo "TOUT_FINI"
