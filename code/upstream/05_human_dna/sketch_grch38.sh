#!/bin/bash
# GRCh38.p14: resumable download, integrity check, sketch (k = 31, scaled = 1000); the big file is deleted afterwards
W=$1; PY=${PYTHON:-python3}                                    # Python with sourmash 4.9.4
cd $W/refdl
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
for i in 1 2 3 4 5 6; do
  curl -s -L -C - --retry 5 --retry-delay 5 --connect-timeout 30 -o grch38.fna.gz "$URL"
  if gzip -t grch38.fna.gz 2>/dev/null; then echo "grch38 download ok ($(stat -c %s grch38.fna.gz) bytes)"; break; fi
  echo "grch38 attempt $i incomplete ($(stat -c %s grch38.fna.gz 2>/dev/null) bytes), resuming"; sleep 5
done
gzip -t grch38.fna.gz && gunzip -c grch38.fna.gz | $PY -I ../sketch_stream.py ../ref_grch38.npy && rm -f grch38.fna.gz
echo GRCH38_DONE
