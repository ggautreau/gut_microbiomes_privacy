#!/bin/bash
# T2T-CHM13v2.0 (RefSeq GCF_009914755.1): resumable download, integrity check, bowtie2 index; the big files are deleted afterwards
W=$1
cd $W/refdl
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
for i in 1 2 3 4 5 6; do
  curl -s -L -C - --retry 5 --retry-delay 5 --connect-timeout 30 -o t2t.fna.gz "$URL"
  if gzip -t t2t.fna.gz 2>/dev/null; then echo "t2t download ok ($(stat -c %s t2t.fna.gz) bytes)"; break; fi
  echo "t2t attempt $i incomplete ($(stat -c %s t2t.fna.gz 2>/dev/null) bytes), resuming"; sleep 5
done
gzip -t t2t.fna.gz || { echo "T2T_DOWNLOAD_FAILED"; exit 1; }
gunzip -c t2t.fna.gz > t2t.fna
grep -c '^>' t2t.fna
mkdir -p $W/idx_t2t
nice -n 10 bowtie2-build --threads 12 t2t.fna $W/idx_t2t/t2t > $W/idx_t2t/build.log 2>&1
if ls $W/idx_t2t/t2t.1.bt2 $W/idx_t2t/t2t.rev.1.bt2 >/dev/null 2>&1; then echo "INDEX_OK"; rm -f t2t.fna t2t.fna.gz; else echo "INDEX_FAILED"; tail -5 $W/idx_t2t/build.log; fi
du -sh $W/idx_t2t
