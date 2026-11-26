#!/usr/bin/env python3
"""For each of the seven pairs (queried run, nearest run of PRJEB41102): streams read 1 of the two public FASTQ files until the first tile ends and compares the clusters by (instrument, lane, tile, x, y).
Nothing is written but the counts. No read is kept on disk.

    python3 -I costea_reads_compare.py INV OUT.json PAIRS.tsv

INV        folder with ena_PRJEB17632.tsv and ena_PRJEB41102.tsv (ENA portal: run_accession ... fastq_ftp)
PAIRS.tsv  one line per pair: subject, queried run, matched run (tab separated)"""
import csv, json, re, subprocess, sys
INV = sys.argv[1]; OUT = sys.argv[2]
PAIRS = [l.strip().split("\t") for l in open(sys.argv[3]) if l.strip()]
ftp = {}
for P in ("PRJEB17632", "PRJEB41102"):
    for r in csv.DictReader(open(f"{INV}/ena_{P}.tsv"), delimiter="\t"):
        ftp[r["run_accession"]] = r["fastq_ftp"].split(";")
HDR = re.compile(r"^@\S+\s+(\S+)")


def key(name):
    name = name.split("#")[0].split("/")[0]
    t = name.split(":")
    if len(t) == 5:                                # instrument:lane:tile:x:y
        return (t[0], t[1], t[2], t[3], t[4])
    if len(t) == 7:                                # instrument:run:flowcell:lane:tile:x:y
        return (t[0], t[3], t[4], t[5], t[6])
    return None


def first_tile(run, maxreads=1500000):
    url = "https://" + ftp[run][0]
    p = subprocess.Popen(f"curl -sS --max-time 400 '{url}' | zcat 2>/dev/null", shell=True, stdout=subprocess.PIPE, text=True, errors="replace")
    d, tile, n, meta = {}, None, 0, None
    while n < maxreads:
        h = p.stdout.readline()
        if not h:
            break
        s = p.stdout.readline().rstrip("\n"); p.stdout.readline(); p.stdout.readline()
        m = HDR.match(h)
        k = key(m.group(1)) if m else None
        if k is None:
            continue
        if tile is None:
            tile = k[:3]; meta = m.group(1)
        if k[:3] != tile:
            break
        d[k[3:]] = s; n += 1
    p.stdout.close(); p.kill()
    return tile, meta, d


res = []
for s, q, nb in PAIRS:
    tq, mq, dq = first_tile(q)
    tn, mn, dn = first_tile(nb)
    common = set(dq) & set(dn)
    same = sum(1 for c in common if dq[c] == dn[c])
    pref = sum(1 for c in common if dq[c] in dn[c] or dn[c] in dq[c])
    lq = sorted(len(v) for v in dq.values()); ln = sorted(len(v) for v in dn.values())
    row = dict(subject=s, same_tile=tq == tn, tile_q=list(tq) if tq else None, tile_n=list(tn) if tn else None, reads_q=len(dq), reads_n=len(dn), common=len(common), identical=same, contained=pref,
               median_len_q=lq[len(lq) // 2] if lq else None, median_len_n=ln[len(ln) // 2] if ln else None)
    res.append(row); print(row, flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
