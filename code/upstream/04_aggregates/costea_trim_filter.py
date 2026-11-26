#!/usr/bin/env python3
"""Is the queried Costea run (PRJEB17632) a trimmed and filtered version of the matched run (PRJEB41102)? For each of the seven pairs (queried run, nearest run of PRJEB41102) this streams read 1 of the two public
FASTQ files until the first tile ends, as costea_reads_compare.py does, and compares the clusters by (instrument, lane, tile, x, y). Nothing is written but aggregates; no read is kept on disk.

    python3 -I costea_trim_filter.py INV OUT.json PAIRS.tsv

INV        folder with ena_PRJEB17632.tsv and ena_PRJEB41102.tsv (ENA portal: run_accession ... fastq_ftp)
PAIRS.tsv  one line per pair: subject, queried run, matched run (tab separated)
On the clusters that both runs have (the common reads) the sequence of the queried run (q) is compared with that of the matched run (n): identical; q a strictly shorter part of n; n a strictly shorter part of q; neither.
Where q is inside n: the number of bases of n that q lacks at the 5' end and at the 3' end, and the mean base quality (Phred+33) of these bases and of the bases that q keeps; whether the kept bases have the same quality
string in both runs. The reads of n that the queried run lacks (clusters of n absent from q) are compared with the common reads by their mean base quality (below Q30 or not) and length."""
import csv, json, re, subprocess, sys, time
import numpy as np

INV, OUT, PAIRS = sys.argv[1:4]
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


def first_tile_once(run, maxreads=1500000):
    url = "https://" + ftp[run][0]
    p = subprocess.Popen(f"curl -sS --max-time 600 '{url}' | zcat 2>/dev/null", shell=True, stdout=subprocess.PIPE, text=True, errors="replace")
    d, tile, n, complete = {}, None, 0, False
    while n < maxreads:
        h = p.stdout.readline()
        if not h:
            break
        s = p.stdout.readline().rstrip("\n"); p.stdout.readline(); q = p.stdout.readline().rstrip("\n")
        m = HDR.match(h)
        k = key(m.group(1)) if m else None
        if k is None:
            continue
        if tile is None:
            tile = k[:3]
        if k[:3] != tile:
            complete = True                        # the next tile begins: the first tile has been read entirely
            break
        d[k[3:]] = (s, q); n += 1
    p.stdout.close(); p.kill()
    return tile, d, complete


def first_tile(run, tries=4):
    for t in range(tries):                          # a stream that stops inside the first tile (network error) is read again
        tile, d, complete = first_tile_once(run)
        if complete and d:
            return tile, d
        time.sleep(5 * (t + 1))
    raise RuntimeError(f"first tile of {run} not read entirely after {tries} tries")


def phred(qs):
    return np.frombuffer(qs.encode("latin-1"), dtype=np.uint8).astype(float) - 33          # Phred+33, as in the ENA files of both projects


def mean_q(qs):
    return float(phred(qs).mean()) if qs else None


def med(x):
    return float(np.median(x)) if len(x) else None


def lengths(d):
    L = [len(v[0]) for v in d.values()]
    return dict(min=int(min(L)), median=float(np.median(L)), max=int(max(L)))


res = []
for line in open(PAIRS):
    if not line.strip():
        continue
    s, q, nb = line.strip().split("\t")
    tq, dq = first_tile(q); tn, dn = first_tile(nb)
    assert tq == tn, (s, tq, tn)
    common = [c for c in dq if c in dn]
    rel = dict(identical=0, queried_inside_matched=0, matched_inside_queried=0, neither=0, ambiguous_position=0)
    c5, c3, q5, q3, q3_kept, same_q = [], [], [], [], [], 0
    for c in common:
        a, qa = dq[c]; b, qb = dn[c]
        if a == b:
            rel["identical"] += 1
        elif len(a) < len(b) and a in b:
            rel["queried_inside_matched"] += 1
            if b.count(a) > 1:
                rel["ambiguous_position"] += 1; continue
            i = b.find(a); j = i + len(a)
            same_q += int(qa == qb[i:j])
            c5.append(i); c3.append(len(b) - j)
            if i:
                q5.append(mean_q(qb[:i]))
            if len(b) - j:
                q3.append(mean_q(qb[j:])); q3_kept.append(mean_q(qb[i:j]))
        elif len(b) < len(a) and b in a:
            rel["matched_inside_queried"] += 1
        else:
            rel["neither"] += 1
    absent = [c for c in dn if c not in dq]
    mq_abs = [mean_q(dn[c][1]) for c in absent]; mq_pres = [mean_q(dn[c][1]) for c in common]
    n_in = rel["queried_inside_matched"] - rel["ambiguous_position"]
    row = dict(subject=s, reads_q=len(dq), reads_n=len(dn), common=len(common), only_in_q=len([c for c in dq if c not in dn]), only_in_n=len(absent), relation=rel,
               same_quality_on_kept_bases=dict(same=same_q, checked=n_in), length_matched=lengths(dn), length_queried=lengths(dq),
               cut_5prime=dict(share_of_common_reads=len([x for x in c5 if x > 0]) / len(common), median_bases=med(c5), median_quality_of_the_bases_cut=med(q5)),
               cut_3prime=dict(share_of_common_reads=len([x for x in c3 if x > 0]) / len(common), median_bases=med([x for x in c3 if x > 0]), median_quality_of_the_tail=med(q3), median_quality_of_the_bases_kept_in_these_reads=med(q3_kept)),
               absent_from_queried_run=dict(reads=len(absent), share_of_matched_reads=len(absent) / len(dn), median_mean_quality=med(mq_abs), share_below_q30=float(np.mean([x < 30 for x in mq_abs]))),
               present_in_both_runs=dict(median_mean_quality=med(mq_pres), share_below_q30=float(np.mean([x < 30 for x in mq_pres]))))
    res.append(row); print(s, row["common"], row["relation"], flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
