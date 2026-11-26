"""Near-human k-mers: human variants are absent from the reference, so exact hash matching cannot see them.
A 31-mer that is not found exactly in T2T-CHM13v2.0 is aligned to it with bowtie2 (end to end, at most 2 mismatches, or a one-base gap).
usage: python3 -I all558_variants.py spike WORKDIR SKDIR PAIRS.json IDX THREADS                 (recall of the alignment on human 31-mers carrying 1 to 3 changes)
       python3 -I all558_variants.py near WORKDIR SKDIR PAIRS.json IDX THREADS OUT_PREFIX [THRESHOLD_PERCENT [SAMPLE_PER_RUN]]   (align the non-exact-human hashes: all of them for the pairs with a human-rich sketch, a sample for the others)
Writes OUT_PREFIX_near.npz (hashes that align, with their mismatches and k-mers), OUT_PREFIX_near_summary.json and OUT_PREFIX_near_runs.csv (one line per run); the messages of bowtie2 go to OUT_PREFIX_bowtie2.err.
The mode spike writes WORKDIR/all558/spike_recall.json; its messages of bowtie2 go to WORKDIR/all558/bowtie2_spike.err."""
import csv, json, os, subprocess, sys, threading, time
import numpy as np

K = 31
SH = (2 * (K - 1 - np.arange(K))).astype(np.uint64)
ACGT = np.frombuffer(b"ACGT", np.uint8)
CTRL = 20000                                                    # unshared hashes sampled per donor pair as a control


def kmers_to_bytes(km):
    codes = ((km[:, None] >> SH[None, :]) & np.uint64(3)).astype(np.uint8)
    out = np.empty((len(km), K + 1), np.uint8)
    out[:, :K] = ACGT[codes]; out[:, K] = 10
    return out.tobytes()


def mutate(km, n, rng):
    codes = ((km[:, None] >> SH[None, :]) & np.uint64(3)).astype(np.uint8)
    for r in range(len(km)):
        for pos in rng.choice(K, n, replace=False):
            codes[r, pos] = (codes[r, pos] + rng.integers(1, 4)) % 4
    return (codes.astype(np.uint64) << SH[None, :]).sum(axis=1).astype(np.uint64)


def align(kmers, idx, threads, errfile):
    """returns {position in kmers: (mismatches, gap opens)} for the k-mers that align"""
    cmd = ["bowtie2", "-x", idx, "-r", "-U", "-", "--end-to-end", "--score-min", "L,0,-0.4", "-N", "1", "-L", "20", "-i", "C,1", "-D", "10", "-R", "1",
           "--no-unal", "--no-hd", "-p", str(threads)]
    with open(errfile, "w") as ef:
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=ef)

        def feed():
            for s in range(0, len(kmers), 2_000_000):
                p.stdin.write(kmers_to_bytes(kmers[s:s + 2_000_000]))
            p.stdin.close()
        t = threading.Thread(target=feed); t.start()
        hits = {}
        for line in p.stdout:
            f = line.split(b"\t")
            xm = go = 0
            for tag in f[11:]:
                if tag.startswith(b"XM:i:"): xm = int(tag[5:])
                elif tag.startswith(b"XO:i:"): go = int(tag[5:])
            hits[int(f[0])] = (xm, go)
        t.join(); rc = p.wait()
    if rc:
        raise RuntimeError(f"bowtie2 exit {rc}, see {errfile}")
    return hits


def load(SK, run):
    z = np.load(f"{SK}/{run}.npz")
    return z["h"], z["k"]


if __name__ == "__main__":
    mode, W, SK, PAIRS, IDX, THREADS = sys.argv[1:7]
    THREADS = int(THREADS)
    pairs = [p for p in json.load(open(PAIRS)) if os.path.exists(f"{SK}/{p['q']}.npz") and os.path.exists(f"{SK}/{p['partner']}.npz")]
    t2t = np.load(f"{W}/ref_t2t.npy"); H = np.union1d(t2t, np.load(f"{W}/ref_grch38.npy"))
    rng = np.random.default_rng(20261024)
    t0 = time.time()
    if mode == "spike":
        pool_h, pool_k = [], []
        for p in pairs[:40]:
            h, k = load(SK, p["q"]); m = np.isin(h, t2t, assume_unique=True)
            pool_h.append(h[m]); pool_k.append(k[m])
        h = np.concatenate(pool_h); k = np.concatenate(pool_k)
        uh, first = np.unique(h, return_index=True); k = k[first]
        k = k[rng.permutation(len(k))[:5000]]
        print(f"{len(k)} k-mers that match T2T exactly (taken from the sketches)")
        out = {}
        for nmut in (0, 1, 2, 3):
            km = k if nmut == 0 else mutate(k, nmut, rng)
            hits = align(km, IDX, THREADS, f"{W}/all558/bowtie2_spike.err")
            xm = np.array([hits[i][0] if i in hits else -1 for i in range(len(km))])
            out[nmut] = dict(aligned=float((xm >= 0).mean()), aligned_with_the_right_count=float((xm == nmut).mean()), aligned_with_at_most_2=float(((xm >= 0) & (xm <= 2)).mean()))
            print(f"changes: {nmut}  aligned {100 * out[nmut]['aligned']:.1f}%  with exactly {nmut} mismatches {100 * out[nmut]['aligned_with_the_right_count']:.1f}%  with <= 2: {100 * out[nmut]['aligned_with_at_most_2']:.1f}%", flush=True)
        json.dump(out, open(f"{W}/all558/spike_recall.json", "w"), indent=1)
        sys.exit(0)
    OUT = sys.argv[7]
    THR = float(sys.argv[8]) if len(sys.argv) > 8 else 1.0          # runs of the pairs with a sketch at least this human (%) are aligned completely
    NS = int(sys.argv[9]) if len(sys.argv) > 9 else 3000            # k-mers sampled from every other run, to measure the background rate
    runs = sorted({r for p in pairs for r in (p["q"], p["partner"])})
    share = {}
    for r in runs:
        h = np.load(f"{SK}/{r}.npz")["h"]; share[r] = 100 * float(np.isin(h, H, assume_unique=True).mean())
    complete = sorted({r for p in pairs if max(share[p["q"]], share[p["partner"]]) >= THR for r in (p["q"], p["partner"])})
    cset = set(complete)
    print(f"{len(runs)} runs; {len(complete)} runs of pairs with a sketch >= {THR}% human are aligned completely, {len(runs) - len(complete)} others are sampled ({NS} k-mers each)", flush=True)
    Hh, Hk, src = [], [], []
    for r in runs:
        h, k = load(SK, r); keep = ~np.isin(h, H, assume_unique=True); h, k = h[keep], k[keep]
        if r not in cset:
            c = rng.permutation(len(h))[:NS]; h, k = h[c], k[c]
        Hh.append(h); Hk.append(k)
    sizes = np.array([len(x) for x in Hh]); Hh = np.concatenate(Hh); Hk = np.concatenate(Hk)
    uh, first = np.unique(Hh, return_index=True); Uk = Hk[first]
    print(f"{len(uh) / 1e6:.2f} M distinct k-mers to align, prepared in {time.time() - t0:.0f} s", flush=True)
    hits = align(Uk, IDX, THREADS, f"{OUT}_bowtie2.err")
    print(f"alignment done in {time.time() - t0:.0f} s: {len(hits)} distinct k-mers aligned", flush=True)
    pos = np.fromiter(hits.keys(), np.int64, len(hits)); order = np.argsort(pos); pos = pos[order]
    xm = np.array([hits[i][0] for i in pos], np.int8); gap = np.array([hits[i][1] for i in pos], np.int8)
    np.savez(f"{OUT}_near.npz", h=uh[pos], xm=xm, gap=gap, k=Uk[pos], complete_runs=np.array(complete), thr=THR)
    # background rate: sampled k-mers of the runs that are not aligned completely
    aligned_h = uh[pos]; rate = []
    off = np.concatenate([[0], np.cumsum(sizes)])
    for i, r in enumerate(runs):
        if r in cset: continue
        c = min(NS, sizes[i]); 
        if c == 0: continue
        hh = Hh[off[i]:off[i + 1]]
        rate.append((r, c, int(np.isin(hh, aligned_h, assume_unique=True).sum())))
    tot_c = sum(x[1] for x in rate); tot_a = sum(x[2] for x in rate)
    S = dict(runs=len(runs), threshold_percent=THR, complete_runs=len(complete), sampled_runs=len(rate), distinct_kmers=int(len(uh)), aligned=int(len(pos)),
             exact=int(((xm == 0) & (gap == 0)).sum()), mismatch1=int(((xm == 1) & (gap == 0)).sum()), mismatch2=int(((xm == 2) & (gap == 0)).sum()), gapped=int((gap > 0).sum()),
             sampled_kmers=int(tot_c), sampled_kmers_aligned=int(tot_a), background_rate_percent=float(100 * tot_a / max(tot_c, 1)), seconds=round(time.time() - t0))
    json.dump(S, open(f"{OUT}_near_summary.json", "w"), indent=1)
    with open(f"{OUT}_near_runs.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["run", "human_share_percent", "aligned_completely", "sampled_kmers", "sampled_aligned"])
        smp = {x[0]: x for x in rate}
        for r in runs:
            w.writerow([r, f"{share[r]:.4f}", int(r in cset), smp[r][1] if r in smp else "", smp[r][2] if r in smp else ""])
    print(json.dumps(S, indent=1))
