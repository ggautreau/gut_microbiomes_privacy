"""Sketch every Logan unitig file needed for the 558 donors and keep the k-mer behind each retained hash.

usage: python3 -I all558_sketch.py selftest RUN REFERENCE.npy      (one run, compared with a sketch made by sketch_stream.py)
       python3 -I all558_sketch.py run OUTDIR RUNS.json NPROC

Sketch = FracMinHash, k = 31, scaled = 1000, canonical k-mers (sourmash 4.9.4), the same sketch as sketch_stream.py; the hash set is identical because
the hashes come from MinHash.seq_to_hashes instead of MinHash.add_sequence. The records of a chunk are joined by N, as in sketch_stream.py.
Output for each run, OUTDIR/sk/<run>.npz:
  h      sorted unique hashes of the sketch (uint64)
  k      canonical 31-mer behind each hash (uint64, 2 bits per base, A=0 C=1 G=2 T=3, first base in the high bits)
  nrec, nbases
Nothing but these files is kept: the unitigs are streamed and discarded."""
import sys, os, json, time, subprocess, warnings, random
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

warnings.simplefilter("ignore")
from sourmash import MinHash

K = 31
MH = MinHash(n=0, ksize=K, scaled=1000)
MAXH = MH.max_hash
CHUNK = 2_000_000                                   # bases per call; the Python list of hashes of a chunk stays below 100 MB
LUT = np.full(256, 255, np.uint8)
for i, ch in enumerate(b"ACGT"):
    LUT[ch] = i
POW = np.uint64(4) ** np.arange(K - 1, -1, -1, dtype=np.uint64)


def process_chunk(seq):
    """hashes retained by the sketch in this chunk, with their canonical 31-mer packed in a uint64"""
    if len(seq) < K:
        return np.empty(0, np.uint64), np.empty(0, np.uint64)
    a = np.array(MH.seq_to_hashes(seq, force=True, bad_kmers_as_zeroes=True), dtype=np.uint64)
    sel = np.nonzero((a != 0) & (a <= MAXH))[0]
    if sel.size == 0:
        return np.empty(0, np.uint64), np.empty(0, np.uint64)
    win = sliding_window_view(np.frombuffer(seq.encode("ascii"), dtype=np.uint8), K)[sel]
    codes = LUT[win]
    assert codes.max() <= 3, "a retained hash comes from a window with a non-ACGT base"
    fwd = (codes.astype(np.uint64) * POW).sum(axis=1)
    rev = ((3 - codes[:, ::-1]).astype(np.uint64) * POW).sum(axis=1)
    return a[sel], np.minimum(fwd, rev)              # A < C < G < T, so the numeric minimum is the lexicographic one used for the canonical k-mer


def decode(v):
    return "".join("ACGT"[(int(v) >> (2 * (K - 1 - i))) & 3] for i in range(K))


def stream_sketch(url):
    p1 = subprocess.Popen(["curl", "-sS", "--fail", "-L", "--connect-timeout", "30", "--speed-limit", "20000", "--speed-time", "60", url],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(["zstd", "-dc"], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close()
    try:
        hs, ks = [], []
        chunk, size, nrec, nbases, rec = [], 0, 0, 0, []

        def push(seq):
            nonlocal chunk, size, nrec, nbases
            if seq:
                chunk.append(seq); size += len(seq) + 1; nrec += 1; nbases += len(seq)
                if size > CHUNK:
                    h, k = process_chunk("N".join(chunk)); hs.append(h); ks.append(k); chunk, size = [], 0

        for line in p2.stdout:
            if line.startswith(b">"):
                push("".join(rec)); rec = []
            else:
                rec.append(line.strip().decode().upper())
        push("".join(rec))
        if chunk:
            h, k = process_chunk("N".join(chunk)); hs.append(h); ks.append(k)
        err2 = p2.stderr.read().decode(errors="replace"); rc2 = p2.wait()
        err1 = p1.stderr.read().decode(errors="replace"); rc1 = p1.wait()
        if rc1 or rc2:
            raise RuntimeError(f"curl rc={rc1} {err1.strip()[:200]} | zstd rc={rc2} {err2.strip()[:200]}")
    except BaseException:
        for p in (p1, p2):
            try: p.kill()
            except Exception: pass
        raise
    h = np.concatenate(hs); k = np.concatenate(ks)
    uh, idx = np.unique(h, return_index=True)
    return uh, k[idx], nrec, nbases


def url_of(run):
    return f"https://s3.amazonaws.com/logan-pub/u/{run}/{run}.unitigs.fa.zst"


def work(a):
    run, outdir = a
    path = f"{outdir}/sk/{run}.npz"
    if os.path.exists(path):
        return run, "cached", {}
    err = ""
    for attempt in range(5):
        t0 = time.time()
        try:
            h, k, nrec, nb = stream_sketch(url_of(run))
            tmp = f"{path}.tmp{os.getpid()}.npz"
            np.savez(tmp, h=h, k=k, nrec=nrec, nbases=nb)
            os.replace(tmp, path)
            return run, "ok", dict(nrec=int(nrec), nbases=int(nb), nhash=int(len(h)), secs=round(time.time() - t0, 1), attempt=attempt + 1)
        except Exception as e:
            err = f"{type(e).__name__}: {str(e)[:300]}"
            time.sleep(15 * (attempt + 1))
    return run, "FAILED", dict(error=err)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "selftest":
        run, ref = sys.argv[2], sys.argv[3]
        t0 = time.time()
        h, k, nrec, nb = stream_sketch(url_of(run))
        R = np.load(ref)
        print(f"{run}: {nrec:,} records, {nb / 1e6:,.1f} Mbases, {len(h):,} hashes in {time.time() - t0:.0f} s ({nb / 1e6 / (time.time() - t0):.1f} Mbases/s)")
        print("identical to the sketch made with add_sequence:", bool(np.array_equal(h, R)), f"({len(R):,} hashes)")
        from sourmash.minhash import hash_murmur
        idx = random.Random(1).sample(range(len(h)), min(5000, len(h)))
        bad = sum(hash_murmur(decode(k[i])) != int(h[i]) for i in idx)
        print(f"k-mer recovery: {len(idx)} hashes recomputed from their decoded 31-mer, {bad} mismatches")
        sys.exit(0 if (np.array_equal(h, R) and bad == 0) else 1)
    elif mode == "run":
        import multiprocessing as mp
        outdir, runs_json, nproc = sys.argv[2], sys.argv[3], int(sys.argv[4])
        os.makedirs(f"{outdir}/sk", exist_ok=True)
        D = json.load(open(runs_json))
        size = {r: s for r, s in D["size"].items() if s}
        runs = sorted(size, key=lambda r: -size[r])      # biggest first: the tail of the job is made of small files
        todo = [(r, outdir) for r in runs]
        print(f"{len(todo)} runs, {sum(size.values()) / 1e9:.1f} GB of unitigs, {nproc} processes", flush=True)
        t0 = time.time(); done = 0; nfail = 0
        with mp.Pool(nproc) as pool, open(f"{outdir}/log.jsonl", "a") as log:
            for run, status, info in pool.imap_unordered(work, todo, chunksize=1):
                done += 1; nfail += status == "FAILED"
                log.write(json.dumps(dict(run=run, status=status, **info)) + "\n"); log.flush()
                if done % 25 == 0 or status == "FAILED":
                    print(f"{done}/{len(todo)} runs, {nfail} failed, {(time.time() - t0) / 60:.1f} min", flush=True)
        print(f"ALL_DONE {done} runs, {nfail} failed, {(time.time() - t0) / 60:.1f} min", flush=True)
