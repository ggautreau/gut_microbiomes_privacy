"""Read a FASTA stream on stdin and write its FracMinHash sketch (k = 31, scaled = 1000, canonical k-mers, as in the index of Faure et al.) to a .npy file.
usage: ... | python3 -I sketch_stream.py OUT.npy"""
import sys, time
import numpy as np
from sourmash import MinHash

out = sys.argv[1]
mh = MinHash(n=0, ksize=31, scaled=1000)
chunk, size, nrec, nbases = [], 0, 0, 0
t0 = time.time()


def flush():
    global chunk, size
    if chunk:
        mh.add_sequence("N".join(chunk), force=True)           # k-mers that span an N are skipped, so records stay separate
        chunk, size = [], 0


rec = []
def push(seq):
    global size, nrec, nbases
    if seq:
        chunk.append(seq); size += len(seq) + 1; nrec += 1; nbases += len(seq)
        if size > 200_000_000:
            flush()


for line in sys.stdin.buffer:
    if line.startswith(b">"):
        push("".join(rec)); rec = []
    else:
        rec.append(line.strip().decode().upper())
push("".join(rec))
flush()
h = np.array(sorted(mh.hashes), dtype=np.uint64)
np.save(out, h)
print(f"{out}: {nrec:,} records, {nbases / 1e6:,.1f} Mbases, {len(h):,} hashes, {time.time() - t0:.0f} s", flush=True)
