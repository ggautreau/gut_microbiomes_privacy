#!/usr/bin/env python3
"""Independent recomputation of the headline numbers of the main text, written from the definitions alone (no scipy, no code of figure1_brief_communication.py or reproduce_numbers.py).

    python3 -I check_independent.py [OUTPUT_DIR]          # run reproduce.sh first: it compares with output/numbers_in_text.json

Input: data/adult_ranks.csv (483 donors). rank_any_narrow = 1 + the number of runs of the HGM set (M = 291,881 human gut metagenomes of the index) that are more similar to the query than the donor's other sample;
an empty value = the other sample is not among the 5,000 nearest runs (a failure at every pool size). k = rank - 1.
A pool of N runs holds the donor's other sample and N - 1 runs drawn without replacement from the HGM set.
  first in the pool:           prod_{j=0}^{N-2} (M - k - j) / (M - j)                                  (0 when M - k < N - 1; ties favour the other sample)
  among the 10 nearest:        sum_{x=0}^{9} C(k, x) C(M - k, N - 1 - x) / C(M, N - 1)                (hypergeometric: at most 9 of the N - 1 runs drawn are among the k runs above the other sample)
R(N) is the mean over the donors. Checked against the numbers of numbers_in_text.json: the share first at N = 2, 100, 1,000, about 1e4, 10,000, 100,000 and M, and among the 10 nearest at about 1e4 and M; the number of donors first,
among the first 10, first with the first sample alone as the other sample, and not among the runs returned; the size of the 30-day subset and its share that is first.
"""
import csv, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "output")
M = 291_881
rows = [r for r in csv.DictReader(open(os.path.join(HERE, "data", "adult_ranks.csv"))) if r["age_class"] in ("adult", "adult_by_design") and r["q_lib"] == "DNA"]
want = json.load(open(os.path.join(OUT, "numbers_in_text.json")))
D = len(rows)
rank = [int(r["rank_any_narrow"]) if r["rank_any_narrow"] != "" else None for r in rows]
ok = True


def check(name, got, ref, tol=1e-6):
    global ok
    good = abs(got - ref) <= tol
    ok = ok and good
    print(f"{'PASS' if good else 'FAIL'}  {name:46s} recomputed {got:12.6f}   numbers_in_text.json {ref:12.6f}")


def log_comb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def p_first(k, N):                                           # the explicit product of the definition
    n1 = N - 1
    if M - k < n1:
        return 0.0
    s = 0.0
    for j in range(n1):
        s += math.log((M - k - j) / (M - j))
    return math.exp(s)


def p_top10(k, N):                                           # the sum of the hypergeometric probabilities of 0 to 9 runs above the other sample among the N - 1 runs drawn
    n1 = N - 1
    t = 0.0
    for x in range(10):
        if x > k or n1 - x > M - k or n1 - x < 0:
            continue
        t += math.exp(log_comb(k, x) + log_comb(M - k, n1 - x) - log_comb(M, n1))
    return t


found = [(r - 1) for r in rank if r is not None]
n1e4 = want["pool"]["N_near_1e4"]
print(f"{D} donors, {len(found)} with the other sample among the 5,000 nearest runs\n")
check("donors first in the HGM set (count)", sum(1 for r in rank if r == 1), want["rank1_narrow"], 0)
check("donors among the first 10 (count)", sum(1 for r in rank if r is not None and r <= 10), want["top10_narrow"], 0)
check("donors first, first sample alone as the other sample (count)", sum(1 for r in rows if r["rank_first_narrow"] == "1"), want["rank1_first_sample"], 0)
check("donors not among the 5,000 nearest (count)", sum(1 for r in rank if r is None), want["absent"], 0)
for label, N, key_first, key_top in (("N = 2", 2, "first_N2", None), ("N = 100", 100, "first_N1e2", None), ("N = 1000", 1000, "first_N1e3", None), (f"N = {n1e4}", n1e4, "first_N1e4", "top10_N1e4"), ("N = 10000", 10_000, "first_N10000", None), ("N = 100000", 100_000, "first_N1e5", None), (f"N = M = {M}", M, "first_NM", "top10_NM")):
    # donors with the same k give the same probability: compute once per distinct k
    cache = {}
    tot = 0.0
    for k in found:
        if k not in cache:
            cache[k] = p_first(k, N)
        tot += cache[k]
    check(f"share first in a pool, {label} (%)", 100 * tot / D, want["pool"][key_first], 1e-6)
    if key_top:
        cache = {}
        tot = 0.0
        for k in found:
            if k not in cache:
                cache[k] = p_top10(k, N)
            tot += cache[k]
        check(f"share among the 10 nearest, {label} (%)", 100 * tot / D, want["pool"][key_top], 1e-6)
far = [r for r in rows if r["far_possible"] == "1"]
check("donors with another sample at least 30 days away", len(far), want["far"]["n"], 0)
check("  of them first (%)", 100 * sum(1 for r in far if r["rank_far_narrow"] == "1") / len(far), want["far"]["far_rank1_narrow"], 1e-9)
print("\nALL INDEPENDENT CHECKS PASSED" if ok else "\nSOME INDEPENDENT CHECKS FAILED")
sys.exit(0 if ok else 1)
