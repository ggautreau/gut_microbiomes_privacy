#!/usr/bin/env python3
"""Extreme-value model of the similarity tail with the donors that the two rules of Supplementary Note 1 remove (last row of Supplementary Table 17). Aggregates only.

    python3 -I evt_removed_donors.py EVT_SIMILARITIES_497.npz ROWS_V4_497.json ROWS_V4.json EVT_SIMILARITIES.npz EVT_SUMMARY.json OUT.json

EVT_SIMILARITIES_497.npz  output of evt_similarities.py for the 497 donors of the cohorts kept (run on rows_v4_497.json and lists_v4_497.npz), in the order of ROWS_V4_497.json
ROWS_V4_497.json          these 497 donors; ROWS_V4.json, the 483 donors of the analysis (exclude_donors.py removes the 14 others)
EVT_SIMILARITIES.npz      output of evt_similarities.py for the 483 donors, the input of evt_model.py
EVT_SUMMARY.json          output of evt_model.py
OUT.json                  the model of Fig. 1B made again on the 497 donors: donors, donors recovered in the HGM set, recovery rate at 1e7 runs and at the population of the world

The model is that of evt_model.py (cells, integrate_sigma, fit and stay are those of that script), with the tail drawn in Fig. 1B and the shape free in sign. Before it writes, this
script checks that the 483 donors of the analysis have the same similarities in the two files, and that the model fitted on them gives the Beta parameters and the recovery
rates of EVT_SUMMARY.json again. It then fits the model on the 497 donors: the distribution of shapes is estimated on all of them, as it is on the donors kept when a cohort
is left out. Numbers are written with 10 decimals."""
import json, os, sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                   # python -I does not put the script's folder on the path
from evt_model import fit, stay

M = 291_881
WORLD = 8e9


if __name__ == "__main__":
    SIM497, ROWS497, ROWS, SIM, EVS, OUT = sys.argv[1:7]
    ev = json.load(open(EVS)); DRAWN = int(ev["drawn"]); LO, HI = ev["ranges"]["free"]; ref = ev["models"]["free"][str(DRAWN)]
    assert ev["hgm_runs"] == M
    d = np.load(SIM497, allow_pickle=True); d0 = np.load(SIM, allow_pickle=True)
    keys497 = [(r["study"], r["subject"]) for r in json.load(open(ROWS497))]; keys = [(r["study"], r["subject"]) for r in json.load(open(ROWS))]
    assert len(set(keys497)) == len(keys497) == len(d["k"]) and len(keys) == len(d0["k"]) == ev["donors"] and set(keys) < set(keys497), "the donors of the analysis are among those of the cohorts kept"
    at = {k_: i for i, k_ in enumerate(keys497)}; idx = np.array([at[k_] for k_ in keys])
    # ---------------------------------------------------------------- the donors of the analysis have the same similarities in the two files
    for name in ("top", "n", "s", "k", "study"):
        assert np.array_equal(d[name][idx], d0[name]), f"`{name}` of the donors of the analysis differs between the two files"
    TOP, n_all, s_all, k_all, study = d["top"].astype(np.float64), d["n"].astype(np.float64), d["s"], d["k"], d["study"]
    D = len(k_all); alive = k_all == 0; kept = np.zeros(D, bool); kept[idx] = True
    # ---------------------------------------------------------------- the model of evt_model.py is found again on the donors of the analysis
    f0 = fit(TOP[idx], n_all[idx], s_all[idx], DRAWN, LO, HI)
    w0 = 100 * float(stay(f0, alive[idx], WORLD - M).mean()); s0 = 100 * float(stay(f0, alive[idx], 1e7 - M).mean())
    assert abs(f0["a"] - ref["a"]) < 1e-9 and abs(f0["b"] - ref["b"]) < 1e-9 and abs(w0 - ref["at_world"]) < 1e-9 and abs(s0 - ref["at_1e7"]) < 1e-9 and int(alive[idx].sum()) == ev["recovered"], "the model of evt_model.py is not found again"
    # ---------------------------------------------------------------- the model on the donors of the cohorts kept, with those that the two rules remove
    integ = {k_: np.empty((D,) + v_.shape[1:], dtype=v_.dtype) for k_, v_ in f0["integ"].items()}
    for k_ in integ:
        integ[k_][idx] = f0["integ"][k_]
    fr = fit(TOP[~kept], n_all[~kept], s_all[~kept], DRAWN, LO, HI)
    for k_ in integ:
        integ[k_][~kept] = fr["integ"][k_]
    f = fit(TOP, n_all, s_all, DRAWN, LO, HI, prof=integ)
    rnd = lambda v: round(float(v), 10)
    out = dict(donors=int(D), donors_of_the_analysis=int(kept.sum()), removed=int((~kept).sum()), recovered=int(alive.sum()), removed_recovered=int(alive[~kept].sum()), hgm_runs=M, tail=DRAWN, range=[LO, HI],
               a=rnd(f["a"]), b=rnd(f["b"]), shape_mean=rnd(LO + (HI - LO) * f["a"] / (f["a"] + f["b"])), at_1e7=rnd(100 * stay(f, alive, 1e7 - M).mean()), at_world=rnd(100 * stay(f, alive, WORLD - M).mean()),
               analysis=dict(a=rnd(f0["a"]), b=rnd(f0["b"]), at_1e7=rnd(s0), at_world=rnd(w0)))
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"{D} donors, {int(alive.sum())} recovered in the HGM set ({int((~kept).sum())} removed by the two rules, {int(alive[~kept].sum())} of them recovered); tail of {DRAWN:,}: {out['at_1e7']:.1f} % at 1e7 and {out['at_world']:.1f} % at 8e9 "
          f"({w0:.1f} % for the {int(kept.sum())} donors of the analysis)")
    print("wrote", OUT)