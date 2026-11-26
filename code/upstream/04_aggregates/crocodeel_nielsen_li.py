#!/usr/bin/env python3
"""CroCoDeEL on the species tables of Nielsen 2014 and Li 2014 merged, to see the events between a Nielsen sample and a Li sample (Supplementary Fig. 4C).

    python3 -I crocodeel_nielsen_li.py OUTDIR

The two cohorts carry the same subject identifiers for 44 patients (Supplementary Note 1) and CroCoDeEL had been run on each cohort alone, so it could not report an event between them. The merged table has one column per sample
(names prefixed N| or L| so that no name is shared) and one row per species of either table (a species missing from one table counts as absent: 0). CroCoDeEL 1.3.1 runs with its default settings,
as for the other cohorts. Writes OUTDIR/nielsen_li.species_relab.tsv and OUTDIR/nielsen_li.conta.tsv.
"""
import os, subprocess, sys

OUT = sys.argv[1]
CRO = "local/crocodeel_valles/matrices"
BIN = os.environ.get("CROCODEEL", "crocodeel")                      # the crocodeel 1.3.1 executable
os.makedirs(OUT, exist_ok=True)
tabs = {}
for tag, st in (("N", "NielsenHB_2014"), ("L", "LiJ_2014")):
    rows = [l.rstrip("\n").split("\t") for l in open(f"{CRO}/{st}.species_relab.tsv") if l.strip()]
    tabs[tag] = (rows[0][1:], {r[0]: r[1:] for r in rows[1:]})
species = sorted(set(tabs["N"][1]) | set(tabs["L"][1]))
cols = [f"N|{s}" for s in tabs["N"][0]] + [f"L|{s}" for s in tabs["L"][0]]
table = os.path.join(OUT, "nielsen_li.species_relab.tsv")
with open(table, "w") as f:
    f.write("species\t" + "\t".join(cols) + "\n")
    for sp in species:
        f.write(sp + "\t" + "\t".join(sum((tabs[t][1].get(sp) or ["0"] * len(tabs[t][0]) for t in ("N", "L")), [])) + "\n")
print(len(cols), "samples,", len(species), "species")
subprocess.run([BIN, "search_conta", "-s", table, "-c", os.path.join(OUT, "nielsen_li.conta.tsv")], check=True)
