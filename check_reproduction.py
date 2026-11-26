#!/usr/bin/env python3
"""Compares what reproduce.sh wrote (output/) with the files used in the manuscript (expected/).

    python3 -I check_reproduction.py [OUTPUT_DIR [EXPECTED_DIR]]

PNG   identical pixels = PASS. A few grey levels of difference on a handful of pixels (another version of matplotlib, FreeType or the C library moves anti-aliasing) = NEAR, with the figures printed.
JSON  the same numbers = PASS (integers and floats compared exactly); floats equal to 3 decimals, or to 1e-9 (relative) for the very large ones, = NEAR (another Python or numpy moves
      the last digits of a fit); anything else = FAIL, with the first differences.
TeX   supplementary.tex must be identical byte for byte.
Exit status 0 when nothing FAILs (PASS and NEAR both count as reproduced).
"""
import json, math, os, sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "output"
EXP = sys.argv[2] if len(sys.argv) > 2 else "expected"
PAIRS = [("figure1.png", "figure1.png"), ("supplementary_figure1.png", "supplementary_figure1.png"), ("supplementary_figure2.png", "supplementary_figure2.png"),
         ("supplementary_figure3.png", "supplementary_figure3.png"), ("supplementary_figure4.png", "supplementary_figure4.png"), ("supplementary_figure5.png", "supplementary_figure5.png"),
         ("supplementary_figure6.png", "supplementary_figure6.png"),
         ("figure1_quadratic.json", "figure1_quadratic.json"), ("figure1_open_set.json", "figure1_open_set.json"), ("supplementary_figure3_values.json", "supplementary_figure3_values.json"),
         ("supplementary_figure4_stats.json", "supplementary_figure4_stats.json"), ("supplementary_figure5_values.json", "supplementary_figure5_values.json"),
         ("supplementary_figure6_values.json", "supplementary_figure6_values.json"),
         ("numbers_in_text.json", "numbers_in_text.json"), ("manuscript/supplementary.tex", "supplementary.tex")]
results = []


ABS_TOL, REL_TOL = 1e-3, 1e-9                                                       # 3 decimals; 1e-9 (relative) for the very large numbers, such as pool sizes of 1e100 and more


def diff_json(a, b, path="", out=None):
    """Returns the differences: a string for each hard one, a tuple for each float equal to 3 decimals (or to REL_TOL relative)."""
    out = [] if out is None else out
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            out.append(f"{path}: keys differ {sorted(set(a) ^ set(b))[:5]}")
        for k in a:
            if k in b:
                diff_json(a[k], b[k], f"{path}/{k}", out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: length {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            diff_json(x, y, f"{path}[{i}]", out)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
        if a != b:
            if isinstance(a, float) or isinstance(b, float):
                if math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL):
                    out.append(("near", path, a, b))
                    return out
            out.append(f"{path}: {a!r} vs {b!r}")
    elif a != b:
        out.append(f"{path}: {a!r} vs {b!r}")
    return out


def compare_png(a, b):
    from PIL import Image
    import numpy as np
    A = np.asarray(Image.open(a).convert("RGBA"), dtype=np.int16); B = np.asarray(Image.open(b).convert("RGBA"), dtype=np.int16)
    if A.shape != B.shape:
        return "FAIL", f"size {A.shape[1]}x{A.shape[0]} vs {B.shape[1]}x{B.shape[0]}"
    d = np.abs(A - B).max(axis=2)
    n = int((d > 0).sum())
    if n == 0:
        return "PASS", f"identical pixels ({A.shape[1]}x{A.shape[0]})"
    frac = n / d.size
    if d.max() <= 24 and frac <= 0.005:
        return "NEAR", f"{n} pixels differ ({100 * frac:.3f}%), largest difference {int(d.max())} of 255"
    return "FAIL", f"{n} pixels differ ({100 * frac:.2f}%), largest difference {int(d.max())} of 255"


for o, e in PAIRS:
    po, pe = os.path.join(OUT, o), os.path.join(EXP, e)
    if not os.path.exists(po):
        results.append((o, "FAIL", "not written")); continue
    if not os.path.exists(pe):
        results.append((o, "FAIL", "no expected file")); continue
    if o.endswith(".png"):
        results.append((o, *compare_png(po, pe)))
    elif o.endswith(".json"):
        d = diff_json(json.load(open(po)), json.load(open(pe)))
        hard = [x for x in d if isinstance(x, str)]
        near = [x for x in d if not isinstance(x, str)]
        if hard:
            results.append((o, "FAIL", "; ".join(hard[:4]) + (f" (+{len(hard) - 4} more)" if len(hard) > 4 else "")))
        elif near:
            results.append((o, "NEAR", f"{len(near)} floats equal to 3 decimals (or to 1e-9 relative), e.g. {near[0][1]}: {near[0][2]!r} vs {near[0][3]!r}"))
        else:
            results.append((o, "PASS", "same numbers"))
    else:
        same = open(po, "rb").read() == open(pe, "rb").read()
        if same:
            results.append((o, "PASS", "identical bytes"))
        else:
            a, b = open(po, encoding="utf-8").read().splitlines(), open(pe, encoding="utf-8").read().splitlines()
            k = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
            results.append((o, "FAIL", f"first difference at line {k + 1}: {a[k][:100] if k < len(a) else '(end)'!r} vs {b[k][:100] if k < len(b) else '(end)'!r}"))

w = max(len(r[0]) for r in results)
for name, status, note in results:
    print(f"{status:5s} {name:{w}s}  {note}")
bad = [r for r in results if r[1] == "FAIL"]
print()
print(f"{len(results) - len(bad)} of {len(results)} files reproduced" + ("" if not bad else f"; {len(bad)} FAILED"))
print("ALL CHECKS PASSED" if not bad else "SOME CHECKS FAILED")
sys.exit(1 if bad else 0)
