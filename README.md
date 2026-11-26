# Reproducing "Adult gut microbiomes single out individuals at population scale"

This reproduction package regenerates every figure, every number of the main text and the Methods that depends on the donors, and the Supplementary Notes and Tables of the Brief Communication, and compares each file with the one used in the manuscript. Subjects and samples carry their identifiers of curatedMetagenomicData, and `data/accessions.tsv` gives the two runs of each donor. `local/` holds the answers of the public service and the public tables from which the upstream steps compute these inputs. Built on 16 November 2026 (483 donors, 7 cohorts).

## Run it

`bash reproduce.sh` (about 40 seconds on a machine with several cores and about a minute on one core, no network). It needs Python 3 with numpy, scipy, matplotlib and Pillow (`requirements.txt`; the versions of the build are in `environment.txt`). The output ends with:

```
SHA256SUMS: every file is as shipped
...
PASS  figure1.png                        identical pixels (4320x3954)
...
15 of 15 files reproduced
ALL CHECKS PASSED
...
ALL INDEPENDENT CHECKS PASSED
```

The seven figures come out with identical pixels, the seven files of numbers with the same numbers, and `supplementary.tex` byte for byte. On another machine a figure can differ by a few grey levels on a few pixels (another FreeType or matplotlib), and a number by its last digits (another Python or numpy moves the optimum of a fit, by about 1e-9 relative); `check_reproduction.py` then says `NEAR` when the numbers are equal to 3 decimals, and the file counts as reproduced.

`python3 -I rerun_upstream.py` (about 5 hours and 20 minutes on 32 cores, of which 4 hours and 55 minutes for `evt_model.py`; about 16 minutes with `--fast`, which leaves out the three slow extrapolation steps; no network) goes one step back. In a scratch copy it re-runs the upstream steps from the files of `local/` and compares what they write with the files of the package. It ends with `35 of 35 steps identical; 1 not run`; the step left out reads a table that is not here (below). As in `check_reproduction.py`, a file of numbers whose floats are equal only to 3 decimals counts as identical, and its step is marked `equal to 3 decimals` in the table.

## What the paper measures

The last stool sample of each of 483 adult donors of 7 longitudinal cohorts of curatedMetagenomicData queries the public Logan sketch index (Faure et al.; 4,791,208 runs), which returns the 5,000 nearest runs. A donor is *recovered* when no run of another person among the 291,881 human gut metagenomes of the index (the *HGM set*) is more similar to the query than the donor's other sample. 454 of 483 donors (94.0%) are recovered; 15 (3.1%) have their other sample outside the runs returned, which counts as a failure. The definitions are those of the Methods of the paper.

## What `reproduce.sh` reproduces

| Item of the paper | Script (`code/`) | Checked against (`expected/`) |
|---|---|---|
| Fig. 1, panels A to H | `figure1_brief_communication.py` | `figure1.png`, `figure1_quadratic.json`, `figure1_open_set.json` |
| Supplementary Fig. 1 (the choice of the cohorts, the donors and their samples) | `supplementary_figure1.py` | `supplementary_figure1.png` |
| Supplementary Fig. 2 (Nielsen and Li) | `supplementary_figure2.py` | `supplementary_figure2.png` |
| Supplementary Fig. 3 (human DNA, the 558 donors with Vincent) | `supplementary_figure3.py` | `supplementary_figure3.png`, `supplementary_figure3_values.json` |
| Supplementary Fig. 4 (CroCoDeEL events) | `supplementary_figure4.py` | `supplementary_figure4.png`, `supplementary_figure4_stats.json` |
| Supplementary Fig. 5 (the three models of Fig. 1B fitted to the measured curve, for five fit windows, and the extreme-value model for its five tails) | `supplementary_figure5.py` | `supplementary_figure5.png`, `supplementary_figure5_values.json` |
| Supplementary Fig. 6 (resampling intervals for the four models of Fig. 1B) | `supplementary_figure6.py` | `supplementary_figure6.png`, `supplementary_figure6_values.json` |
| Numbers of the main text and the Methods | `reproduce_numbers.py` | `numbers_in_text.json` |
| Supplementary Notes 1 to 4 and Tables 1 to 18 | `build_supplementary.py` | `supplementary.tex` |

`check_independent.py` then recomputes the headline numbers from `data/adult_ranks.csv` and the definitions alone, without scipy and without the code of the figures.

**What the checks prove, and what they do not.** `reproduce.sh` first checks `SHA256SUMS` (`SKIP_SHA=1` skips it for inputs changed on purpose), then recomputes every output from the inputs; no script reads `expected/` and no headline value is written into a script. The comparison is a regression test: this code, applied to these tables, writes the files of the manuscript. It does not recompute the ranks or the aggregates of `aggregates/`, which are typeset and drawn. `rerun_upstream.py` recomputes them from the files kept in `local/`, except `aggregates/human_dna_pool_summary.json`, a copy of `local/humancheck/all558/pool_final_summary.json` (`all558_pool.py`, not re-run). Neither queries the service again.

## Contents

```
reproduce.sh            the whole run, then the two checks
check_reproduction.py   compares output/ with expected/ (pixels of the PNG, numbers of the JSON, bytes of the TeX)
check_independent.py    recomputes the headline numbers from the definitions
rerun_upstream.py       re-runs the upstream steps from local/ and compares what they write with local/, data/, aggregates/ and reference/
code/                   the scripts of the figures and tables
code/upstream/          the scripts that turn public data into the tables below (README.md there)
local/                  what the upstream scripts read: public tables as downloaded, the answers of the service, the tables of each step
data/                   the inputs of reproduce.sh: the per-donor tables and accessions.tsv
aggregates/             outputs of the upstream steps: counts, medians and arrays, without names except the subject identifiers of the seven Costea donors of Supplementary Note 1 (costea_same_people.json)
reference/              the 172,075 runs labelled human gut or faecal, and the HGM set (291,881 runs)
expected/               the files of the manuscript
fonts/                  Liberation Sans, and the solid icons of Font Awesome 5 Free for the DNA symbol of Fig. 1A, with their licences
provenance/             MANIFEST_upstream.tsv (the files of local/), VERIFIED_LOCALLY.md (the upstream steps re-run)
requirements.txt, environment.txt, SHA256SUMS
```

* `data/adult_ranks.csv`: one row per donor, in the order of `data/donors.json`. `rank_any_narrow` is the rank of the other sample among the runs of the HGM set (empty when it is not among the 5,000 nearest runs); the other columns give the cohort, the days between the samples, the size of the queried run, the rank of Fig. 1E (`rank_first_narrow`) and that of the 30-day subset of the Methods (`rank_far_narrow`).
* `data/open_set_scores.csv`: one row per donor in the same order, for Fig. 1H. `partner` is the similarity of the other sample (0 when it is not among the 5,000 nearest runs), `best_other` that of the nearest run of other people in the HGM set.
* `data/donors.json`, `data/donors_left_out.json` (the donors of the 4 cohorts left out, for the samples studied of Supplementary Note 1), `data/donors_excluded.json` (the 14 donors of the cohorts kept that the rules of contamination and re-deposit remove, Supplementary Note 1), `data/query_gbp.json`, `data/crocodeel/` and `data/human_dna/`: the inputs of the supplementary figures and tables.
* `data/accessions.tsv`: the subject, the queried run and the run of the other sample of each of the 558 donors of the human-DNA check (the 497 donors of the 7 cohorts kept, 483 of them analysed, and the 61 of Vincent), with their BioProjects.

## What is not here

One table that the upstream steps read is not in the package. `hmp2_metadata_2018-08-20.csv` is the public metadata table of the HMP consortium (https://ibdmdb.org); `local/hmp_sibling_links.json`, which `hmp_sibling_links.py` writes from it, is here, and `provenance/MANIFEST_upstream.tsv` gives the sha256 of the table. The sketches of the human-DNA check (3.3 GB) are not here either, nor is the hypervector database of the service (18 GB unpacked; `code/upstream/03_ranks/evt_similarities.py` gives its address), from which `evt_similarities.py` computed `local/evt_similarities.npz`, the input of the extreme-value model of the similarity tail.

## From public data to these tables

`code/upstream/README.md` lists the upstream steps in order, with their commands. `rerun_upstream.py` re-runs 36 of them from `local/`, one of them only with the table named above, and they write the files of the package exactly; the steps that query the service or the ENA, run CroCoDeEL, sketch the runs of the human-DNA check or read the hypervector database are not repeated. `provenance/MANIFEST_upstream.tsv` gives the size and sha256 of each file of `local/`.

## Environment

Python 3.11.13, numpy 1.26.4, scipy 1.14.0, matplotlib 3.10.8, Pillow 11.0.0 (`environment.txt`). The random draws have fixed seeds; another numpy release could change the bands and the jitter of the figures, not the point values. Liberation Sans is shipped in `fonts/` (GNU GPL v2 with the font exception for Regular and Italic, SIL OFL 1.1 for Bold and Bold Italic).

## Licence, citation and contact

The code, the tables computed here and the text are released under the GNU General Public License, version 3 or any later version (`LICENSE`; `GPL-3.0-or-later`), unless a file says otherwise. The fonts keep their own licences. The tables of `local/` taken from curatedMetagenomicData (`sm.rda`, `cmd_meta/cmd_meta.tsv` and `crocodeel_valles/matrices/`) remain under the licence of that package, Artistic-2.0; `code/upstream/01_cohort_tables/extract_all.R` shows how the tables of `matrices/` are cut from it. The other files of `local/` copied from the ENA, the NCBI and the similarity service, the run accessions and the public metadata keep the terms of use of their sources. `.gitignore` keeps `output/` and the files that `prepare_local.sh` writes out of git.

**Repository.** https://github.com/ggautreau/gut_microbiomes_privacy.

**Citation.** The paper above (in preparation); metadata in `CITATION.cff`.

**Contact.** Guillaume Gautreau (corresponding author), guillaume.gautreau@inrae.fr.
