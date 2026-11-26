"""Per-subject rows for Fig. 1: the rank of the donor's other sample under the definitions below, and the age class.

    python3 -I rows_with_partner.py DATA_DIR OUT_JSON [--hmp LINKS.json]

DATA_DIR holds pairs.json, pairs_hmp2.json, runs_hmp2.json, runs_hmp2_dates.json, biosample_runs.json,
curated_sample_biosample.json, logan_index_ids.txt.gz, logan_gut_wgs_narrow_ids.txt.gz (the HGM set), results.jsonl, results_pass2.jsonl and matrices/<study>.metadata.tsv; logan_gut_wgs_broad_ids.txt.gz and mehta_sibling_runs.json are read when present.

For each answered query (the last-visit sample of one subject) the hits are classified by identity (study, subject, sample):
  `partner_any` : highest similarity over every OTHER sample of the same donor                       -> rank_any   (the headline definition)
  `partner_v1`  : similarity of the FIRST-visit sample only (the sample furthest in time)             -> rank_v1
  `partner_far` : highest similarity over other samples of the donor collected >= 30 days from the query -> rank_far
  In each variant the donor's remaining samples are neither the other sample nor runs of other people; runs of other people are every hit that is not the donor
  (other donors, other studies, unmapped runs). rank = 1 + number of runs of other people with a strictly higher similarity; None = not in the 5000.
"""
import csv, collections, datetime, gzip, json, os, sys
import numpy as np

DATA = sys.argv[1]
OUT = sys.argv[2]
LINKS = json.load(open(sys.argv[sys.argv.index('--hmp') + 1])) if '--hmp' in sys.argv else {}   # option --hmp LINKS.json (hmp_sibling_links.py): run -> [study, subject], a run of
                                                                    # another deposit that public records give to the donor, neither its other sample nor another person
DUP = 0.99
runs = lambda s: [a for a in s.split(';') if a not in ('', 'NA')]
P = json.load(open(f'{DATA}/pairs.json')) + [p for p in json.load(open(f'{DATA}/pairs_hmp2.json')) if p['study'] == 'HMP_2019_ibdmdb']
idx = set(gzip.open(f'{DATA}/logan_index_ids.txt.gz', 'rt').read().split('\n'))
# the HGM set competes (build_stock.py), and rank_any counts every run of the index; broad: any gut / feces / human metagenome label
NARROW = set(gzip.open(f'{DATA}/logan_gut_wgs_narrow_ids.txt.gz', 'rt').read().split('\n'))
BROAD = set(gzip.open(f'{DATA}/logan_gut_wgs_broad_ids.txt.gz', 'rt').read().split('\n')) if os.path.exists(f'{DATA}/logan_gut_wgs_broad_ids.txt.gz') else set()
subset = {(p['study'], p['subject']): p for p in P
          if runs(p['v1_acc']) and runs(p['v2_acc']) and any(a in idx for a in runs(p['v1_acc'])) and any(a in idx for a in runs(p['v2_acc']))}
res = {}
for fn in ('results.jsonl', 'results_pass2.jsonl'):
    path = f'{DATA}/{fn}'
    if os.path.exists(path):
        for l in open(path):
            try:
                q = json.loads(l)
            except Exception:
                continue
            if 'results' in q and (q['study'], q['subject']) not in res:
                res[(q['study'], q['subject'])] = q

# identity of every run: (study, subject, sample)
allruns = json.load(open(f'{DATA}/biosample_runs.json'))
cs = json.load(open(f'{DATA}/curated_sample_biosample.json'))
MAT = f'{DATA}/matrices'
meta_subj, day_of, ages = {}, {}, collections.defaultdict(list)
studies = {s for s, _ in subset}
for study in studies:
    for r in csv.DictReader(open(f'{MAT}/{study}.metadata.tsv'), delimiter='\t'):
        if r['body_site'] != 'stool' or not r['subject_id']:
            continue
        meta_subj[(study, r['sample_id'])] = r['subject_id']
        try:
            day_of[(study, r['sample_id'])] = float(r['days_from_first_collection'])
        except ValueError:
            pass
        try:
            ages[(study, r['subject_id'])].append(float(r['age']))
        except ValueError:
            pass
bs2s = {}
for k, bs in cs.items():
    st, sid = k.split('|', 1)
    if (st, sid) in meta_subj:
        bs2s[bs] = (st, meta_subj[(st, sid)], sid)
acc = {run: bs2s[v['biosample']] for run, v in allruns.items() if v['biosample'] in bs2s}
ibd_runs = {run: tuple(v) for run, v in json.load(open(f'{DATA}/runs_hmp2.json')).items() if v[0] == 'HMP_2019_ibdmdb'}
acc.update(ibd_runs)
# runs of the Mehta project that the curated BioSamples do not name (the DNA and the RNA library of one stool sample sit in two BioSamples), identified by their ENA alias (mehta_siblings.py)
if os.path.exists(f'{DATA}/mehta_sibling_runs.json'):
    acc.update({run: tuple(v) for run, v in json.load(open(f'{DATA}/mehta_sibling_runs.json')).items() if run not in acc})
# day of each sample: curated days_from_first_collection; HMP IBD collection dates (days from that donor's first date)
ibd_dates = json.load(open(f'{DATA}/runs_hmp2_dates.json'))
first_date = {}
for run, (st, subj, bs) in ibd_runs.items():
    if run in ibd_dates:
        d = datetime.date.fromisoformat(ibd_dates[run]).toordinal()
        first_date[subj] = min(first_date.get(subj, d), d)
for run, (st, subj, bs) in ibd_runs.items():
    if run in ibd_dates:
        day_of[(st, bs)] = datetime.date.fromisoformat(ibd_dates[run]).toordinal() - first_date[subj]


def adult_class(study, subj):
    a = ages.get((study, subj), [])
    if a:
        return 'adult' if min(a) >= 18 else 'under18_or_mixed'
    return 'adult_by_design' if study == 'MehtaRS_2018' else 'unknown'


indexed_samples = collections.defaultdict(set)          # (study, subject) -> samples having at least one indexed run
for run, m in acc.items():
    if run in idx:
        indexed_samples[(m[0], m[1])].add(m[2])
rows = []
for key, p in subset.items():
    q = res.get(key)
    if q is None:
        continue
    study, subj = key
    v1s, v2s = p['v1_sample'], p['v2_sample']
    d2 = day_of.get((study, v2s))
    imp, imp_acc, same = [], [], []         # similarities and accessions of the runs of other people ; (similarity, sample) of the donor's other samples
    same_acc = []                           # accession of each entry of `same`
    rep = 0.0
    for a, sc in q['results'].items():
        m = acc.get(a)
        if m and m[0] == study and m[2] == v2s and a != q['run']:
            rep = max(rep, sc)
        if sc >= DUP:
            continue
        if LINKS and LINKS.get(a) == [study, subj]:
            continue
        if m and m[0] == study and m[1] == subj:
            if m[2] != v2s:
                same.append((sc, m[2]))
                same_acc.append(a)
        else:
            imp.append(sc)
            imp_acc.append(a)
    imp = np.array(imp)
    in_narrow = np.array([a in NARROW for a in imp_acc], dtype=bool)
    in_broad = np.array([a in BROAD for a in imp_acc], dtype=bool)

    def rank_for(partner, mask=None):
        if partner <= 0:
            return None
        return int((imp > partner).sum()) + 1 if mask is None else int(((imp > partner) & mask).sum()) + 1

    p_any = max([s for s, _ in same], default=0.0)
    # search restricted to a labelled population S: the donor's other sample must itself be in S, and only S competes
    def in_set_rank(S, mask):
        partner = max([sc for (sc, _), a in zip(same, same_acc) if a in S], default=0.0)
        return int((imp > partner) .__and__(mask).sum()) + 1 if partner > 0 else None, partner > 0
    p_v1 = max([s for s, x in same if x == v1s], default=0.0)
    eligible = {x for x in indexed_samples[(study, subj)] if x != v2s and d2 is not None and day_of.get((study, x)) is not None and abs(day_of[(study, x)] - d2) >= 30}
    p_far = max([s for s, x in same if x in eligible], default=0.0)
    best_any = max(same, default=(0.0, None))[1]
    if len(imp):
        ti = int(np.argmax(imp)); tm = acc.get(imp_acc[ti])
        top_kind = 'unmapped' if tm is None else ('same_cohort_other_person' if tm[0] == study else 'other_cohort')
    else:
        top_kind = None
    best_far = max([(s, x) for s, x in same if x in eligible], default=(0.0, None))[1]
    rows.append(dict(v1s=v1s, v2s=v2s, best_any=best_any, best_far=best_far, top_kind=top_kind, eligible=sorted(eligible), study=study, subject=subj, rank_any=rank_for(p_any), rank_v1=rank_for(p_v1), rank_far=rank_for(p_far),
                     rank_any_narrow=rank_for(p_any, in_narrow), rank_in_narrow=in_set_rank(NARROW, in_narrow)[0], rank_in_broad=in_set_rank(BROAD, in_broad)[0],
                     partner_in_narrow=in_set_rank(NARROW, in_narrow)[1], partner_in_broad=in_set_rank(BROAD, in_broad)[1], rank_any_broad=rank_for(p_any, in_broad), rank_v1_narrow=rank_for(p_v1, in_narrow), rank_far_narrow=rank_for(p_far, in_narrow),
                     far_possible=bool(eligible), n_far_eligible=len(eligible),
                     partner_any=p_any, partner_v1=p_v1, partner_far=p_far, rep=rep,
                     max_imp=float(imp.max()) if len(imp) else 0.0, days=float(p['days']),
                     age_class=adult_class(study, subj), age_min=min(ages[(study, subj)]) if ages.get((study, subj)) else None,
                     n_same_samples_hit=len({x for _, x in same}),
                     q_lib='RNA' if '.rna' in allruns.get(q['run'], {}).get('alias', '').lower() else 'DNA'))
json.dump(rows, open(OUT, 'w'))
n_sub = len(subset)
print(f'subset {n_sub}; answered {len(rows)}; rows written to {OUT}')
