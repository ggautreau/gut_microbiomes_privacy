"""Per-donor lists of the similarities of the runs of other people in the HGM set, for the thinning test (extrapolation_tables.py) and panel H of Figure 1 (make_open_set_scores.py).

    python3 -I extract_lists.py DATA_DIR ROWS_JSON OUT_NPZ [--hmp LINKS.json]

Same identities and exclusions as rows_with_partner.py (`partner_any`, runs of other people = every hit that is not the donor,
quasi-exact re-deposits >= 0.99 removed).  For each donor of ROWS_JSON that is an adult with a DNA query, stores
  `partner`   highest similarity over the donor's other samples (0 = other sample not in the 5,000 returned)
  `imp_nar`   similarities of the runs of other people that belong to the HGM set (build_stock.py), descending, float64 (the service returns 5 decimals; float32 would
              turn exact ties with the other sample into spurious differences, so no float32 anywhere)
and checks that 1 + #(`imp_nar` > `partner`) equals rank_any_narrow of ROWS_JSON.
"""
import csv, collections, datetime, gzip, json, os, sys
import numpy as np

DATA, ROWS, OUT = sys.argv[1:4]
LINKS = json.load(open(sys.argv[sys.argv.index('--hmp') + 1])) if '--hmp' in sys.argv else {}   # option --hmp LINKS.json (hmp_sibling_links.py), as in rows_with_partner.py
DUP = 0.99
runs = lambda s: [a for a in s.split(';') if a not in ('', 'NA')]
P = json.load(open(f'{DATA}/pairs.json')) + [p for p in json.load(open(f'{DATA}/pairs_hmp2.json')) if p['study'] == 'HMP_2019_ibdmdb']
idx = set(gzip.open(f'{DATA}/logan_index_ids.txt.gz', 'rt').read().split('\n'))
NARROW = set(gzip.open(f'{DATA}/logan_gut_wgs_narrow_ids.txt.gz', 'rt').read().split('\n'))     # the HGM set (build_stock.py)
subset = {(p['study'], p['subject']): p for p in P
          if runs(p['v1_acc']) and runs(p['v2_acc']) and any(a in idx for a in runs(p['v1_acc'])) and any(a in idx for a in runs(p['v2_acc']))}
res = {}
for fn in ('results.jsonl', 'results_pass2.jsonl'):
    for l in open(f'{DATA}/{fn}'):
        try:
            q = json.loads(l)
        except Exception:
            continue
        if 'results' in q and (q['study'], q['subject']) not in res:
            res[(q['study'], q['subject'])] = q
allruns = json.load(open(f'{DATA}/biosample_runs.json'))
cs = json.load(open(f'{DATA}/curated_sample_biosample.json'))
MAT = f'{DATA}/matrices'
meta_subj = {}
for study in {s for s, _ in subset}:
    for r in csv.DictReader(open(f'{MAT}/{study}.metadata.tsv'), delimiter='\t'):
        if r['body_site'] != 'stool' or not r['subject_id']:
            continue
        meta_subj[(study, r['sample_id'])] = r['subject_id']
bs2s = {}
for k, bs in cs.items():
    st, sid = k.split('|', 1)
    if (st, sid) in meta_subj:
        bs2s[bs] = (st, meta_subj[(st, sid)], sid)
acc = {run: bs2s[v['biosample']] for run, v in allruns.items() if v['biosample'] in bs2s}
acc.update({run: tuple(v) for run, v in json.load(open(f'{DATA}/runs_hmp2.json')).items() if v[0] == 'HMP_2019_ibdmdb'})
# runs of the Mehta project that the curated BioSamples do not name (the DNA and the RNA library of one stool sample sit in two BioSamples), identified by their ENA alias (mehta_siblings.py)
if os.path.exists(f'{DATA}/mehta_sibling_runs.json'):
    acc.update({run: tuple(v) for run, v in json.load(open(f'{DATA}/mehta_sibling_runs.json')).items() if run not in acc})

rows = {(r['study'], r['subject']): r for r in json.load(open(ROWS))}
keep = [k for k, r in rows.items() if r['age_class'] in ('adult', 'adult_by_design') and r['q_lib'] == 'DNA']
out = dict(study=[], subject=[], partner=[], n_list=[], floor=[], k_rows=[], imp_nar=[])
bad = 0
for key in keep:
    p, q = subset[key], res[key]
    study, subj = key
    v2s = p['v2_sample']
    imp, imp_acc, same = [], [], []
    for a, sc in q['results'].items():
        if sc >= DUP:
            continue
        if LINKS and LINKS.get(a) == [study, subj]:
            continue
        m = acc.get(a)
        if m and m[0] == study and m[1] == subj:
            if m[2] != v2s:
                same.append(sc)
        else:
            imp.append(sc); imp_acc.append(a)
    imp = np.array(imp); nar = np.array([a in NARROW for a in imp_acc], dtype=bool)
    partner = max(same, default=0.0)
    k_calc = int((imp[nar] > partner).sum()) if partner > 0 else -1
    r = rows[key]
    if (r['rank_any_narrow'] is None) != (partner <= 0) or (partner > 0 and k_calc != r['rank_any_narrow'] - 1) or abs(partner - r['partner_any']) > 1e-9:
        bad += 1
    out['study'].append(study); out['subject'].append(subj); out['partner'].append(partner)
    out['n_list'].append(len(q['results'])); out['floor'].append(min(q['results'].values()))
    out['k_rows'].append(-1 if r['rank_any_narrow'] is None else r['rank_any_narrow'] - 1)
    out['imp_nar'].append(np.sort(imp[nar])[::-1].astype(np.float64))
print('donors', len(keep), '| mismatches with ROWS_JSON:', bad)
np.savez_compressed(OUT, study=np.array(out['study']), subject=np.array(out['subject']), partner=np.array(out['partner']),
                    n_list=np.array(out['n_list']), floor=np.array(out['floor']), k_rows=np.array(out['k_rows']),
                    imp_nar=np.array(out['imp_nar'], dtype=object))
