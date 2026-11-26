"""HMP IBD identity, version 2. A sample is a participant-visit (alias <participant>C<n>, with an optional library suffix
_MGX, _MVX or _TR_MGX that only says how the same stool was sequenced). Every WGS metagenomic run of PRJNA398089 is mapped,
not only the 1000 whose alias has no suffix. Writes runs_hmp2.json (run -> [study, participant, visit key]),
runs_hmp2_dates.json (run -> collection date) and pairs_hmp2.json with the visit keys of v1 and v2."""
import csv, collections, json, re
rows = [r for r in csv.DictReader(open('ibdmdb_runs_ena.tsv'), delimiter='\t') if r['library_strategy'] == 'WGS' and r['library_source'] == 'METAGENOMIC']
pat = re.compile(r'^([A-Z]\d{4})C(\d+)(?:_(MGX|MVX|TR_MGX))?$')
runs, dates, keydates, suffix_n = {}, {}, collections.defaultdict(set), collections.Counter()
for r in rows:
    m = pat.match(r['sample_alias'])
    if not m:
        continue
    part, visit, suf = m.groups()
    key = f'{part}C{visit}'
    runs[r['run_accession']] = ['HMP_2019_ibdmdb', part, key]
    dates[r['run_accession']] = r['collection_date']
    keydates[key].add(r['collection_date'])
    suffix_n[suf or 'plain'] += 1
print(len(rows), 'WGS metagenomic runs;', len(runs), 'mapped;', dict(suffix_n))
conflict = {k: v for k, v in keydates.items() if len(v) > 1}
print('visit keys:', len(keydates), '| keys whose runs carry different collection dates:', len(conflict), list(conflict.items())[:3])
old = json.load(open('pairs_hmp2.json'))
new = []
for p in old:
    if p['study'] != 'HMP_2019_ibdmdb':
        new.append(p); continue
    p = dict(p)
    p['v1_sample'] = runs[p['v1_acc']][2] if p['v1_acc'] in runs else p['v1_sample']
    p['v2_sample'] = runs[p['v2_acc']][2] if p['v2_acc'] in runs else p['v2_sample']
    new.append(p)
json.dump(new, open('pairs_hmp2.json', 'w'))
json.dump(runs, open('runs_hmp2.json', 'w'))
json.dump(dates, open('runs_hmp2_dates.json', 'w'))
print('pairs updated:', sum(1 for p in new if p['study'] == 'HMP_2019_ibdmdb'))
