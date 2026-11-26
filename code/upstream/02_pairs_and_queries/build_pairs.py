"""Rebuild, for the 14 adult cohorts of adult_cohorts.txt (Supplementary Table 2), the same first/last stool pair per subject as cohort_age.py,
with the run accession of each sample (curatedMetagenomicData NCBI_accession)."""
import csv,collections,json,re
SRC='local/crocodeel_valles/matrices'
adults=[l.strip() for l in open('adult_cohorts.txt') if l.strip()]        # the 14 cohorts of curatedMetagenomicData labelled as adults (cohort_selection.py, Supplementary Table 2)
day=lambda r:float(r.get('days_from_first_collection') or 0) if re.fullmatch(r'-?[\d.]+',r.get('days_from_first_collection') or '0') else 0.0
pairs=[]
for study in adults:
    md=[r for r in csv.DictReader(open(f'{SRC}/{study}.metadata.tsv'),delimiter='\t') if r['body_site']=='stool' and r['subject_id']]
    by=collections.defaultdict(list)
    for r in md: by[r['subject_id']].append(r)
    for subj,rows in by.items():
        if len(rows)<2: continue
        rows.sort(key=day)
        if rows[0]['sample_id']==rows[-1]['sample_id']: continue
        a,b=rows[0],rows[-1]
        pairs.append({'study':study,'subject':subj,'v1_sample':a['sample_id'],'v2_sample':b['sample_id'],
                      'v1_acc':a['NCBI_accession'],'v2_acc':b['NCBI_accession'],'days':day(b)-day(a),
                      'v1_reads':a['number_reads'],'v2_reads':b['number_reads'],'n_samples':len(rows)})
json.dump(pairs,open('pairs.json','w'),indent=0)
print(len(pairs),'subject pairs in',len(adults),'adult cohorts (expected 949 if all subjects have two distinct stool samples)')
c=collections.Counter(p['study'] for p in pairs); print(dict(c))
pref=collections.Counter(re.sub(r'\d+','#',p[k].split(';')[0]) for p in pairs for k in('v1_acc','v2_acc'))
print('accession prefixes:',pref)
print('multi-run samples:',sum(1 for p in pairs for k in('v1_acc','v2_acc') if ';' in p[k]))
print('empty accessions:',sum(1 for p in pairs for k in('v1_acc','v2_acc') if p[k] in('','NA')))
