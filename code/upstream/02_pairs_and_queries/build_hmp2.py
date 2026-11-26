"""HMP IBD (PRJNA398089) and HMP 2019 T2D (PRJNA281230): first/last stool metagenome per participant, from the deposits themselves.
Participant and date come from the ENA/BioSample record (sample_alias, collection_date), not from curatedMetagenomicData.
For HMP IBD the sample-level link to the curated profiles is made later, by read counts (hmp_curated_ids.py, build_hmp2_v3.py)."""
import urllib.request,urllib.parse,csv,io,re,json,collections,datetime,gzip
idx=set(gzip.open('logan_index_ids.txt.gz','rt').read().split('\n'))
def tsv(url): return list(csv.DictReader(io.StringIO(urllib.request.urlopen(url,timeout=180).read().decode()),delimiter='\t'))
F='run_accession,sample_accession,sample_alias,collection_date,library_strategy,library_source,read_count,first_public'
ibd=tsv(f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA398089&result=read_run&fields={F}&format=tsv&limit=0')
t2d=tsv(f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA281230&result=read_run&fields={F}&format=tsv&limit=0')
print('PRJNA398089 runs',len(ibd),'| PRJNA281230 runs',len(t2d))
def build(rows,study,subj_of,keep):
    by=collections.defaultdict(list)
    for r in rows:
        if r['library_strategy']!='WGS' or r['library_source']!='METAGENOMIC' or not keep(r): continue
        s=subj_of(r)
        if not s or not re.fullmatch(r'\d{4}-\d\d-\d\d',r['collection_date'] or ''): continue
        by[s].append(r)
    runs={}          # run -> (study, subject, sample=biosample)
    pairs=[]
    for s,rs in by.items():
        for r in rs: runs[r['run_accession']]=(study,s,r['sample_accession'])
        dates=sorted({r['collection_date'] for r in rs})
        if len(dates)<2: continue
        d1,d2=dates[0],dates[-1]
        a=[r for r in rs if r['collection_date']==d1]; b=[r for r in rs if r['collection_date']==d2]
        # choose the deepest run on each date that is indexed in the service
        pick=lambda L:max([r for r in L if r['run_accession'] in idx] or [None],key=lambda r:int(r['read_count'] or 0) if r else -1)
        a,b=pick(a),pick(b)
        pairs.append(dict(study=study,subject=s,v1_acc=a['run_accession'] if a else 'NA',v2_acc=b['run_accession'] if b else 'NA',
            v1_sample=a['sample_accession'] if a else None,v2_sample=b['sample_accession'] if b else None,
            days=(datetime.date.fromisoformat(d2)-datetime.date.fromisoformat(d1)).days,v1_reads=a['read_count'] if a else '0',v2_reads=b['read_count'] if b else '0',n_samples=len(rs)))
    return pairs,runs
ibd_pairs,ibd_runs=build(ibd,'HMP_2019_ibdmdb',lambda r:(re.fullmatch(r'([A-Z]\d{4})C\d+',r['sample_alias']) or [None,None])[1] if re.fullmatch(r'([A-Z]\d{4})C\d+',r['sample_alias']) else None,lambda r:True)
t2d_pairs,t2d_runs=build(t2d,'HMP_2019_t2d',lambda r:(re.search(r'_(Z[A-Z0-9]{6})-',r['sample_alias']) or [None,None])[1] if re.search(r'_(Z[A-Z0-9]{6})-',r['sample_alias']) else None,lambda r:'_ST_' in r['sample_alias'])
for n,p,r in (('HMP IBD',ibd_pairs,ibd_runs),('HMP 2019 T2D',t2d_pairs,t2d_runs)):
    ok=[x for x in p if x['v1_acc']!='NA' and x['v2_acc']!='NA']
    print(f'{n}: participants with >=2 collection dates {len(p)}; both ends indexed in the service {len(ok)}; runs mapped {len(r)}; median gap days {sorted(x["days"] for x in ok)[len(ok)//2] if ok else None}')
json.dump(ibd_pairs+t2d_pairs,open('pairs_hmp2.json','w')); json.dump({**ibd_runs,**t2d_runs},open('runs_hmp2.json','w'))
