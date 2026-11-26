"""For every curated stool sample of the subset studies: curated run(s) -> BioSample (ENA), then BioSample -> ALL its runs.
Output biosample_runs.json: {run: {biosample, study, alias, strategy}} for every run of every curated BioSample, and
curated_sample_biosample.json: {(study|sample_id): biosample}."""
import csv,json,gzip,urllib.request,urllib.parse,io,time,collections
SRC='local/crocodeel_valles/matrices'
P=json.load(open('pairs.json')); studies=sorted({p['study'] for p in P if p['v1_acc'] not in('NA','') })
runs=lambda s:[a for a in s.split(';') if a not in('','NA')]
cur=[]   # (study, subject, sample_id, run)
for st in studies:
    for r in csv.DictReader(open(f'{SRC}/{st}.metadata.tsv'),delimiter='\t'):
        if r['body_site']=='stool' and r['subject_id']:
            for a in runs(r['NCBI_accession']): cur.append((st,r['subject_id'],r['sample_id'],a))
print(len(cur),'curated runs in',len(studies),'studies',flush=True)
def ena(field,values,fields):
    out=[]
    for i in range(0,len(values),60):
        q=' OR '.join(f'{field}="{v}"' for v in values[i:i+60])
        url='https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&limit=0&fields='+fields+'&query='+urllib.parse.quote(q)
        for att in range(4):
            try: txt=urllib.request.urlopen(url,timeout=120).read().decode(); break
            except Exception as e: time.sleep(4*(att+1)); txt=None
        if txt is None: print('fail batch',i); continue
        out+=list(csv.DictReader(io.StringIO(txt),delimiter='\t'))
        time.sleep(0.3)
    return out
r1=ena('run_accession',sorted({c[3] for c in cur}),'run_accession,sample_accession,study_accession')
run2bs={r['run_accession']:r['sample_accession'] for r in r1}
print('curated runs resolved to a BioSample:',len(run2bs),'/',len({c[3] for c in cur}),flush=True)
bss=sorted(set(run2bs.values()))
print(len(bss),'distinct BioSamples',flush=True)
r2=ena('sample_accession',bss,'run_accession,sample_accession,study_accession,sample_alias,library_strategy')
allruns={r['run_accession']:{'biosample':r['sample_accession'],'study':r['study_accession'],'alias':r['sample_alias'][:40],'strategy':r['library_strategy']} for r in r2}
print('all runs of those BioSamples:',len(allruns),' (curated runs:',len(run2bs),')')
extra=[a for a in allruns if a not in run2bs]
print('runs of the same BioSamples that curatedMetagenomicData does not list:',len(extra))
cs={f'{st}|{sid}':run2bs.get(a) for st,sub,sid,a in cur if run2bs.get(a)}
json.dump(allruns,open('biosample_runs.json','w')); json.dump(cs,open('curated_sample_biosample.json','w'))
coll=collections.defaultdict(set)
for st,sub,sid,a in cur:
    b=run2bs.get(a)
    if b: coll[b].add((st,sub,sid))
print('BioSamples shared by >1 curated sample:',sum(1 for v in coll.values() if len(v)>1))
