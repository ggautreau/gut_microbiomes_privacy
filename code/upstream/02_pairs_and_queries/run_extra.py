"""Second pass: (1) curated pairs whose first v2 run is not indexed but another run of the sample is; (2) HMP IBD pairs."""
import json,gzip,os,time,urllib.request
idx=set(gzip.open('logan_index_ids.txt.gz','rt').read().split('\n'))
runs=lambda s:[a for a in s.split(';') if a not in('','NA')]
todo=[]
for p in json.load(open('pairs.json')):
    r1,r2=runs(p['v1_acc']),runs(p['v2_acc'])
    if not r1 or not r2 or not any(a in idx for a in r1) or not any(a in idx for a in r2): continue
    if r2[0] in idx: continue                      # first run is indexed: already queried in pass 1
    todo.append((p['study'],p['subject'],[a for a in r2 if a in idx][0]))
n_cur=len(todo)
for p in json.load(open('pairs_hmp2.json')):
    if p['study']=='HMP_2019_ibdmdb' and p['v1_acc']!='NA' and p['v2_acc']!='NA': todo.append((p['study'],p['subject'],p['v2_acc']))
done=set()
if os.path.exists('results_pass2.jsonl'):
    for l in open('results_pass2.jsonl'): done.add(json.loads(l)['run'])
print(n_cur,'curated re-queries +',len(todo)-n_cur,'HMP IBD queries;',len(done),'already done',flush=True)
out=open('results_pass2.jsonl','a'); n=0
for study,subj,run in todo:
    if run in done: continue
    body=json.dumps({'sample_id':run,'sync':True,'topk':5000}).encode()
    req=urllib.request.Request('https://similarity.logan-search.org/search',body,{'Content-Type':'application/json'})
    d=None
    for att in range(3):
        try: d=json.load(urllib.request.urlopen(req,timeout=600)); break
        except Exception as e: err=str(e); time.sleep(5*(att+1))
    if not d or d.get('status')!='completed': out.write(json.dumps({'run':run,'study':study,'subject':subj,'error':(d or {}).get('error') or err})+'\n')
    else: out.write(json.dumps({'run':run,'study':study,'subject':subj,'results':{k:round(v['score'],5) for k,v in d['results'].items()}})+'\n')
    out.flush(); n+=1; time.sleep(1)
    if n%20==0: print(n,'done',flush=True)
print('finished',n,flush=True)
