"""Query the similarity service of Faure et al. with the last-visit sample of each adult subject (first run of the sample),
top 5000 neighbours. Resumable: one JSON line per query in results.jsonl."""
import json,os,time,urllib.request
P=json.load(open('pairs.json'))
todo=[]
for p in P:
    if p['v1_acc'] in('','NA') or p['v2_acc'] in('','NA'): continue
    todo.append((p['study'],p['subject'],p['v2_acc'].split(';')[0]))
done=set()
if os.path.exists('results.jsonl'):
    for l in open('results.jsonl'): done.add(json.loads(l)['run'])
print(len(todo),'queries;',len(done),'already done',flush=True)
out=open('results.jsonl','a'); n=0; fails=0
for study,subj,run in todo:
    if run in done: continue
    body=json.dumps({'sample_id':run,'sync':True,'topk':5000}).encode()
    req=urllib.request.Request('https://similarity.logan-search.org/search',body,{'Content-Type':'application/json'})
    d=None
    for att in range(3):
        try: d=json.load(urllib.request.urlopen(req,timeout=600)); break
        except Exception as e:
            err=str(e); time.sleep(5*(att+1))
    if not d or d.get('status')!='completed':
        fails+=1; out.write(json.dumps({'run':run,'study':study,'subject':subj,'error':(d or {}).get('error') or err})+'\n'); out.flush()
    else:
        out.write(json.dumps({'run':run,'study':study,'subject':subj,'results':{k:round(v['score'],5) for k,v in d['results'].items()}})+'\n'); out.flush()
    n+=1; time.sleep(1)
    if n%25==0: print(n,'done, failures',fails,flush=True)
print('finished',n,'failures',fails,flush=True)
