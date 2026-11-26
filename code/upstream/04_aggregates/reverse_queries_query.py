import json,os,time,urllib.request
Q=json.load(open('local/nielsen_li_validation/reverse/prep.json'))['q']
done=set()
if os.path.exists('local/nielsen_li_validation/reverse/answers.jsonl'):
    for l in open('local/nielsen_li_validation/reverse/answers.jsonl'):
        j=json.loads(l)
        if 'results' in j: done.add(j['run'])
out=open('local/nielsen_li_validation/reverse/answers.jsonl','a');consec=0;n=0
for x in Q:
    run=x['run']
    if run in done: continue
    body=json.dumps({'sample_id':run,'sync':True,'topk':5000}).encode()
    req=urllib.request.Request('https://similarity.logan-search.org/search',body,{'Content-Type':'application/json'})
    d=None;err=''
    for att in range(3):
        try: d=json.load(urllib.request.urlopen(req,timeout=600));
        except Exception as e: err=str(e);d=None
        if d and d.get('status')=='completed': break
        if d: err=str(d.get('error') or d.get('status'))
        time.sleep(5*(att+1))
    if d and d.get('status')=='completed':
        consec=0
        out.write(json.dumps({'run':run,'donor':x['donor'],'results':{k:round(v['score'],5) for k,v in d['results'].items()}})+'\n')
    else:
        consec+=1
        out.write(json.dumps({'run':run,'donor':x['donor'],'error':err})+'\n')
    out.flush();n+=1
    if consec>5: print('STOP consecutive failures');break
    time.sleep(1)
print('finished',n,flush=True)
