"""Download the list of sample ids indexed by similarity.logan-search.org (GET /samples, 10000 per page)."""
import json,urllib.request,time,gzip
ids=[]; off=0; total=None
while total is None or off<total:
    for att in range(4):
        try:
            d=json.load(urllib.request.urlopen(f'https://similarity.logan-search.org/samples?offset={off}&size=10000',timeout=120)); break
        except Exception as e: time.sleep(3*(att+1)); d=None
    if d is None: print('failed at',off); break
    total=d['total']; ids+=d['ids']; off+=10000
    if (off//10000)%50==0: print(off,'/',total,flush=True)
    time.sleep(0.2)
gzip.open('logan_index_ids.txt.gz','wt').write('\n'.join(ids))
print('saved',len(ids),'of',total)
