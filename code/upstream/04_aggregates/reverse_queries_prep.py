import csv,json,gzip
M='local/crocodeel_valles/matrices/'
S='local/nielsen_li_validation/shared_donors_runs.tsv'                  # runs of the patients that carry the same identifier in Nielsen and Li (ENA)
idx=set(l.strip() for l in gzip.open('local/adultes_data_v2/logan_index_ids.txt.gz','rt'))
bc={}
for r in csv.DictReader(open('local/nielsen_li_validation/ena_PRJEB5224.tsv'),delimiter='\t'): bc[r['run_accession']]=int(r['base_count'])        # by name: the table of the ENA and the short table of the same name both have this column
meta={}
for c in ('NielsenHB_2014','LiJ_2014'):
    meta[c]=list(csv.DictReader(open(M+c+'.metadata.tsv'),delimiter='\t'))
shared=list(csv.DictReader(open(S),delimiter='\t'))
# donor -> subject id via run -> metadata sample
run2samp={}
for c in meta:
    for m in meta[c]:
        for a in m['NCBI_accession'].split(';'): run2samp[(c,a)]=m
subj={}
for x in shared:
    m=run2samp.get((x['cohort'],x['run']))
    if m: subj.setdefault(x['donor'],set()).add(m['subject_id'])
bad={d:s for d,s in subj.items() if len(s)!=1}
print('donors',len(subj),'ambiguous',bad)
subj={d:next(iter(s)) for d,s in subj.items()}
# all runs per subject (both cohorts, all samples), nielsen stool runs
allruns={};nstool={}
for c in meta:
    for m in meta[c]:
        s=m['subject_id']
        if s not in subj.values(): continue
        rs=set(m['NCBI_accession'].split(';'))
        allruns.setdefault(s,set()).update(rs)
        if c=='NielsenHB_2014' and m['body_site']=='stool': nstool.setdefault(s,set()).update(rs)
# also shared file runs
for x in shared:
    if x['donor'] in subj: allruns[subj[x['donor']]].add(x['run'])
# Li queries
q=[]
for d in sorted(subj):
    lir=[x for x in shared if x['donor']==d and x['cohort']=='LiJ_2014']
    for samp in sorted({x['sample_id'] for x in lir}):
        rs=[x['run'] for x in lir if x['sample_id']==samp]
        cand=[r for r in rs if r in idx]
        if not cand: q.append(dict(donor=d,sample=samp,run=None,runs=rs));continue
        best=max(cand,key=lambda r:bc[r])
        q.append(dict(donor=d,sample=samp,run=best,runs=rs,nruns=len(rs),inidx=len(cand),bases=bc[best],subject=subj[d]))
json.dump(dict(q=q,subj=subj,allruns={k:sorted(v) for k,v in allruns.items()},nstool={k:sorted(v) for k,v in nstool.items()}),open('local/nielsen_li_validation/reverse/prep.json','w'))
print(len(q),sum(1 for x in q if x['run']),[ (x['nruns'],x['inidx']) for x in q][:10])
print('no stool nielsen',[d for d in subj if subj[d] not in nstool])
