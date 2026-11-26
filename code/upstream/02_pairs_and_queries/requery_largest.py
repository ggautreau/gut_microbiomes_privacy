"""Re-query the public similarity service of Faure et al. with the largest indexed run of the query sample, for the donors whose queried run
(the first run listed by curatedMetagenomicData) was not that run. One accession per request, top 5000 neighbours. Resumable."""
import json, os, sys, time, urllib.request
SP = os.environ.get("WORK", "local/work")                                  # folder of requery_todo.json (the donors whose queried run was not the largest indexed run of the sample)
todo = json.load(open(f"{SP}/requery_todo.json"))
out_path = f"{SP}/results_largest.jsonl"
done = set()
if os.path.exists(out_path):
    for l in open(out_path):
        try:
            q = json.loads(l)
        except Exception:
            continue
        if "results" in q:
            done.add(q["run"])
print(len(todo), "queries;", len(done), "already done", flush=True)
out = open(out_path, "a"); n = fails = 0; t0 = time.time()
for t in todo:
    run = t["new"]
    if run in done:
        continue
    body = json.dumps({"sample_id": run, "sync": True, "topk": 5000}).encode()
    req = urllib.request.Request("https://similarity.logan-search.org/search", body, {"Content-Type": "application/json"})
    d = None; err = ""
    for att in range(3):
        try:
            d = json.load(urllib.request.urlopen(req, timeout=600)); break
        except Exception as e:
            err = str(e); time.sleep(5 * (att + 1))
    if not d or d.get("status") != "completed":
        fails += 1
        out.write(json.dumps({"run": run, "study": t["study"], "subject": t["subject"], "error": (d or {}).get("error") or err}) + "\n"); out.flush()
    else:
        out.write(json.dumps({"run": run, "study": t["study"], "subject": t["subject"], "results": {k: round(v["score"], 5) for k, v in d["results"].items()}}) + "\n"); out.flush()
    n += 1; time.sleep(1)
    if n % 10 == 0:
        print(n, "done, failures", fails, "elapsed %.0f s" % (time.time() - t0), flush=True)
print("finished", n, "failures", fails, "elapsed %.0f s" % (time.time() - t0), flush=True)
