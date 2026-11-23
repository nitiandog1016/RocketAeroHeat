"""归档原生小记录及远场/壁节点明细。大场留在本机，清单给出其来源及哈希。"""
import csv
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np

from run_shared import ROOT,HERE,read,sha
from thermistorWireSurface import residual


def as_csv(names,rows):
    f=io.StringIO();w=csv.writer(f,lineterminator="\n");w.writerow(names);w.writerows(rows)
    return f.getvalue().encode("ascii")


def pack(out):
    base=ROOT/"outputs/v329_validation";run=base/"run01";coldUi=base/"coldUi01"
    cold=Path(read(coldUi/"ui.json")["job"]);payload={};onlyLocal={}
    def add(path):
        key=path.relative_to(ROOT).as_posix();payload[key]=path.read_bytes()
    def native(work):
        for d in sorted(work.iterdir()):
            if d.is_dir():
                for f in d.iterdir():
                    if f.name in ("debugger_process.json","debugger.log","case.cfg","record.json","stages.csv","far_pairs.csv","history.csv","flow_far.csv"):add(f)
                    if f.suffix==".bin" or f.name.startswith("restart_") and f.suffix==".csv":
                        onlyLocal[f.relative_to(ROOT).as_posix()]=sha(f)
            elif d.suffix==".json":add(d)
    native(run);native(cold/"run")
    for d in (base/"ui01",coldUi):
        for name in ("ui.json","debugger_process.json"):add(d/name)
    for f in (base/"build01").glob("*.json"):add(f)
    add(cold/"completed.json");add(base/"verified_endpoints.json")
    for name,rel in (("early","releaseEarly"),("late","releaseLate")):
        obs=run/name;n=np.fromfile(obs/"nodes.bin",dtype="<f8").reshape(-1,16);n=n[n[:,0].argsort()]
        far=np.loadtxt(obs/"flow_far.csv",delimiter=",");ids=far[:,0].astype(int)
        payload[name+"_far.csv"]=as_csv(("id","old","flowMass","scalarMass","scalarAdv","divergence"),zip(ids,n[ids,1],far[:,1],n[ids,11],n[ids,12],n[ids,13]))
        wall=np.flatnonzero(n[:,4]);raw=residual(run/rel/"restart.csv")
        payload[name+"_wall.csv"]=as_csv(("id","releaseRHS","fullCorrection","limitedCorrection"),zip(wall,raw[wall],n[wall,9],n[wall,10]))
    a=read(run/"thermal.json");b=read(cold/"run/thermal.json")
    with ZipFile(base/"surface01/surface.zip") as za, ZipFile(cold/"surface/surface.zip") as zb:
        get=lambda z:np.loadtxt(io.BytesIO(z.read("wall.csv")),delimiter=",",skiprows=1)
        diff=np.max(abs(get(za)-get(zb)),axis=0)
    assert max(diff)<1e-8 and b["passed"] and b["updates"]==a["updates"]
    payload["repeat.json"]=json.dumps(dict(gate=b["passed"],frame=b["early"]["sourceFrame"],maxAbsoluteColumnDifferences=diff.tolist())).encode()
    manifest=dict(included={k:hashlib.sha256(v).hexdigest() for k,v in payload.items()},local=onlyLocal)
    payload["manifest.json"]=json.dumps(manifest,indent=2).encode()
    with ZipFile(out,"x",compression=ZIP_DEFLATED) as z:
        for name,data in payload.items():z.writestr(name,data)
    with ZipFile(out) as z:assert all(z.read(k)==v for k,v in payload.items())


if __name__=="__main__":pack(HERE/"evidence/receipts.zip")
