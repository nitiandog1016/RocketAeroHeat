"""随小版本保存复核用的边界表和求解记录，大场登记位置与哈希。"""
import csv
import hashlib
import io
import json
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np

from check_grid import HERE,ROOT,sha
from thermistorWireSurface import residual


def as_csv(header,rows):
    data=io.StringIO();w=csv.writer(data,lineterminator="\n");w.writerow(header);w.writerows(rows)
    return data.getvalue().encode("ascii")


if __name__=="__main__":
    base=ROOT/"outputs/v330_validation";work=base/"run01";payload={};large={}
    files=[base/"grid.json",base/"audit.json",base/"prepared01/prepared.json",base/"ui01/ui.json",base/"ui01/debugger_process.json",
           base/"export01/completed.json",*list((base/"build01").glob("*.json"))]
    files.extend(work.glob("*.json"))
    for d in work.iterdir():
        if not d.is_dir():continue
        for p in d.iterdir():
            if p.name in ("debugger_process.json","debugger.log","record.json","case.cfg","history.csv","stages.csv","far_pairs.csv","flow_far.csv"):files.append(p)
            if p.suffix==".bin" or p.name.startswith("restart_") and p.suffix==".csv":large[p.relative_to(base).as_posix()]=sha(p)
    for p in files:payload[p.relative_to(base).as_posix()]=p.read_bytes()
    for k,rel in (("early","releaseEarly"),("late","releaseLate")):
        n=np.fromfile(work/k/"nodes.bin","<f8").reshape(-1,16);n=n[n[:,0].argsort()]
        flow=np.loadtxt(work/k/"flow_far.csv",delimiter=",");fi=flow[:,0].astype(int);w=np.flatnonzero(n[:,4])
        payload[k+"_far.csv"]=as_csv(("id","old","flow","scalar"),zip(fi,n[fi,1],flow[:,1],n[fi,11]))
        payload[k+"_wall.csv"]=as_csv(("id","release","full","limited"),zip(w,residual(work/rel/"restart.csv")[w],n[w,9],n[w,10]))
    manifest=dict(included={k:hashlib.sha256(v).hexdigest() for k,v in payload.items()},largeLocal=large)
    payload["manifest.json"]=json.dumps(manifest,indent=2).encode()
    out=HERE/"evidence";out.mkdir(exist_ok=True)
    with ZipFile(out/"receipts.zip","x",compression=ZIP_DEFLATED) as z:
        for name,data in payload.items():z.writestr(name,data)
    with ZipFile(out/"receipts.zip") as z:assert all(z.read(k)==v for k,v in payload.items())
