"""归档真实历史、失败回执和两个端点边界小表，保留大场的本地身份。"""
import csv
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

from run import BASE, HERE, ROOT, read, sha
from thermistorWireSurface import residual


def csv_bytes(header, rows):
    buffer=io.StringIO()
    writer=csv.writer(buffer,lineterminator="\n");writer.writerow(header);writer.writerows(rows)
    return buffer.getvalue().encode("ascii")


def pack(out):
    folders=[BASE/"run01"/f"steady{i:02d}" for i in range(1,4)]+[BASE/"repeat03"]
    test=BASE/"thermal01"
    files=[BASE/"build01"/"affine.log",test/"samples.json",BASE/"ui01"/"ui.json",BASE/"ui01"/"debugger_process.json"]
    large=[]
    for build in (BASE/"build01",BASE/"readBuild01"):
        files.extend(build.glob("*.json"));large.extend(build.glob("*.exe"))
    for folder in folders+[p for p in test.iterdir() if p.is_dir()]:
        files.extend(p for p in folder.iterdir() if p.name in
                     ("case.cfg","history.csv","stages.csv","debugger_process.json","stage.json","sample.json","debugger.log"))
        large.extend(folder.glob("restart_*.csv"));large.extend(folder.glob("*.bin"))
    def key(path):return path.relative_to(ROOT).as_posix()
    data={key(p):p.read_bytes() for p in files}
    for name,release in (("early","releaseEarly"),("late","releaseLate")):
        folder=test/name
        n=np.fromfile(folder/"nodes.bin",dtype="<f8").reshape(-1,16);n=n[n[:,0].argsort()]
        far=np.loadtxt(folder/"flow_far.csv",delimiter=",");ids=far[:,0].astype(int)
        data[name+"_far.csv"]=csv_bytes(("id","theta","scalarMass","flowFDSMass","scalarOut","divergenceSource"),
            zip(ids,n[ids,1],n[ids,11],far[:,1],n[ids,12],n[ids,13]))
        w=np.flatnonzero(n[:,4]);rhs=residual(test/release/"restart.csv")
        data[name+"_wall.csv"]=csv_bytes(("id","releasedRHS","rawWallSpatial","fullCorrection","limitedCorrection"),
            zip(w,rhs[w],n[w,8],n[w,9],n[w,10]))
    hashes={k:hashlib.sha256(v).hexdigest() for k,v in data.items()}
    manifest={"included":hashes,"largeLocal":{key(p):sha(p) for p in large}}
    data["manifest.json"]=(json.dumps(manifest,indent=2)+"\n").encode()
    out.parent.mkdir(parents=True,exist_ok=True)
    with ZipFile(out,"x",compression=ZIP_DEFLATED) as archive:
        for name,value in data.items():archive.writestr(name,value)
    with ZipFile(out) as archive:
        assert set(archive.namelist())==set(data) and all(archive.read(k)==v for k,v in data.items())


if __name__=="__main__":pack(HERE/"evidence/receipts.zip")
