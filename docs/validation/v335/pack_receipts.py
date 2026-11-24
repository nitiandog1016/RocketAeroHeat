"""为 20D 银球保存可独立积分的热量证据、原生回执和大场指纹。"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile,ZIP_DEFLATED

from outer_grid import HERE, ROOT, read, sha
sys.path.insert(0,str(HERE.parent/"v333"))
from pack import check
from compare_domains import CASE, compare
import numpy as np


def pack(work,dest):
    large={};included={}
    with ZipFile(dest,"x",compression=ZIP_DEFLATED,compresslevel=9) as z:
        def add(name,raw):
            assert name not in included
            included[name]=hashlib.sha256(raw).hexdigest()
            z.writestr(name,raw)
        paths=[work/"prepared/prepared.json",work/"prepared/extension.json",work/"prepare_process.json",
               work/"time/job.json",work/"completed.json",work/"steady_complete.json"]
        for driver in ("driver","time_driver"):paths.extend((work/driver).glob("*.json"))
        paths.extend((work/"steady").glob("*.json"))
        for base in (work/"time",work/"steady"):
            for folder in base.rglob("debugger_process.json"):
                for f in folder.parent.iterdir():
                    if not f.is_file():continue
                    if f.suffix==".json" or f.name in ("case.cfg","debugger.log","wall_budget.npz","stages.csv","far_pairs.csv") or f.name.startswith("history"):
                        paths.append(f)
                    elif f.suffix in (".bin",".vtk") or f.name.startswith(("seed_","restart_")):
                        large[f.relative_to(work).as_posix()]=sha(f)
        for f in paths:add(f.relative_to(work).as_posix(),f.read_bytes())
        for branch,step in (("full",1),("half",1),("half",2)):
            d=work/"time"/branch/f"step{step}";e=read(d/"energy.json")
            old=np.loadtxt(d/f"seed_{step-1:05d}.csv",skiprows=1,delimiter=",",usecols=9)
            new=np.loadtxt(d/f"restart_{step:05d}.csv",skiprows=1,delimiter=",",usecols=9)
            n=np.fromfile(d/"nodes.bin","<f8").reshape(-1,18);order=np.argsort(n[:,0]);n=n[order]
            np.testing.assert_array_equal(n[:,0],np.arange(len(new)))
            before=np.fromfile(d/"before_time.bin","<f8")[order];after=np.fromfile(d/"after_time.bin","<f8")[order]
            starts=np.arange(0,len(new),4096)
            row=np.c_[starts,np.minimum(starts+4096,len(new)),np.add.reduceat(n[:,5]*(new-old)*e["scale"]/e["dt"],starts),
                      np.add.reduceat((after-before)*e["scale"],starts)]
            f=io.StringIO();np.savetxt(f,row,delimiter=",",fmt=["%d","%d","%.17g","%.17g"],comments="",
                                       header="firstID,endExclusive,storageFromFramesW,storageFromNativeW")
            add(f"storage/{branch}{step}.csv",f.getvalue().encode("ascii"))
        z.writestr("manifest.json",json.dumps(dict(files=included,largeLocal=large),indent=2))


def verify_published(evidence):
    check(evidence/"receipts.zip",prefix="time",nodes=269778)
    actual=compare(CASE/"outer20")
    # 重算可能使用不同的 BLAS 归约次序；逐项比较数值，字符串和原面编号保留精确比较。
    def equal(a,b):
        if isinstance(a,dict):
            assert a.keys()==b.keys()
            for k in a:equal(a[k],b[k])
        elif isinstance(a,list):
            assert len(a)==len(b)
            for x,y in zip(a,b):equal(x,y)
        elif isinstance(a,float):
            np.testing.assert_allclose(a,b,rtol=1e-12,atol=1e-12)
        else:assert a==b
    equal(actual,read(evidence/"comparison.json"))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--work",type=Path);args=p.parse_args()
    evidence=HERE/"evidence"
    if args.work:pack(args.work.resolve(),evidence/"receipts.zip")
    verify_published(evidence);print("Checked wall reaction, far transport, storage and both domain comparisons.")
