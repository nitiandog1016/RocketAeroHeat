"""从 v314 同一接受态等步推进原 6 K 工况，检验 12 K 边界读数的迭代位置效应。"""
import argparse
import csv
import os
from pathlib import Path
import re
import sys
import uuid

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path[:0]=[str(ROOT/"scripts"),str(HERE.parent/"v311")]
from thermistorBead import read,sha,write_json,solve
from thermistorReaction import terminal
from run_steady import frame


def run(base,out):
    source=terminal(base);s=read(base/"surface/summary.json");k=read(ROOT/"cases/thermistor_b/sidewindAngular96.json")
    assert source["outputSHA256"]["restart.csv"]==s["sourceFrameSHA256"]==s["energy"]["inputFrameSHA256"]
    assert all(sha(base/n)==k[h] for n,h in (("mesh.su2","meshSHA256"),("nodes.npz","nodeMapSHA256"),("body.json","bodySHA256"),("case.json","caseSHA256")))
    assert max(s["residualWindow"].values())<-8 and abs(s["energy"]["relativeDifference"])<.005
    out.mkdir(parents=True,exist_ok=False)
    names=("mesh.su2","nodes.npz","case.json","body.json")
    for n in names:os.link(base/n,out/n)
    os.link(base/"restart.csv",out/"initial_flow.csv")
    cfg=(base/"case.cfg").read_text("ascii");assert cfg.count("OUTPUT_WRT_FREQ= (100, 100)")==1
    cfg,n=re.subn(r"^ITER= \d+$","ITER= 40",cfg,flags=re.M);assert n==1
    (out/"case.cfg").write_text(cfg.replace("OUTPUT_WRT_FREQ= (100, 100)","OUTPUT_WRT_FREQ= (20, 20)"),encoding="ascii")
    locks=[HERE/"PLAN.md",HERE/"MATCHED.md",Path(__file__).resolve(),HERE.parent/"v311/run_steady.py",
           base/"engine_process.json",base/"surface/summary.json",ROOT/"scripts/thermistorCrossflowSurface.py",*(out/n for n in names)]
    write_json(out/"engine_request.json",dict(workDirectory=str(out),configFile="case.cfg",meshFile="mesh.su2",threads=4,
               maximumSeconds=240,jobToken=str(uuid.uuid4()),sourceFiles=[p.relative_to(ROOT).as_posix() for p in locks]))
    solve(out,memory_bytes=4*1024**3,output_bytes=1024**3);rec=terminal(out)
    with np.load(out/"nodes.npz") as z:xyz=z["xyz"];wet=np.unique(z["wallMap"])
    with (out/"history.csv").open(encoding="ascii") as f:
        header=[k.strip() for k in next(csv.reader(f,skipinitialspace=True))];history=np.loadtxt(f,delimiter=",",ndmin=2)
    np.testing.assert_array_equal(history[:,header.index("Inner_Iter")],np.arange(40))
    assert np.isfinite(history).all() and rec["progress"]["Inner_Iter"]==39
    frames={}
    for i in (20,39):
        n=f"restart_{i:06d}.csv";frames[n]=frame(out/n,xyz,wet)
        assert frames[n]["sha256"]==rec["outputSHA256"][n]
        surface=f"surface_{i:06d}.vtk";assert sha(out/surface)==rec["outputSHA256"][surface]
    assert frames[n]["sha256"]==rec["outputSHA256"]["restart.csv"]
    window={k:float(history[-20:,header.index(k)].max()) for k in ("rms[P]","rms[U]","rms[V]","rms[W]","rms[h]","rms[rho*Y_0]")}
    write_json(out/"field_audit.json",dict(steps=40,cumulativeSteps=40,sourceFrame=s["sourceFrameSHA256"],
               initialSHA256=sha(out/"initial_flow.csv"),lastFrame=frames[n]["sha256"],frames=frames,last20Residual=window,allResidualsBelowMinus8=max(window.values())<-8))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("base",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();out=a.out.resolve();out.relative_to(ROOT/"outputs");run(a.base.resolve(),out)
