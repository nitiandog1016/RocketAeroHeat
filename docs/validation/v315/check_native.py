"""以不同原生强温差边界验算 B 珠 12 K 测点，沿用 96 档流场与保守壁面读数。"""
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
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha,solve,write_json
from thermistorReaction import terminal


def run(base,out,parent=None):
    key=read(ROOT/"cases/thermistor_b/sidewindAngular96.json");src=terminal(base)
    c=read(base/"case.json");s=read(base/"surface/summary.json")
    assert c["angularCells"]==96 and c["wallK"]==299.15 and c["airK"]==293.15 and c["velocityMps"]==[0,-3,0]
    assert src["outputSHA256"]["restart.csv"]==s["sourceFrameSHA256"]==s["energy"]["inputFrameSHA256"]
    assert max(s["residualWindow"].values())<-8 and max(s["plus20Window"].values())<-8
    assert s["relativeChange20"]<.002 and abs(s["energy"]["relativeDifference"])<.005
    assert sha(base/"surface/wall.csv")==s["wallSHA256"] and sha(base/"case.json")==key["caseSHA256"]
    for n,k in (("mesh.su2","meshSHA256"),("nodes.npz","nodeMapSHA256"),("body.json","bodySHA256")):
        assert sha(base/n)==s[k]==key[k]
    out.mkdir(parents=True,exist_ok=False)
    for n in ("mesh.su2","nodes.npz","body.json"):os.link(base/n,out/n)
    c.update(wallK=305.15,scalarReferenceDeltaK=24,sourceSixKFrame=s["sourceFrameSHA256"])
    write_json(out/"case.json",c)
    with np.load(out/"nodes.npz") as z:xyz=z["xyz"];wet=np.unique(z["wallMap"])
    def check_frame(path):
        data=np.loadtxt(path,delimiter=",",skiprows=1)
        np.testing.assert_array_equal(data[:,0],np.arange(len(xyz)))
        np.testing.assert_allclose(data[:,1:4],xyz,rtol=0,atol=1e-12)
        assert np.isfinite(data).all() and np.all(data[wet,9]==.5) and np.max(abs(data[wet,5:8]))<1e-10
        assert data[:,9].min()>=-1e-12 and data[:,9].max()<=.5+1e-12
        pressure=c["airPa"]+c["rho"]*9*data[:,4];assert pressure.min()>0
        return dict(sha256=sha(path),theta=[data[:,9].min(),data[:,9].max()],pressurePa=[pressure.min(),pressure.max()])
    previous=0
    locks=[HERE/"PLAN.md",Path(__file__).resolve(),base/"engine_process.json",base/"surface/summary.json",base/"case.cfg",
           ROOT/"scripts/thermistorCrossflowSurface.py",ROOT/"scripts/thermistorReaction.py",ROOT/"scripts/thermistorBead.py"]
    if parent:
        a=read(parent/"field_audit.json");r=terminal(parent)
        assert a["lastFrame"]==r["outputSHA256"]["restart.csv"] and sha(parent/"case.json")==sha(out/"case.json")
        assert sha(parent/"mesh.su2")==key["meshSHA256"] and not a["allResidualsBelowMinus8"]
        os.link(parent/"restart.csv",out/"initial_flow.csv");previous=a["cumulativeSteps"];steps=300
        locks.extend((parent/"field_audit.json",parent/"engine_process.json"))
    else:
        initial=np.loadtxt(base/"restart.csv",delimiter=",",skiprows=1)
        assert np.all(initial[wet,9]==1)
        initial[:,9]*=.5
        initial[:,8]=(1-298.15)+.5*(initial[:,8]-(1-298.15))
        with (base/"initial_flow.csv").open(encoding="ascii") as f:header=f.readline().strip()
        np.savetxt(out/"initial_flow.csv",initial,delimiter=",",header=header,comments="",fmt=["%d"]+["%.17g"]*9)
        steps=40
    cfg=(base/"case.cfg").read_text("ascii")
    for old,new in (("(wall, VALUE, 1.0)","(wall, VALUE, 0.5)"),("MARKER_ISOTHERMAL= (wall, 2)","MARKER_ISOTHERMAL= (wall, 1.5)"),
                    ("OUTPUT_WRT_FREQ= (100, 100)","OUTPUT_WRT_FREQ= (20, 20)")):
        assert cfg.count(old)==1;cfg=cfg.replace(old,new)
    cfg,count=re.subn(r"^ITER= \d+$",f"ITER= {steps}",cfg,flags=re.M);assert count==1
    (out/"case.cfg").write_text(cfg,encoding="ascii");initial_info=check_frame(out/"initial_flow.csv")
    names=("mesh.su2","case.json","nodes.npz","body.json")
    locks.extend(out/n for n in names)
    write_json(out/"engine_request.json",dict(workDirectory=str(out),configFile="case.cfg",meshFile="mesh.su2",
               maximumSeconds=600 if parent else 240,threads=4,jobToken=str(uuid.uuid4()),sourceFiles=[p.relative_to(ROOT).as_posix() for p in locks]))
    solve(out,memory_bytes=4*1024**3,output_bytes=(1536 if parent else 1024)*1024**2)
    rec=terminal(out)
    with (out/"history.csv").open(encoding="ascii") as f:
        hdr=[k.strip() for k in next(csv.reader(f,skipinitialspace=True))]
        rows=np.loadtxt(f,delimiter=",",ndmin=2)
    assert len(rows)>=20 and np.isfinite(rows).all()
    np.testing.assert_array_equal(rows[:,hdr.index("Inner_Iter")],np.arange(len(rows)))
    last=len(rows)-1;assert last==rec["progress"]["Inner_Iter"]
    frames={}
    for i in sorted(set(range(20,last+1,20))|{last}):
        name=f"restart_{i:06d}.csv";frames[name]=check_frame(out/name)
        assert frames[name]["sha256"]==rec["outputSHA256"][name]
        surface=f"surface_{i:06d}.vtk";assert sha(out/surface)==rec["outputSHA256"][surface]
    assert frames[name]["sha256"]==rec["outputSHA256"]["restart.csv"]
    window={k:float(rows[-20:,hdr.index(k)].max()) for k in ("rms[P]","rms[U]","rms[V]","rms[W]","rms[h]","rms[rho*Y_0]")}
    write_json(out/"field_audit.json",dict(steps=last+1,cumulativeSteps=previous+last+1,sourceFrame=s["sourceFrameSHA256"],
               lastFrame=frames[name]["sha256"],initial=initial_info,frames=frames,last20Residual=window,allResidualsBelowMinus8=max(window.values())<-8))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("base",type=Path);p.add_argument("out",type=Path);p.add_argument("--parent",type=Path)
    a=p.parse_args();out=a.out.resolve();out.relative_to(ROOT/"outputs");run(a.base.resolve(),out,a.parent.resolve() if a.parent else None)
