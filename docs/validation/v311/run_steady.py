"""从已审核的全周侧吹输入，或同一工况的正常保存帧，进行有限稳态推进。"""
import argparse
import os
from pathlib import Path
import sys
import uuid

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import ENGINE, ENGINE_SHA, read, sha, solve, write_json
from thermistorReaction import terminal


def frame(path,xyz,wall):
    """检查官方输出的压力、无滑移速度及真实强温差场。"""
    state=np.loadtxt(path,delimiter=",",skiprows=1)
    assert state.shape[0]==len(xyz) and state.shape[1]>=10 and np.isfinite(state).all()
    np.testing.assert_array_equal(state[:,0],np.arange(len(xyz)))
    np.testing.assert_allclose(state[:,1:4],xyz,rtol=0,atol=1e-12)
    theta=state[:,9];pressure=101325+1.18*9*state[:,4]
    assert pressure.min()>0 and theta.min()>=-1e-12 and theta.max()<=1+1e-12
    assert np.all(theta[wall]==1) and np.max(abs(state[wall,5:8]))<1e-10
    return dict(sha256=sha(path),pressurePa=[float(pressure.min()),float(pressure.max())],
                thetaRange=[float(theta.min()),float(theta.max())])


def run(cold,out,steps,parent,threads):
    info=read(cold/"mesh_audit.json")
    assert info["coldStateChecked"] and all(sha(cold/k)==h for k,h in info["files"].items())
    assert sha(ENGINE)==ENGINE_SHA and info["cells"]==278688
    start=cold/"initial_flow.csv";prior=0
    if parent:
        rec=terminal(parent);previous=read(parent/"field_audit.json")
        assert previous["meshAuditSHA256"]==sha(cold/"mesh_audit.json")
        start=parent/"restart.csv";assert sha(start)==previous["lastFrame"]==rec["outputSHA256"]["restart.csv"]
        prior=previous["cumulativeSteps"]
    out.mkdir(parents=True,exist_ok=False)
    for name in ("mesh.su2","body.json","nodes.npz","case.json"):
        os.link(cold/name,out/name)
    os.link(start,out/"initial_flow.csv")
    cfg=(cold/"case.cfg").read_text("ascii");assert cfg.count("ITER= 7000\n")==1
    (out/"case.cfg").write_text(cfg.replace("ITER= 7000\n",f"ITER= {steps}\n"),encoding="ascii")
    files=[HERE/"PLAN.md",HERE/"audit_mesh.py",Path(__file__).resolve(),ROOT/"scripts/thermistorBead.py",
           ROOT/"scripts/thermistorLeadCrossflow.py",cold/"mesh_audit.json",out/"case.json",out/"nodes.npz",out/"body.json"]
    if parent:files.extend((HERE/"STAGES.md",parent/"engine_process.json",parent/"field_audit.json"))
    write_json(out/"engine_request.json",dict(workDirectory=str(out),configFile="case.cfg",meshFile="mesh.su2",
        threads=threads,maximumSeconds=720,jobToken=str(uuid.uuid4()),sourceFiles=[p.relative_to(ROOT).as_posix() for p in files]))
    solve(out,memory_bytes=3*1024**3,output_bytes=1536*1024**2)
    receipt=terminal(out)
    with (out/"history.csv").open(encoding="ascii") as f:header=[s.strip(' "\n') for s in f.readline().split(",")]
    history=np.loadtxt(out/"history.csv",delimiter=",",skiprows=1,ndmin=2)
    iteration=history[:,header.index("Inner_Iter")].astype(int)
    assert np.array_equal(iteration,np.arange(len(iteration))) and np.isfinite(history).all()
    last=int(iteration[-1]);assert receipt["progress"]["Inner_Iter"]==last
    with np.load(out/"nodes.npz") as z:xyz=z["xyz"];wall=np.unique(z["wallMap"])
    frames={}
    for i in sorted(set([*range(100,last+1,100),last])):
        file=f"restart_{i:06d}.csv";record=frame(out/file,xyz,wall)
        surface=f"surface_{i:06d}.vtk"
        assert record["sha256"]==receipt["outputSHA256"][file] and sha(out/surface)==receipt["outputSHA256"][surface]
        record["surfaceSHA256"]=receipt["outputSHA256"][surface];frames[file]=record
    assert sha(out/"restart.csv")==frames[f"restart_{last:06d}.csv"]["sha256"]
    residual={k:float(max(history[-20:,header.index(k)])) for k in ("rms[P]","rms[U]","rms[V]","rms[W]","rms[h]","rms[rho*Y_0]")}
    write_json(out/"field_audit.json",dict(meshAuditSHA256=sha(cold/"mesh_audit.json"),steps=last+1,cumulativeSteps=prior+last+1,
        frames=frames,lastFrame=sha(out/"restart.csv"),last20Residual=residual,initialSHA256=sha(start),
        allResidualsBelowMinus8=len(history)>=20 and max(residual.values())<-8,seconds=receipt["elapsedSeconds"]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("cold",type=Path);p.add_argument("out",type=Path)
    p.add_argument("--steps",type=int,required=True);p.add_argument("--parent",type=Path);p.add_argument("--threads",type=int,default=4)
    a=p.parse_args();out=a.out.resolve();out.relative_to(ROOT/"outputs")
    assert 1<=a.steps<=7000 and 1<=a.threads<=8
    run(a.cold.resolve(),out,a.steps,a.parent.resolve() if a.parent else None,a.threads)
