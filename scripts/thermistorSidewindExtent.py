"""对带线 B 型热敏珠侧吹，在完整 20D 原网外接球壳，运行 30D 全周外域对照。"""
import argparse
import csv
import os
from pathlib import Path
import shutil
import sys
import uuid

import numpy as np

from thermistorBead import ENGINE, ENGINE_SHA, ROOT, read, sha, solve, write_json
from thermistorOuterDomain import append_shell, audit, load_mesh, write_mesh
from thermistorReaction import terminal

HERE=ROOT/"docs/validation/v312"
sys.path.insert(0,str(HERE.parent/"v311"))
from run_steady import frame


def prepare(base,seed,out):
    """锁定已交付的网格和终帧，生成新半径下的同一试件、初猜和证据。"""
    key=read(ROOT/"cases/thermistor_b/sidewind3.json")
    files={"mesh.su2":"meshSHA256","nodes.npz":"nodeMapSHA256","case.json":"caseSHA256",
           "body.json":"bodySHA256","case.cfg":"configSHA256","initial_flow.csv":"initialSHA256"}
    assert all(sha(base/k)==key[v] for k,v in files.items())
    rec=terminal(seed)
    # 安装在另一台机器时，接受该机已通过全周收支检查的同工况帧，并登记其真实指纹。
    seed_hash=sha(seed/"restart.csv");checked=read(seed/"field_audit.json");wall=read(seed/"surface/summary.json")
    assert seed_hash==rec["outputSHA256"]["restart.csv"]==checked["lastFrame"]==wall["sourceFrameSHA256"]
    assert checked["allResidualsBelowMinus8"] and max(wall["plus20Window"].values())<-8
    assert wall["energy"]["inputFrameSHA256"]==seed_hash and abs(wall["energy"]["relativeDifference"])<.005
    assert wall["relativeChange20"]<.002 and sha(seed/"surface/wall.csv")==wall["wallSHA256"]
    assert sha(seed/"mesh.su2")==key["meshSHA256"] and sha(seed/"case.json")==key["caseSHA256"]
    c=read(base/"case.json")
    assert c["velocityMps"]==[0,-3,0] and c["sectorDeg"]==360 and c["outerRadiusOverD"]==20
    out.mkdir(parents=True,exist_ok=False)
    cells,coords,marks=load_mesh(base/"mesh.su2",("far","wall"))
    with np.load(base/"nodes.npz") as src:
        np.testing.assert_array_equal(coords,src["xyz"])
        full,xyz,boundary,rays=append_shell(cells,coords,marks,src["ids"])
        np.savez_compressed(out/"nodes.npz",xyz=xyz,ids=rays,wallMap=src["wallMap"],firstMap=src["firstMap"],area=src["area"])
    write_mesh(out/"mesh.su2",full,xyz,boundary)
    # 壁面、工况和离散沿用原本定义，唯一几何改动为外域半径。
    for name in ("case.cfg","body.json"):shutil.copyfile(base/name,out/name)
    c["outerRadiusOverD"]=30;write_json(out/"case.json",c)
    shell=audit(base,out,("far","wall"))
    np.testing.assert_array_equal(rays[0],rays[-1])
    far=np.unique([i for face in marks["far"] for i in face])
    prev=np.loadtxt(seed/"restart.csv",skiprows=1,delimiter=",")
    np.testing.assert_array_equal(prev[:,0],np.arange(len(coords)))
    np.testing.assert_allclose(prev[:,1:4],coords,rtol=0,atol=1e-13)
    ambient=np.array([0,0,-1,0,1-298.15,0])
    state=np.zeros((len(xyz),10));state[:len(prev)]=prev
    # SU2 保存坐标可能有末位打印舍入。初猜的物理列直接抄录，坐标绑定成品网格。
    state[:,0]=np.arange(len(xyz));state[:,1:4]=xyz
    for layer in range(1,6):
        pos=len(prev)+(layer-1)*len(far)+np.arange(len(far))
        state[pos,4:]=ambient+(1-layer/5)*(prev[far,4:]-ambient)
        np.testing.assert_allclose(xyz[pos],coords[far]*(1+layer/10),atol=1e-13,rtol=0)
    with (base/"initial_flow.csv").open(encoding="ascii") as f:header=f.readline().strip()
    np.savetxt(out/"initial_flow.csv",state,delimiter=",",fmt=["%d"]+["%.17g"]*9,comments="",header=header)
    restored=np.loadtxt(out/"initial_flow.csv",skiprows=1,delimiter=",")
    np.testing.assert_array_equal(restored,state)
    np.testing.assert_array_equal(restored[:len(prev),4:],prev[:,4:])
    assert np.all(restored[-len(far):,4:]==ambient) and restored[:,9].min()>=-1e-12 and restored[:,9].max()<=1+1e-12
    assert len(full)==len(cells)+5*len(marks["far"]) and len(xyz)==len(coords)+5*len(far)
    write_json(out/"mesh_audit.json",dict(shell=shell,baseMeshSHA256=key["meshSHA256"],
               seedSHA256=seed_hash,files={name:sha(out/name) for name in files},
               nodes=len(xyz),cells=len(full),physicalFaces={k:len(v) for k,v in boundary.items()},
               initialState="20D 侧吹终场，外壳沿半径过渡至自由流",axisAndSeamChecked=True))


def run(prepared,out,steps,parent=None):
    """新域从已确认的初猜开始，续算只接受同域的真实终帧。"""
    info=read(prepared/"mesh_audit.json");digest=sha(prepared/"mesh_audit.json")
    assert info["baseMeshSHA256"]==read(ROOT/"cases/thermistor_b/sidewind3.json")["meshSHA256"]
    assert info["axisAndSeamChecked"] and all(sha(prepared/k)==h for k,h in info["files"].items())
    assert sha(ENGINE)==ENGINE_SHA and 1<=steps<=7000
    initial=prepared/"initial_flow.csv";prior=0
    locks=[HERE/"PLAN.md",HERE/"NATIVE.md",Path(__file__).resolve(),ROOT/"scripts/thermistorOuterDomain.py",
           ROOT/"scripts/thermistorBead.py",HERE.parent/"v311/run_steady.py",prepared/"mesh_audit.json"]
    if parent:
        old=terminal(parent);a=read(parent/"field_audit.json");initial=parent/"restart.csv"
        assert a["meshAuditSHA256"]==digest and a["lastFrame"]==sha(initial)==old["outputSHA256"]["restart.csv"]
        prior=a["cumulativeSteps"];locks.extend((parent/"engine_process.json",parent/"field_audit.json"))
    out.mkdir(parents=True,exist_ok=False)
    for name in ("mesh.su2","body.json","nodes.npz","case.json"):
        os.link(prepared/name,out/name);locks.append(out/name)
    os.link(initial,out/"initial_flow.csv")
    cfg=(prepared/"case.cfg").read_text("ascii");assert cfg.count("ITER= 7000\n")==1
    (out/"case.cfg").write_text(cfg.replace("ITER= 7000\n",f"ITER= {steps}\n"),encoding="ascii")
    write_json(out/"engine_request.json",dict(workDirectory=str(out),configFile="case.cfg",meshFile="mesh.su2",
        threads=4,maximumSeconds=720,jobToken=str(uuid.uuid4()),sourceFiles=[p.relative_to(ROOT).as_posix() for p in locks]))
    solve(out,memory_bytes=3*1024**3,output_bytes=1536*1024**2)
    receipt=terminal(out)
    with (out/"history.csv").open(encoding="ascii") as f:
        header=[word.strip() for word in next(csv.reader(f,skipinitialspace=True))]
        history=np.loadtxt(f,delimiter=",",ndmin=2)
    it=history[:,header.index("Inner_Iter")]
    np.testing.assert_array_equal(it,np.arange(len(it)));assert np.isfinite(history).all()
    last=int(it[-1]);assert last==receipt["progress"]["Inner_Iter"]
    with np.load(out/"nodes.npz") as z:xyz=z["xyz"];wall=np.unique(z["wallMap"])
    frames={}
    for i in sorted(set([*range(100,last+1,100),last])):
        name=f"restart_{i:06d}.csv";f=frame(out/name,xyz,wall);surface=f"surface_{i:06d}.vtk"
        assert f["sha256"]==receipt["outputSHA256"][name] and sha(out/surface)==receipt["outputSHA256"][surface]
        f["surfaceSHA256"]=receipt["outputSHA256"][surface];frames[name]=f
    end=sha(out/"restart.csv");assert frames[f"restart_{last:06d}.csv"]["sha256"]==end
    residual={k:float(history[-20:,header.index(k)].max()) for k in ("rms[P]","rms[U]","rms[V]","rms[W]","rms[h]","rms[rho*Y_0]")}
    write_json(out/"field_audit.json",dict(meshAuditSHA256=digest,initialSHA256=sha(initial),lastFrame=end,frames=frames,
        cumulativeSteps=prior+last+1,steps=last+1,last20Residual=residual,allResidualsBelowMinus8=len(it)>=20 and max(residual.values())<-8))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest="action",required=True)
    a=sub.add_parser("prepare");a.add_argument("base",type=Path);a.add_argument("seed",type=Path);a.add_argument("output",type=Path)
    a=sub.add_parser("run");a.add_argument("prepared",type=Path);a.add_argument("output",type=Path)
    a.add_argument("--steps",required=True,type=int);a.add_argument("--parent",type=Path)
    a=p.parse_args();out=a.output.resolve();out.relative_to(ROOT/"outputs")
    if a.action=="prepare":prepare(a.base.resolve(),a.seed.resolve(),out)
    else:run(a.prepared.resolve(),out,a.steps,a.parent.resolve() if a.parent else None)
