"""核验同速扇区基流，在冻结全周网格上应用相同扰动并比较两次真实启动。"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PREVIOUS = HERE.parent/"v287"
sys.path.insert(0, str(ROOT/"scripts"))
sys.path.insert(0, str(PREVIOUS))
from audit_step import RMS, audit, frame
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe
from thermistorMatchedBase import CONTRACT, ENGINE, ENGINE_SHA, FULL, FULL_CONFIG_SHA, PARENT, SOURCE, matched_config


def read(path):
    return json.loads(path.read_bytes())


def save(path, record):
    with path.open("x",encoding="utf8") as f:
        json.dump(record,f,ensure_ascii=False,indent=2,allow_nan=False)
        f.write("\n")


def terminal(folder):
    """核对实际启动及退出、预算、进程树和绑定输入，而后才读取保存帧。"""
    rec=read(folder/"debugger_process.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0
    assert rec["stopReason"] is None and rec["supervisorError"] is None and not rec["survivingProcessIDs"]
    assert rec["inputsStable"] and all(rec["inputsStable"].values())
    for path,digest in rec["inputSHA256"].items():
        assert sha(Path(path))==digest,path
    return rec


def base_state(folder):
    rec=terminal(folder)
    assert rec["command"]==[str(ENGINE),"case.cfg","--threads","1"] and sha(ENGINE)==ENGINE_SHA
    assert (folder/"case.cfg").read_text("ascii")==matched_config()
    for key in ("mesh.su2","body.json","nodes.npz"):
        assert sha(folder/key)==SOURCE[key]
    assert sha(folder/"initial_flow.csv")==SOURCE["restart.csv"]
    log=(folder/"debugger.log").read_text(encoding="utf8",errors="replace")
    assert "Exit Success (SU2_CFD)" in log
    hist=folder/"history.csv"
    with hist.open(encoding="ascii") as f:
        keys=[v.strip() for v in next(csv.reader(f,skipinitialspace=True))]
        rows=np.loadtxt(f,delimiter=",",ndmin=2)
    assert len(rows)>=20 and np.isfinite(rows).all()
    np.testing.assert_array_equal(rows[:,keys.index("Inner_Iter")],np.arange(len(rows)))
    assert np.all(rows[:,keys.index("Time_Iter")]==0)
    index=len(rows)-1
    endpoint=folder/f"restart_{index:06d}.csv"
    assert endpoint.name in log and sha(endpoint)==sha(folder/"restart.csv")
    with np.load(folder/"nodes.npz") as z:
        xyz,ids=z["xyz"],z["ids"]
    values=frame(endpoint,xyz)
    assert sha(endpoint)!=SOURCE["restart.csv"]
    w=ids[:,:65].ravel()
    win=rows[-20:,[keys.index(key) for key in RMS]].max(axis=0)
    theta=values[:,9]
    p=101325+1.18*6**2*values[:,4]
    wall_error=float(np.max(abs(theta[w]-1)))
    no_slip=float(np.max(abs(values[w,5:8])))
    bounded=theta.min()>=-1e-12 and theta.max()<=1+1e-12
    accepted=win.max()<-8 and bounded and p.min()>0 and wall_error<=1e-12 and no_slip<=1e-12
    # 五条射线各自的解可能受角向离散影响，量化使用第零射线所舍弃的差异。
    drift={key:float(np.max(np.ptp(values[ids,col],axis=0))) for key,col in
           (("pressureOverRhoU2",4),("axialVelocityOverU",5),("temperatureScalar",9))}
    report=dict(acceptedForInitialization=bool(accepted),lastSavedIndex=index,innerRows=len(rows),
                last20MaxLogRMS=dict(zip(RMS,win)),temperatureScalarRange=[float(theta.min()),float(theta.max())],
                wallScalarError=wall_error,wallVelocityMaxOverU=no_slip,minPressurePa=float(p.min()),
                angularPeakToPeak=drift,frameSHA256=sha(endpoint),frameName=endpoint.name,
                historySHA256=sha(hist),configSHA256=sha(folder/"case.cfg"),
                elapsedSeconds=rec["elapsedSeconds"],peakProcessTreeRSSBytes=rec["peakFamilyRSSBytes"])
    return values, ids[0], report


def write_seed(path, base, mother, xyz, source, phi, perturbation):
    """全周节点按共同母线点取值；恢复径向速度向量，扰动只加到轴向速度。"""
    val=base[mother[source]]
    a=np.column_stack((np.arange(len(xyz)),xyz,val[:,4],val[:,5]+perturbation,
                       val[:,6]*np.cos(phi),val[:,6]*np.sin(phi),val[:,8:10]))
    with (PARENT/"restart.csv").open(encoding="ascii") as f:
        header=f.readline().strip()
    np.savetxt(path,a,delimiter=",",header=header,comments="",fmt=["%d"]+["%.17g"]*9)


def prepare(base_folder,out):
    """先回放 3 m/s 投影并逐字节核对，再在同一位置以新基流替换六个场量。"""
    base,mother,report=base_state(base_folder)
    save(out/"base_check.json",report)
    if not report["acceptedForInitialization"]:
        raise ValueError("同速基流未满足启动条件")
    frozen=read(PREVIOUS/"evidence/local_assets.json")["files"]
    for key in ("case.cfg","body.json","mesh.su2","nodes.npz","seed_00000.csv"):
        assert sha(FULL/key)==frozen["first02/"+key]["sha256"],key
    with np.load(FULL/"nodes.npz") as grid:
        xyz,ids=grid["xyz"],grid["ids"]
    assert np.array_equal(ids[0],ids[-1])
    unique,first=np.unique(ids[:-1].ravel(),return_index=True)
    np.testing.assert_array_equal(unique,np.arange(len(xyz)))
    ray,source=np.divmod(first,ids.shape[1])
    phi=np.linspace(0,2*np.pi,len(ids))[ray]
    xy=xyz[ids[0],:2]
    np.testing.assert_allclose(base[mother,1:3],xy,rtol=0,atol=1e-13)
    x,r=xy[source].T
    ellipse=((x+1.7)/.8)**2+((r-.55)/.45)**2
    perturbation=1e-3*np.maximum(1-ellipse,0.)**2*np.cos(phi)
    with np.load(PARENT/"nodes.npz") as native:
        old=frame(PARENT/"restart.csv",native["xyz"])
    replay=out/"reference_replay.csv"
    write_seed(replay,old,mother,xyz,source,phi,perturbation)
    assert sha(replay)==frozen["first02/seed_00000.csv"]["sha256"]
    seed=out/"seed_00000.csv"
    write_seed(seed,base,mother,xyz,source,phi,perturbation)
    geom=read(FULL/"geometry.json")
    assert geom["peakAxialPerturbationOverU"]==float(np.max(abs(perturbation)))
    geom["initialGuessFrameSHA256"]=report["frameSHA256"]
    save(out/"geometry.json",geom)
    save(out/"transfer.json",dict(referenceReplaySHA256=sha(replay),newSeedSHA256=sha(seed),
         matchedBaseSHA256=report["frameSHA256"],fullMeshSHA256=frozen["first02/mesh.su2"]["sha256"],
         frozenConfigSHA256=FULL_CONFIG_SHA,nodes=len(xyz),perturbationMaxOverU=float(np.max(abs(perturbation)))))


def supervised(stage,source,out):
    """映射和数值推进各用独有目录、实名进程与事前资源预算。"""
    out.mkdir(parents=True,exist_ok=False)
    files=[CONTRACT,Path(__file__).resolve(),ROOT/"scripts/thermistorMatchedBase.py",
           ROOT/"scripts/thermistorFullRing.py",ROOT/"scripts/runSU2NativeDebugger.py",
           PREVIOUS/"audit_step.py",PREVIOUS/"evidence/local_assets.json"]
    if stage=="prepare":
        files.extend(source/k for k in ("debugger_process.json","restart.csv","history.csv","case.cfg","nodes.npz"))
        files.extend(FULL/k for k in ("geometry.json","nodes.npz","mesh.su2","seed_00000.csv","case.cfg"))
        command=[sys.executable,str(Path(__file__).resolve()),"transfer",str(source),str(out)]
        seconds,limit=90,80
    else:
        terminal(source)
        receipt=read(source/"transfer.json")
        assert sha(ENGINE)==ENGINE_SHA and sha(source/"seed_00000.csv")==receipt["newSeedSHA256"]
        assert sha(FULL/"case.cfg")==FULL_CONFIG_SHA and sha(FULL/"mesh.su2")==receipt["fullMeshSHA256"]
        for path in [source/"geometry.json",source/"seed_00000.csv",FULL/"nodes.npz",FULL/"mesh.su2",
                     FULL/"body.json",FULL/"case.cfg"]:
            os.link(path,out/path.name)
            files.append(out/path.name)
        files.extend([ENGINE,source/"transfer.json",source/"debugger_process.json",source/"base_check.json"])
        command=[str(ENGINE),"case.cfg","--threads","1"]
        seconds,limit=360,180
    result=observe(command,out,{str(p):sha(p) for p in files},seconds,limit*1024**2)
    return result["debuggerExitCode"]==0 and result["stopReason"] is None


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=("prepare","transfer","solve","audit"))
    parser.add_argument("source",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    tree=ROOT/"outputs/v288_validation"
    source,out=args.source.resolve(),args.output.resolve()
    source.relative_to(tree)
    out.relative_to(tree)
    if args.stage=="transfer":
        prepare(source,out)
    elif args.stage=="audit":
        save(out,audit(source))
    else:
        raise SystemExit(0 if supervised(args.stage,source,out) else 1)
