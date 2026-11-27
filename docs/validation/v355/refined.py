"""单独加密固定低速风洞的热网，复用已合格的同一组全周动量终帧。"""
import argparse
import json
import os
from pathlib import Path
import sys

for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[name]="1"

import numpy as np

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
PREVIOUS=DOC.parent/"v354"
sys.path[:0]=[str(PREVIOUS),str(ROOT/"scripts")]
from flux_transfer import read
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe
from simpleLauncherAngularHeat import SECTORS,domain,heat_rows,steady_fields,wall_results
from simpleLauncherTemperature import evaluate


def save(path,record):
    path.write_bytes((json.dumps(record,allow_nan=False,indent=2,ensure_ascii=False)+"\n").encode())


def wall_resolution(path,diameter):
    """由实际贴壁三角形量取边长和相对该边的法向高度。"""
    with np.load(path) as z:
        p=z["xy"];tri=z["triangles"];wet=z["wall"]
    boundary={tuple(sorted((i,j))):k for k,(i,j) in enumerate(zip(wet[:-1],wet[1:]))}
    length=np.linalg.norm(np.diff(p[wet],axis=0),axis=1)*diameter
    height=np.full(len(length),np.nan)
    for side in range(3):
        a,b,c=tri[:,side],tri[:,(side+1)%3],tri[:,(side+2)%3]
        wall_pair=np.isin(a,wet)&np.isin(b,wet)
        for i,j,k in zip(a[wall_pair],b[wall_pair],c[wall_pair]):
            n=boundary.get(tuple(sorted((i,j))))
            if n is not None:
                assert np.isnan(height[n])
                height[n]=abs(np.cross(p[j]-p[i],p[k]-p[i]))*diameter**2/length[n]
    assert np.isfinite(height).all() and height.min()>0
    return dict(wallEdges=len(length),meridianNodes=len(p),triangles=len(tri),
                medianWallEdgeMm=1000*float(np.median(length)),medianFirstHeightMm=1000*float(np.median(height)))


def prepare(upstream,work):
    ref=read(DOC.parent/"v352/baseline.json")["coarse"];coarse=upstream/"heat_coarse"
    assert sha(upstream/"coarse/restart.csv")==ref["flowSHA256"] and sha(coarse/"temperature.npz")==ref["heatMeshSHA256"]
    native=read(upstream/"coarse/condition.json");assert native["grid"]=="coarse" and native["extentD"]==14
    fine=work/"steady";initial=work/"initial"
    evaluate(upstream/"coarse",fine,.002)
    seed=read(fine/"summary.json")
    assert seed["heatSpacingD"]==.002 and seed["sourceRestartSHA256"]==ref["flowSHA256"]
    fine_mesh=sha(fine/"temperature.npz")
    model=domain(upstream,"coarse",dict(flowSHA256=ref["flowSHA256"],heatMeshSHA256=fine_mesh),fine)
    wall=np.ones((1,SECTORS,len(model["wet"])))
    theta=steady_fields(model,wall);reactions,balance=heat_rows(model,theta[0])
    flux,integrals=wall_results(model,reactions)
    with np.load(fine/"temperature.npz") as z:diff=np.max(abs(theta[0]-z["theta"]))
    wall_table=np.loadtxt(fine/"wall.csv",delimiter=",",skiprows=1)
    area=wall_table[:,1]
    rms=np.sqrt((area@(flux-wall_table[:,3])**2)/(area@wall_table[:,3]**2))
    resolutions={k:wall_resolution(p/"temperature.npz",model["diameter"]) for k,p in (("coarse",coarse),("fine",fine))}
    gates=dict(refined=all(resolutions["fine"][k]<resolutions["coarse"][k] for k in ("medianWallEdgeMm","medianFirstHeightMm")),
        bounds=bool(theta.min()>=-1e-12 and theta.max()<=1+1e-12),referenceTemperature=bool(diff<1e-9),
        residual=bool(balance["relativeFreeRowResidual"]<1e-8),balance=bool(balance["energyImbalance"]<.005),
        baselineWall=bool(rms<.02),heatStability=bool(seed["heatStability"]<.002))
    initial.mkdir()
    np.savez_compressed(initial/"fields.npz",theta=theta)
    result=dict(accepted=all(gates.values()),gates=gates,sourceFlowSHA256=ref["flowSHA256"],sourceHeatMeshSHA256=fine_mesh,
                fieldsSHA256=sha(initial/"fields.npz"),resolutions=resolutions,heatStability=seed["heatStability"],
                thermalDOF=SECTORS*(len(model["xy"])-len(model["axis"]))+len(model["axis"]),
                referenceThetaDifference=float(diff),wallFluxRMS=float(rms),**integrals,**balance)
    save(initial/"results.json",result)
    assert result["accepted"],result
    return result


def launch_step(upstream,momentum,prepared,work,name):
    heat=work/"heat";dest=heat/name
    step=read(momentum/name/"step.json");frame=momentum/name/f"restart_{step['step']:05d}.csv"
    assert step["accepted"] and sha(frame)==step["frameSHA256"]
    index=ROOT/"cases/simple_launcher/wake/flow_index.json"
    assert all(step[k]==v for k,v in read(index)[name].items())
    dest.mkdir(parents=True,exist_ok=False)
    initial=work/"initial/fields.npz";fine=work/"steady"
    paths=[DOC/"PLAN.md",PREVIOUS/"heat_step.py",PREVIOUS/"flux_transfer.py",ROOT/"scripts/simpleLauncherAngularHeat.py",
           index,frame,momentum/name/"step.json",prepared/"prepared.json",prepared/"control.csv",prepared/"stations.npz",
           initial,initial.parent/"results.json",fine/"temperature.npz",fine/"summary.json",upstream/"coarse/restart.csv"]
    if step["step"]==2:paths.append(heat/f"{step['branch']}_half_1/thermal.npz")
    cmd=[sys.executable,str(PREVIOUS/"heat_step.py"),str(upstream),str(momentum),str(prepared),str(heat),name,
         "--initial-field",str(initial),"--heat-folder",str(fine)]
    rec=observe(cmd,dest,{str(p):sha(p) for p in paths},240,300*1024**2)
    if rec["debuggerExitCode"] or rec["stopReason"] or rec["survivingProcessIDs"]:
        raise RuntimeError(f"{name} 未通过，详见 {dest/'debugger.log'}")
    result=read(dest/"heat.json")
    assert result["accepted"] and result["heatMeshSHA256"]==sha(fine/"temperature.npz")
    assert result["sourceFrameSHA256"]==sha(frame) and sha(dest/"thermal.npz")==result["fieldSHA256"]
    print(name,result["outwardW"],rec["elapsedSeconds"],rec["peakFamilyRSSBytes"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("upstream","momentum","prepared","work"):parser.add_argument(key,type=Path)
    parser.add_argument("stage",choices=("prepare","pilot","remaining"))
    args=parser.parse_args();up,mom,prep,work=[getattr(args,k).resolve() for k in ("upstream","momentum","prepared","work")]
    work.relative_to(ROOT/"outputs")
    if args.stage=="prepare":prepare(up,work)
    elif args.stage=="pilot":launch_step(up,mom,prep,work,"control_full_1")
    else:
        pilot=work/"heat/control_full_1"
        record=read(pilot/"debugger_process.json");heat=read(pilot/"heat.json")
        assert record["status"]=="terminal" and record["debuggerExitCode"]==0 and all(record["inputsStable"].values())
        assert not (record["stopReason"] or record["survivingProcessIDs"] or record["supervisorError"])
        assert record["peakFamilyRSSBytes"]<4*1024**3 and heat["accepted"]
        assert heat["fieldSHA256"]==sha(pilot/"thermal.npz") and heat["wallSHA256"]==sha(pilot/"wall.csv")
        assert heat["heatMeshSHA256"]==sha(work/"steady/temperature.npz")
        for n in ("disturbed_full_1","control_half_1","control_half_2","disturbed_half_1","disturbed_half_2"):
            launch_step(up,mom,prep,work,n)
