"""逐节点检查实际保存场及强标量温界，在同迭代已有帧上分离压力、速度与温差变化。"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from transport_case import BASE, HERE, ROOT, candidate_config, read, sha
from thermistorBead import D, DT, P0, RHO, T0, U


def frame(folder, path, xyz, walls, near):
    raw=np.loadtxt(path,delimiter=",",skiprows=1)
    assert raw.shape==(len(xyz),10) and np.isfinite(raw).all()
    np.testing.assert_array_equal(raw[:,0],np.arange(len(xyz)))
    np.testing.assert_allclose(raw[:,1:4],xyz,atol=1e-13,rtol=0)
    theta=raw[:,9]
    def point(j):
        return dict(id=int(j),xyzM=(xyz[j]*D).tolist(),theta=float(theta[j]),temperatureK=float(T0+DT*theta[j]),
                    atWall=bool(j in walls),atFirstLayer=bool(j in near))
    low=theta < -1e-12
    high=theta > 1+1e-12
    bad=low|high
    result=dict(file=path.name,sha256=sha(path),nodes=len(theta),min=point(np.argmin(theta)),max=point(np.argmax(theta)),
                below=int(low.sum()),above=int(high.sum()),belowAtWall=int(low[walls].sum()),aboveAtWall=int(high[walls].sum()),
                belowAtFirstLayer=int(low[near].sum()),aboveAtFirstLayer=int(high[near].sum()),
                outOfBoundsXRangeM=(D*np.asarray([xyz[bad,0].min(),xyz[bad,0].max()])).tolist() if bad.any() else [],
                minAbsolutePressurePa=float(P0+RHO*U**2*raw[:,4].min()),wallSpeedMaxMps=float(U*np.linalg.norm(raw[walls,5:8],axis=1).max()),
                maxWallScalarError=float(np.max(abs(theta[walls]-1))))
    base_path=BASE/path.name
    if base_path.is_file():
        base_rec=read(BASE/"engine_process.json")
        digests=base_rec["outputSHA256"] if path.name!="initial_flow.csv" else {path.name:read(folder/"inputs.json")["unchanged"][path.name]}
        assert sha(base_path)==digests[path.name]
        other=np.loadtxt(base_path,delimiter=",",skiprows=1)
        np.testing.assert_array_equal(other[:,:4],raw[:,:4])
        result["sameIteration"]=dict(baseSHA256=sha(base_path),maxAbsPressureChangePa=float(RHO*U**2*np.max(abs(raw[:,4]-other[:,4]))),
            maxVelocityVectorChangeMps=float(U*np.linalg.norm(raw[:,5:8]-other[:,5:8],axis=1).max()),
            maxTemperatureChangeK=float(DT*np.max(abs(theta-other[:,9]))),baseThetaRange=[float(other[:,9].min()),float(other[:,9].max())])
    return result


def check(folder, complete=True):
    assert (folder/"case.cfg").read_text("ascii")==candidate_config()
    spec=read(folder/"inputs.json")
    assert all(sha(folder/k)==v for k,v in spec["unchanged"].items())
    rec=read(folder/"engine_process.json")
    with (folder/"history.csv").open() as stream:
        lines=list(csv.reader(stream,skipinitialspace=True))
    labels=[s.strip() for s in lines[0]];history=[r for r in lines[1:] if len(r)==len(labels)]
    last=int(history[-1][labels.index("Inner_Iter")])
    with np.load(folder/"nodes.npz") as z:
        xyz=z["xyz"];walls=np.unique(z["wallMap"]);near=np.unique(z["firstMap"])
    frames=[]
    for path in sorted(folder.glob("restart_*.csv")):
        step=int(path.stem.split("_")[-1]);vtk=folder/f"surface_{step:06d}.vtk"
        if not complete and (step>=last or not vtk.is_file()):continue
        assert vtk.is_file()
        if complete:
            assert rec["status"]=="terminal" and rec["outputSHA256"][path.name]==sha(path) and rec["outputSHA256"][vtk.name]==sha(vtk)
        record=frame(folder,path,xyz,walls,near);record.update(step=step,surfaceSHA256=sha(vtk));frames.append(record)
    if complete and rec["exitCode"]==0 and not rec["cancelled"]:
        assert [f["step"] for f in frames]==[*range(500,last,500),last]
        assert sha(folder/"restart.csv")==frames[-1]["sha256"] and sha(folder/"surface.vtk")==frames[-1]["surfaceSHA256"]
    keys=["rms["+s+"]" for s in ("P","U","V","W","h","rho*Y_0")]
    assert len(history)>=20 and [int(r[labels.index("Inner_Iter")]) for r in history[-20:]]==list(range(last-19,last+1))
    residual={k:max(float(r[labels.index(k)]) for r in history[-20:]) for k in keys}
    return dict(initial=frame(folder,folder/"initial_flow.csv",xyz,walls,near),frames=frames,
                lastHistoryIteration=last,residualWindow=residual,residualPass=max(residual.values()) < -8,
                checkpointsWithinBounds=all(f["above"]+f["below"]==0 for f in frames) if frames else None,
                nativeStatus=rec["status"],exitCode=rec.get("exitCode"),cancelled=rec.get("cancelled"),terminalAudit=complete)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("folder",type=Path);p.add_argument("--live",action="store_true")
    p.add_argument("--out",type=Path);a=p.parse_args();folder=a.folder.resolve();folder.relative_to(ROOT/"outputs")
    result=check(folder,not a.live)
    if a.out:
        dest=a.out.resolve();dest.relative_to(ROOT/"outputs")
        with dest.open("x",encoding="utf8") as f:json.dump(result,f,indent=2,allow_nan=False)
    print([(f["step"],f["min"]["theta"],f["max"]["theta"],f["below"],f["above"]) for f in result["frames"]])
