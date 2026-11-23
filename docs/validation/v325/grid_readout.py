"""从实际求解输入及原生回执提取可随工作台发布的两个网格剖面。"""
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from shortSupportMesh import BODY_HASH, PHASE, read, sha
from thermistorOuterDomain import load_mesh
from thermistorMeasuredCold import RESIDUALS


def meridian(mesh):
    cells,xyz,marks=load_mesh(mesh,("wall","far"))
    r=xyz[:,0]*np.cos(PHASE)-xyz[:,2]*np.sin(PHASE)
    off=xyz[:,0]*np.sin(PHASE)+xyz[:,2]*np.cos(PHASE)
    ids=np.flatnonzero((abs(off)<1e-10)&(r>=-1e-12))
    local=np.full(len(xyz),-1,int);local[ids]=np.arange(len(ids))
    p=np.c_[xyz[ids,1],r[ids]]
    mask=local>=0
    polygons=set()
    for cell in cells:
        on=tuple(sorted(local[v] for v in cell if mask[v]))
        if len(on)>=3:polygons.add(on)
    edges=set()
    for poly in polygons:
        nodes=np.array(poly);v=p[nodes];y,z=(v-v.mean(0)).T
        loop=nodes[np.argsort(np.arctan2(z,y))]
        for a,b in zip(loop,np.roll(loop,1)):edges.add((min(a,b),max(a,b)))
    return dict(pointsD=p.tolist(),edges=(np.array(sorted(edges))+1).tolist(),pointIDs=ids.tolist(),
                planeDegrees=float(np.degrees(PHASE)),sectionCells=len(polygons))


def scan(prep,run):
    s=read(prep/"prepared.json");qualified=read(run/"qualified.json")
    assert qualified["prepared"]==sha(prep/"prepared.json") and sha(prep/"mesh.su2")==s["files"]["mesh.su2"]
    with np.load(prep/"mapping.npz") as m:
        wet=m["wet"];xyz=m["xyz"];theta=m["theta"];tree=cKDTree(xyz[wet])
    stages=[];source=s["files"]["initial.csv"]
    for dest in sorted(run.glob("steady[0-9][0-9]")):
        r=read(dest/"debugger_process.json")
        assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
        assert not(r["survivingProcessIDs"] or r["stopReason"] or r["supervisorError"])
        assert r["inputSHA256"][str(dest/"initial_flow.csv")]==sha(dest/"initial_flow.csv")==source
        assert "Exit Success (SU2_CFD)" in (dest/"debugger.log").read_text("ascii",errors="replace")
        with (dest/"history.csv").open() as stream:
            keys=[k.strip() for k in next(csv.reader(stream,skipinitialspace=True))]
            history=np.loadtxt(stream,delimiter=",",ndmin=2)
        np.testing.assert_array_equal(history[:,keys.index("Inner_Iter")],np.arange(len(history)))
        frames=[]
        for i in sorted({len(history)-1,*range(100,len(history),100)}):
            path=dest/f"restart_{i:06d}.csv"
            a=np.loadtxt(path,skiprows=1,delimiter=",",usecols=range(10))
            np.testing.assert_allclose(a[:,1:4],xyz,rtol=0,atol=1e-12)
            np.testing.assert_array_equal(a[:,0],np.arange(len(xyz)));assert np.isfinite(a).all()
            np.testing.assert_allclose(a[wet,9],theta,atol=2e-14,rtol=0)
            assert abs(a[wet,5:8]).max()<1e-10
            t=a[:,9];lower=t< -1e-12;upper=t>theta.max()+1e-12;bad=lower|upper;worst=int(t.argmin())
            distances=tree.query(xyz[bad])[0] if bad.any() else np.array([])
            frames.append(dict(index=i,sha256=sha(path),minScalar=t.min(),maxScalar=t.max(),
                               below=int(lower.sum()),above=int(upper.sum()),worstPointID=worst,
                               worstXYZD=xyz[worst].tolist(),badNearestWallNodeD=[distances.min(),distances.max()] if bad.any() else [],
                               minAbsolutePa=float(np.min(s["condition"]["airPa"]+s["condition"]["rho"]*s["condition"]["speedMps"]**2*a[:,4]))))
        source=frames[-1]["sha256"];assert sha(dest/"restart.csv")==source
        window={k:history[-20:,keys.index(k)].max() for k in RESIDUALS}
        stages.append(dict(stage=dest.name,iterations=len(history),window=window,frames=frames,seconds=r["elapsedSeconds"],
                           maxRSSBytes=r["peakFamilyRSSBytes"],nativeReceiptSHA256=sha(dest/"debugger_process.json")))
    return dict(meshSHA256=s["files"]["mesh.su2"],transferSHA256=s["files"]["transfer.npz"],geometry=s["geometry"],
                stages=stages,section=meridian(prep/"mesh.su2"),equatorK=s["equatorK"],qualification=qualified["manufacturedPowerErrorsW"])


if __name__=="__main__":
    base=ROOT/"outputs/v325_validation"
    all_runs={name:scan((base/prep/"prepared").resolve(),(base/run).resolve())
              for name,prep,run in [("hybrid","prep02","run01"),("structured","prep03","run02")]}
    transfer=Path(sys.argv[2])
    with np.load(base/"prep02/prepared/mapping.npz") as m,np.load(base/"prep02/prepared/transfer.npz") as t,transfer.open("xb") as stream:
        wet=m["wet"]
        np.savez_compressed(stream,**{k:t[k] for k in t.files},nativeCorners=np.searchsorted(wet,m["corners"]),
                            nativeWallXYZD=m["xyz"][wet],firstXYZD=m["xyz"][m["first"]],nativeArea=m["nativeArea"],nativeRegion=m["group"])
    out=dict(bodySHA256=BODY_HASH,caseSHA256=sha(ROOT/"cases/short_sphere/test65.json"),**all_runs,
             publicTransferSHA256=sha(transfer))
    with Path(sys.argv[1]).open("xb") as f:f.write(json.dumps(out,separators=(",",":"),allow_nan=False).encode("utf8"))
