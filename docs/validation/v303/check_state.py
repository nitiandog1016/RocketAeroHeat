"""从实际冷态、周期帧和终点逐点检查温界、压力和壁条件。"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import D,P0,RHO,U,read,sha
from thermistorReaction import terminal
from beadSurfaceResult import thermal_frame

BASE=ROOT/"outputs/v300_validation/wired01"
NEW=ROOT/"outputs/v303_validation/outer01"


def save(path,value):
    with path.open("x",encoding="utf8") as f:json.dump(value,f,indent=2,ensure_ascii=False,allow_nan=False)


def audit(folder,live=False):
    locks=read(folder/"inputs.json");assert all(sha(folder/k)==v for k,v in locks["inputs"].items())
    with (folder/"history.csv").open() as f:
        lines=list(csv.reader(f,skipinitialspace=True))
    col=[name.strip() for name in lines[0]].index("Inner_Iter")
    last=int([r for r in lines[1:] if len(r)==len(lines[0])][-1][col])
    rec=None if live else terminal(folder)
    with np.load(folder/"nodes.npz") as z:xyz=z["xyz"];wall=np.unique(z["wallMap"]);near=np.unique(z["firstMap"])
    frames=[]
    paths=[folder/"initial_flow.csv",*sorted(folder.glob("restart_*.csv"))]
    for path in paths:
        step=None if path==paths[0] else int(path.stem.split("_")[-1])
        if live and step is not None and step>=last:continue
        raw=np.loadtxt(path,delimiter=",",skiprows=1)
        assert raw.shape==(len(xyz),10) and np.isfinite(raw).all()
        np.testing.assert_array_equal(raw[:,0],np.arange(len(xyz)))
        np.testing.assert_allclose(raw[:,1:4],xyz,rtol=0,atol=1e-13)
        theta=raw[:,9];j=int(np.argmin(theta));outside=(theta < -1e-12)|(theta > 1+1e-12)
        f=dict(file=path.name,step=step,sha256=sha(path),thetaRange=[float(theta.min()),float(theta.max())],
               minimumNode=j,minimumXYZM=(D*xyz[j]).tolist(),below=int((theta< -1e-12).sum()),above=int((theta>1+1e-12).sum()),
               outsideInnerDomain=int(outside[:16581].sum()),outsideAtWall=int(outside[wall].sum()),outsideAtFirstLayer=int(outside[near].sum()),
               minimumPressurePa=float(P0+RHO*U**2*raw[:,4].min()),maxWallSpeedMps=float(U*np.linalg.norm(raw[wall,5:8],axis=1).max()),
               maxAbsWallThetaError=float(np.max(abs(theta[wall]-1))))
        if step is not None:
            vtk=folder/f"surface_{step:06d}.vtk";f["surfaceSHA256"]=sha(vtk)
            if rec:assert rec["outputSHA256"][vtk.name]==f["surfaceSHA256"] and rec["outputSHA256"][path.name]==f["sha256"]
        frames.append(f)
    if rec:
        assert [f["step"] for f in frames]==[None,*range(500,last,500),last]
        assert frames[-1]["sha256"]==sha(folder/"restart.csv") and frames[-1]["surfaceSHA256"]==sha(folder/"surface.vtk")
    good=all(f["below"]==f["above"]==0 and f["minimumPressurePa"]>0 and f["maxAbsWallThetaError"]==0 and
             f["maxWallSpeedMps"]<1e-10*U for f in frames)
    return dict(directory=folder.relative_to(ROOT).as_posix(),frames=frames,fieldChecksPassed=good,lastHistoryIteration=last,
                residualWindow=thermal_frame(folder)[1] if rec else None,nativeSeconds=rec["elapsedSeconds"] if rec else None)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--live",action="store_true");p.add_argument("--out",type=Path)
    a=p.parse_args();result=audit(NEW,a.live)
    if a.out:save(a.out,result)
    print([(f["step"],f["thetaRange"],f["below"],f["above"]) for f in result["frames"]])
    if not result["fieldChecksPassed"]:sys.exit(1)
