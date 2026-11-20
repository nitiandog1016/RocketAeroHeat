"""对终场检查对称面速度、近壁网格尺寸、原表面闭合性和极值位置，登记本机证据。"""
import csv
import gzip
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import D, read, sha
from thermistorBeadSector import RAYS


def write(path, data):
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf8")


def audit():
    base=ROOT/"outputs/v282_validation";out=Path(__file__).resolve().parent/"evidence"
    processes={};originals={};assets={}
    for path in sorted(base.rglob("engine_process.json")):
        folder=path.parent;rec=read(path)
        assert rec["status"]=="terminal" and not rec.get("supervisorError") and rec["exitCode"]==0
        assert not rec["timedOut"] and not rec["cancelled"] and all(rec["inputsStable"].values())
        name=folder.relative_to(base).as_posix()
        originals[name]=rec
        with (folder/"history.csv").open() as f:
            table=csv.DictReader(f,skipinitialspace=True)
            table.fieldnames=[k.strip() for k in table.fieldnames];history=list(table)
        processes[name]=dict(elapsedSeconds=rec["elapsedSeconds"],lastInnerIteration=int(history[-1]["Inner_Iter"]),
            resources=read(folder/"resources.json"),receiptSHA256=sha(path),terminalSHA256=sha(folder/"restart.csv"),
            savedRestartFiles=sorted(p.name for p in folder.glob("restart_*.csv")))
        assert processes[name]["terminalSHA256"]==rec["outputSHA256"]["restart.csv"]
    (out/"native_receipts.json.gz").write_bytes(gzip.compress(json.dumps(originals,ensure_ascii=False).encode("utf8")))
    for name in ["source01/tsubouchi1959.pdf","source01/p815.png","source01/figure7.png","scalar_observer01/energy.json",
                 "scalar_observer01/observed/original_scalar.csv","scalar_observer01/observed/scalar_stages.csv",
                 "ui02/export/workbench/index.html","ui02/export/workbench/case-data.js","ui02/wall_pair.csv"]:
        file=base/name
        assets[name]=dict(bytes=file.stat().st_size,sha256=sha(file))
    body=read(base/"strong_outer01/body.json")
    vertices,faces=np.asarray(body["V"]),np.asarray(body["F"])-1
    edges=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]])
    _,inv,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    assert np.all(counts==2) and np.all(np.bincount(inv,weights=np.sign(edges[:,1]-edges[:,0]))==0)
    grids={}
    for g in ("coarse","fine","outer"):
        folder=base/f"strong_{g}01";case=read(folder/"case.json")
        with np.load(folder/"nodes.npz") as grid:
            xyz,ids=grid["xyz"],grid["ids"]
        m=case["meridianIntervals"]+1;wall=ids[:,:m]
        sep=D*np.linalg.norm(xyz[ids[:,m:2*m]]-xyz[wall],axis=2)
        native=np.loadtxt(folder/"restart.csv",delimiter=",",skiprows=1)
        symmetry=max(np.max(abs(native[ids[k],7]*np.cos(RAYS[k])-native[ids[k],6]*np.sin(RAYS[k]))) for k in (0,-1))
        assert symmetry<1e-8 and np.max(abs(native[wall,5:8]))<1e-8
        table=np.genfromtxt(base/f"result_{g}01/wall.csv",delimiter=",",names=True)
        imin=int(np.argmin(table["heatFluxWm2"]));imax=int(np.argmax(table["heatFluxWm2"]))
        grids[g]=dict(firstNodeDistanceMicrons=[sep.min()*1e6,sep.max()*1e6],maxSymmetryNormalVelocityOverU=symmetry,
                      maxWallSpeedComponentOverU=np.max(abs(native[wall,5:8])),
                      mostNegativeFluxFaceID=imin,faceCentroidM=vertices[faces[imin]].mean(axis=0).tolist(),
                      leastNegativeFluxFaceID=imax,leastNegativeFluxCentroidM=vertices[faces[imax]].mean(axis=0).tolist())
    write(out/"local_assets.json",dict(localRoot="outputs/v282_validation",assets=assets,processes=processes,
                                      grids=grids,closedEdges=len(counts),originalFaces=len(faces)))


if __name__=="__main__":
    audit()
