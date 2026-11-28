"""0.05 m/s 风洞在原子午网内加点，延拓同一批初场，并按登记预算运行独立的动量时间步。"""
import argparse
import json
import os
from pathlib import Path
import sys

for env in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
    os.environ[env]="1"
import numpy as np

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
sys.path[:0]=[str(DOC.parent/"v352"),str(DOC.parent/"v354"),str(ROOT/"scripts")]
from cm011TunnelSector import mesh
from prepare_wake import write_full_grid
from runEmbeddedSU2 import sha
from run_flow import run as flow_step
from simpleLauncherTunnel import geometry

OLD=ROOT/"outputs/v352_validation/prepared01"
UP=ROOT/"outputs/launcher_6639759c-6b80-41b7-bb07-55437bcd3a4a"
SEED=ROOT/"outputs/v355_validation/run01"
BASE=ROOT/"outputs/v357_validation"


def read(p):
    return json.loads(p.read_bytes())


def save(p,data):
    p.write_bytes((json.dumps(data,ensure_ascii=False,allow_nan=False,indent=2)+"\n").encode())


def indices(width):
    """新增点均在原三角形棱边上。双向加密的四边形中心取原 0–2 对角线中点。"""
    layers=np.r_[np.arange(32)/2,np.arange(16,37)]
    cols=np.arange(2*width-1)/2
    row,col=np.meshgrid(layers,cols,indexing="ij")
    a=row.astype(int)*width+col.astype(int)
    b=np.ceil(row).astype(int)*width+np.ceil(col).astype(int)
    return a,b,layers,cols


def resolution(p,d):
    direction=np.diff(p[0],axis=0)
    edges=np.linalg.norm(direction,axis=1)
    first=p[1]-p[0]
    h=abs(np.cross(direction,(first[:-1]+first[1:])/2))/edges
    return dict(wallEdgeMedianMm=float(np.median(edges)*d*1000),
                firstNormalMedianMm=float(np.median(h)*d*1000),
                firstRayMedianMm=float(np.median(np.linalg.norm(first,axis=1))*d*1000))


def prepare(dest):
    dest.mkdir(parents=True,exist_ok=False)
    definition,tri,d,wall=geometry()
    coarse=read(OLD/"prepared.json")
    assert coarse["files"]==read(DOC.parent/"v352/evidence.json")["geometry"]["files"]
    for f,value in coarse["files"].items():assert sha(OLD/f)==value
    frame=UP/"coarse/restart.csv"
    assert sha(frame)==coarse["parentFrameSHA256"]==read(DOC.parent/"v352/baseline.json")["coarse"]["flowSHA256"]
    frozen=read(DOC.parent/"v355/evidence.json")["initial"]
    for path,key in ((SEED/"steady/temperature.npz","sourceHeatMeshSHA256"),(SEED/"initial/fields.npz","fieldsSHA256")):
        assert sha(path)==frozen[key]
    old=np.loadtxt(frame,delimiter=",",skiprows=1);width=len(wall)
    a,b,layers,cols=indices(width)
    extended=(old[a]+old[b])/2
    extended[:,:,0]=np.arange(a.size).reshape(a.shape)
    fine=extended[:,:,1:3];coarse_xy=old[:,1:3].reshape(37,width,2)
    np.testing.assert_array_equal(fine[layers==layers.astype(int)][:,::2],coarse_xy)
    np.testing.assert_array_equal(fine[-1,1::2],(coarse_xy[-1,:-1]+coarse_xy[-1,1:])/2)
    sector=mesh(np.array([0.,2*np.pi/48]),layers=len(layers)-1,wall=fine[0],extent=14,diameter=d,meridian_nodes=fine)
    old_shape=old[:,2]==0
    lookup=np.empty((48,len(old)),int)
    lookup[:,old_shape]=np.arange(old_shape.sum())
    lookup[:,~old_shape]=old_shape.sum()+np.arange(48*(~old_shape).sum()).reshape(48,-1)
    control=np.loadtxt(OLD/"control.csv",delimiter=",",skiprows=1,usecols=(5,6,7))
    disturbed=np.loadtxt(OLD/"disturbed.csv",delimiter=",",skiprows=1,usecols=(5,6,7))
    diff=disturbed-control
    def extend_perturbation(ids):
        delta=np.empty((ids.max()+1,3))
        for j,line in enumerate(ids):
            delta[line]=.5*(diff[lookup[j,a.ravel()]]+diff[lookup[j,b.ravel()]])
        return delta
    out=write_full_grid(dest,definition,tri,d,fine[0],sector,extended.reshape(-1,old.shape[1]),frame,
                        extend_perturbation,[Path(__file__).resolve(),DOC/"PLAN.md"])
    assert out["nodes"]==547066 and out["cells"]==539136
    assert abs(out["volumeD3"]/coarse["volumeD3"]-1)<1e-12 and abs(out["wallAreaM2"]-coarse["wallAreaM2"])<1e-12
    assert abs(out["peakPerturbationOverU"]-.01)<1e-14
    # 在所有共有三维节点上核对实际写出的两支种子，各自保留原节点顺序。
    axis=fine[:,:,1].ravel()==0
    ids=np.empty((48,a.size),int);ids[:,axis]=np.arange(axis.sum())
    ids[:,~axis]=axis.sum()+np.arange(48*(~axis).sum()).reshape(48,-1)
    old_to_new=np.empty(len(control),int)
    old_to_new[lookup]=ids.reshape(48,len(layers),len(cols))[:,layers==layers.astype(int),::2].reshape(48,-1)
    for label in ("control","disturbed"):
        prev=np.loadtxt(OLD/(label+".csv"),delimiter=",",skiprows=1)
        current=np.loadtxt(dest/(label+".csv"),delimiter=",",skiprows=1)
        np.testing.assert_allclose(current[old_to_new,1:],prev[:,1:],rtol=0,atol=1e-14)
    np.savez_compressed(dest/"meridian.npz",xy=fine,station=cols,oldToNew=old_to_new)
    result=dict(nodes=out["nodes"],cells=out["cells"],wallFaces=out["wallFaces"],newFiles=out["files"],
                oldFiles=coarse["files"],meridianSHA256=sha(dest/"meridian.npz"),oldResolution=resolution(coarse_xy,d),
                refinedResolution=resolution(fine,d),savedInitial=frozen["fieldsSHA256"],heatMesh=frozen["sourceHeatMeshSHA256"],
                volumeRelativeChange=out["volumeD3"]/coarse["volumeD3"]-1,oldNodesChecked=len(old_to_new),
                maximumPatchAreaRelativeDifference=out["maxPatchRelativeDifference"],
                predictedPeakWithMarginBytes=1212968960*out["nodes"]/coarse["nodes"]*1.15)
    assert result["predictedPeakWithMarginBytes"]<4*1024**3
    save(dest/"refinement.json",result)
    return result


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("root",type=Path);ap.add_argument("stage",choices=("prepare","first","second"))
    args=ap.parse_args();work=args.root.resolve();work.relative_to(BASE)
    dest=work/"prepared"
    if args.stage=="prepare":
        prepare(dest)
    else:
        proof=read(dest/"refinement.json")
        assert sha(dest/"meridian.npz")==proof["meridianSHA256"] and read(dest/"prepared.json")["files"]==proof["newFiles"]
        n=1 if args.stage=="first" else 2
        if n==2:assert read(work/"flow/disturbed_half_1/step.json")["accepted"]
        sys.exit(0 if flow_step(work/"flow","disturbed","half",n,dest,720,600,DOC/"PLAN.md") else 2)
