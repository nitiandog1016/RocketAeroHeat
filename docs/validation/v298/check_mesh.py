"""以两个实际 SU2 文件的原坐标和单元连接证明远区相同，登记新工况的启动输入。"""
from pathlib import Path
import json
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import A,D,read,sha


def cells(directory):
    with (directory/"mesh.su2").open() as f:
        assert next(f).strip()=="NDIME= 3"
        rows=[list(map(int,next(f).split())) for _ in range(int(next(f).split("=")[1]))]
        xyz=np.array([list(map(float,next(f).split())) for _ in range(int(next(f).split("=")[1]))])
    distance=np.hypot(abs(xyz[:,0])-A/D,np.hypot(xyz[:,1],xyz[:,2]))
    points=list(map(tuple,xyz))
    inside=set();outside=set()
    for kind,*indices in rows:
        key=(kind,tuple(points[i] for i in indices))
        (outside if distance[indices].min()>=.45 else inside).add(key)
    assert len(inside)+len(outside)==len(rows)
    return inside,outside


def shared_faces(wired,smooth):
    def triangles(body,region):
        v=np.array(body["V"]);f=np.array(body["F"])-1;g=np.array(body["regionID"])
        # 保留有向三角面，直接配对原面真实顶点。
        return {tuple(map(tuple,v[f[i]])):int(i) for i in np.flatnonzero(g==region)}
    left=triangles(wired,4);right=triangles(smooth,2)
    assert len(left)==len(right)==6144 and left.keys()==right.keys()
    return np.array([(left[k],right[k]) for k in left])


def run():
    wired=ROOT/"outputs/v297_validation/fine01"
    smooth=ROOT/"outputs/v298_validation/smooth01"
    old=read(ROOT/"docs/validation/v297/evidence/fine_inputs.json")
    for name in ("body.json","mesh.su2","nodes.npz","case.cfg","case.json"):
        assert sha(wired/name)==old["inputs"][name]
    assert (wired/"case.cfg").read_bytes()==(smooth/"case.cfg").read_bytes()
    cw,cs=read(wired/"case.json"),read(smooth/"case.json")
    for key in ("sourceAreaM2","diameterM","speedMps","airPa","airK","wallK","rho","mu","cp","k",
                "sectorDeg","angularCells","firstLayerM","outerRadiusOverD"):
        assert cw[key]==cs[key]
    wi,wo=cells(wired);si,so=cells(smooth)
    # 判据按旧单元定义。锯齿状切口内新生成的三角形可以有三个顶点都在半径之外。
    assert wo.issubset(so) and len(wo)>0
    common=shared_faces(read(wired/"body.json"),read(smooth/"body.json"))
    with np.load(wired/"nodes.npz") as w, np.load(smooth/"nodes.npz") as s:
        assert np.array_equal(w["area"][common[:,0]],s["area"][common[:,1]])
        a=float(s["area"][common[:,1]].sum())
    inputs={k:sha(smooth/k) for k in ("mesh.su2","nodes.npz","body.json","initial_flow.csv","case.cfg",
                                      "case.json","engine_request.json","geometry.json","debugger_process.json")}
    prep=read(smooth/"debugger_process.json")
    assert prep["status"]=="terminal" and prep["debuggerExitCode"]==0 and all(prep["inputsStable"].values())
    data=dict(path=smooth.relative_to(ROOT).as_posix(),inputs=inputs,geometry=read(smooth/"geometry.json"),
              outerCellsIdentical=len(wo),newPatchCellsWithAllVerticesBeyondRadius=len(so-wo),
              identicalCellsTotal=len(wo)+len(wi&si),
              wiredCells=len(wi)+len(wo),smoothCells=len(si)+len(so),exactCommonFaces=len(common),commonAreaM2=a,
              optionsByteIdentical=True,threads=1,seconds=720)
    destination=HERE/"evidence";destination.mkdir(exist_ok=True)
    with (destination/"inputs.json").open("x",encoding="utf8") as f:json.dump(data,f,indent=2)


if __name__=="__main__":run()
