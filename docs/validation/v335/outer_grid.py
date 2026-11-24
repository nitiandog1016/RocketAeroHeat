"""从实际 14D 全周网向外接三层。共同区保留原始浮点坐标、连通与五区壁载荷身份。"""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys

for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT/"scripts"), str(HERE.parent/"v329")]
from shortSupportMesh import audit_volume
from thermistorOuterDomain import load_mesh, write_mesh
from thermistorMeasuredCold import HEADER
from shortSphere65 import CASE, read, sha, write_json

BASE = ROOT/"outputs/v332_validation"
PREP = BASE/"prepared01"
PREP_KEY = "22ef3bf27abf1a093fae6bfb03145fac1b9f401229f11a3d12042be19cd9a04c"
FRAME = "c23c9ca1beb024105c600642f3fb016c9b7077599c8cecff227eb6d3ab0dbe63"
RADII = (16.,18.,20.)


def prepare(dest):
    dest.mkdir(parents=True,exist_ok=False)
    assert sha(PREP/"prepared.json")==PREP_KEY
    spec=read(PREP/"prepared.json")
    assert sha(CASE)==spec["caseSHA256"] and all(sha(PREP/k)==v for k,v in spec["files"].items())
    conv=read(BASE/"run01/converged.json");stage=BASE/"run01"/conv["stage"]
    frame=stage/"restart_000199.csv";rec=read(stage/"debugger_process.json");heat=read(BASE/"run01/thermal.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not(rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
    assert conv["frame"]==sha(frame)==FRAME==heat["early"]["sourceFrame"] and heat["passed"]
    assert heat["meshSHA256"]==spec["files"]["mesh.su2"]
    cells,xyz,marks=load_mesh(PREP/"mesh.su2",("far","wall"))
    original_cells=len(cells);n=len(xyz);far=marks["far"]
    nodes=np.unique([i for face in far for i in face])
    np.testing.assert_allclose(np.linalg.norm(xyz[nodes],axis=1),14,rtol=0,atol=1e-11)
    sets={w:np.array([f for f in far if len(f)==w],dtype=int) for w in (3,4)}
    edge_pairs={tuple(sorted((a,f[(j+1)%len(f)]))) for f in far for j,a in enumerate(f)}
    assert len(edge_pairs)==len(nodes)+len(far)-2
    prefix=xyz;blocks=[xyz];lower=nodes
    for k,r in enumerate(RADII):
        upper=np.arange(len(nodes))+n+k*len(nodes)
        blocks.append(prefix[nodes]*(r/14))
        for faces in sets.values():
            local=np.searchsorted(nodes,faces)
            cells.extend(map(tuple,np.column_stack((lower[local],upper[local])).tolist()))
        lower=upper
    outer=[tuple(lower[np.searchsorted(nodes,f)]) for f in far]
    xyz=np.vstack(blocks)
    write_mesh(dest/"mesh.su2",cells,xyz,dict(wall=marks["wall"],far=outer))
    # 独立重读文件；所有共享区节点、单元和湿壁必须逐值一致，新最外层逐片对应旧外边界。
    reread,vertices,tags=load_mesh(dest/"mesh.su2",("wall","far"))
    assert reread==cells and tags["wall"]==marks["wall"] and tags["far"]==outer
    np.testing.assert_array_equal(vertices,xyz)
    assert len(vertices)==269778 and len(cells)==406272
    del reread,cells,vertices
    stats=audit_volume(dest,RADII[-1])
    with np.load(PREP/"mapping.npz") as m:
        np.testing.assert_array_equal(m["xyz"],prefix)
        fields={k:m[k] for k in m.files if k!="xyz"}
    np.savez_compressed(dest/"mapping.npz",xyz=xyz,**fields)
    for file in ("profile.dat","release.dat","transfer.npz"):
        shutil.copyfile(PREP/file,dest/file);assert sha(dest/file)==spec["files"][file]
    initial=np.loadtxt(frame,skiprows=1,delimiter=",",usecols=range(10));np.testing.assert_array_equal(initial[:,1:4],prefix)
    np.testing.assert_array_equal(initial[:,0],np.arange(n))
    shell=np.zeros((len(xyz)-n,10));shell[:,0]=np.arange(n,len(xyz));shell[:,1:4]=xyz[n:];shell[:,5]=-1; shell[:,8]=1-298.15
    np.savetxt(dest/"initial.csv",np.vstack((initial,shell)),delimiter=",",header=HEADER,comments="",fmt=["%d"]+["%.17g"]*9)
    saved=np.loadtxt(dest/"initial.csv",delimiter=",",skiprows=1)
    np.testing.assert_array_equal(saved[:n],initial);np.testing.assert_array_equal(saved[n:],shell)
    geom={**spec["geometry"],**stats,"outerRadiusD":RADII[-1]}
    filekeys={k:sha(dest/k) for k in spec["files"]}
    write_json(dest/"prepared.json",{**spec,"geometry":geom,"files":filekeys})
    result=dict(basePrepared=PREP_KEY,prepared=sha(dest/"prepared.json"),baseMesh=spec["files"]["mesh.su2"],mesh=filekeys["mesh.su2"],
                initial=filekeys["initial.csv"],guessFrame=FRAME,baseCells=original_cells,baseNodes=n,baseVolumeD3=spec["geometry"]["volumeD3"],
                shellRadiiD=RADII,interfaceFaces=len(far),interfaceNodes=len(nodes),interfaceEdges=len(edge_pairs),
                addedEdges=len(RADII)*(len(edge_pairs)+len(nodes)),actualShellRadialPitchM=2*spec["condition"]["diameterM"],
                geometry=geom,condition=spec["condition"])
    write_json(dest/"extension.json",result)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path);a=parser.parse_args()
    dest=a.output.resolve();dest.relative_to(ROOT/"outputs")
    print(json.dumps(prepare(dest),ensure_ascii=False))
