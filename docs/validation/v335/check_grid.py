"""从两份磁盘体网独立核对共同区、外接单元、面积与体积。"""
import argparse
import json
from pathlib import Path
import numpy as np
from outer_grid import BASE,HERE,PREP,PREP_KEY,FRAME,read,sha,load_mesh


def enclosed_volume(xyz,faces):
    value=0.
    for size in (3,4):
        p=xyz[np.array([f for f in faces if len(f)==size],int)]
        for i in range(1,size-1):
            value+=np.einsum("ij,ij->i",p[:,0],np.cross(p[:,i],p[:,i+1])).sum()/6
    return float(value)


def check(new):
    assert sha(PREP/"prepared.json")==PREP_KEY
    key=read(new/"prepared.json");s=key["geometry"]
    a,x,m=load_mesh(PREP/"mesh.su2",("wall","far"));b,y,n=load_mesh(new/"mesh.su2",("wall","far"))
    assert a==b[:len(a)] and m["wall"]==n["wall"]
    np.testing.assert_array_equal(x,y[:len(x)])
    for k in ("profile.dat","release.dat","transfer.npz"):assert sha(PREP/k)==sha(new/k)==key["files"][k]
    with np.load(PREP/"mapping.npz") as orig,np.load(new/"mapping.npz") as other:
        for k in orig.files:np.testing.assert_array_equal(orig[k],other[k][:len(x)] if k=="xyz" else other[k])
    previous=set(m["far"]);first=len(x);pitch=[];layers=[]
    for start in range(len(a),len(b),len(previous)):
        cells=b[start:start+len(previous)]
        lower={c[:len(c)//2] for c in cells};upper={c[len(c)//2:] for c in cells}
        assert lower==previous and len(upper)==len(cells);points=np.unique([i for f in upper for i in f])
        np.testing.assert_array_equal(points,np.arange(first,first+len(points)));first+=len(points)
        for size in (6,8):
            v=y[np.asarray([c for c in cells if len(c)==size])]
            segments=v[:,size//2:]-v[:,:size//2];pitch.extend([np.linalg.norm(segments,axis=2).min(),np.linalg.norm(segments,axis=2).max()])
        layers.append(float(np.linalg.norm(y[points],axis=1).mean()))
        previous=upper
    assert previous==set(n["far"]);np.testing.assert_allclose(layers,[16,18,20],atol=1e-11,rtol=0)
    np.testing.assert_allclose(pitch,2.,atol=1e-11,rtol=0)
    delta_volume=enclosed_volume(y,n["far"])-enclosed_volume(x,m["far"])
    old=read(PREP/"prepared.json")["geometry"];np.testing.assert_allclose(s["volumeD3"]-old["volumeD3"],delta_volume,rtol=1e-12)
    seed=np.loadtxt(new/"initial.csv",skiprows=1,delimiter=",")
    orig_frame=BASE/"run01/steady05/restart_000199.csv";assert sha(orig_frame)==FRAME
    np.testing.assert_array_equal(seed[:len(x)],np.loadtxt(orig_frame,skiprows=1,delimiter=",",usecols=range(10)))
    return dict(commonCells=len(a),commonNodes=len(x),newCells=len(b)-len(a),newNodes=len(y)-len(x),wetFaces=len(m["wall"]),
                addedShellVolumeD3=delta_volume,shellRadialPitchD=[float(min(pitch)),float(max(pitch))],radiiD=layers,
                outerRadiusM=20*key["condition"]["diameterM"],minCellVolumeD3=s["minCellVolumeD3"],minFaceConeD3=s["minConeVolumeD3"],
                binding=sha(new/"prepared.json"),mesh=sha(new/"mesh.su2"),initialFrame=sha(new/"initial.csv"),baselineFrame=FRAME)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("prepared",type=Path);p.add_argument("receipt",type=Path);a=p.parse_args()
    proof=check(a.prepared.resolve())
    with a.receipt.open("x",encoding="utf8") as f:json.dump(proof,f,indent=2);f.write("\n")
