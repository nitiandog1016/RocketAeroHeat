"""从两套已算光滑体网的真实对称面边读取法向穿越间距，定位后续离散对照。"""
import json
from pathlib import Path

import numpy as np

from check_mesh import HERE,ROOT
from thermistorBead import A,B,D,sha


def boundary_grid(folder):
    with (folder/"mesh.su2").open("r") as f:
        assert next(f).strip()=="NDIME= 3"
        for _ in range(int(next(f).split("=")[1])):next(f)
        count=int(next(f).split("=")[1])
        xy=np.array([list(map(float,next(f).split()[:2])) for _ in range(count)])
        edges=set()
        for _ in range(int(next(f).split("=")[1])):
            name=next(f).split("=")[1].strip();n=int(next(f).split("=")[1])
            for _ in range(n):
                line=next(f)
                if name=="sym0":
                    kind,*pts=map(int,line.split())
                    assert len(pts)=={5:3,9:4}[kind]
                    edges.update(tuple(sorted((i,j))) for i,j in zip(pts,pts[1:]+pts[:1]))
    with np.load(folder/"nodes.npz") as z:
        ids=z["ids"]
    m=json.loads((folder/"case.json").read_bytes())["meridianIntervals"]+1
    return xy,np.array(sorted(edges)),xy[ids[0,:m]]


def crossings(folder):
    xy,edges,profile=boundary_grid(folder)
    a=xy[edges[:,0]];e=xy[edges[:,1]]-a
    angles=np.arctan2(profile[:,1]/B,-profile[:,0]/A)
    results=[]
    for degrees in (25,60,85,120,155):
        angle=np.pi-np.deg2rad(degrees)
        point=np.array([np.interp(angle,angles,profile[:,j]) for j in range(2)])
        normal=np.array([-np.cos(angle)/A,np.sin(angle)/B]);normal/=np.linalg.norm(normal)
        den=normal[0]*e[:,1]-normal[1]*e[:,0]
        use=abs(den)>1e-14
        ap=a[use]-point
        t=np.cross(ap,e[use])/den[use];u=np.cross(ap,normal)/den[use]
        hits=np.sort(t[(t>1e-9)&(u>=-1e-10)&(u<=1+1e-10)])
        hits=hits[np.r_[True,np.diff(hits)>1e-9]]
        assert hits[-1]>.25
        after=hits[hits>.015+1e-7][:3]
        line=np.r_[.015,hits[(hits>.015+1e-7)&(hits<.1-1e-7)],.1]
        results.append(dict(degFromUpstream=degrees,intervalsBetween015DAnd100D=len(line)-1,
                            maxCrossingIntervalOverD=float(np.diff(line).max()),nextThreeCrossingsOverD=after.tolist()))
    return results


def run():
    folder=ROOT/"outputs/v298_validation/comparison01"
    table=np.genfromtxt(str(folder/"historical.csv"),names=True,delimiter=",")
    old_folder=ROOT/"outputs/v283_validation/result01"
    old=np.genfromtxt(str(old_folder/"wall.csv"),names=True,delimiter=",")
    from thermistorBead import body
    geometry=body();v=np.asarray(geometry["V"]);tri=v[np.asarray(geometry["F"])-1]
    x=tri[:,:,0].mean(axis=1)/A
    area=old["areaM2"];delta=table["pairedHeatFluxWm2"]-old["heatFluxWm2"]
    bins=[]
    for lo,hi in zip((-1,-.5,0,.5),(-.5,0,.5,1)):
        keep=(x>=lo)&(x<hi)
        bins.append(dict(xOverA=[lo,hi],areaM2=float(area[keep].sum()),
                         extraOutwardW=float(-area[keep]@delta[keep]),faces=int(keep.sum())))
    s=json.loads((folder/"comparison.json").read_bytes())
    assert abs(sum(r["extraOutwardW"] for r in bins)-s["methodDeltaW"])<1e-13
    old_grid=ROOT/"outputs/v283_validation/strong01";new_grid=ROOT/"outputs/v298_validation/smooth01"
    report=dict(smoothPowerChangeByAxialBand=bins,oldStructured=crossings(old_grid),newMixed=crossings(new_grid),
                oldGridSHA256=sha(old_grid/"mesh.su2"),newGridSHA256=sha(new_grid/"mesh.su2"))
    with (HERE/"evidence/spatial.json").open("x",encoding="utf8") as f:json.dump(report,f,indent=2)


if __name__=="__main__":run()
