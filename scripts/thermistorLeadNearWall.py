"""用原生切应力和首层流速定位反向流动，并在同一原面上核对首层壁面单位距离。"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from readCM011TunnelSector import vtk_tokens
from thermistorBead import D,RHO,U,MU,ROOT,read,sha
from thermistorLeadMesh import meridian,original,SECTORS


def assess(folder):
    rec=read(folder/"engine_process.json");vtk=folder/"surface.vtk"
    assert rec["outputSHA256"]["surface.vtk"]==sha(vtk)
    saved=folder/"restart.csv"
    assert rec["outputSHA256"]["restart.csv"]==sha(saved)
    surf,field=vtk_tokens(vtk)
    with np.load(folder/"nodes.npz") as z:
        xyz,wall,near=z["xyz"],z["wallMap"],z["firstMap"]
    all_nodes=np.unique(wall)
    dist,order=cKDTree(surf).query(xyz[all_nodes])
    assert len(np.unique(order))==len(surf)==len(all_nodes) and dist.max()<2e-6
    vtk_ids=order[np.searchsorted(all_nodes,wall)]
    cf=field["Skin_Friction_Coefficient"][vtk_ids]
    flow=np.loadtxt(saved,delimiter=",",skiprows=1)
    assert np.allclose(flow[:,1:4],xyz,rtol=0,atol=1e-13)
    smooth=read(folder/"case.json")["caseKind"]=="thermistor-B-paired-smooth"
    profile,_=meridian(smooth)
    body,links,area=original(smooth)
    assert body==read(folder/"body.json")
    strip=links[:,:,1].min(axis=1);phi=links[:,:,0]*2*np.pi/SECTORS
    direction=-np.diff(profile,axis=0)
    direction/=np.linalg.norm(direction,axis=1)[:,None]
    tx,tr=direction[strip].T
    tangent=np.stack(np.broadcast_arrays(tx[:,None],tr[:,None]*np.cos(phi),tr[:,None]*np.sin(phi)),axis=-1)
    shear=np.sum(cf*tangent,axis=2).mean(axis=1)
    speed=np.sum(flow[near,5:8]*tangent,axis=2).mean(axis=1)
    # 用原生 Cf 对应的摩擦速度与首层节点到本面弦平面的投影距离估算 y+。
    pts=xyz[wall];normal=np.cross(pts[:,1]-pts[:,0],pts[:,2]-pts[:,0])
    normal/=np.linalg.norm(normal,axis=1)[:,None]
    y=D*np.abs(np.einsum("fij,fj->fi",xyz[near]-pts,normal)).mean(axis=1)
    wall_units=y*RHO*U/MU*np.sqrt(np.linalg.norm(cf,axis=2).mean(axis=1)/2)
    assert np.isfinite(wall_units).all() and np.all(y>0)
    reverse=(shear<0)&(speed<0)
    raw=np.genfromtxt(str(folder/"surface/wall.csv"),names=True,delimiter=",")
    assert np.array_equal(raw["areaM2"],area) and np.array_equal(raw["regionID"],body["regionID"])
    tri_x=(np.array(body["V"])[np.array(body["F"])-1])[:,:,0]
    strips=np.unique(strip[reverse])
    summary=dict(sourceFrameSHA256=sha(saved),surfaceVTKSHA256=sha(vtk),
                 nativeYPlusRange=[float(field["Y_Plus"].min()),float(field["Y_Plus"].max())],
                 estimatedFirstLayerYPlusRange=[float(wall_units.min()),float(wall_units.max())],
                 signDisagreementFaces=int(((shear<0)!=(speed<0)).sum()),reverseFaces=int(reverse.sum()),
                 reverseAreaM2=float(area[reverse].sum()),reverseReleaseW=float(-np.sum(raw["heatFluxWm2"][reverse]*area[reverse])),
                 reverseXRangeM=[float(tri_x[reverse].min()),float(tri_x[reverse].max())] if reverse.any() else [],
                 reverseStripsZeroBased=strips.tolist(),reverseRegionNames=[body["regionNames"][i-1] for i in np.unique(raw["regionID"][reverse].astype(int))])
    return summary,np.column_stack((raw["faceID"]+1,strip,shear,speed,wall_units,reverse))


COLUMNS="face_id,meridian_strip_0,downstream_Cf,downstream_first_U_over_Uinf,estimated_first_yplus,reverse_both"


def write_faces(target,rows):
    np.savetxt(target,rows,delimiter=",",fmt=["%d","%d","%.17g","%.17g","%.17g","%d"],comments="",header=COLUMNS)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("folder",type=Path)
    cli.add_argument("--out",type=Path)
    a=cli.parse_args();folder=a.folder.resolve();folder.relative_to(ROOT/"outputs")
    dest=a.out.resolve() if a.out else folder/"surface"
    dest.relative_to(ROOT/"outputs")
    summary,rows=assess(folder)
    target=dest/"flow.csv"
    with target.open("x",encoding="ascii",newline="") as out:write_faces(out,rows)
    summary["tableSHA256"]=sha(target)
    with (dest/"flow.json").open("x",encoding="utf8") as out:json.dump(summary,out,indent=2,ensure_ascii=False)
