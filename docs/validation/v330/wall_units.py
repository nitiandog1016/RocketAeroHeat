"""区分层流 VTK 的 Y_Plus 字段和由实际切应力、首层法向距离估算的壁面单位。"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from check_grid import COARSE,HERE,ROOT,sha
from readCM011TunnelSector import vtk_tokens


def measure(stage,prepared):
    record=json.loads((stage/"record.json").read_bytes())
    frame=stage/"restart_000199.csv";vtk=stage/"surface_000199.vtk"
    assert record["end"]==sha(frame) and sha(vtk)==sha(stage/"surface.vtk")
    points,fields=vtk_tokens(vtk);flow=np.loadtxt(frame,delimiter=",",skiprows=1,usecols=(0,4))
    s=json.loads((prepared/"prepared.json").read_bytes())
    with np.load(prepared/"mapping.npz") as m:
        xyz,c,wet,first=(m[k] for k in ("xyz","corners","wet","first"))
    dist,order=cKDTree(points).query(xyz[wet]);assert len(set(order))==len(points)==len(wet) and max(dist)<1e-6
    np.testing.assert_allclose(fields["Pressure"][order].ravel(),flow[wet,1],atol=1e-7,rtol=7e-7)
    ix=np.searchsorted(wet,c);p=xyz[c]
    n=np.cross(p[:,2]-p[:,0],p[:,3]-p[:,1]);n/=np.linalg.norm(n,axis=1)[:,None]
    normalDistance=np.einsum("ijk,ik->ij",xyz[first[ix]]-p,n)
    cf=fields["Skin_Friction_Coefficient"][order[ix]];tan=cf-n[:,None,:]*np.einsum("ijk,ik->ij",cf,n)[:,:,None]
    y=s["condition"]["re"]*normalDistance*np.sqrt(np.linalg.norm(tan,axis=2)/2)
    assert np.isfinite(y).all() and normalDistance.min()>0
    return dict(frame=record["end"],vtkSHA256=sha(vtk),matchedNodes=len(wet),coordinateToleranceObservedD=dist.max(),
                laminarVTKYPlusRange=[float(fields["Y_Plus"].min()),float(fields["Y_Plus"].max())],
                estimatedFirstYPlusCornerRange=[y.min(),y.max()],minFirstM=normalDistance.min()*s["condition"]["diameterM"],
                maxFirstM=normalDistance.max()*s["condition"]["diameterM"])


if __name__=="__main__":
    b=ROOT/"outputs/v330_validation"
    report=dict(coarse=measure(ROOT/"outputs/v329_validation/run01/steady03",COARSE),fine=measure(b/"run01/steady03",b/"prepared01"))
    with (HERE/"evidence/wall_units.json").open("x",encoding="utf8",newline="\n") as f:
        json.dump(report,f,indent=2,allow_nan=False);f.write("\n")
    print(report["coarse"]["estimatedFirstYPlusCornerRange"],report["fine"]["estimatedFirstYPlusCornerRange"])
