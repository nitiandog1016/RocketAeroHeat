"""独立重读全周带线热敏珠体网、真实原面与冷态，核对侧吹的方向和参考通量。"""
import argparse
from itertools import combinations
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import D, RHO, read, sha, write_json
from thermistorLeadMesh import SIDES
from thermistorOuterDomain import load_mesh


def face_keys(faces):
    """将有向多边形绕至最小节点，反向取得同键而保留相反符号。"""
    n=faces.shape[1];j=np.arange(n)[None,:]
    rot=np.take_along_axis(faces,(np.argmin(faces,axis=1)[:,None]+j)%n,axis=1)
    forward=rot[:,1]<rot[:,-1]
    keys=np.full((len(faces),4),-1,np.int32)
    keys[:,:n]=np.where(forward[:,None],rot,rot[:,(-j[0])%n])
    return keys,forward


def normals(p):
    v=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])/2
    if p.shape[1]==4:v+=np.cross(p[:,2]-p[:,0],p[:,3]-p[:,0])/2
    return v


def audit(folder):
    saved=read(folder/"inputs.json")
    assert all(sha(folder/k)==h for k,h in saved["files"].items())
    base=Path(saved["base"]);assert all(sha(base/k)==h for k,h in saved["source"].items())
    cells,xyz,marks=load_mesh(folder/"mesh.su2",("wall","far"))
    spec=read(folder/"case.json");body=read(folder/"body.json")
    assert spec["velocityMps"]==[0,-3,0] and spec["wallK"]-spec["airK"]==6
    assert len(cells)==saved["cells"] and len(xyz)==saved["nodes"]
    keys=[];signs=[];exterior=[];volume=0.;lowest=float("inf");cone_min=lowest
    def record(faces,physical=False):
        key,sign=face_keys(faces);keys.append(key)
        signs.append(~sign if physical else sign)
        exterior.append(np.full(len(key),physical))
    for width,sides in SIDES.items():
        group=np.array([c for c in cells if len(c)==width],np.int32)
        if not len(group):continue
        center=xyz[group].mean(axis=1);vol=np.zeros(len(group))
        for side in sides:
            faces=group[:,side];p=xyz[faces];v=normals(p)
            dv=(v*(p.mean(axis=1)-center)).sum(axis=1)/3
            assert np.min(dv)>0
            vol+=dv;cone_min=min(cone_min,float(np.min(dv)));record(faces)
        volume+=vol.sum();lowest=min(lowest,float(np.min(vol)))
    fluxes=[]
    for name,faces in marks.items():
        for width in (3,4):
            group=np.array([f for f in faces if len(f)==width],np.int32)
            if not len(group):continue
            record(group,True)
            if name=="far":
                p=xyz[group];v=normals(p)
                assert np.min((v*p.mean(axis=1)).sum(axis=1))>0
                np.testing.assert_allclose(np.linalg.norm(p,axis=2),20,rtol=0,atol=1e-11)
                fluxes.extend(v@np.asarray(spec["velocityMps"]))
    key=np.vstack(keys);sign=np.concatenate(signs);boundary=np.concatenate(exterior)
    order=np.lexsort(key.T[::-1])
    a,b=order[::2],order[1::2]
    assert np.array_equal(key[a],key[b]) and np.all(np.any(key[a[1:]]!=key[a[:-1]],axis=1))
    assert np.all(sign[a]!=sign[b]) and not np.any(boundary[a]&boundary[b])
    assert np.isclose(volume,12*read(base/"geometry.json")["volumeD3"],rtol=2e-13)

    # 直接按三角形所在的物理壁片归属核对面积和方向。首层沿其本地方位取样。
    with np.load(folder/"nodes.npz") as z:
        rays=z["ids"];wall=z["wallMap"];first=z["firstMap"];area=z["area"]
        assert np.array_equal(z["xyz"],xyz)
    assert len(rays)==49 and np.array_equal(rays[0],rays[-1])
    axis=rays[0]==rays[1];assert np.all(rays[:,axis]==rays[0,axis]) and axis.sum()==saved["axisNodes"]
    assert len(np.unique(rays[:-1,~axis]))==48*sum(~axis)
    np.testing.assert_allclose(xyz[wall]*D,np.asarray(body["V"])[np.asarray(body["F"])-1],rtol=0,atol=1e-16)
    with np.load(base/"nodes.npz") as parent:
        step=parent["xyz"][parent["firstMap"]]-parent["xyz"][parent["wallMap"]]
        np.testing.assert_allclose(np.linalg.norm(xyz[first]-xyz[wall],axis=2),np.linalg.norm(step,axis=2),rtol=0,atol=2e-15)
    width=spec["meridianIntervals"]+1
    assert np.array_equal(np.unique(rays[:,:width]),np.unique(wall))
    first_at=np.full(len(xyz),-1);first_at[rays[:,:width]]=rays[:,width:2*width]
    np.testing.assert_array_equal(first_at[wall],first)
    patches=marks["wall"];own={};normal=[]
    for i,f in enumerate(patches):
        normal.append(normals(xyz[np.array([f])])[0])
        for tri in combinations(f,3):
            k=tuple(sorted(tri));assert k not in own;own[k]=i
    normal=np.array(normal);owner=np.array([own[tuple(sorted(t))] for t in wall])
    tri_normal=normals(xyz[wall])*D**2
    assert np.all(np.sum(normal[owner]*tri_normal,axis=1)<0)
    np.testing.assert_allclose(np.linalg.norm(tri_normal,axis=1),area,rtol=1e-12,atol=0)
    np.testing.assert_allclose(np.bincount(owner,weights=area),np.linalg.norm(normal,axis=1)*D**2,rtol=1e-12)
    assert len(area)==10272 and abs(sum(area)-3.19665061256752e-6)<1e-17

    state=np.loadtxt(folder/"initial_flow.csv",delimiter=",",skiprows=1);cfg=(folder/"case.cfg").read_text("ascii")
    wet=np.isin(np.arange(len(xyz)),np.unique(wall))
    np.testing.assert_array_equal(state[:,0],np.arange(len(xyz)))
    np.testing.assert_array_equal(state[:,1:4],xyz)
    assert state.shape==(len(xyz),10) and np.all(state[:,[4,5,7]]==0)
    assert np.all(state[:,6]==np.where(wet,0,-1)) and np.all(state[:,9]==wet)
    assert np.all(state[:,8]==np.where(wet,2,1)-298.15)
    assert "INC_VELOCITY_INIT= (0, -1, 0)" in cfg and "MARKER_FAR= (far)" in cfg and "MARKER_SYM" not in cfg
    flux=np.array(fluxes)*D**2
    entering=flux<0;leaving=flux>0
    assert entering.sum()==leaving.sum() and np.isclose(-sum(flux[entering]),sum(flux[leaving]),rtol=1e-13)
    return dict(nodes=len(xyz),cells=len(cells),physicalFaces={k:len(v) for k,v in marks.items()},
                seamInternalFaces=saved["sectionFaces"]*12,sharedAxisNodes=int(axis.sum()),
                interiorFaces=len(a)-sum(map(len,marks.values())),volumeD3=float(volume),minCellVolumeD3=lowest,
                minPositiveFaceConeD3=cone_min,faces=len(area),wallAreaM2=float(area.sum()),
                groups=np.bincount(np.asarray(body["regionID"],int))[1:].tolist(),speedMps=spec["velocityMps"],
                inletReferenceFaces=int(entering.sum()),outletReferenceFaces=int(leaving.sum()),
                uniformInflowKgs=float(-RHO*flux[entering].sum()),uniformOutflowKgs=float(RHO*flux[leaving].sum()),
                files=saved["files"],baseFiles=saved["source"],coldStateChecked=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("folder",type=Path)
    args=parser.parse_args();folder=args.folder.resolve()
    path=folder/"mesh_audit.json";assert not path.exists();write_json(path,audit(folder))
