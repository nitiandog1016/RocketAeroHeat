"""从固定风洞试件的已审核母线和扇区构建无对称边界的全周对照。"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from cm011TunnelSector import mesh, HEX_FACES, PRISM_FACES
from runEmbeddedSU2 import sha
from simpleLauncherTunnel import CASE, geometry


def surface_area(points):
    vector = np.cross(points[:, 1]-points[:, 0], points[:, -1]-points[:, 0])
    if points.shape[1] == 4:
        vector += np.cross(points[:, 2]-points[:, 1], points[:, 3]-points[:, 1])
    return np.linalg.norm(vector, axis=1)/2


def prepare(out, flow_run=None):
    definition, tri, diameter, wall = geometry()
    if flow_run is None:
        saved = json.loads((ROOT/"docs/validation/v350/evidence.json").read_bytes())["fresh"]["runDirectory"]
        flow_run = ROOT/"outputs"/Path(saved.replace("\\","/")).name
    parent = flow_run/"coarse"
    frame = parent/"restart.csv"
    record = json.loads((parent/"engine_process.json").read_bytes())
    assert record["status"] == "terminal" and record["exitCode"] == 0 and all(record["inputsStable"].values())
    assert sha(frame) == record["outputSHA256"]["restart.csv"] == "720e4ce803e78c984b3015b390b1f95bc95d70dde1676216e161ec07692cb8d6"
    condition = json.loads((CASE/"tunnel.json").read_bytes())
    assert condition["velocityMps"] == [-.05,0,0] and condition["surfaceSHA256"] == definition["fileSHA256"]
    out.relative_to(ROOT/"outputs")
    out.mkdir(parents=True,exist_ok=False)
    sectors = 48
    xy, sector_xyz, _, _, rows, cells, markers, single = mesh(
        np.array([0.,2*np.pi/sectors]), 1,36,wall,14,diameter)
    state = np.loadtxt(frame,delimiter=",",skiprows=1)
    np.testing.assert_allclose(xy,state[:,1:3],rtol=0,atol=1e-13)
    np.testing.assert_array_equal(state[:,0],np.arange(len(xy)))
    axis = rows[0] == rows[1]
    assert np.array_equal(axis,xy[:,1] == 0)
    off = np.flatnonzero(~axis)
    ids = np.empty((sectors,len(xy)),dtype=int)
    ids[:,axis] = np.arange(axis.sum())
    ids[:,off] = axis.sum()+np.arange(sectors*len(off)).reshape(sectors,-1)
    count = int(ids.max()+1)
    xyz = np.empty((count,3)); vel = np.empty_like(xyz)
    p = np.empty(count); h = np.empty(count)
    maps = np.empty((sectors,len(sector_xyz)),dtype=int)
    for j in range(sectors):
        c,s = np.cos(2*np.pi*j/sectors),np.sin(2*np.pi*j/sectors)
        idx = ids[j]
        xyz[idx] = np.column_stack((xy[:,0],c*xy[:,1],s*xy[:,1]))
        vel[idx] = np.column_stack((state[:,4],c*state[:,5],s*state[:,5]))
        p[idx],h[idx] = state[:,3],state[:,6]
        maps[j,rows[0]] = idx
        maps[j,rows[1]] = ids[(j+1)%sectors]
    # 验证全部接缝使用同一个节点；绕轴刚性旋转保留源单元绕向和体积。
    for j, index in enumerate(maps):
        c,s = np.cos(2*np.pi*j/sectors),np.sin(2*np.pi*j/sectors)
        rotation = np.array([[1,0,0],[0,c,-s],[0,s,c]])
        np.testing.assert_allclose(xyz[index],sector_xyz@rotation.T,atol=4e-14,rtol=0)
    def tiles(faces, width):
        part = np.array([f for f in faces if len(f)==width],int)
        return maps[:,part].reshape(-1,width)
    volumes = {13:tiles(cells,6),12:tiles(cells,8)}
    bounds = {name:{5:tiles(markers[name],3),9:tiles(markers[name],4)} for name in ("wall","far")}
    # 直接枚举整网的各个面。内面必须有两个相邻单元，剩余面必须由壁面及远场完整覆盖。
    for width, tag in [(3,5),(4,9)]:
        faces = np.vstack([block[:,list(f)] for kind,block in volumes.items()
            for f in (PRISM_FACES if kind==13 else HEX_FACES) if len(f)==width])
        unique,counts = np.unique(np.sort(faces,axis=1),axis=0,return_counts=True)
        assert np.all(counts<=2)
        expected = np.vstack([group[tag] for group in bounds.values()])
        np.testing.assert_array_equal(unique[counts==1],np.unique(np.sort(expected,axis=1),axis=0))
        assert (counts==1).sum() == len(expected)
    # 按同一轴向带和周向扇片对照原 STL。底面径向细分的面积先加回对应扇片。
    def bands(polys):
        x = polys[:,:,0]
        mid = polys.mean(axis=1)
        sector = np.floor(np.mod(np.arctan2(mid[:,2],mid[:,1]),2*np.pi)*sectors/(2*np.pi))
        return np.column_stack((np.rint(x.min(axis=1)*1e7),np.rint(x.max(axis=1)*1e7),sector)).astype(int)
    ref,ri = np.unique(bands(tri),axis=0,return_inverse=True)
    parts = [xyz[idx]*diameter for idx in bounds["wall"].values()]
    obs,oi = np.unique(np.vstack([bands(b) for b in parts]),axis=0,return_inverse=True)
    np.testing.assert_array_equal(ref,obs)
    ar = np.bincount(ri,weights=surface_area(tri))
    av = np.bincount(oi,weights=np.concatenate([surface_area(b) for b in parts]))
    np.testing.assert_allclose(av,ar,rtol=1e-5,atol=1e-12)
    with (out/"mesh.su2").open("x",encoding="ascii") as stream:
        def put(kind,block):
            np.savetxt(stream,block,fmt=str(kind)+" "+" ".join(["%d"]*block.shape[1]))
        stream.write(f"NDIME= 3\nNELEM= {sum(map(len,volumes.values()))}\n")
        for kind,block in volumes.items(): put(kind,block)
        stream.write(f"NPOIN= {count}\n")
        np.savetxt(stream,xyz,fmt="%.17g")
        stream.write("NMARK= 2\n")
        for name,group in bounds.items():
            stream.write(f"MARKER_TAG= {name}\nMARKER_ELEMS= {sum(map(len,group.values()))}\n")
            for kind,block in group.items(): put(kind,block)
    # 局部矢势 A_y 的旋度。扰动在壁面、远界及其邻域严格为零。
    x,y,z = xyz.T; dx=x+.45; a=.35; b=.25
    weight = np.maximum(1-(dx/a)**2-(y*y+z*z)/b**2,0.)**2
    change = np.column_stack((6*z/b**2*weight,np.zeros(count),-6*dx/a**2*weight))
    change *= .01/np.linalg.norm(change,axis=1).max()
    assert all(not change[idx].any() for group in bounds.values() for idx in group.values())
    header='"PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy"'
    for name,v in [("control",vel),("disturbed",vel+change)]:
        np.savetxt(out/(name+".csv"),np.column_stack((np.arange(count),xyz,p,v,h)),fmt=["%d"]+["%.17g"]*8,
                   delimiter=",",header=header,comments="")
    base_ids = ids[:,:len(wall)][:,wall[:,0]==0]
    np.savez_compressed(out/"stations.npz",base=base_ids,wall=ids[:,:len(wall)],xyz=xyz)
    paths = [out/name for name in ("mesh.su2","control.csv","disturbed.csv","stations.npz")]
    result=dict(nodes=count,cells=sum(map(len,volumes.values())),markers=list(bounds),sectors=sectors,
                wallFaces=sum(map(len,bounds["wall"].values())),originalFaces=len(tri),originalAreaM2=ar.sum(),
                wallAreaM2=av.sum(),checkedPatches=len(ar),maxPatchRelativeDifference=float(max(abs(av/ar-1))),
                minCellVolumeD3=single["minCellVolumeD3"],volumeD3=sectors*single["volumeD3"],
                diameterM=diameter,baseDiameterM=2*wall[wall[:,0]==0,1].max()*diameter,
                perturbedNodes=int(np.count_nonzero(np.linalg.norm(change,axis=1))),
                peakPerturbationOverU=float(np.linalg.norm(change,axis=1).max()),parentFrameSHA256=sha(frame),
                sources={str(q.relative_to(ROOT)):sha(q) for q in [CASE/"tunnel.json",CASE/definition["file"],
                    Path(__file__).resolve(),ROOT/"scripts/cm011Tunnel.py",ROOT/"scripts/cm011TunnelSector.py"]},
                files={q.name:sha(q) for q in paths})
    (out/"prepared.json").write_text(json.dumps(result,indent=2,allow_nan=False),encoding="utf8")
    print(count,"nodes;",len(ar),"original surface patches checked",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--flow-run",type=Path,help="从默认工作台冷启动复算产生的完整风洞作业目录")
    args=parser.parse_args()
    prepare(args.output.resolve(),args.flow_run)
