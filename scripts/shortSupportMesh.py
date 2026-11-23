"""为 Short 低速双支管银球生成全周近壁体网；计算轴为支管轴，文件坐标为项目 XYZ。"""
import argparse
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np

from shortSphere65 import CASE, ROOT, condition, kelvin, read, sha, write_json
from thermistorLeadMesh import SIDES, outer_mesh, rotate, sweep, wall_layers
from thermistorMeasuredCold import HEADER
from thermistorOuterDomain import load_mesh, write_mesh

BODY = ROOT/"cases/short_sphere/supported65.zip"
BODY_HASH = "d1f30c067083efb7731d1a9de74eed60064dc590b0277be2e7d59183a7d06778"
SECTORS = 64
PHASE = np.pi/SECTORS
COUNTS = (3,20,48,20,3)
REGIONS = (5,3,1,2,4)


def reference():
    assert sha(BODY)==BODY_HASH
    with ZipFile(BODY) as z:
        data=json.loads(z.read("geometry.json"))
        info=json.loads(z.read("summary.json"))
    assert info["caseSHA256"]==sha(CASE) and info["extensionD"]==1 and info["polygonSides"]==SECTORS
    return data,info


def section(radius, end, structured, split_near=False):
    angle=np.arcsin(2*radius)
    t=np.linspace(np.pi-angle,angle,COUNTS[2]+1)
    sphere=.5*np.c_[np.cos(t),np.sin(t)]
    joint=sphere[-1,0]
    pipe=np.c_[np.linspace(joint,end,COUNTS[3]+1),np.full(COUNTS[3]+1,radius)]
    cap=np.c_[np.full(COUNTS[4]+1,end),np.linspace(radius,0,COUNTS[4]+1)]
    mirror=np.array([-1,1])
    wall=np.vstack([cap[:0:-1]*mirror,pipe[:0:-1]*mirror,sphere[:-1],pipe[:-1],cap])
    layers=wall_layers(wall,.002)
    if split_near:
        # 对每个既有间隔插入中面，原五层及混合外区接口的位置由原数组直接保留。
        nested=np.empty((2*len(layers)-1,*layers.shape[1:]))
        nested[::2]=layers
        nested[1::2]=(layers[:-1]+layers[1:])/2
        layers=nested
    if structured:
        # 球身中段沿近径向伸展，端盖和直管各留连续的外圈角度，避免极区挤在同一方向。
        alpha=np.linspace(np.pi,0,len(wall))
        nominal=np.arctan2(wall[:,1],wall[:,0])
        blend=np.minimum(1.,abs(wall[:,0])/.5)**3
        beta=nominal+blend*(alpha-nominal)
        assert np.all(np.diff(beta)<0)
        outer=14*np.c_[np.cos(beta),np.sin(beta)];outer[[0,-1],1]=0
        intervals=30
        from scipy.optimize import brentq
        ratio=brentq(lambda a:.012*sum(a**np.arange(intervals))-14,1,2)
        t=np.cumsum(ratio**np.arange(intervals));t/=t[-1]
        far_layers=layers[-1]+t[:,None,None]*(outer-layers[-1])
        lattice=np.vstack((layers,far_layers))
        ids=np.arange(lattice.shape[0]*len(wall)).reshape(-1,len(wall))
        cells=np.stack((ids[:-1,:-1],ids[:-1,1:],ids[1:,1:],ids[1:,:-1]),axis=-1).reshape(-1,4)
        far=np.c_[ids[-1,1:],ids[-1,:-1]]
        return wall,lattice.reshape(-1,2),cells.tolist(),far
    background,tri,far,interface,active=outer_mesh(layers[-1],14.)
    extra=np.setdiff1d(active,interface)
    idx=np.full(len(background),-1,dtype=int)
    idx[interface]=np.arange(len(wall))+layers[1:].size//2
    idx[extra]=layers.size//2+np.arange(len(extra))
    np.testing.assert_allclose(background[interface],layers[-1],rtol=0,atol=1e-12)
    grid=np.arange(layers.size//2).reshape(layers.shape[:2])
    quads=np.stack((grid[:-1,:-1],grid[:-1,1:],grid[1:,1:],grid[1:,:-1]),axis=-1).reshape(-1,4)
    xy=np.vstack((layers.reshape(-1,2),background[extra]))
    return wall,xy,[*quads.tolist(),*idx[tri].tolist()],idx[far]


def audit_volume(folder):
    """独立重读 SU2 文件；无须依赖生成器的单元序号判定内外邻接。"""
    sys.path.insert(0,str(ROOT/"docs/validation/v311"))
    from audit_mesh import face_keys, normals
    cells,xyz,marks=load_mesh(folder/"mesh.su2",("far","wall"))
    keys=[];signs=[];volume=0.;lowest=np.inf;cone=np.inf
    for width,sides in SIDES.items():
        group=np.array([c for c in cells if len(c)==width],dtype=np.int32)
        if not len(group):continue
        center=xyz[group].mean(axis=1)
        v=np.zeros(len(group))
        for ids in sides:
            faces=group[:,ids];p=xyz[faces]
            dv=np.einsum("ij,ij->i",normals(p),p.mean(1)-center)/3
            assert dv.min()>0,(width,ids,dv.min())
            cone=min(cone,dv.min());v+=dv
            k,s=face_keys(faces);keys.append(k);signs.append(s)
        volume+=v.sum();lowest=min(lowest,v.min())
    wall_and_far=[]
    for name,faces in marks.items():
        for width in (3,4):
            group=np.array([c for c in faces if len(c)==width],dtype=np.int32)
            if not len(group):continue
            k,s=face_keys(group);keys.append(k);signs.append(~s);wall_and_far.append(k)
            if name=="far":
                np.testing.assert_allclose(np.linalg.norm(xyz[group],axis=2),14,rtol=0,atol=1e-10)
                assert np.min(np.einsum("ij,ij->i",normals(xyz[group]),xyz[group].mean(1)))>0
    exterior=np.vstack(wall_and_far)
    assert len(np.unique(exterior,axis=0))==len(exterior)
    k=np.vstack(keys);s=np.concatenate(signs);order=np.lexsort(k.T[::-1]);a,b=order[::2],order[1::2]
    assert np.array_equal(k[a],k[b]) and np.all(s[a]!=s[b]) and np.all((k[a[1:]]!=k[a[:-1]]).any(axis=1))
    return dict(nodes=len(xyz),cells=len(cells),wallFaces=len(marks["wall"]),farFaces=len(marks["far"]),
                minCellVolumeD3=lowest,minConeVolumeD3=cone,volumeD3=volume)


def prepare(out,structured=False,split_near=False,split_sectors=False):
    data,info=reference();case=read(CASE);s=condition(case);D=s["diameterM"]
    out.mkdir(parents=True,exist_ok=False)
    wall,xy,cells2,far=section(info["tubeDiameterM"]/(2*D),info["endYM"]/D,structured,split_near)
    if split_sectors and (structured or not split_near):
        raise ValueError("这项周向对照以十层混合网为基线")
    sectors=SECTORS*(2 if split_sectors else 1)
    rays=PHASE+np.arange(sectors+1)*2*np.pi/sectors
    rotated,ids=rotate(xy,rays,True)
    if split_sectors:
        # 湿壁与全部近壁层使用旧弦的几何中点。14D 远场新增射线仍落在原半径上。
        radial=np.linalg.norm(xy,axis=1)
        blend=np.clip((radial-2.)/12.,0.,1.)[:,None]
        for k in range(1,sectors,2):
            mid=(rotated[ids[k-1]]+rotated[ids[k+1]])/2
            rotated[ids[k]]=mid+blend*(rotated[ids[k]]-mid)
    xyz=rotated[:,[1,0,2]].copy();xyz[:,2]*=-1
    cells=[sweep(c,lo,hi) for lo,hi in zip(ids,ids[1:]) for c in cells2]
    width=len(wall)
    def surface(edges):
        a,b=edges.T
        return np.stack([ids[:-1,a],ids[1:,a],ids[1:,b],ids[:-1,b]],axis=-1).reshape(-1,4)
    corners=surface(np.c_[np.arange(width-1),np.arange(1,width)])
    def native_faces(array,invert=False):
        # 在旋转轴合并两重顶点；流体侧的内边界法向指向固体。
        result=[]
        for row in array.tolist():
            face=list(dict.fromkeys(row))
            result.append(face[::-1] if invert else face)
        return result
    # Gmsh 外圆弧由正轴经赤道到负轴，其母线方向与物体母线相反。
    far_faces=native_faces(surface(far),True)
    boundary=dict(wall=native_faces(corners,True),far=far_faces)
    write_mesh(out/"mesh.su2",cells,xyz,boundary)
    del cells
    wet=np.unique(corners)
    lookup=np.full(len(xyz),-1,int);lookup[ids[:,:width]]=ids[:,width:2*width]
    first=lookup[wet];assert np.all(first>=0)
    # 检查每个壁片角点的实际法向首层距离，含弯头与端盖。
    p=xyz[corners];vec=np.cross(p[:,2]-p[:,0],p[:,3]-p[:,1])/2
    area=np.linalg.norm(vec,axis=1);normals=vec/area[:,None]
    heights=np.einsum("ijk,ik->ij",xyz[lookup[corners]]-p,normals)
    divisor=2 if split_near else 1
    assert heights.min()>.00195/divisor and heights.max()<.00205/divisor
    group=np.tile(np.repeat(REGIONS,COUNTS),sectors)
    # 接头共享球面温界，管壁其余排与端盖给定相同赤道温度。
    tags=ids[:,:width];sphere_nodes=np.unique(tags[:,23:72])
    eq=np.interp(90,case["angleDeg"],kelvin(case["surfaceF"]))
    wall_t=np.full(len(wet),eq)
    angles=np.degrees(np.arctan2(np.linalg.norm(xyz[sphere_nodes,1:],axis=1),xyz[sphere_nodes,0]))
    wall_t[np.searchsorted(wet,sphere_nodes)]=np.interp(angles,case["angleDeg"],kelvin(case["surfaceF"]))
    theta=(wall_t-s["airK"])/s["scalarUnitK"]
    for name,v in (("profile",theta),("release",np.zeros(len(theta)))):
        np.savetxt(out/(name+".dat"),np.c_[wet,1+theta,v],fmt=["%d","%.17g","%.17g"])
    state=np.zeros((len(xyz),6));state[:,1]=-1;state[wet,1]=0;state[:,4]=1-298.15
    state[wet,4]+=theta;state[wet,5]=theta
    np.savetxt(out/"initial.csv",np.c_[np.arange(len(xyz)),xyz,state],delimiter=",",comments="",header=HEADER,fmt=["%d"]+["%.17g"]*9)
    np.savez_compressed(out/"mapping.npz",xyz=xyz,wet=wet,theta=theta,first=first,corners=corners,
                        nativeArea=area*D**2,group=group)
    seal_prepared(out,split_near,split_sectors)


def seal_prepared(out,split_near=False,split_sectors=False):
    data,info=reference();case=read(CASE);s=condition(case);D=s["diameterM"]
    stats=audit_volume(out)
    with np.load(out/"mapping.npz") as z:
        xyz=z["xyz"];c=z["corners"];wet=z["wet"];area=z["nativeArea"];group=z["group"]
        p=xyz[c];q=xyz[z["first"][np.searchsorted(wet,c)]]
    normal=np.cross(p[:,2]-p[:,0],p[:,3]-p[:,1])
    height=np.einsum("ijk,ik->ij",q-p,normal/np.linalg.norm(normal,axis=1)[:,None])
    divisor=2 if split_near else 1
    assert height.min()>.00195/divisor and height.max()<.00205/divisor
    stats.update(minFirstM=height.min()*D,maxFirstM=height.max()*D,outerRadiusD=14,
                 nativeAreasM2=np.bincount(group,area)[1:].tolist())
    if split_near:stats["nearWallIntervals"]=10
    if split_sectors:
        stats["angularSectors"]=2*SECTORS
        # 球形远场两个轴端共用节点，其余节点均来自真实远场标记。
        _,_,marks=load_mesh(out/"mesh.su2",("far","wall"))
        stats["farNodes"]=len({i for face in marks["far"] for i in face})
    from shortSupportTransfer import transfer
    stats["transfer"]=transfer(out,data,s,stats)
    files=["mesh.su2","initial.csv","profile.dat","release.dat","mapping.npz","transfer.npz"]
    write_json(out/"prepared.json",dict(geometry=stats,caseSHA256=sha(CASE),supportedSHA256=BODY_HASH,
               condition=s,equatorK=float(np.interp(90,case["angleDeg"],kelvin(case["surfaceF"]))),
               files={n:sha(out/n) for n in files}))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--outer",choices=("hybrid","structured"),default="hybrid")
    parser.add_argument("--split-near",action="store_true",help="将原五个近壁间隔全部对半")
    parser.add_argument("--split-sectors",action="store_true",help="在十层近壁网的每个旧扇区内增加中间射线")
    args=parser.parse_args()
    out=args.output.resolve();out.relative_to(ROOT/"outputs");prepare(out,args.outer=="structured",args.split_near,args.split_sectors)
