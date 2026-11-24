"""在 Short Test 65 的同一物理终帧上截取七个极角的湿壁周向分布。"""
import csv
import io
import json
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

import numpy as np
from scipy.spatial import cKDTree

from shortSphere65 import CASE, ROOT, conductivity, read, sha
from map_surface import corner_weights

BODY=ROOT/"cases/short_sphere/supported65.zip"
FIELDS=("angleDeg","phiDeg","kind","loDeg","hiDeg","status","faceID","faceCount","patchID",
        "xM","yM","zM","qInWm2","wallK","nu","nuLower","nuUpper")


def unit(v):
    return v/np.linalg.norm(v,axis=-1,keepdims=True)


def direction(angle, phi):
    a,b=np.deg2rad([angle,phi])
    return np.array([np.cos(a),np.sin(a)*np.sin(b),np.sin(a)*np.cos(b)])


def cut_angles(values):
    values=np.sort(np.mod(values,360.))
    values=values[(values>1e-8)&(values<360-1e-8)]
    return np.r_[0.,values[np.r_[True,np.diff(values)>1e-8]] if len(values) else [],360.]


class Wall:
    """对实际三角面作从试件中心出发的射线求交；共享边保留全部命中面。"""
    def __init__(self,xyz,tri,groups,parents):
        self.xyz=xyz;self.tri=tri;self.groups=groups;self.parents=parents
        vertices=xyz[tri]
        self.inverse=np.linalg.inv(vertices.transpose(0,2,1))
        dirs=unit(vertices);centers=unit(dirs.sum(axis=1))
        self.tree=cKDTree(centers)
        self.radius=np.linalg.norm(dirs-centers[:,None,:],axis=2).max()+1e-9
        assert self.radius<.4
        edges=np.unique(np.sort(np.vstack([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]]),axis=1),axis=0)
        self.u,self.v=unit(xyz[edges]).transpose(1,0,2)
        self.normals=unit(np.cross(self.u,self.v))

    def crossings(self,angle):
        """锥面与有限边相交的全部方位；共面边贡献端点，开区间内命中集合恒定。"""
        a=np.deg2rad(angle);sn=np.sin(a);cs=np.cos(a)
        n=self.normals
        A=sn*n[:,2];B=sn*n[:,1];C=cs*n[:,0]
        size=np.hypot(A,B)
        ids=np.flatnonzero((size>1e-13)&(abs(C)<=size+1e-13))
        base=np.arctan2(B[ids],A[ids])
        offset=np.arccos(np.clip(-C[ids]/size[ids],-1.,1.))
        phi=np.r_[base+offset,base-offset];idx=np.tile(ids,2)
        rays=np.column_stack((np.full(len(phi),cs),sn*np.sin(phi),sn*np.cos(phi)))
        on_arc=(np.einsum("ij,ij->i",np.cross(self.u[idx],rays),n[idx])>=-1e-12)&(
            np.einsum("ij,ij->i",np.cross(rays,self.v[idx]),n[idx])>=-1e-12)
        aligned=(size<=1e-13)&(abs(C)<=1e-13)
        ends=np.vstack((self.u[aligned],self.v[aligned]))
        return np.degrees(np.r_[phi[on_arc],np.arctan2(ends[:,1],ends[:,2])])

    def hit(self,ray):
        candidates=np.asarray(self.tree.query_ball_point(ray,self.radius),dtype=int)
        coeff=self.inverse[candidates]@ray
        total=coeff.sum(axis=1)
        valid=(total>0)&(coeff.min(axis=1)>=-2e-10*np.abs(total))
        ids=candidates[valid];den=total[valid]
        if not len(ids):return "missing",ids,None,None
        if np.ptp(1/den)>2e-8:return "multiple",ids,None,None
        tags=set(self.groups[ids])
        status="seam" if len(tags)>1 else ("sphere" if tags=={1} else ("tube" if next(iter(tags)) in (2,3) else "cap"))
        return status,ids,coeff[valid]/den[:,None],ray/den[0]


def trig_bounds(v, angle, lo, hi):
    """求 v·d 的端点及区间内驻点值，用于给线性分式插值作包围。"""
    a=np.deg2rad(angle);left,right=np.deg2rad([lo,hi])
    A=np.sin(a)*v[2];B=np.sin(a)*v[1];C=np.cos(a)*v[0]
    points=[left,right];peak=np.arctan2(B,A)
    for k in range(-1,4):
        t=peak+k*np.pi
        if left<t<right:points.append(t)
    vals=C+A*np.cos(points)+B*np.sin(points)
    return min(vals),max(vals)


def field_bounds(values,inverse,angle,lo,hi):
    # 节点凸包与弦面上的分式插值界取交；中心化降低近常温场的消去误差。
    base=values.mean()
    num=trig_bounds((values-base)@inverse,angle,lo,hi)
    den=trig_bounds(inverse.sum(axis=0),angle,lo,hi)
    assert den[0]>0
    quot=[x/y for x in num for y in den]
    pad=2e-10*max(1.,max(abs(values)))
    return max(min(values),base+min(quot))-pad,min(max(values),base+max(quot))+pad


def sample(native,body,ray,q,tw,condition,angle,phi,lo,hi,kind):
    ns,ids,w,point=native.hit(ray);rs,faces,_,_=body.hit(ray)
    owners=np.unique(body.parents[faces])
    status=ns if ns==rs else ("seam" if "seam" in (ns,rs) else "mesh-disagreement")
    row=dict(angleDeg=angle,phiDeg=phi,kind=kind,loDeg=lo,hiDeg=hi,status=status,
             faceID=int(min(owners)) if len(owners) else 0,faceCount=len(owners),
             patchID=int(min(native.parents[ids])) if len(ids) else 0)
    if status!="sphere":return row
    flux=np.sum(w*q[native.tri[ids]],axis=1);wall=np.sum(w*tw[native.tri[ids]],axis=1)
    assert np.ptp(flux)<1e-5 and np.ptp(wall)<1e-7
    D=condition["diameterM"];air=condition["airK"]
    def nu(f,T):return -f*D/((T-air)*conductivity(T))
    f=flux[0];T=wall[0];value=nu(f,T);lower=upper=value
    if kind=="interval":
        i=ids[0];tri=native.tri[i]
        qb=field_bounds(q[tri],native.inverse[i],angle,lo,hi)
        tb=field_bounds(tw[tri],native.inverse[i],angle,lo,hi)
        assert tb[0]>air
        bounds=[nu(x,y) for x in qb for y in tb]
        lower=min(bounds);upper=max(bounds)
        assert lower<=value<=upper
    row.update(zip(("xM","yM","zM"),point*D))
    row.update(qInWm2=f,wallK=T,nu=value,nuLower=lower,nuUpper=upper)
    return row


def ranges(rows,key):
    available=[r for r in rows if r["status"]=="sphere"]
    if not available:return None
    def info(r):return dict(value=r[key],phiDeg=r["phiDeg"],faceID=r["faceID"])
    return dict(minimum=info(min(available,key=lambda r:r[key])),maximum=info(max(available,key=lambda r:r[key])))


def compare_rear(front,back):
    valid=[];lower=[];upper=[];degrees=0.
    front=[r for r in front if r["kind"]!="inside-gap-boundary"]
    back=[r for r in back if r["kind"]!="inside-gap-boundary"]
    assert len(front)==len(back)
    for a,b in zip(front,back):
        assert (a["phiDeg"],a["kind"],a["hiDeg"])==(b["phiDeg"],b["kind"],b["hiDeg"])
        if a["status"]!= "sphere" or b["status"]!="sphere":continue
        delta=a["nu"]-b["nu"]
        valid.append(dict(phiDeg=a["phiDeg"],deltaNu=delta,face119=a["faceID"],face166=b["faceID"]))
        lower.append(a["nuLower"]-b["nuUpper"]);upper.append(a["nuUpper"]-b["nuLower"])
        degrees+=a["hiDeg"]-a["loDeg"]
    if not valid:return dict(status="no-common-direction",validDegrees=0)
    best=min(valid,key=lambda r:r["deltaNu"])
    status="observed-reversal" if best["deltaNu"]<0 else ("positive-throughout" if min(lower)>0 else "unresolved")
    return dict(status=status,validDegrees=degrees,minimumObserved=best,maximumObserved=max(valid,key=lambda r:r["deltaNu"]),
                lowerBound=min(lower),upperBound=max(upper))


def compute(prepared,work,out,*,branch="half"):
    assert branch in ("full","half")
    s=read(out/branch/"summary.json");step=s["energy"]["step"]
    stage=work/branch/f"step{step}";frame=stage/f"restart_{step:05d}.csv"
    endpoint=read(stage/"energy.json");p=read(prepared/"prepared.json")
    c=s["condition"];angles=read(CASE)["angleDeg"];D=c["diameterM"]
    assert s["branch"]==branch and p["condition"]==c and p["geometry"]["outerRadiusD"]==14
    assert s["energy"]==endpoint and sha(frame)==s["frameSHA256"]==endpoint["frame"]
    assert s["caseSHA256"]==sha(CASE) and s["bodySHA256"]==sha(BODY)
    assert s["meshSHA256"]==p["files"]["mesh.su2"] and sha(prepared/"mapping.npz")==p["files"]["mapping.npz"]
    assert sha(stage/"wall_budget.npz")==endpoint["files"]["wall_budget.npz"]
    with np.load(prepared/"mapping.npz") as m:
        xyz,wet,theta,corners,area,group=[m[k] for k in ("xyz","wet","theta","corners","nativeArea","group")]
    with np.load(stage/"wall_budget.npz") as b:
        np.testing.assert_array_equal(b["wet"],wet);power=b["nodeInW"]
    assert abs(-power.sum()-endpoint["wallW"])<1e-10
    weights=corner_weights(corners,xyz)*area[:,None]
    faces=np.searchsorted(wet,corners)
    q=power/np.bincount(faces.ravel(),weights=weights.ravel(),minlength=len(wet))
    tw=c["airK"]+c["scalarUnitK"]*theta
    tri=np.concatenate((faces[:,[0,1,2]],faces[:,[0,2,3]]))
    nonzero=(tri[:,0]!=tri[:,1])&(tri[:,0]!=tri[:,2])&(tri[:,1]!=tri[:,2])
    native=Wall(xyz[wet],tri[nonzero],np.tile(group,2)[nonzero],np.tile(np.arange(len(group))+1,2)[nonzero])
    with ZipFile(BODY) as z:g=json.loads(z.read("geometry.json"))
    body=Wall(np.asarray(g["V"])/D,np.array(g["F"])-1,np.array(g["regionID"]),np.arange(len(g["F"]))+1)
    grids={a:cut_angles(np.r_[np.arange(360),native.crossings(a),body.crossings(a)]) for a in angles if a}
    a,b=119.8,166.7
    grids[a]=grids[b]=cut_angles(np.r_[grids[a],grids[b]])
    stations=np.genfromtxt(out/branch/"stations.csv",names=True,delimiter=",")
    assert sha(out/branch/"stations.csv")==s["files"]["stations.csv"]
    curves={};stats=[]
    for angle,st in zip(angles,stations):
        assert angle==st["angleDeg"];rows=[]
        def add(phi,lo,hi,kind):
            rows.append(sample(native,body,direction(angle,phi),q,tw,c,angle,phi,lo,hi,kind))
        if angle==0:add(0.,0.,0.,"pole")
        else:
            for lo,hi in zip(grids[angle][:-1],grids[angle][1:]):
                add(lo,lo,lo,"edge");add((lo+hi)/2,lo,hi,"interval")
            # 接缝本身保持缺项。另在暴露侧取实际内点，使靠近管根的陡梯度进入实测范围。
            edges={r["phiDeg"]:r["status"] for r in rows if r["kind"]=="edge"}
            for r in list(rows):
                if r["kind"]!="interval" or r["status"]!="sphere":continue
                lo,hi=r["loDeg"],r["hiDeg"];inset=min(1e-5,(hi-lo)/4)
                for bound,near in ((lo,lo+inset),(hi,hi-inset)):
                    if edges[bound%360.]!="sphere":add(near,near,near,"inside-gap-boundary")
            rows.sort(key=lambda r:r["phiDeg"])
        r=rows[0];assert abs(r["nu"]-st["computedNu"])<1e-8 and abs(r["qInWm2"]-st["fluxInWm2"])<1e-6
        available=[v for v in rows if v["status"]=="sphere"]
        spans={}
        for v in rows:
            spans[v["status"]]=spans.get(v["status"],0.)+v["hiDeg"]-v["loDeg"]
        curves[angle]=rows
        stats.append(dict(angleDeg=angle,experimentNu=st["reportedNu"],plusZNu=r["nu"],spanDeg=spans,
                          nu=ranges(rows,"nu"),qInWm2=ranges(rows,"qInWm2"),
                          nuLower=min(v["nuLower"] for v in available),nuUpper=max(v["nuUpper"] for v in available)))
    summary=dict(schema="short65-azimuth-v1",caseSHA256=s["caseSHA256"],bodySHA256=s["bodySHA256"],
                 meshSHA256=s["meshSHA256"],frameSHA256=s["frameSHA256"],reactionSHA256=s["releaseSHA256"],
                 wallBudgetSHA256=endpoint["files"]["wall_budget.npz"],wallSHA256=s["files"]["wall.csv"],
                 elapsedS=s["elapsedSeconds"],dtS=s["dtSeconds"],stations=stats,rear=compare_rear(curves[a],curves[b]))
    text=io.StringIO();writer=csv.DictWriter(text,fieldnames=FIELDS,lineterminator="\n")
    writer.writeheader()
    for rows in curves.values():writer.writerows(rows)
    path=out/("stationAzimuth.zip" if branch=="half" else "stationAzimuthFull.zip")
    with ZipFile(path,"x") as z:
        for name,content in (("summary.json",json.dumps(summary,allow_nan=False)),("points.csv",text.getvalue())):
            z.writestr(ZipInfo(name),content.encode("utf8"),compress_type=ZIP_DEFLATED,compresslevel=9)
    return sha(path)
