"""由实际释放行、边修正及体温重建七个周向场，以独立射线求交和整圈包围审核。"""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import sys
from zipfile import ZipFile

for name in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):os.environ[name]="1"
import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];DATA=ROOT/"cases/short_sphere"
sys.path.insert(0,str(HERE.parent/"v333"))
from verify_time import normals_and_weights,residual,sha


def normalize(v):
    return v/np.linalg.norm(v,axis=-1,keepdims=True)


class Triangles:
    def __init__(self,xyz,tri,regions,ids):
        self.tri=tri;self.r=regions;self.ids=ids
        p=xyz[tri];a,b,c=p.transpose(1,0,2)
        n=np.cross(b-a,c-a);n2=(n*n).sum(axis=1)
        self.n=n;self.dist=(a*n).sum(axis=1)
        self.u=np.cross(c-a,n)/n2[:,None];self.v=np.cross(n,b-a)/n2[:,None]
        self.u0=(self.u*a).sum(axis=1);self.v0=(self.v*a).sum(axis=1)
        d=normalize(p);center=normalize(d.sum(axis=1))
        self.polar=np.arccos(np.clip(center[:,0],-1.,1.))
        self.cap=np.arccos(np.clip(np.min((d*center[:,None]).sum(axis=2),axis=1),-1.,1.))
        assert max(self.cap)<np.pi/2
        h=normalize(np.cross(d,np.roll(d,-1,axis=1)))
        self.inward=h*np.sign((h*np.roll(d,-2,axis=1)).sum(axis=2))[:,:,None]

    def band(self,angle):
        return np.flatnonzero(abs(self.polar-np.deg2rad(angle))<=self.cap+1e-10)

    def hit(self,ix,d):
        with np.errstate(divide="ignore",invalid="ignore"):
            t=self.dist[ix]/(self.n[ix]@d)
            u=t*(self.u[ix]@d)-self.u0[ix];v=t*(self.v[ix]@d)-self.v0[ix]
            good=(t>0)&(u>=-2.1e-10)&(v>=-2.1e-10)&(1-u-v>=-2.1e-10)
        faces=ix[good]
        if not len(faces):return "missing",faces,None,None
        if np.ptp(t[good])>2e-8:return "multiple",faces,None,None
        tags=set(self.r[faces])
        status="seam" if len(tags)>1 else {1:"sphere",2:"tube",3:"tube",4:"cap",5:"cap"}[next(iter(tags))]
        return status,faces,np.column_stack((1-u[good]-v[good],u[good],v[good])),t[good][0]*d


def ray(angle,phi):
    a,p=np.deg2rad([angle,phi]);return np.array([np.cos(a),np.sin(a)*np.sin(p),np.sin(a)*np.cos(p)])


def nu(q,t,c):
    k=.0241*(t/273.15)**1.5*(273.15+194)/(t+194)
    return -q*c["diameterM"]/(k*(t-c["airK"]))


def continuous_order(wall,body,q,tw,c):
    # 独立证明：对均匀 0.25° 方位窗，用三条球面半空间和最远可能方向偏移选出所有可能命中片。
    # 对每片热流与温度取顶点包围。该证明完全独立于生产算法的边交点、驻点和角度分区。
    windows=np.arange(0,360,.25);bands={}
    for a in (119.8,166.7):
        ix=wall.band(a);j=body.band(a)
        nodes=wall.tri[ix];wt=tw[nodes];outward=-q[nodes]
        assert np.min(outward)>0 and wt.min()>c["airK"]
        lower=nu(-outward.min(axis=1),wt.max(axis=1),c);upper=nu(-outward.max(axis=1),wt.min(axis=1),c)
        radius=np.sin(np.deg2rad(a))*np.deg2rad(.25)/2+1e-10
        low=[];high=[]
        for left in windows:
            d=ray(a,left+.125)
            possible=(wall.inward[ix]@d>=-radius).all(axis=1)
            original=(body.inward[j]@d>=-radius).all(axis=1)
            assert possible.any() and original.any() and np.all(wall.r[ix[possible]]==1) and np.all(body.r[j[original]]==1)
            low.append(lower[possible].min());high.append(upper[possible].max())
        bands[a]=(np.array(low),np.array(high))
    lower=bands[119.8][0]-bands[166.7][1]
    return dict(windows=len(windows),minCommonDifferenceBound=float(min(lower)),
                lower119=float(bands[119.8][0].min()),upper166=float(bands[166.7][1].max()))


def raw_sources():
    app=ROOT/"outputs/v336_validation/clean01"
    job=app/"outputs/short65_time_6dea1a68-073d-4129-867d-3772b804b131"
    stage=job/"run/half/step2";e=json.loads((stage/"energy.json").read_bytes())
    assert e["frame"]==sha(stage/"restart_00002.csv") and sha(stage/"edges.bin")==e["files"]["edges.bin"]
    assert sha(stage/"release/restart_00002.csv")==json.loads((stage/"release/native.json").read_bytes())["frame"]
    with np.load(DATA/"supportedAngularTransfer65.npz") as g:
        xyz,corners,area=(g[k] for k in ("nativeWallXYZD","nativeCorners","nativeArea"))
    prep=app/json.loads((job/"construction.json").read_bytes())["prepared"]
    with np.load(prep/"mapping.npz") as m:r=m["group"]
    state=np.loadtxt(stage/"restart_00002.csv",delimiter=",",skiprows=1,usecols=(1,2,3,9))
    dist,wet=cKDTree(state[:,:3]).query(xyz);assert max(dist)<1e-12
    edges=np.fromfile(stage/"edges.bin","<f8").reshape(-1,6);i,j=edges[:,:2].astype(int).T
    delta=edges[:,4]*(1-edges[:,5])
    power=(residual(stage/"release/restart_00002.csv")+np.bincount(i,delta,minlength=len(state))-np.bincount(j,delta,minlength=len(state)))[wet]*e["scale"]
    _,weights=normals_and_weights(xyz,corners)
    dual=np.bincount(corners.ravel(),(area[:,None]*weights).ravel(),minlength=len(wet))
    assert abs(power.sum()+e["wallW"])<1e-10
    tri=np.concatenate([corners[:,[0,1,2]],corners[:,[2,3,0]]]);keep=np.any(np.cross(xyz[tri[:,1]]-xyz[tri[:,0]],xyz[tri[:,2]]-xyz[tri[:,0]])!=0,axis=1)
    wall=Triangles(xyz,tri[keep],np.tile(r,2)[keep],np.tile(np.arange(len(r))+1,2)[keep])
    pair=json.loads((DATA/"supportedTime65.json").read_bytes());c=pair["condition"]
    with ZipFile(DATA/"supported65.zip") as z:g=json.loads(z.read("geometry.json"))
    f=np.array(g["F"])-1
    body=Triangles(np.array(g["V"])/c["diameterM"],f,np.array(g["regionID"]),np.arange(len(f))+1)
    return wall,body,power/dual,c["airK"]+c["scalarUnitK"]*state[wet,3],c


def verify(file,volumes):
    with ZipFile(file) as z:
        summary=json.loads(z.read("summary.json"));rows=list(csv.DictReader(io.StringIO(z.read("points.csv").decode())))
    pair=json.loads((DATA/"supportedTime65.json").read_bytes())
    assert summary["frameSHA256"]==pair["halfFrame"] and summary["meshSHA256"]==pair["meshSHA256"]
    assert summary["caseSHA256"]==sha(DATA/"test65.json") and summary["bodySHA256"]==sha(DATA/"supported65.zip")
    assert summary["elapsedS"]==pair["elapsedS"] and len(summary["stations"])==7
    if volumes:wall,body,q,tw,c=raw_sources()
    gaps=[];max_error=0.;tested=0;rear={}
    for station,expected in zip(summary["stations"],pair["halfNu"]):
        angle=station["angleDeg"];band=[r for r in rows if float(r["angleDeg"])==angle]
        assert abs(float(band[0]["nu"])-expected)<1e-8
        intervals=[r for r in band if r["kind"]=="interval"]
        if angle:
            ends=[float(r["hiDeg"]) for r in intervals]
            assert [float(r["loDeg"]) for r in intervals]==[0.]+ends[:-1] and ends[-1]==360.
        else:assert len(band)==1 and band[0]["kind"]=="pole"
        valid=[r for r in band if r["status"]=="sphere"]
        for k,field in (("nu","nu"),("qInWm2","qInWm2")):
            for agg,method in (("minimum",min),("maximum",max)):
                r=method(valid,key=lambda r:float(r[field]));ref=station[k][agg]
                assert float(r[field])==ref["value"] and int(r["faceID"])==ref["faceID"] and float(r["phiDeg"])==ref["phiDeg"]
        if angle==90:
            for r in intervals:
                if r["status"]=="sphere":continue
                a,b=float(r["loDeg"]),float(r["hiDeg"])
                if gaps and gaps[-1][1]==a:gaps[-1][1]=b
                else:gaps.append([a,b])
        if angle in (119.8,166.7):rear[angle]=band
        if not volumes:continue
        wi,bi=wall.band(angle),body.band(angle)
        for r in band:
            phis=[float(r["phiDeg"])]
            if r["kind"]=="interval":
                a,b=float(r["loDeg"]),float(r["hiDeg"]);phis.extend([a+.2113*(b-a),a+.7887*(b-a)])
            for j,p in enumerate(phis):
                direction=ray(angle,p)
                s,faces,w,pos=wall.hit(wi,direction);g,index,_,_=body.hit(bi,direction)
                status=s if s==g else ("seam" if "seam" in (s,g) else "mesh-disagreement")
                assert status==r["status"],(angle,p,status,r["status"])
                if j==0 and len(index):assert int(r["faceID"])==body.ids[index].min() and int(r["faceCount"])==len(set(body.ids[index]))
                if status!="sphere":
                    assert r["nu"]==r["qInWm2"]=="";continue
                f=np.sum(w*q[wall.tri[faces]],axis=1);T=np.sum(w*tw[wall.tri[faces]],axis=1);v=nu(f,T,c)
                if j==0:
                    max_error=max(max_error,max(abs(f-float(r["qInWm2"]))))
                    assert max(abs(v-float(r["nu"])))<1e-8 and max(abs(T-float(r["wallK"])))<1e-8
                    np.testing.assert_allclose(pos*c["diameterM"],[float(r[k]) for k in ("xM","yM","zM")],rtol=0,atol=1e-11)
                else:
                    assert v.min()>=float(r["nuLower"])-1e-8 and v.max()<=float(r["nuUpper"])+1e-8
                    tested+=1
    front,back=rear.values();assert len(front)==len(back)
    assert all(a["phiDeg"]==b["phiDeg"] and a["status"]==b["status"]=="sphere" for a,b in zip(front,back))
    min_bound=min(float(a["nuLower"])-float(b["nuUpper"]) for a,b in zip(front,back))
    assert min_bound==summary["rear"]["lowerBound"] and min_bound>0 and summary["rear"]["status"]=="positive-throughout"
    proof=dict(rows=len(rows),equatorGapDegrees=gaps,productionDifferenceLower=min_bound)
    if volumes:
        independent=continuous_order(wall,body,q,tw,c)
        assert independent["minCommonDifferenceBound"]>0 and max_error<1e-6
        proof.update(maxRawHeatDifferenceWm2=max_error,interiorChecks=tested,independentCoverage=independent)
    return proof


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--volumes",action="store_true");p.add_argument("--output",type=Path)
    a=p.parse_args();proof=verify(DATA/"stationAzimuth65.zip",a.volumes)
    if a.output:
        with a.output.open("x",encoding="utf8",newline="\n") as f:json.dump(proof,f,indent=2);f.write("\n")
    print(json.dumps(proof))
