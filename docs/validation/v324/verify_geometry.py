"""从原平面三角形重新裁剪覆盖面积，并独立核对支管接缝、完整边界和七站。"""
import io
import json
import math
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def close(a,b,tol):
    np.testing.assert_allclose(a,b,atol=tol,rtol=0)


def clipped_area(points,polygon):
    # 在三维中依次裁掉各管侧平面以外的部分，独立求原三角面被圆管多边形遮盖的面积。
    for p,q in zip(polygon,np.roll(polygon,-1,axis=0)):
        edge=q-p
        signed=(points[:,2]-p[1])*edge[0]-(points[:,0]-p[0])*edge[1]
        if (signed<0).all():
            return 0.
        if (signed>=0).all():
            continue
        out=[]
        for j,end in enumerate(points):
            a,b=signed[j-1],signed[j]
            if (a>=0)!=(b>=0):
                out.append(points[j-1]+a/(a-b)*(end-points[j-1]))
            if b>=0:
                out.append(end)
        points=np.asarray(out)
    return np.linalg.norm(np.cross(points-points[0],np.roll(points,-1,axis=0)-points[0]).sum(0))/2


def read_bundle(path):
    with ZipFile(path) as z:
        assert z.testzip() is None
        geo=json.loads(z.read("geometry.json")); summary=json.loads(z.read("summary.json"))
        table=lambda name:np.genfromtxt(io.BytesIO(z.read(name)),delimiter=",",names=True)
        parents,stations=table("parents.csv"),table("stations.csv")
        binary=z.read("supported_m.stl")
    vertices=np.array(geo["V"]);faces=np.array(geo["F"])-1
    assert int.from_bytes(binary[80:84],"little")==len(faces) and len(binary)==84+50*len(faces)
    stl=np.frombuffer(binary,offset=84,dtype=[("normal","<f4",(3,)),("v","<f4",(3,3)),("attribute","<u2")])
    np.testing.assert_array_equal(stl["v"],vertices[faces].astype("<f4"))
    return geo,summary,parents,stations,vertices,faces


def check(path,body,extension):
    g,s,p,stations,v,f=read_bundle(path)
    parent=np.array(g["parentID"]);region=np.array(g["regionID"]);tri=v[f]
    vec=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(vec,axis=1)/2
    close(area,g["areaM2"],1e-18)
    assert area.min()>0 and np.unique(np.sort(f,axis=1),axis=0).shape==f.shape
    assert np.array_equal(np.unique(f),np.arange(len(v))) and np.isfinite(v).all()
    edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
    pairs,idx,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    assert np.all(count==2) and np.all(np.bincount(idx,np.sign(edges[:,1]-edges[:,0]))==0)
    adjacent=np.tile(region,3)[np.argsort(idx)].reshape(-1,2)
    graph=coo_matrix((np.ones(len(pairs)),pairs.T),shape=(len(v),len(v)))
    assert connected_components(graph,directed=False)[0]==1 and len(v)-len(pairs)+len(f)==2
    signed_volume=np.einsum("ij,ij->",tri[:,0],vec)/6
    assert signed_volume>0
    D=.0127;R=D/2;tube_r=.096*.0254/2;n=64
    phi=(np.arange(n)+.5)*2*np.pi/n
    poly=tube_r*np.c_[np.cos(phi),np.sin(phi)]
    diam_end=D*(.5+extension)
    assert s["extensionD"]==extension
    ref_v=np.array(body["V"]);ref_f=np.array(body["F"])-1;ref=ref_v[ref_f]
    old_normal=np.cross(ref[:,1]-ref[:,0],ref[:,2]-ref[:,0])
    original=np.linalg.norm(old_normal,axis=1)/2
    covered=np.zeros(len(ref))
    for i,t in enumerate(ref):
        xz=t[:,[0,2]]
        if np.all(xz.min(0)<=tube_r) and np.all(xz.max(0)>=-tube_r):
            covered[i]=clipped_area(t,poly)
    sphere=(region==1);pid=parent[sphere]-1
    assert np.all(parent[~sphere]==0) and np.all((pid>=0)&(pid<len(ref)))
    # 每个外露子三角形与父三角形共面、同向，三个新顶点都在父面的重心范围内。
    pr=ref[pid];normal=old_normal[pid];nn=np.sum(normal**2,axis=1)
    assert np.all(np.sum(vec[sphere]*normal,axis=1)>0)
    for t in tri[sphere].transpose(1,0,2):
        e=t-pr[:,0]
        close(np.sum(e*normal,axis=1)/np.sqrt(nn),0,2e-14)
        b=np.sum(np.cross(e,pr[:,2]-pr[:,0])*normal,axis=1)/nn
        c=np.sum(np.cross(pr[:,1]-pr[:,0],e)*normal,axis=1)/nn
        assert min(b.min(),c.min(),(1-b-c).min())>=-1e-10
    exposed=np.bincount(pid,area[sphere],minlength=len(ref))
    close(p["originalM2"],original,1e-18);close(p["exposedM2"],exposed,1e-16)
    close(p["coveredM2"],covered,1e-16)
    np.testing.assert_array_equal(p["parentID"],np.arange(1,len(ref)+1))
    close(original,exposed+covered,1e-16)
    region_sum=np.array([math.fsum(area[region==k]) for k in range(1,6)])
    close(s["regionAreaM2"],region_sum,1e-16)
    close([sum(region_sum[:3]),sum(region_sum[3:])],[s["solverWallAreaM2"],s["artificialCutAreaM2"]],1e-16)
    perimeter=2*n*tube_r*np.sin(np.pi/n)
    expected_cap=.5*n*tube_r**2*np.sin(2*np.pi/n)
    joins=[]
    for sign,wall,cap,key in [(1,2,4,"joinPositive"),(-1,3,5,"joinNegative")]:
        seam=np.array(g[key])-1
        np.testing.assert_array_equal(seam[:,1],np.roll(seam[:,0],-1))
        assert len(np.unique(seam[:,0]))==len(seam)
        boundary=pairs[np.all(np.sort(adjacent,axis=1)==[1,wall],axis=1)]
        assert set(map(tuple,boundary))==set(map(tuple,np.sort(seam,axis=1)))
        corners=v[seam]; ds=np.linalg.norm((corners[:,1]-corners[:,0])[:,[0,2]],axis=1)
        height=diam_end-sign*corners[:,:,1].mean(1)
        close(ds.sum(),perimeter,1e-14);close(ds@height,area[region==wall].sum(),1e-16)
        close(area[region==cap].sum(),expected_cap,1e-16)
        close(tri[region==cap,:,1],sign*diam_end,1e-14)
        assert np.all(sign*vec[region==cap,1]>0)
        close(vec[region==wall,1],0,1e-16)
        assert np.all(np.sum(vec[region==wall][:,[0,2]]*tri[region==wall,0][:,[0,2]],axis=1)>0)
        joins.append(len(seam))
    close(stations["angleDeg"],[0,41.3,60.2,75.8,90,119.8,166.7],1e-12)
    old_meridian=ref_v[(abs(ref_v[:,1])<1e-15)&(ref_v[:,2]>=-1e-15)&(abs(np.linalg.norm(ref_v,axis=1)-R)<1e-14)]
    theta=np.arctan2(old_meridian[:,2],old_meridian[:,0]);order=np.argsort(theta)
    assert len(theta)==49
    for row in stations:
        i=int(row["faceID"])-1
        assert region[i]==1 and parent[i]==row["parentID"] and exposed[parent[i]-1]>0
        b=np.array([row["w1"],row["w2"],row["w3"]])
        close(b.sum(),1.,1e-13);assert min(b)>-1e-10
        pos=np.array([row["xM"],row["yM"],row["zM"]]);close(b@tri[i],pos,1e-14)
        angle=np.radians(row["angleDeg"]);j=np.clip(np.searchsorted(theta[order],angle),1,48)
        l,h=theta[order][j-1:j+1];distance=R*np.cos((h-l)/2)/np.cos(angle-(l+h)/2)
        close(pos,distance*np.array([np.cos(angle),0,np.sin(angle)]),1e-14)
    lengths=np.linalg.norm(np.roll(tri,-1,axis=1)-tri,axis=2)
    unit=(np.roll(tri,-1,axis=1)-tri)/lengths[:,:,None]
    cos=np.sum(-np.roll(unit,1,axis=1)*unit,axis=2)
    angles=np.degrees(np.arccos(np.clip(cos,-1,1))).min(1)
    thin=angles<5
    return dict(faces=len(f),vertices=len(v),volumeM3=signed_volume,joinEdgeCounts=joins,originalAreaM2=original.sum(),
                exposedSphereM2=exposed.sum(),coveredSphereM2=covered.sum(),regionAreasM2=region_sum.tolist(),
                maxParentErrorM2=max(abs(original-covered-exposed)),analyticCoveredM2=4*np.pi*R*(R-np.sqrt(R**2-tube_r**2)),
                unaffected=int((p["coveredM2"]==0).sum()),fullyCovered=int((p["exposedM2"]==0).sum()),
                partiallyCovered=int(((p["exposedM2"]>0)&(p["coveredM2"]>0)).sum()),
                minimumAngleDeg=angles.min(),minAngleFaceID=int(np.argmin(angles)+1),belowFiveDegrees=int(thin.sum()),
                belowFiveByRegion=np.bincount(region[thin],minlength=6)[1:].tolist(),thinAreaM2=area[thin].sum(),
                minEdgeM=lengths.min(),physicalPressureHeatCoverageM2=s["pressureHeatSolvedM2"],stationParentIDs=stations["parentID"].astype(int).tolist())


if __name__=="__main__":
    nominal,short,ui,dest=map(Path,sys.argv[1:5])
    resources={}
    for path in (nominal,short,ui):
        rec=json.loads((path/"debugger_process.json").read_bytes())
        assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
        assert not (rec["stopReason"] or rec["survivingProcessIDs"] or rec["supervisorError"])
        resources[path.name]={"seconds":rec["elapsedSeconds"],"peakRSSBytes":rec["peakFamilyRSSBytes"]}
    with ZipFile(ROOT/"cases/short_sphere/surface65.zip") as z:body=json.loads(z.read("body.json"))
    a=check(nominal/"supported.zip",body,1)
    b=check(short/"supported.zip",body,.5)
    close(a["regionAreasM2"][0],b["regionAreasM2"][0],1e-16)
    close(a["regionAreasM2"][3:],b["regionAreasM2"][3:],1e-16)
    close(np.subtract(a["regionAreasM2"][1:3],b["regionAreasM2"][1:3]),.00635*64*.096*.0254*np.sin(np.pi/64),1e-16)
    with ZipFile(nominal/"supported.zip") as z, ZipFile(ui/"supported.zip") as other:
        assert set(z.namelist())==set(other.namelist())
        for name in z.namelist():assert z.read(name)==other.read(name),name
    proof=json.loads((ui/"ui.json").read_bytes())
    assert proof["faces"]==a["faces"] and proof["pixels"]>0 and proof["steadyReturned"]
    with dest.open("x",encoding="utf8") as f:
        json.dump(dict(nominal=a,shorter=b,ui=proof,resources=resources),f,indent=2);f.write("\n")
