"""独立复核银球近壁诊断。--public 只需随版本交付的小表，完整核对另读原始场。"""
import argparse
import csv
from collections import defaultdict
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/"cases/short_sphere"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(archive,name):
    with ZipFile(archive) as z:
        return np.genfromtxt(io.BytesIO(z.read(name)),delimiter=",",names=True)


def close(a,b):
    np.testing.assert_allclose(a,b,atol=4e-13,rtol=2e-12,equal_nan=True)


def unit(x):
    return x/np.linalg.norm(x,axis=-1,keepdims=True)


def normals_and_weights(coords,corners):
    # 将壁片划成顶点、下一顶点及质心组成的扇形，半个扇形归给各相邻顶点。
    normals=np.empty((len(corners),3));weights=np.zeros(corners.shape)
    for i,ids in enumerate(corners):
        slots=sorted(np.unique(ids,return_index=True)[1]);p=coords[ids[slots]]
        offset=p-p.mean(axis=0);fan=np.cross(offset,np.roll(offset,-1,axis=0))
        normals[i]=unit(fan.sum(axis=0))
        area=fan@normals[i];weights[i,slots]=.5*(area+np.roll(area,1))/sum(area)
    return normals, weights


def cross_zone(faces,region):
    incidence=defaultdict(set)
    for f,r in zip(faces,region):
        for i in f: incidence[i].add(r)
    boundary=[i for i,r in incidence.items() if len(r)>1]
    return np.any(np.isin(faces,boundary),axis=1)


def verify_ledger(section,ids,area,q,status):
    for code,key in ((-1,"reverse"),(0,"unknown"),(1,"forward")):
        local=ids[status[ids]==code];r=section[key];power=-area[local]*q[local]
        assert r["faces"]==len(local)
        close([r["areaM2"],r["outW"],r["releaseW"],r["absorbedW"]],
              [sum(area[local]),sum(power),sum(np.maximum(power,0)),sum(np.maximum(-power,0))])


def raw_check(k,report,m,n,tau_dir,defined,rows,folder=None,sectors=64,frame_suffix="000199"):
    if folder is None:
        folder=ROOT/"outputs"/("v329_validation" if k=="coarse" else "v330_validation")/"run01/steady03"
    frame=folder/("restart_"+frame_suffix+".csv");surface=folder/("surface_"+frame_suffix+".vtk")
    assert digest(frame)==report["frame"] and digest(surface)==report["vtk"]
    state=np.loadtxt(frame,skiprows=1,delimiter=",",usecols=range(8))
    lines=surface.read_text("ascii").splitlines()
    def vtk_array(header):
        i=next(i for i,line in enumerate(lines) if line.startswith(header))
        return np.fromstring(lines[i+1],sep=" ").reshape(-1,3)
    wall=m["nativeWallXYZD"];tree=cKDTree(state[:,1:4])
    def identify(points):
        d,i=tree.query(points);assert max(d)<1e-12 and len(set(i))==len(points)
        return i
    wet=identify(wall);first=m["firstXYZD"]
    if k=="fine":first=(first+wall)/2
    first=identify(first)
    d,order=cKDTree(wall).query(vtk_array("POINTS "))
    assert max(d)<1e-6 and len(set(order))==len(wall)
    u=4.06*.3048;rho=1.1015571271965614;scales=[rho*u*u/2,u]
    fields=[vtk_array("VECTORS Skin_Friction_Coefficient ")[np.argsort(order)]*scales[0],state[first,5:8]*u]
    c=m["nativeCorners"];sn,weights=normals_and_weights(wall,c)
    target,source=m["target"],m["source"]
    transfer=csr_matrix((m["targetFraction"],(target,source)),shape=(len(n),len(sn)))
    result={};signs=[];mixed=[];small=[]
    for i,(signed,mag) in enumerate((("ShearPa","ShearMagnitudePa"),("SpeedMps","SpeedMagnitudeMps"))):
        nodal=fields[i][c];averaged=np.einsum("ij,ijk->ik",weights,nodal)
        source_v=np.cross(sn,np.cross(averaged,sn))
        projected=np.cross(n,np.cross(transfer@source_v,n))
        actual=(projected*tau_dir).sum(axis=1);actual[~defined]=np.nan
        size=np.linalg.norm(projected,axis=1)
        close(actual,rows[k+signed]);close(size,rows[k+mag])
        result[signed]=float(np.nanmax(abs(actual-rows[k+signed])))
        result[mag]=float(np.max(abs(size-rows[k+mag])))
        amp=np.sum(weights*np.linalg.norm(nodal,axis=2),axis=1);eps=64*np.finfo(float).eps
        tol=(1e-6 if i==0 else eps)*amp+eps*scales[i]
        local=np.einsum("ij,ij->i",source_v[source],tau_dir[target])
        conflicting=np.intersect1d(target[local>tol[source]],target[local< -tol[source]])
        mixed.append(np.isin(np.arange(len(n)),conflicting)&defined)
        threshold=transfer@tol
        small.append((abs(actual)<=threshold)&defined)
        signs.append(np.where(abs(actual)>threshold,np.sign(actual),0))
    zero=~defined
    axis=m["axis"];boundary=cross_zone(m["f"],m["region"])
    near_boundary=np.unique(target[cross_zone(c,m["nativeRegion"])[source]])
    masks=[axis,zero,boundary,np.isin(np.arange(len(n)),near_boundary),*small,signs[0]*signs[1]<0,*mixed]
    flags=sum(b.astype(int)*2**j for j,b in enumerate(masks))
    np.testing.assert_array_equal(flags,rows[k+"Flags"])
    np.testing.assert_array_equal(np.where(flags,0,signs[0]),rows[k+"Status"])
    # 球面测线用解析的增角切向和射线/弦的二维行列式，重新定位原节点及插值权重。
    ring=np.flatnonzero(abs(wall[:,1])<1e-12)
    phi=np.arctan2(wall[ring,2],wall[ring,0]);ring=ring[np.argsort(phi)];phi=np.sort(phi)
    assert len(ring)==sectors
    for r in report["stations"]:
        a=np.deg2rad(r["angleDeg"]);j=np.searchsorted(phi,a)
        ids=ring[[j-1,j]]
        tangent=np.array([-np.sin(a),0.,np.cos(a)]);dist=wall[ids]@tangent;weight=dist[0]/(dist[0]-dist[1])
        close(weight,r["rightWeight"]);np.testing.assert_array_equal(wet[ids],[r["leftNode"],r["rightNode"]])
        if a==0:
            assert r["status"]==0 and r["shearPa"] is None and r["speedMps"] is None
            continue
        expected=[np.array([1-weight,weight])@(field[ids]@tangent) for field in fields]
        close(expected,[r["shearPa"],r["speedMps"]])
        assert np.sign(expected[0])==np.sign(expected[1])==r["status"]
    positive=wall[ring,2]>1e-12
    half=ring[positive];phi=phi[positive];t=np.column_stack((-np.sin(phi),np.zeros_like(phi),np.cos(phi)))
    signs=np.sign(np.stack([np.einsum("ij,ij->i",field[half],t) for field in fields]))
    code=np.where(signs[0]==signs[1],signs[0],0)
    i=np.flatnonzero(np.diff(code));trans=report["transitions"]
    close(np.rad2deg(np.column_stack((phi[i],phi[i+1]))),[[e["left"]["angleDeg"],e["right"]["angleDeg"]] for e in trans])
    return result


def verify(folder,public):
    info=json.loads((folder/"supportedDirection65.json").read_bytes())
    archive=folder/"supportedDirection65.zip";assert digest(archive)==info["archiveSHA256"]
    w=table(archive,"direction.csv");ids=np.arange(len(w));np.testing.assert_array_equal(w["faceID"],ids+1)
    with ZipFile(DATA/"supported65.zip") as z:body=json.loads(z.read("geometry.json"))
    with np.load(DATA/"supportedTransfer65.npz") as data:m=dict(data)
    f=np.array(body["F"])-1;m["f"]=f;p=np.array(body["V"])[f];area=m["targetArea"]
    n=unit(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]))
    t=np.cross(n,np.cross(np.array([-1.,0.,0.]),n));length=np.linalg.norm(t,axis=1)
    defined=length>1e-10;t[defined]/=length[defined,None]
    m["axis"]=np.any(np.sum(p[:,:,1:]**2,axis=2)<1e-24,axis=1)&(m["region"]==1)
    np.testing.assert_array_equal(w["regionID"],m["region"]);close(area.sum(),info["areaM2"])
    edge_faces=defaultdict(list)
    for i,tri in enumerate(f):
        for a,b in zip(tri,np.roll(tri,1)): edge_faces[tuple(sorted((a,b)))].append(i)
    adj=[set() for _ in w]
    for a,b in edge_faces.values(): adj[a].add(b);adj[b].add(a)
    result=dict(faces=len(ids),publicOnly=public,archiveSHA256=info["archiveSHA256"])
    for k,name in (("coarse","supportedSurface65.zip"),("fine","supportedFine65.zip")):
        report=info[k];q=table(DATA/name,"wall.csv")["fluxInWm2"];status=w[k+"Status"]
        verify_ledger(report["total"],ids,area,q,status)
        for zone,r in report["regions"].items():
            verify_ledger(r,ids[m["region"]==int(zone)],area,q,status)
        def neighbourhood(r,core):
            around=set().union(*(adj[i] for i in core))-set(core)
            verify_ledger(r["core"],np.array(core),area,q,status)
            verify_ledger(r["neighbours"],np.array(sorted(around)),area,q,status)
        for r in [report["peak"],*report["stations"]]:neighbourhood(r,[r["faceID"]-1])
        neighbourhood(report["nearPeak"],ids[-q>=.99*max(-q)])
        remaining=set(ids[status==-1]);components={}
        while remaining:
            first=remaining.pop();stack=[first];group={first}
            while stack:
                batch=adj[stack.pop()]&remaining;remaining-=batch;stack.extend(batch);group.update(batch)
            components[min(group)+1]=np.array(sorted(group))
        assert len(components)==len(report["components"])
        for r in report["components"]:
            members=components[r["firstFaceID"]];assert len(members)==r["faces"]
            close([sum(area[members]),-area[members]@q[members]],[r["areaM2"],r["outW"]])
        if not public:result[k+"MaxAbsError"]=raw_check(k,report,m,n,t,defined,w)
    c=w["coarseStatus"]==-1;f=w["fineStatus"]==-1;o=info["overlap"]
    close([sum(area[c&f]),sum(area[c|f])],[o["bothM2"],o["eitherM2"]]);assert sum(c&f)==o["bothFaces"]
    return result


def verify_page(folder,ui):
    proof=json.loads((ui/"ui.json").read_bytes());g=proof["grid"]
    page=json.loads((ui/"workbench/case-data.js").read_text("utf8").split("=",1)[1].strip().rstrip(";"))
    rows=table(folder/"supportedDirection65.zip","direction.csv")
    info=json.loads((folder/"supportedDirection65.json").read_bytes())
    fields={f["key"]:f for f in page["fields"]};prefix="support"+g.title()
    assert page["defaultField"]==prefix+"Status" and proof["frame"]==info[g]["frame"]
    assert page["evidenceIdentity"]["sourceFrame"]["shearVTKSHA256"]==info[g]["vtk"]
    for suffix in ("ShearPa","SpeedMps","ShearMagnitudePa","SpeedMagnitudeMps","Status"):
        field=fields[prefix+suffix];expected=rows[g+suffix].copy()
        if suffix=="Status":expected[expected==0]=np.nan
        close(np.array(field["values"],float),expected)
        np.testing.assert_array_equal(field["status"],np.isfinite(expected))
        assert all(field["missingReason"][i] for i in np.flatnonzero(~np.isfinite(expected)))
    cover=page["wallCoverage"]
    assert cover["jointComputedFaceCount"]==len(rows) and cover["fullSurfaceComputed"]
    close(cover["wholeSurfaceHeatPowerW"],-sum(v["outW"] for v in info[g]["total"].values()))
    with (ui/"shear_section.csv").open(encoding="utf-8-sig",newline="") as file:
        section=list(csv.DictReader(file))
    assert section
    for row in section:
        assert row["field_key"]==prefix+"ShearPa" and row["unit"]=="Pa"
        close(float(row["value"]),rows[g+"ShearPa"][int(row["face_id"])-1])
    return dict(grid=g,sectionRows=len(section),faces=len(rows),fieldsChecked=5)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("folder",type=Path,nargs="?",default=DATA)
    parser.add_argument("output",type=Path,nargs="?");parser.add_argument("--public",action="store_true")
    parser.add_argument("--ui",type=Path,help="核对指定工作台导出，不读取大型原生场")
    args=parser.parse_args();result=verify_page(args.folder,args.ui) if args.ui else verify(args.folder,args.public)
    if args.output:
        with args.output.open("x",encoding="ascii",newline="\n") as out:json.dump(result,out,indent=2);out.write("\n")
    print(result)
