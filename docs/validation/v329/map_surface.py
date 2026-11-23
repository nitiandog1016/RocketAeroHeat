"""把同一守恒端点的节点反力和压力经原生壁片交叠分配至名义双支管五区。"""
import argparse
import csv
import io
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np

from run_shared import HERE, ROOT, PREP, CASE, read, sha, budget, table
from shortSphere65 import conductivity
from thermistorWireSurface import residual

DATA=ROOT/"cases/short_sphere"


def corner_weights(corners,xyz):
    # 用有向面积在面法向上的投影划分中点对偶面积；轴上三角片只计三个不同节点。
    weights=np.zeros(corners.shape)
    for k,ids in enumerate(corners):
        _,positions=np.unique(ids,return_index=True);positions.sort()
        p=xyz[ids[positions]];center=p.mean(axis=0)
        forward=(p+np.roll(p,-1,axis=0))*.5
        previous=np.roll(forward,1,axis=0)
        sector=.5*(np.cross(forward-p,center-p)+np.cross(center-p,previous-p))
        area=sector.sum(axis=0);share=sector@area/(area@area)
        assert np.all(share>0) and abs(share.sum()-1)<1e-13
        weights[k,positions]=share
    return weights


def measure_stations(xyz,wet,sphere,q,tw,physical,case):
    # 测线在 Y=0 的真实 64 点环上。各测站的射线与环的直弦相交，插值节点量。
    ring=sphere[np.abs(xyz[sphere,1])<1e-10]
    assert len(ring)==64
    phi=np.mod(np.arctan2(xyz[ring,2],xyz[ring,0]),2*np.pi);order=np.argsort(phi)
    ring=ring[order];phi=phi[order];rows=[]
    for deg,reference in zip(case["angleDeg"],case["localNuSurface"]):
        a=np.deg2rad(deg);j=np.searchsorted(phi,a)%len(phi);i=(j-1)%len(phi)
        g,h=ring[i],ring[j]
        tangent=np.array([-np.sin(a),0,np.cos(a)])
        t=-(xyz[g]@tangent)/((xyz[h]-xyz[g])@tangent)
        assert -1e-12<=t<=1+1e-12
        ix=np.searchsorted(wet,[g,h]);w=np.array([1-t,t])
        temp=w@tw[ix];qin=w@q[ix]
        nu=-qin*physical["diameterM"]/((temp-physical["airK"])*conductivity(temp))
        rows.append([deg,reference,nu,nu-reference,nu/reference-1,qin,temp,int(g),int(h),t])
    return np.array(rows)


def mapped(work,out,*,prepared=PREP,expected_mesh=None):
    thermal=read(work/"thermal.json")
    assert thermal["passed"] and thermal["thermalChange"]<.002 and thermal["controlError"]<2e-12
    config=read(prepared/"prepared.json");s=config["condition"];case=read(CASE)
    grid=read(DATA/"supportedGrid65.json")
    assert sha(DATA/"supported65.zip")==grid["bodySHA256"] and sha(DATA/"supportedTransfer65.npz")==grid["publicTransferSHA256"]
    mesh_key=grid["hybrid"]["meshSHA256"] if expected_mesh is None else expected_mesh
    assert sha(prepared/"mesh.su2")==mesh_key==thermal["meshSHA256"] and grid["caseSHA256"]==thermal["caseSHA256"]
    with ZipFile(DATA/"supported65.zip") as z:
        body=json.loads(z.read("geometry.json"))
        reference=list(csv.DictReader(io.StringIO(z.read("stations.csv").decode("ascii"))))
    with np.load(prepared/"mapping.npz") as mesh:
        xyz,wet,first,corners,sourceArea,group,theta=(mesh[k] for k in ("xyz","wet","first","corners","nativeArea","group","theta"))
    for name,release in (("early","releaseEarly"),("late","releaseLate")):
        actual=budget(work/name,work/release,wet,thermal["scale"])
        assert actual==thermal[name] and abs(actual["physicalGapFraction"])<.005
    out.mkdir(parents=True,exist_ok=False)
    with np.load(prepared/"transfer.npz") as m:
        target,source,wTarget,wSource,area,regions=(m[k] for k in ("target","source","targetFraction","sourceFraction","targetArea","region"))
    assert sha(prepared/"mapping.npz")==config["files"]["mapping.npz"] and sha(prepared/"transfer.npz")==config["files"]["transfer.npz"]
    with np.load(DATA/"supportedTransfer65.npz") as public:
        np.testing.assert_array_equal(public["nativeWallXYZD"],xyz[wet])
        np.testing.assert_array_equal(public["nativeCorners"],np.searchsorted(wet,corners))
        for name,values in (("target",target),("source",source),("targetFraction",wTarget),("sourceFraction",wSource)):
            np.testing.assert_array_equal(public[name],values)
    np.testing.assert_array_equal(area,body["areaM2"]);np.testing.assert_array_equal(regions,body["regionID"])
    np.testing.assert_array_equal(regions[target],group[source])
    v=np.asarray(body["V"]);f=np.array(body["F"])-1;tri=v[f]
    np.testing.assert_allclose(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2,area,rtol=1e-12,atol=1e-19)
    state=table(work/"early","initial_flow.csv");frame=sha(work/"early/initial_flow.csv")
    assert frame==thermal["early"]["sourceFrame"]
    n=np.fromfile(work/"early/nodes.bin",dtype="<f8").reshape(-1,16);n=n[n[:,0].argsort()]
    np.testing.assert_array_equal(n[:,1],state[:,9])
    raw=residual(work/"releaseEarly/restart.csv")
    # 正号入壁；端点原生空间行的正壁反力对应热壁向空气放热。
    nodeW=(raw[wet]+n[wet,9]-n[wet,10])*thermal["scale"]
    assert abs(nodeW.sum()+thermal["early"]["wallW"])<1e-10
    ix=np.searchsorted(wet,corners);np.testing.assert_array_equal(wet[ix],corners)
    weights=corner_weights(corners,xyz);shares=weights*sourceArea[:,None]
    lump=np.bincount(ix.ravel(),weights=shares.ravel(),minlength=len(wet));assert np.min(lump)>0
    nodeQ=nodeW/lump;sourceW=(shares*nodeQ[ix]).sum(axis=1)
    outW=np.bincount(target,wSource*sourceW[source],minlength=len(area))
    def intensive(nodal):
        mean=(weights*nodal[ix]).sum(axis=1)
        return np.bincount(target,wTarget*mean[source],minlength=len(area))
    gauge=intensive(s["rho"]*s["speedMps"]**2*state[wet,4])
    tw=s["airK"]+s["scalarUnitK"]*theta
    near=intensive(s["airK"]+s["scalarUnitK"]*state[first,9]);wall=intensive(tw)
    q=outW/area;p=s["airPa"]+gauge
    assert np.isfinite(np.stack([q,p,near,wall])).all() and np.min(p)>0
    regions_out=[]
    for region in range(1,6):
        mask=regions==region;power=float(-outW[mask].sum());src=float(-sourceW[group==region].sum())
        assert abs(power-src)<1e-10
        ids=np.flatnonzero(mask)
        regions_out.append(dict(id=region,faces=int(mask.sum()),areaM2=float(area[mask].sum()),outW=power,nativeOutW=src,
                minFluxIn=float(q[mask].min()),maxFluxIn=float(q[mask].max()),peakOutFace=int(ids[np.argmin(q[mask])]+1),
                minPa=float(p[mask].min()),maxPa=float(p[mask].max())))
    assert abs(sum(r["outW"] for r in regions_out)-thermal["early"]["wallW"])<1e-10
    stations=measure_stations(xyz,wet,np.unique(corners[group==1]),nodeQ,tw,s,case)
    np.testing.assert_array_equal(stations[:,0],[float(r["angleDeg"]) for r in reference])
    station_faces=np.array([int(r["faceID"]) for r in reference])
    assert np.all(regions[station_faces-1]==1)
    stations=np.column_stack([stations,station_faces,q[station_faces-1]])
    np.savetxt(out/"stations.csv",stations,delimiter=",",header="angleDeg,reportedNu,computedNu,deltaNu,relativeError,fluxInWm2,wallK,leftNode,rightNode,rightWeight,originalFaceID,faceMeanFluxInWm2",
               comments="",fmt=["%.17g"]*7+["%d","%d","%.17g","%d","%.17g"])
    np.savetxt(out/"wall.csv",np.column_stack([np.arange(1,len(area)+1),regions,area,p,gauge,q,wall,near]),
               delimiter=",",header="faceID,regionID,areaM2,absolutePa,gaugePa,fluxInWm2,wallInputK,firstK",
               comments="",fmt=["%d","%d"]+["%.17g"]*6)
    silver=regions_out[0];equiv=silver["outW"]*np.pi*s["diameterM"]**2/silver["areaM2"]
    observed=case["observedBtuS"]*1055.05585262
    summary=dict(thermal=thermal,condition=s,caseSHA256=sha(CASE),bodySHA256=grid["bodySHA256"],transferSHA256=grid["publicTransferSHA256"],
                 frameSHA256=frame,releaseSHA256=sha(work/"releaseEarly/restart.csv"),faces=len(area),areaM2=area.sum(),regions=regions_out,
                 minGaugePa=float(gauge.min()),maxGaugePa=float(gauge.max()),peakOutFace=int(q.argmin()+1),peakOutWm2=float(-q.min()),
                 solvedAreaM2=area.sum(),jointExperimentalAreaM2=0.,silverEquivalentFullSphereW=equiv,observedFullSphereW=observed,
                 equivalentRelativeError=equiv/observed-1,stationRMSE=float(np.sqrt(np.mean(stations[:,3]**2))),
                 stationMaxAbs=float(np.max(abs(stations[:,3]))),files={name:sha(out/name) for name in ("wall.csv","stations.csv")})
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False),encoding="utf8")
    with ZipFile(out/"surface.zip","x",compression=ZIP_DEFLATED,compresslevel=9) as z:
        for name in ("summary.json","wall.csv","stations.csv"):z.write(out/name,name)
    return summary


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("work",type=Path);parser.add_argument("out",type=Path)
    a=parser.parse_args();r=mapped(a.work.resolve(),a.out.resolve());print(r["faces"],r["areaM2"],r["thermal"]["early"]["wallW"])
