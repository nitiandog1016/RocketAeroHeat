"""从发布交叠、原生释放反力和第一层坐标独立复算全部面及科研出口。"""
import argparse
import csv
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/"v328"))
from verify import terminal,table
from run import sha, read
from thermistorWireSurface import residual
from shortSphere65 import conductivity, CASE


def check(work,ui):
    data=ROOT/"cases/short_sphere"
    with ZipFile(data/"supportedSurface65.zip") as z:
        surface=np.genfromtxt(io.BytesIO(z.read("wall.csv")),names=True,delimiter=",")
        s=json.loads(z.read("summary.json"));stations=np.genfromtxt(io.BytesIO(z.read("stations.csv")),names=True,delimiter=",")
    with ZipFile(data/"supported65.zip") as z:g=json.loads(z.read("geometry.json"))
    t=s["thermal"];prod=read(work/"converged.json");start=sha(work/"steady01/initial_flow.csv")
    iterations=0;pairsChecked=0
    for idx in range(len(prod["trend"])):
        d=work/f"steady{idx+1:02d}";terminal(d);rec=read(d/"record.json")
        assert rec["start"]==start and sha(d/"initial_flow.csv")==start and sha(d/"history.csv")==rec["history"]
        for name,h in rec["frames"].items():assert sha(d/name)==h
        start=sha(d/"restart.csv");assert start==rec["end"]
        history=table(d,"history.csv");pairs=np.loadtxt(d/"far_pairs.csv",delimiter=",",ndmin=2)
        trace=np.genfromtxt(d/"stages.csv",delimiter=",",names=True,dtype=None,encoding="ascii")
        np.testing.assert_array_equal(pairs[:,0],np.arange(len(history)))
        assert np.all(pairs[:,1]==1346) and np.all(pairs[:,4]<=1e-15)
        assert len(trace)==400 and np.all(trace["below"]==0) and np.all(trace["above"]==0)
        iterations+=len(history);pairsChecked+=int(pairs[:,1].sum())
    assert iterations==t["updates"] and max(history[-20:,3:9].max(axis=0))<-8 and max(trace["rawRMS"][1::2][-20:])<1e-8
    assert start==prod["frame"]==s["frameSHA256"]==sha(work/"early/initial_flow.csv")
    assert t["passed"] and t["thermalChange"]<.002 and all(abs(t[k]["physicalGapFraction"])<.005 for k in ("early","late"))
    saved=table(work/"early","initial_flow.csv");cond=s["condition"]
    scale=cond["rho"]*cond["cp"]*cond["speedMps"]*cond["diameterM"]**2*cond["scalarUnitK"]
    state=np.fromfile(work/"early/nodes.bin",dtype="<f8").reshape(-1,16);state=state[state[:,0].argsort()]
    edge=np.fromfile(work/"early/edges.bin",dtype="<f8").reshape(-1,6)
    # 直接重建有向内部边对两端的反力。独立释放读取壁行，额外逐边修正只用同一快照。
    i,j=edge[:,:2].astype(int).T
    adjustment=edge[:,4]*(edge[:,5]-1)
    correction=np.bincount(i,adjustment,minlength=len(saved))-np.bincount(j,adjustment,minlength=len(saved))
    release=residual(work/"releaseEarly/restart.csv")
    with np.load(data/"supportedTransfer65.npz") as transfer:
        required=("target","source","targetFraction","sourceFraction","targetArea","region","nativeRegion","nativeCorners",
                  "nativeArea","nativeWallXYZD","firstXYZD")
        target,source,wt,ws,area,group,srcGroup,corners,srcArea,wallXYZ,firstXYZ=[transfer[k] for k in required]
    tree=cKDTree(saved[:,1:4]);dw,wall=tree.query(wallXYZ);df,first=tree.query(firstXYZ)
    assert max(dw.max(),df.max())<1e-12 and len(np.unique(wall))==len(wall)==5954 and not np.intersect1d(wall,first).size
    load=scale*(release[wall]-correction[wall])
    patchIDs=[];vertexIDs=[];pieceAreas=[]
    for k,v in enumerate(corners):
        ids=list(dict.fromkeys(v.tolist()));tri=wallXYZ[ids];center=tri.mean(axis=0)
        vectors=np.cross(tri-center,np.roll(tri,-1,axis=0)-center)
        normal=vectors.sum(axis=0);wedges=vectors@normal/(normal@normal)
        shares=.5*(wedges+np.roll(wedges,1))*srcArea[k]
        assert min(shares)>0
        patchIDs.extend([k]*len(ids));vertexIDs.extend(ids);pieceAreas.extend(shares)
    B=coo_matrix((pieceAreas,(patchIDs,vertexIDs)),shape=(len(srcArea),len(wall))).tocsr()
    lump=np.asarray(B.sum(axis=0))[0];np.testing.assert_allclose(np.asarray(B.sum(axis=1))[:,0],srcArea,atol=1e-19)
    watts=B@(load/lump);q=(coo_matrix((ws,(target,source)),shape=(len(area),len(watts))).tocsr()@watts)/area
    np.testing.assert_array_equal(group[target],srcGroup[source])
    I=coo_matrix((wt,(target,source)),shape=(len(area),len(watts))).tocsr()
    def field(values):return I@((B@values)/srcArea)
    gauge=field(cond["rho"]*cond["speedMps"]**2*saved[wall,4]);p=cond["airPa"]+gauge
    tw=field(cond["airK"]+cond["scalarUnitK"]*saved[wall,9]);near=field(cond["airK"]+cond["scalarUnitK"]*saved[first,9])
    expect={"absolutePa":p,"gaugePa":gauge,"wallInputK":tw,"firstK":near,"fluxInWm2":q}
    maxErrors={k:float(np.max(abs(surface[k]-v))) for k,v in expect.items()}
    assert max(maxErrors.values())<1e-8
    np.testing.assert_array_equal(surface["faceID"],np.arange(len(area))+1)
    np.testing.assert_array_equal(surface["regionID"],group)
    np.testing.assert_allclose(np.bincount(group,-q*area)[1:],np.bincount(srcGroup,-watts)[1:],atol=1e-11,rtol=0)
    assert abs(q@area+t["early"]["wallW"])<1e-10 and q.argmin()+1==s["peakOutFace"]
    # 在第二个端点独立映射所有原面的压力与热流，避免只用总功率掩盖局部漂移。
    lateEdges=np.fromfile(work/"late/edges.bin",dtype="<f8").reshape(-1,6)
    le=lateEdges[:,4]*(lateEdges[:,5]-1);li,lj=lateEdges[:,:2].astype(int).T
    lc=np.bincount(li,le,minlength=len(saved))-np.bincount(lj,le,minlength=len(saved))
    lateLoad=scale*(residual(work/"releaseLate/restart.csv")[wall]-lc[wall])
    latePatch=B@(lateLoad/lump)
    lateQ=(coo_matrix((ws,(target,source)),shape=(len(area),len(watts))).tocsr()@latePatch)/area
    lateState=table(work/"late","initial_flow.csv")
    lateGauge=field(cond["rho"]*cond["speedMps"]**2*lateState[wall,4])
    assert abs(lateQ@area+t["late"]["wallW"])<1e-10
    # 以实际测线附近两个原生节点核对角度、弦上位置及测站热流。
    case=read(CASE);np.testing.assert_array_equal(stations["reportedNu"],case["localNuSurface"])
    for station in stations:
        ids=np.searchsorted(wall,[int(station["leftNode"]),int(station["rightNode"])]);a=station["rightWeight"]
        xyz=(1-a)*wallXYZ[ids[0]]+a*wallXYZ[ids[1]]
        assert abs(xyz[1])<1e-12 and abs(np.rad2deg(np.arctan2(xyz[2],xyz[0]))-station["angleDeg"])<1e-10
        qst=(1-a)*load[ids[0]]/lump[ids[0]]+a*load[ids[1]]/lump[ids[1]]
        nu=-qst*cond["diameterM"]/((station["wallK"]-cond["airK"])*conductivity(station["wallK"]))
        assert abs(nu-station["computedNu"])<1e-10
    # 检查 WebGL 与 CSV 的每一行、未解温度、研究资格以及截面所指原面。
    page=json.loads((ui/"workbench/case-data.js").read_text("utf8").split("=",1)[1].strip().rstrip(";"))
    np.testing.assert_array_equal(np.asarray(page["faces"])+1,g["F"])
    fields={f["key"]:f for f in page["fields"]}
    for key,values in (("heatFlux",q/1000),("pressure",p/1000),("pressureGauge",gauge),("wallInput",tw),("nearWallCellTemperature",near)):
        f=fields[key];np.testing.assert_allclose(f["values"],values,atol=1e-9,rtol=0);assert f["methodNames"]
    for key in ("wallSolved","edgeTemperature","recoveryTemperature"):
        assert all(v is None for v in fields[key]["values"]) and not any(fields[key]["status"])
    with (ui/"wall.csv").open(encoding="utf-8-sig",newline="") as f:
        count=0
        for count,row in enumerate(csv.DictReader(f),1):
            k=count-1;assert int(row["face_id"])==count and float(row["area_m2"])>0 and row["joint_verified"] in ("0","false")
            assert abs(float(row["pressure"])*1000-p[k])<1e-8 and abs(float(row["heat_flux"])*1000-q[k])<1e-8
            assert row["pressure_method"] and row["heat_flux_method"]
        assert count==len(area)
    with (ui/"section.csv").open(encoding="utf-8-sig",newline="") as f:cut=list(csv.DictReader(f))
    assert cut and all(abs(float(r["value"])-gauge[int(r["face_id"])-1])<1e-8 for r in cut)
    return dict(faces=len(q),nodes=len(state),edges=len(edge),pairsChecked=pairsChecked,updates=iterations,sectionRows=len(cut),
                maxIndependentAbsoluteErrors=maxErrors,summedOutW=-q@area,peakOutFace=int(q.argmin()+1),webFields=list(fields),
                plus20MaxPressureChangePa=float(np.max(abs(lateGauge-gauge))),plus20MaxHeatFluxChangeWm2=float(np.max(abs(lateQ-q))),
                plus20PeakOutFace=int(lateQ.argmin()+1))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("run",type=Path);p.add_argument("ui",type=Path);p.add_argument("output",type=Path)
    opts=p.parse_args();result=check(opts.run.resolve(),opts.ui.resolve())
    with opts.output.open("x",encoding="utf8",newline="\n") as stream:json.dump(result,stream,indent=2);stream.write("\n")
    print(result)
