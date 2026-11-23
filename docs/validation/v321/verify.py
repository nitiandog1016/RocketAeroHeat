"""用原生边界面、反力、首层坐标及原论文表值，独立审核工作台和完整球面收支。"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[3];DATA=ROOT/"cases/short_sphere"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for part in iter(lambda:f.read(1024**2),b""):h.update(part)
    return h.hexdigest()
def close(a,b,tol):np.testing.assert_allclose(a,b,rtol=0,atol=tol)


def wall_faces(path):
    lines=iter(path.read_text("ascii").splitlines())
    for line in lines:
        if line.strip()=="MARKER_TAG= wall":
            count=int(next(lines).split("=")[1]);return [list(map(int,next(lines).split()))[1:] for _ in range(count)]
    raise ValueError("未找到球壁")


def check(work,ui):
    prep=work/"prepared";c=read(DATA/"test65.json");rec=read(prep/"prepared.json");s=rec["condition"]
    e=read(work/"endpoints.json");stage=work/e["stage"];d=read(ui.parent/"ui.json")
    assert d[ui.name]["frame"]==digest(stage/"restart.csv")==e["frame"] and d[ui.name]["glPixels"]>0
    assert digest(DATA/"test65.json")==e["case"] and all(digest(prep/k)==h for k,h in rec["files"].items())
    # 原生全程：检查初态接续、连续残差、所有周期帧及监督终态。
    src=digest(prep/"initial.csv");phases=[];weakRange=[float("inf"),-float("inf")]
    for dest in sorted(work.glob("steady[0-9][0-9]")):
        seal=read(dest/"native.json");assert src==seal["input"]==digest(dest/"initial_flow.csv")
        for name,h in seal["frames"].items():
            assert digest(dest/name)==h
            weak=np.loadtxt(dest/name,skiprows=1,delimiter=",",usecols=8)
            temp=s["airK"]+s["scalarUnitK"]*(weak+297.15);assert np.min(temp)>0
            weakRange=[min(weakRange[0],float(temp.min())),max(weakRange[1],float(temp.max()))]
        src=digest(dest/"restart.csv");assert src==seal["frame"];phases.append(dest)
    phases.extend(stage/n for n in ["reaction","plus20","reaction20","control","observer"])
    for dest in phases:
        seal=read(dest/"native.json");r=read(dest/"debugger_process.json")
        assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
        assert not(r["stopReason"] or r["survivingProcessIDs"] or r["supervisorError"])
        assert seal["frame"]==digest(dest/"restart.csv") and seal["history"]==digest(dest/"history.csv")
        if dest==stage or dest.name=="plus20":
            with (dest/"history.csv").open() as f:
                header=[x.strip() for x in next(csv.reader(f,skipinitialspace=True))];hist=np.loadtxt(f,delimiter=",",ndmin=2)
            assert len(hist)>=20
            for key,value in seal["window"].items():assert hist[-20:,header.index(key)].max()==value and value<-8
    state=np.loadtxt(stage/"restart.csv",delimiter=",",skiprows=1)
    pts=state[:,1:4]*s["diameterM"];wall=wall_faces(prep/"mesh.su2");wet=np.unique([i for face in wall for i in face])
    # 由原生面几何重建面积和节点分摊：四边形中心温度由四角点决定。
    lump=np.zeros(len(pts));pieces=[];corner_ids=[]
    for face in wall:
        xyz=pts[face];center=xyz.mean(axis=0)
        if len(face)==3:
            area=np.linalg.norm(np.cross(xyz[1]-xyz[0],xyz[2]-xyz[0]))/2
            lump[face]+=area/3;pieces.append(xyz.mean(axis=0));corner_ids.append((face,[]))
        else:
            ar=np.linalg.norm(np.cross(xyz-center,np.roll(xyz,-1,axis=0)-center),axis=1)/2
            lump[face]+=(ar+np.roll(ar,1))/3+ar.sum()/12
            for j in range(4):
                pair=[face[j],face[(j+1)%4]];pieces.append((pts[pair].sum(axis=0)+center)/3);corner_ids.append((pair,face))
    with (stage/"reaction/restart.csv").open() as f:
        keys=next(csv.reader(f));rhs=np.loadtxt(f,delimiter=",",usecols=keys.index("Residual_Species_0"))
    assert np.all(lump[wet]>0)
    scale=s["rho"]*s["speedMps"]*s["cp"]*s["diameterM"]**2*s["scalarUnitK"]
    vals=np.zeros((len(pts),4));vals[wet,0]=rhs[wet]*scale/lump[wet]
    vals[wet,1]=state[wet,4]*s["rho"]*s["speedMps"]**2
    vals[wet,2]=s["airK"]+s["scalarUnitK"]*state[wet,9]
    target=pts[wet]*(.502/.5);dist,nearest=cKDTree(pts).query(target)
    assert max(dist)<1e-14 and len(np.unique(nearest))==len(wet)
    vals[wet,3]=s["airK"]+s["scalarUnitK"]*state[nearest,9]
    expected=np.array([(vals[pair].sum(axis=0)+(vals[quad].mean(axis=0) if quad else 0))/3 for pair,quad in corner_ids])
    raw=(ui/"workbench/case-data.js").read_text("utf8").strip();page=json.loads(raw[len("window.ROCKET_CASE="):-1])
    tri=np.asarray(page["vertices"])[np.asarray(page["faces"])];dist,idx=cKDTree(pieces).query(tri.mean(axis=1));assert max(dist)<1e-14
    assert len(np.unique(idx))==len(idx)==4464
    expected=expected[idx];fields={f["key"]:f for f in page["fields"]}
    for j,key in enumerate(("heatFlux","pressureGauge","wallInput","nearWallCellTemperature")):
        close(fields[key]["values"],expected[:,j]/(1000 if j==0 else 1),1e-9)
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    with (ui/"wall.csv").open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
    assert len(rows)==len(expected)
    for i,r in enumerate(rows):
        assert int(r["face_id"])==i+1 and r["steady_frame_sha256"]==e["frame"] and r["experimental_doi"]==c["doi"]
        assert r["joint_verified"] in ("0","false") and r["pressure_method"] and r["heat_flux_method"]
        assert abs(float(r["area_m2"])-area[i])<1e-18 and abs(float(r["heat_flux"])*1000-expected[i,0])<1e-7
        assert abs(float(r["pressure"])*1000-(s["airPa"]+expected[i,1]))<1e-8
        assert abs(float(r["wall_input_K"])-expected[i,2])<1e-9 and abs(float(r["first_fluid_K"])-expected[i,3])<1e-9
    power=-expected[:,0]@area;close(power,-rhs[wet].sum()*scale,1e-11)
    # 按空间位置找 +Z 子午线上的壁节点，分别核对原始观测、温差、Nu 与误差。
    angles=np.rad2deg(np.arctan2(np.linalg.norm(pts[:,1:],axis=1),pts[:,0]))
    line=wet[(np.abs(pts[wet,1])<1e-15)&(pts[wet,2]>=-1e-15)];line=line[np.argsort(angles[line])];assert len(line)==49
    stations=np.asarray(c["angleDeg"]);q=np.interp(stations,angles[line],vals[line,0]);t=np.interp(stations,angles[line],vals[line,2])
    k=.0241*(t/273.15)**1.5*467.15/(t+194);nu=-q*.0127/(k*(t-(100.12-32)*5/9-273.15))
    with (ui/"stations.csv").open(encoding="utf-8-sig",newline="") as f:table=list(csv.DictReader(f))
    close([float(r["reportedNu"]) for r in table],c["localNuSurface"],1e-12);close([float(r["computedNu"]) for r in table],nu,1e-10)
    close([float(r["relativeError"]) for r in table],nu/np.asarray(c["localNuSurface"])-1,1e-11)
    assert all(r["frameSHA256"]==e["frame"] for r in table)
    x=.00137
    with (ui/"section.csv").open(encoding="utf-8-sig",newline="") as f:cut=list(csv.DictReader(f))
    cutids=np.array([int(r["face_id"])-1 for r in cut])
    np.testing.assert_array_equal(np.sort(cutids),np.flatnonzero((tri[:,:,0].min(1)<x)&(tri[:,:,0].max(1)>x)))
    close([float(r["value"]) for r in cut],expected[cutids,1],1e-9)
    return dict(faces=len(idx),sectionRows=len(cut),wholeW=float(power),analyticAreaM2=np.pi*.0127**2,gridAreaM2=float(area.sum()),
                computedNu=nu.tolist(),reportedNu=c["localNuSurface"],errorPercent=(100*(nu/c["localNuSurface"]-1)).tolist(),
                weakEnergyDiagnosticRangeK=weakRange,frame=e["frame"],source=d[ui.name],physicalSteps=0,
                stages=[dict(name=p.relative_to(work).as_posix(),**{k:read(p/"native.json")[k] for k in ("steps","seconds","peakRSS","window")}) for p in phases])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("ui",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();u=read(a.ui/"ui.json");terminal=read(a.ui/"debugger_process.json")
    assert terminal["status"]=="terminal" and terminal["debuggerExitCode"]==0 and all(terminal["inputsStable"].values())
    assert not(terminal["stopReason"] or terminal["survivingProcessIDs"] or terminal["supervisorError"])
    result={mode:check(ROOT/"outputs/v321_validation/run01" if mode=="archive" else Path(u[mode]["job"]),a.ui/mode) for mode in ("archive","solve")}
    close(result["archive"]["computedNu"],result["solve"]["computedNu"],5e-8)
    result["uiAndNewSolveSeconds"]=terminal["elapsedSeconds"];result["uiFamilyPeakBytes"]=terminal["peakFamilyRSSBytes"]
    with a.output.open("x",encoding="utf8") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
