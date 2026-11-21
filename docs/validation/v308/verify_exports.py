"""独立读取实际 CSV 和三维页，复算七区放热，与 299.15 K 独立冷态逐行比较。"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def read(path):
    return json.loads(path.read_text(encoding="utf8"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def archive(name):
    path = ROOT/"cases/thermistor_b"/(name+".zip")
    assert sha(path)==read(path.with_suffix(".json"))["archiveSHA256"]
    with zipfile.ZipFile(path) as z:
        raw=z.read("wall.csv"); summary=json.loads(z.read("summary.json"))
        assert hashlib.sha256(raw).hexdigest()==summary["wallSHA256"]
        return np.genfromtxt(io.BytesIO(raw),delimiter=",",names=True), summary


def check_exports(folder,level):
    folder=folder.resolve()
    ui=read(folder/"ui.json");run=read(folder/"debugger_process.json")
    assert run["status"]=="terminal" and run["debuggerExitCode"]==0 and all(run["inputsStable"].values())
    assert not (run["survivingProcessIDs"] or run["supervisorError"] or run["stopReason"])
    assert ui["faces"]==10272 and ui["level"]==level and ui["webglError"]==0 and ui["visibleSamples"]>0
    raw=(folder/"workbench/case-data.js").read_text("utf8").strip()
    assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
    page=json.loads(raw[len("window.ROCKET_CASE="):-1]);fields={f["key"]:f for f in page["fields"]}
    assert page["resultCondition"]["wallUniformValueK"]==299.15 and page["resultCondition"]["wallInputType"]=="uniform"
    native,ref=archive(("leadsMidbody3","leadsThird3")[level-1])
    pair=rows(folder/"wall.csv");first=rows(folder/"first.csv");section=rows(folder/"section.csv")
    assert read(folder/"first.json")["thermalConversion"]==ui["thermal"]==page["evidenceIdentity"]["sourceFrame"]["thermalConversion"]
    assert len(pair)==len(first)==len(native)==10272 and ui["source"]["sha256"]==ref["sourceFrameSHA256"]
    assert np.all(native["faceID"]==np.arange(len(native)))
    for i,(p,t) in enumerate(zip(pair,first),1):
        assert int(p["face_id"])==int(t["face_id"])==i and p["region"]==t["region"]==page["region"][i-1]
        assert float(p["reference_wall_K"])==305.15 and float(p["temperature_difference_ratio"])==.5
        assert p["reference_frame_sha256"]==ref["sourceFrameSHA256"] and "299.15 K" in p["condition"]
        assert p["joint_verified"]=="0" and p["pressure_availability"]==p["heat_flux_availability"]=="research"
        assert p["pressure_unit"]=="kPa" and p["heat_flux_unit"]=="kW/m²" and "参照温差" in p["heat_flux_method"]
        assert t["unit"]=="K" and t["field_key"]=="nearWallCellTemperature" and t["availability_code"]=="2"
        assert "299.15 K" in t["scope"] and ref["sourceFrameSHA256"] in t["scope"]
    p=np.array([float(r["pressure"])*1000 for r in pair])
    q=np.array([float(r["heat_flux"])*1000 for r in pair])
    first_k=np.array([float(r["value"]) for r in first])
    area=np.array([float(r["area_m2"]) for r in pair])
    def close(actual,expected,tol):
        np.testing.assert_allclose(actual,expected,rtol=0,atol=tol)
    close(area,native["areaM2"],1e-18);close(area.sum(),3.19665061256752e-6,1e-17)
    close(p,native["pressurePa"],1e-9);close(q,.5*native["heatFluxWm2"],1e-9)
    close(first_k,293.15+.5*(native["nearWallCellTemperatureK"]-293.15),1e-10)
    for key,vals,status in (("pressure",p/1000,1),("heatFlux",q/1000,1),("nearWallCellTemperature",first_k,1),("wallInput",299.15,2)):
        close(fields[key]["values"],vals,1e-10);assert all(code==status for code in fields[key]["status"])
    assert not any(fields["wallSolved"]["status"]) and section
    for r in section:
        assert float(r["section_x_m"])==.000172 and r["field_key"]=="heatFlux" and "299.15 K" in r["scope"]
        assert abs(float(r["value"])*1000-q[int(r["face_id"])-1])<1e-9
    groups=native["regionID"].astype(int)
    region_w=[float(-area[groups==j]@q[groups==j]) for j in range(1,8)]
    close(region_w,ui["regionW"],1e-14);close(sum(region_w),ui["totalOutwardW"],1e-14)
    old,_=archive("leads3")
    if level==1:
        other,_=archive("outer30");delta=.5*other["heatFluxWm2"]-q
        close(fields["outerDeltaHeat"]["values"],delta,1e-9)
        close(fields["outerDeltaFirst"]["values"],293.15+.5*(other["nearWallCellTemperatureK"]-293.15)-first_k,1e-10)
        close(fields["leadGridDeltaHeat"]["values"],q-.5*old["heatFluxWm2"],1e-9)
    else:
        medium,_=archive("leadsMidbody3")
        close(fields["gridHeat32"]["values"],q-.5*medium["heatFluxWm2"],1e-9)
        close(fields["gridHeat21"]["values"],.5*(medium["heatFluxWm2"]-old["heatFluxWm2"]),1e-9)
        assert all(g["wallK"]==299.15 for g in page["midbodyGridComparison"]["grids"])
    assert "bandRelease" not in fields and int(np.argmin(q))+1==ui["peakFaceID"]
    stats=dict(totalOutwardW=sum(region_w),regions=region_w,peakFaceID=ui["peakFaceID"],pressurePa=[min(p),max(p)],
               intoWallWm2=[min(q),max(q)],firstLayerK=[min(first_k),max(first_k)],peakFamilyRSSBytes=run["peakFamilyRSSBytes"],
               seconds=run["elapsedSeconds"],cutRows=len(section),pageDirectory=folder.relative_to(ROOT).as_posix(),
               conversion=ui["thermal"],powerSolvedWallK=ui["powerSolvedWallK"],otherWallsK=ui["otherWallsK"],
               jointExperimentVerifiedM2=0,source=ui["source"]["sha256"],visibleSamples=ui["visibleSamples"],
               files={k:sha(folder/k) for k in ("ui.json","debugger_process.json","first.csv","first.json","wall.csv","section.csv","workbench/case-data.js")})
    return stats,area,groups,p,q,first_k


def verify(medium,fine,out):
    m,area,groups,p,q,t=check_exports(medium,1)
    f,a2,g2,p2,q2,t2=check_exports(fine,2)
    np.testing.assert_array_equal(groups,g2);np.testing.assert_array_equal(area,a2)
    half=read(HERE.parent/"v302/evidence/half.json");root=ROOT/half["directory"]
    terminal=read(root/"engine_process.json");direct=ROOT/"outputs/v302_validation/native01/direct_wall.csv"
    assert terminal["status"]=="terminal" and terminal["exitCode"]==0 and not (terminal["timedOut"] or terminal["cancelled"])
    assert all(terminal["inputsStable"].values()) and terminal["progress"]["Inner_Iter"]==half["lastIteration"]
    assert sha(root/"restart.csv")==terminal["outputSHA256"]["restart.csv"]==half["energy"]["inputFrameSHA256"]
    assert sha(direct)==half["wallSHA256"] and read(root/"case.json")["wallK"]==299.15
    assert max(half["residualWindow"].values())<=-8 and abs(half["energy"]["relativeDifference"])<.005 and half["relativeHeatChange20"]<.002
    independent=np.genfromtxt(direct,names=True,delimiter=",")
    np.testing.assert_array_equal(independent["faceID"],np.arange(len(q)))
    np.testing.assert_array_equal(independent["regionID"],groups)
    dq=q-independent["heatFluxWm2"]
    err=dict(maxPressurePa=max(abs(p-independent["pressurePa"])),maxFluxWm2=max(abs(dq)),
             maxRelativeFlux=max(abs(dq/independent["heatFluxWm2"])),maxFirstLayerK=max(abs(t-independent["firstLayerK"])),
             relativePower=abs(m["totalOutwardW"]/half["wallReleaseW"]-1))
    assert err["maxPressurePa"]<.001 and err["maxFirstLayerK"]<.001 and max(err["maxRelativeFlux"],err["relativePower"])<.0005
    out.mkdir(exist_ok=False,parents=True)
    # 逐面表按工作台编号，保存实际导出读数和独立冷态，所有热流正号入壁。
    table=io.StringIO(newline="");writer=csv.writer(table)
    writer.writerow(["face_id","region_id","area_m2","given_wall_K","medium_p_Pa","medium_q_W_m2","medium_first_K",
                     "fine_p_Pa","fine_q_W_m2","fine_first_K","direct_medium_p_Pa","direct_medium_q_W_m2","direct_medium_first_K"])
    for i in range(len(q)):
        writer.writerow([i+1,groups[i],area[i],299.15,p[i],q[i],t[i],p2[i],q2[i],t2[i],
                         independent["pressurePa"][i],independent["heatFluxWm2"][i],independent["firstLayerK"][i]])
    with zipfile.ZipFile(out/"wall299K.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("faces.csv",table.getvalue().encode("utf-8-sig"))
    result=dict(faces=len(q),areaM2=float(area.sum()),regionNames=[r["name"] for r in half["regions"]],
                medium=m,fine=f,independentMedium=dict(frame=half["energy"]["inputFrameSHA256"],table=sha(direct),errors=err,
                powerW=half["wallReleaseW"],lastIteration=half["lastIteration"]),tableSHA256=sha(out/"wall299K.zip"))
    (out/"verification.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf8")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("medium","fine","output"):parser.add_argument(key,type=Path)
    args=parser.parse_args()
    verify(args.medium,args.fine,args.output)
