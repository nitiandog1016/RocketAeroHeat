"""按发布的 12 K 原面表独立核查 4 mW 的离线云图、CSV、截面和实验资格。"""
import csv
import argparse
import io
import json
import shutil
import zipfile

import numpy as np

from half_case import HERE,ROOT,read,sha
from audit_half import save


def csv_rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))


def main(output_name):
    root=ROOT/"outputs/v302_validation";ui=root/"wired_ui01";ref_ui=root/"smooth_ui01";native=root/"native01"
    proof=read(ui/"ui.json");regression=read(ref_ui/"result.json")
    assert sha(native/"direct_wall.csv")==read(native/"half.json")["wallSHA256"]
    assert regression["boundaryInvalidationChecked"] and regression["workbench"]["faces"]==5216
    runs={}
    for folder in (ui,ref_ui):
        receipt=read(folder/"debugger_process.json")
        assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and all(receipt["inputsStable"].values())
        assert not (receipt["supervisorError"] or receipt["stopReason"] or receipt["survivingProcessIDs"])
        runs[folder.name]=dict(seconds=receipt["elapsedSeconds"],peakFamilyRSSBytes=receipt["peakFamilyRSSBytes"])
    archive=ROOT/"cases/thermistor_b/leadsMidbody3.zip";assert sha(archive)==read(archive.with_suffix(".json"))["archiveSHA256"]
    with zipfile.ZipFile(archive) as z:
        base=np.genfromtxt(io.BytesIO(z.read("wall.csv")),names=True,delimiter=",");summary=json.loads(z.read("summary.json"))
    q0=.004;t=proof["thermal"];factor=q0/summary["totalOutwardW"];wall=293.15+12*factor
    assert t["netHeatingPowerW"]==q0 and t["scale"]==factor and t["wallTemperatureK"]==wall
    assert t["nativeFrameSHA256"]==summary["sourceFrameSHA256"] and proof["faces"]==len(base)==10272
    q=factor*base["heatFluxWm2"];near=293.15+factor*(base["nearWallCellTemperatureK"]-293.15)
    exported=csv_rows(ui/"wall.csv");temperatures=csv_rows(ui/"temperature.csv")
    assert len(exported)==len(temperatures)==len(q)
    for j,(p,tr) in enumerate(zip(exported,temperatures)):
        assert int(p["face_id"])==int(tr["face_id"])==j+1 and abs(float(p["area_m2"])-base["areaM2"][j])<1e-18
        assert p["pressure_unit"]=="kPa" and p["heat_flux_unit"]=="kW/m²" and p["joint_verified"]=="0"
        assert p["condition_sha256"]==proof["conditionSHA256"] and abs(float(p["pressure"])*1000-base["pressurePa"][j])<1e-9
        assert abs(float(p["heat_flux"])*1000-q[j])<1e-9 and float(tr["value"])==wall
        assert tr["field_key"]=="wallSolved" and tr["availability"]=="研究计算" and tr["qualification"]=="未验证"
    raw=(ui/"workbench/case-data.js").read_text("utf8").strip()
    assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
    data=json.loads(raw[len("window.ROCKET_CASE="):-1]);fields={f["key"]:f for f in data["fields"]}
    assert not any(fields["wallInput"]["status"]) and all(v is None for v in fields["wallInput"]["values"])
    assert data["resultCondition"]["wallInputType"]=="netHeatingPower" and data["resultCondition"]["solvedWallTemperatureK"]==wall
    for k,expected in (("heatFlux",q/1000),("pressure",base["pressurePa"]/1000),("nearWallCellTemperature",near),("wallSolved",np.full(len(q),wall))):
        np.testing.assert_allclose(fields[k]["values"],expected,rtol=0,atol=1e-10)
        assert all(code==1 for code in fields[k]["status"])
    cuts={}
    for key in ("heatFlux","wallSolved"):
        cut=csv_rows(ui/(key+"_section.csv"))
        assert cut
        for row in cut:
            assert row["field_key"]==key and float(row["section_x_m"])==proof["sectionX"]
            assert abs(float(row["value"])-fields[key]["values"][int(row["face_id"])-1])<1e-10
        cuts[key]=len(cut)
    regions=base["regionID"].astype(int)
    region_w=[float(-base["areaM2"][regions==j]@q[regions==j]) for j in range(1,8)]
    np.testing.assert_allclose(region_w,proof["regionW"],atol=1e-15,rtol=0)
    assert abs(sum(region_w)-q0)<=1e-13 and proof["webglError"]==0
    out=root/output_name;out.mkdir(exist_ok=False)
    # 精简原面表直接取自已核对的工作台数值；原始大场单独留在 E 盘。
    field_table=io.StringIO(newline="")
    writer=csv.writer(field_table)
    writer.writerow(("face_id","region","area_m2","wall_pressure_Pa","into_wall_heat_flux_W_m2","first_layer_fluid_K","solved_isothermal_wall_K"))
    for j,name in enumerate(data["region"]):
        writer.writerow((j+1,name,data["faceArea"][j],1000*fields["pressure"]["values"][j],1000*fields["heatFlux"]["values"][j],
                         fields["nearWallCellTemperature"]["values"][j],fields["wallSolved"]["values"][j]))
    with zipfile.ZipFile(out/"power4mW.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("faces.csv",field_table.getvalue().encode("utf-8-sig"))
        info=dict(condition=data["resultCondition"],thermal=t,source=data["evidenceIdentity"],qualifiedJointAreaM2=0,
                  fields={k:fields[k]["methodNames"] for k in ("pressure","heatFlux","wallSolved","nearWallCellTemperature")})
        z.writestr("source.json",json.dumps(info,ensure_ascii=False,indent=2).encode("utf8"))
    save(out/"workbench.json",dict(ui=proof,cutRows=cuts,runs=runs,allRowsChecked=len(q),smoothRegression=regression["workbench"],
        pressurePa=[float(base["pressurePa"].min()),float(base["pressurePa"].max())],intoWallWm2=[float(q.min()),float(q.max())],
        firstLayerK=[float(near.min()),float(near.max())],peakFaceID=int(np.argmin(q))+1,
        beadAndJointsW=sum(region_w[2:5]),leadSidesAndEndsW=sum(region_w[:2]+region_w[5:]),
        pageDataSHA256=sha(ui/"workbench/case-data.js"),tableArchiveSHA256=sha(out/"power4mW.zip"),
        exports={n:sha(ui/n) for n in ("wall.csv","temperature.csv","heatFlux_section.csv","wallSolved_section.csv")}))
    shutil.copyfile(native/"receipts.zip",out/"receipts.zip")
    with zipfile.ZipFile(out/"receipts.zip","a",compression=zipfile.ZIP_DEFLATED) as z:
        for folder in (ui,ref_ui):z.write(folder/"debugger_process.json",folder.name+"/debugger_process.json")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output_name",nargs="?",default="delivery01");a=p.parse_args()
    if not a.output_name.isidentifier():p.error("请用一个新的目录名")
    main(a.output_name)
