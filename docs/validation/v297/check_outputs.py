"""核对实际侧栏新算、归档与最终导出。以独立热流表为准逐面比较单位、编号与未知值。"""
import csv
import gzip
import io
import json
from pathlib import Path
import sys
import zipfile

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha


def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))


def check():
    area=ROOT/"outputs/v297_validation";records={}
    with zipfile.ZipFile(ROOT/"cases/thermistor_b/leads3.zip") as z:
        raw=z.read("wall.csv")
        wall=np.genfromtxt(io.BytesIO(raw),delimiter=",",names=True)
        labels=json.loads(z.read("body.json"))["regionNames"]
    flow=np.genfromtxt(io.BytesIO(gzip.decompress((HERE/"evidence/nearwall_faces.csv.gz").read_bytes())),delimiter=",",names=True)
    for folder in (area/"ui01",area/"ui02",area/"ui03"):
        rec=read(folder/"debugger_process.json");ui=read(folder/"ui.json")
        assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and not rec["stopReason"]
        assert not rec["survivingProcessIDs"] and not rec["supervisorError"] and all(rec["inputsStable"].values())
        assert ui["webglError"]==0 and ui["faces"]==len(wall) and len(ui["regionPeaks"])==7
        if ui["mode"]=="solve":
            job=Path(ui["runDirectory"]);actual=job/"surface/wall.csv"
            assert actual.read_bytes()==raw and ui["sourceFrame"]["sha256"]==sha(job/"restart.csv")
        page=Path(ui["page"]);datafile=page.parent/"case-data.js"
        content=datafile.read_text("utf8").strip()
        assert content.startswith("window.ROCKET_CASE=") and content.endswith(";")
        d=json.loads(content[len("window.ROCKET_CASE="):-1])
        f={field["key"]:field for field in d["fields"]}
        errors={}
        for key,column,scale in (("pressure","pressurePa",1000),("pressureGauge","pressureGaugePa",1),
                                 ("heatFlux","heatFluxWm2",1000),("nearWallCellTemperature","nearWallCellTemperatureK",1)):
            assert f[key]["status"]==[1]*len(wall)
            errors[key]=float(np.max(abs(np.array(f[key]["values"])*scale-wall[column])))
            assert errors[key]<1e-9
        for k in ("wallSolved","edgeTemperature","recoveryTemperature"):
            assert f[k]["values"]==[None]*len(wall) and f[k]["status"]==[0]*len(wall)
        assert f["wallInput"]["values"]==[305.15]*len(wall) and f["wallInput"]["status"]==[2]*len(wall)
        exported=rows(folder/"pairs.csv")
        assert len(exported)==len(wall)
        assert all(int(row["face_id"])==i+1 and row["region"]==labels[int(wall["regionID"][i])-1]
                   and row["joint_verified"]=="0" and row["pressure_method"] and row["heat_flux_method"]
                   and row["pressure_unit"]=="kPa" and row["heat_flux_unit"]=="kW/m²"
                   and abs(float(row["pressure"])*1000-wall["pressurePa"][i])<1e-9
                   and abs(float(row["heat_flux"])*1000-wall["heatFluxWm2"][i])<1e-9 for i,row in enumerate(exported))
        if "reverseFaces" in ui:
            assert ui["reverseFaces"]==int(flow["reverse_both"].sum())
            if folder.name=="ui03":assert ui["freshDiagnosticsCompared"]
            for k,col in (("leadReverse","reverse_both"),("leadDownstreamCf","downstream_Cf"),("leadFirstYPlus","estimated_first_yplus")):
                np.testing.assert_allclose(f[k]["values"],flow[col],rtol=0,atol=1e-14)
            exported=rows(folder/"reverse.csv")
            assert len(exported)==len(wall) and all(int(row["face_id"])==i+1 and int(row["value"])==flow["reverse_both"][i]
                                                      for i,row in enumerate(exported))
        records[folder.name]=dict(workbench=ui,seconds=rec["elapsedSeconds"],peakRSSBytes=rec["peakFamilyRSSBytes"],
                                  pageDataSHA256=sha(datafile),pairedCSVSHA256=sha(folder/"pairs.csv"),maxSIError=errors)
    with (HERE/"evidence/workbench.json").open("x",encoding="utf8") as f:json.dump(records,f,ensure_ascii=False,indent=2)


if __name__=="__main__":check()
