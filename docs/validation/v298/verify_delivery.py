"""逐行核对两种工作台的实际导出，归集此次计算的小型原生记录。"""
import csv
import argparse
import io
import json
from pathlib import Path
import shutil
import zipfile

import numpy as np

from check_mesh import HERE,ROOT
from thermistorBead import read,sha


def audit(proof):
    output=ROOT/"outputs/v298_validation"
    job=output/"smooth01";comp=output/"comparison01"
    package=ROOT/"cases/thermistor_b/paired3.zip"
    assert sha(package)==sha(comp/"paired3.zip")==read(package.with_suffix(".json"))["archiveSHA256"]
    with zipfile.ZipFile(package) as z:
        def table(name):return np.genfromtxt(io.BytesIO(z.read(name)),names=True,delimiter=",")
        ref=table("wall.csv");shared=table("shared.csv");history=table("historical.csv")
        assert z.testzip() is None and z.read("wall.csv")==(job/"surface/wall.csv").read_bytes()
    with zipfile.ZipFile(ROOT/"cases/thermistor_b/leads3.zip") as z:
        wired=np.genfromtxt(io.BytesIO(z.read("wall.csv")),names=True,delimiter=",")
    with zipfile.ZipFile(ROOT/"cases/thermistor_b/surface_7p5.zip") as z:
        original=np.genfromtxt(io.BytesIO(z.read("wall.csv")),names=True,delimiter=",")
    vis={}
    with zipfile.ZipFile(proof/"receipts.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for path in job.rglob("*.json"):
            if path.name not in ("body.json",):z.write(path,"native/"+path.relative_to(job).as_posix())
        for path in (job/"case.cfg",job/"history.csv",job/"engine_stdout.log"):
            z.write(path,"native/"+path.name)
        for mode,folder in (("wired","wired_ui03"),("historical","smooth_ui03")):
            direct=output/folder;ui=read(direct/"ui.json");rec=read(direct/"debugger_process.json")
            assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and not rec["supervisorError"] and not rec["stopReason"]
            assert not rec["survivingProcessIDs"] and all(rec["inputsStable"].values()) and ui["webglError"]==0
            z.write(direct/"debugger_process.json",folder+"/debugger_process.json")
            script=direct/"workbench/case-data.js";raw=script.read_text("utf8").strip()
            assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
            viewer=json.loads(raw[len("window.ROCKET_CASE="):-1])
            fields={r["key"]:r for r in viewer["fields"]}
            if mode=="wired":
                i=shared["wiredID"].astype(int)-1;j=shared["smoothID"].astype(int)-1
                delta=np.full(len(wired),np.nan)
                delta[i]=wired["heatFluxWm2"][i]-ref["heatFluxWm2"][j]
            else:
                delta=history["pairedHeatFluxWm2"]-original["heatFluxWm2"]
                assert ui["changedTemperatureExcluded"]
            mask=np.isfinite(delta)
            indices=np.flatnonzero(mask)
            assert ui["buttonPeakFaceID"]==1+indices[np.argmax(abs(delta[mask]))]
            view=np.array(fields["pairedDeltaHeat"]["values"],float)
            np.testing.assert_allclose(view,delta,rtol=0,atol=1e-9,equal_nan=True)
            assert np.array_equal(np.array(fields["pairedDeltaHeat"]["status"])>0,mask)
            with (direct/"faces.csv").open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
            assert len(rows)==len(delta)
            for k,row in enumerate(rows):
                assert int(row["face_id"])==k+1 and row["field_key"]=="pairedDeltaHeat" and row["unit"]=="W/m²"
                if mask[k]:assert abs(float(row["value"])-delta[k])<1e-9
                else:assert row["value"]=="" and int(row["availability_code"])==0 and "光滑控制没有" in row["scope"]
            with (direct/"section.csv").open(encoding="utf-8-sig",newline="") as f:section=list(csv.DictReader(f))
            assert len(section)==ui["sectionRows"] and all(abs(float(r["value"])-delta[int(r["face_id"])-1])<1e-9 for r in section)
            ranges=viewer["comparisonColorRanges"]
            assert ranges["pairedSmoothHeat"]==ranges["heatFlux"] and ranges["pairedReverseChange"]==[-1,1]
            vis[mode]=dict(workbench=ui,seconds=rec["elapsedSeconds"],maxRSSBytes=rec["peakFamilyRSSBytes"],
                           pairedHeatColorRangeKWm2=ranges["heatFlux"],allRowsCompared=len(rows),
                           exportSHA256=sha(direct/"faces.csv"),pageDataSHA256=sha(script))
    for source,name in ((comp/"comparison.json","comparison.json"),(job/"surface/summary.json","surface.json")):
        with source.open("rb") as src,(proof/name).open("xb") as dst:shutil.copyfileobj(src,dst)
    with (proof/"workbench.json").open("x",encoding="utf8") as f:json.dump(vis,f,ensure_ascii=False,indent=2)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__);cli.add_argument("out",type=Path)
    dest=cli.parse_args().out.resolve();dest.relative_to(ROOT/"outputs/v298_validation")
    dest.mkdir(exist_ok=False);audit(dest)
