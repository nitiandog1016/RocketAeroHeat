"""核对侧栏、三维页、整表及截面与三档网独立原面结果的物理值及来源。"""
import argparse
import csv
import io
import json
from pathlib import Path
import zipfile

import numpy as np

from audit_third_mesh import FINE, ROOT, read, sha, store


def equal(left, right, tol=1e-9):
    np.testing.assert_allclose(left, right, rtol=0, atol=tol)


def verify(ui, reference, out):
    rec=read(ui/"debugger_process.json");proof=read(ui/"ui.json");source=read(reference/"comparison.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"]) and proof["zeroPowerAndInvalidation"]
    assert proof["secondGridFrame"]==source["grids"][1]["frameSHA256"]
    report={}
    with zipfile.ZipFile(reference/"face_comparison.zip") as z:
        assert z.testzip() is None and json.loads(z.read("sources.json"))==source
        for mode in ("givenWall","givenPower"):
            folder=ui/mode;p=proof[mode];comp=source[mode]
            assert p["webglError"]==0 and p["visibleSamples"]>0 and p["faces"]==source["faces"]==10272
            t=np.genfromtxt(io.BytesIO(z.read(mode+".csv")),delimiter=",",names=True)
            raw=(folder/"workbench/case-data.js").read_text("utf8").strip()
            assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
            data=json.loads(raw[len("window.ROCKET_CASE="):-1]);fields={f["key"]:f for f in data["fields"]}
            side=read(folder/"faces.json");info=side["midbodyGridComparison"]
            assert info==p["summary"]==data["midbodyGridComparison"] and info["available"] and info["mode"]==mode
            assert (info["boundaryValue"],info["boundaryUnit"])==((305.15,"K") if mode=="givenWall" else (.004,"W"))
            assert side["rows"]==10272 and side["csvSHA256"]==sha(folder/"faces.csv")
            assert data["evidenceIdentity"]["sourceFrame"]["midbodyGrids"]==info["grids"]
            for j,(r,s) in enumerate(zip(info["grids"],source["grids"])):
                assert all(r[k]==s[k] for k in ("frameSHA256","meshSHA256","reactionSHA256","wallSHA256"))
                equal(r["conductanceWK"],s["conductanceWK"],1e-14);equal(r["wallK"],comp["wallK"][j],1e-12)
                equal(r["totalOutwardW"],comp["totalOutwardW"][j],1e-13)
            assert len(info["grids"])==3 and len(info["regions"])==7
            for r,s in zip(info["regions"],comp["regions"]):
                assert (r["name"],r["faces"])==(s["name"],s["faces"])
                equal(r["outwardW"],s["outwardW"],1e-13)
            q1,q2,q3=(t[f"q{j}_into_W_m2"] for j in range(1,4));delta=q3-q2
            for key,values in (("gridHeat1",q1/1000),("gridHeat2",q2/1000),("heatFlux",q3/1000),
                               ("gridHeat21",q2-q1),("gridHeat32",delta),("pressure",t["p3_Pa"]/1000),
                               ("nearWallCellTemperature",t["first3_K"])):
                equal(fields[key]["values"],values);assert all(v==1 for v in fields[key]["status"])
            tw="wallInput" if mode=="givenWall" else "wallSolved"
            equal(fields[tw]["values"],t["Tw3_K"])
            assert set(fields[tw]["status"])==({2} if mode=="givenWall" else {1})
            assert set(fields["wallSolved" if mode=="givenWall" else "wallInput"]["status"])=={0}
            assert p["probeFace"]==int(np.argmax(abs(delta)))+1
            for name in ("faces","section"):
                with (folder/(name+".csv")).open(encoding="utf-8-sig",newline="") as stream:
                    rows=list(csv.DictReader(stream))
                assert len(rows)==(10272 if name=="faces" else p["cutRows"])
                for j,row in enumerate(rows):
                    i=int(row["face_id"])-1
                    assert 0<=i<len(t) and row["field_key"]=="gridHeat32" and row["unit"]=="W/m²"
                    assert abs(float(row["value"])-delta[i])<1e-9
                    if name=="faces":assert i==j and abs(float(row["area_m2"])-t["area_m2"][i])<1e-18
                    else:assert float(row["section_x_m"])==p["sectionX"]
            report[mode]=dict(ui=p,facesSHA256=sha(folder/"faces.csv"),sectionSHA256=sha(folder/"section.csv"),
                              sidecarSHA256=sha(folder/"faces.json"),pageDataSHA256=sha(folder/"workbench/case-data.js"))
    out.mkdir(parents=True,exist_ok=False)
    store(out/"workbench.json",dict(modes=report,zeroPowerAndInvalidation=True,elapsedSeconds=rec["elapsedSeconds"],
                                   peakFamilyRSSBytes=rec["peakFamilyRSSBytes"],faceComparisonSHA256=sha(reference/"face_comparison.zip")))
    with zipfile.ZipFile(out/"receipts.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for path in FINE.rglob("*.json"):
            if path.name!="body.json":z.write(path,"native/"+path.relative_to(FINE).as_posix())
        for name in ("case.cfg","history.csv","engine_stdout.log"):z.write(FINE/name,"native/"+name)
        for name in ("debugger_process.json","ui.json", "givenWall/faces.json","givenPower/faces.json"):
            z.write(ui/name,"ui/"+name)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("ui","reference","output"):p.add_argument(name,type=Path)
    a=p.parse_args();paths=[v.resolve() for v in (a.ui,a.reference,a.output)]
    for path in paths:path.relative_to(ROOT/"outputs")
    verify(*paths)
