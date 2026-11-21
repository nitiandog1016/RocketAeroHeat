"""独立按原面编号检查工作台导出，归集可审查的小型分带证据。"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import zipfile

import numpy as np

from common_layer import BASE, GRIDS, HERE, ROOT, read, save, sha


def verify(ui,budget,geom,build,samples,out):
    proof=read(ui/"ui.json");rec=read(ui/"debugger_process.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not(rec["survivingProcessIDs"] or rec["supervisorError"] or rec["stopReason"])
    assert proof["powerCaseSeparated"] and proof["samples"]>0 and proof["webglError"]==0
    file=ROOT/"cases/thermistor_b/layerBands.json";pkg=read(file);native=read(budget/"budget.json")
    assert read(budget/"layerBands.json")==pkg and native["layerBandsSHA256"]==sha(budget/"layerBands.json")==sha(file)
    assert pkg["geometryAuditSHA256"]==sha(geom/"geometry.json") and pkg["bands"]==native["bands"]
    for s,g,path in zip(pkg["sourceGrids"],native["grids"],samples):
        assert sha(path/"sample.json")==s["sampleSHA256"] and sha(path/"edges.csv")==s["edgeFileSHA256"]
        assert all(s[k]==g[k] for k in s)
    raw=(ui/"workbench/case-data.js").read_text("utf8").strip()
    assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
    page=json.loads(raw[len("window.ROCKET_CASE="):-1]);side=read(ui/"faces.json")
    details={k:v for k,v in pkg.items() if k!="faceBand"};details["packageSHA256"]=sha(file)
    assert page["layerBudget"]==side["layerBudget"]==details
    assert page["evidenceIdentity"]["sourceFrame"]["layerBudget"]=={"packageSHA256":sha(file),"sourceGrids":pkg["sourceGrids"]}
    fields={f["key"]:f for f in page["fields"]}
    with np.load(geom/"grid3.npz") as mapping:
        ids=mapping["faceBand"]
    with np.load(GRIDS[1]/"nodes.npz") as mapping:
        area=mapping["area"]
    np.testing.assert_array_equal(ids+1,pkg["faceBand"])
    bands=native["bands"];assert len(ids)==side["rows"]==proof["faces"]==10272 and len(bands)==proof["bands"]==10
    q=[]
    for folder in GRIDS:
        wall=np.genfromtxt(folder/"surface/wall.csv",names=True,delimiter=",")
        q.append(wall["heatFluxWm2"])
    for key,target in (("heatFlux",q[1]/1000),("pressure",wall["pressurePa"]/1000),
                       ("nearWallCellTemperature",wall["nearWallCellTemperatureK"]),("gridHeat32",q[1]-q[0])):
        np.testing.assert_allclose(fields[key]["values"],target,rtol=0,atol=1e-9)
        assert set(fields[key]["status"])=={1}
    direct=np.bincount(ids,area*(q[0]-q[1]),minlength=len(bands))
    np.testing.assert_allclose(direct,[b["delta"]["faceOutW"] for b in bands],rtol=0,atol=1e-13)
    assert ids[proof["faceSelectedByButton"]-1]==np.argmax(abs(direct))
    targets={"bandRelease":direct,"bandConvection":[b["delta"]["outerConvW"] for b in bands],
             "bandDiffusion":[b["delta"]["outerDiffW"] for b in bands],"bandSide":[b["delta"]["sideConvW"]+b["delta"]["sideDiffW"] for b in bands]}
    for key,powers in targets.items():
        f=fields[key];assert f["unit"]=="W/m²" and set(f["status"])=={1}
        density=np.asarray(powers)/np.bincount(ids,area)
        np.testing.assert_allclose(f["values"],density[ids],rtol=0,atol=1e-8)
        np.testing.assert_allclose(np.bincount(ids,area*np.asarray(f["values"])),powers,rtol=0,atol=1e-13)
    f=fields["bandRelease"];assert side["csvSHA256"]==sha(ui/"faces.csv")
    for name in ("faces","section"):
        with (ui/(name+".csv")).open(encoding="utf-8-sig",newline="") as stream:rows=list(csv.DictReader(stream))
        assert len(rows)==(len(ids) if name=="faces" else proof["sectionRows"])
        for j,row in enumerate(rows):
            i=int(row["face_id"])-1;assert 0<=i<len(ids) and row["field_key"]=="bandRelease" and row["unit"]=="W/m²"
            assert abs(float(row["value"])-f["values"][i])<1e-9
            if name=="faces":
                assert j==i and abs(float(row["area_m2"])-area[i])<1e-18
                assert bands[ids[i]]["name"] in row["scope"] and "按带积分" in row["method"]
            else:assert float(row["section_x_m"])==proof["sectionX"] and ids[i]==5
    out.mkdir(exist_ok=False)
    save(out/"workbench.json",dict(**proof,elapsedSeconds=rec["elapsedSeconds"],packageSHA256=sha(file),
         pageDataSHA256=sha(ui/"workbench/case-data.js"),facesSHA256=sha(ui/"faces.csv"),sectionSHA256=sha(ui/"section.csv")))
    with zipfile.ZipFile(out/"receipts.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for p in build.glob("*.json"):z.write(p,"build/"+p.name)
        for j,path in enumerate(samples,2):
            for n in ("sample.json","step_process.json","case.cfg","step.log"):
                z.write(path/n,f"grid{j}/"+n)
        for n in ("debugger_process.json","ui.json","faces.json","section.csv"):
            z.write(ui/n,"ui/"+n)
    for p in (geom/"geometry.json",budget/"budget.json"):shutil.copyfile(p,out/p.name)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("ui","budget","geometry","build","grid2","grid3","output"):p.add_argument(name,type=Path)
    args=p.parse_args();paths=[v.resolve() for v in vars(args).values()]
    for path in paths:path.relative_to(BASE)
    verify(*paths[:4],paths[4:6],paths[6])
