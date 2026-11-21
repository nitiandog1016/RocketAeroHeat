"""逐行核对真实导出、冻结原面归档、末帧及保存序列；汇集轻量原生回执。"""
import csv
import json
import zipfile

import numpy as np

from check_grid import HERE, NEW, OLD, ROOT, read, save, sha
from compare_meshes import table
from thermistorReaction import terminal


def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as src:
        return list(csv.DictReader(src))


def audit():
    proof=HERE/"evidence"
    cases=ROOT/"cases/thermistor_b"
    for key, folder in (("midbody3",NEW),("paired3",OLD)):
        archive=cases/(key+".zip")
        assert sha(archive)==read(archive.with_suffix(".json"))["archiveSHA256"]
        with zipfile.ZipFile(archive) as z:
            assert z.testzip() is None and z.read("body.json")== (folder/"body.json").read_bytes()
            for name in ("wall.csv","flow.csv","summary.json"):
                assert z.read(name)==(folder/"surface"/name).read_bytes()
    fresh=table((NEW/"surface/wall.csv").read_bytes())
    base=table((OLD/"surface/wall.csv").read_bytes())
    delta=fresh["heatFluxWm2"]-base["heatFluxWm2"]
    path=NEW.parent/"midbody_ui01"
    ui=read(path/"ui.json")
    raw=(path/"workbench/case-data.js").read_text("utf8").strip()
    assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
    data=json.loads(raw[len("window.ROCKET_CASE="):-1])
    fields={f["key"]:f for f in data["fields"]}
    assert (ui["sourceFrame"]==sha(NEW/"restart.csv") and ui["meshSHA256"]==sha(NEW/"mesh.su2")
            and ui["faces"]==len(fresh) and ui["peakDeltaFaceID"]==int(np.argmax(abs(delta)))+1)
    for key, column in (("heatFlux",fresh["heatFluxWm2"]/1000),("nearWallCellTemperature",fresh["nearWallCellTemperatureK"]),
                        ("pressure",fresh["pressurePa"]/1000),("pairedSmoothHeat",base["heatFluxWm2"]/1000),
                        ("pairedDeltaPressure",fresh["pressurePa"]-base["pressurePa"]),("pairedDeltaHeat",delta)):
        np.testing.assert_allclose(fields[key]["values"],column,atol=1e-9,rtol=0)
        assert all(i==1 for i in fields[key]["status"])
    assert not any(fields["wallSolved"]["status"])
    assert all(i==2 for i in fields["wallInput"]["status"])
    for name in ("faces","section"):
        exported=rows(path/(name+".csv"))
        assert len(exported)==(len(delta) if name=="faces" else ui["sectionRows"])
        for j,r in enumerate(exported):
            face=int(r["face_id"])-1
            assert 0<=face<len(delta) and r["field_key"]=="pairedDeltaHeat" and r["unit"]=="W/m²"
            assert abs(float(r["value"])-delta[face])<1e-9 and r["qualification"]=="未验证"
            if name=="faces":
                assert face==j and abs(float(r["area_m2"])-fresh["areaM2"][j])<1e-18
            else:
                assert float(r["section_x_m"])==ui["sectionX"]
    wired=NEW.parent/"wired_ui01"
    wired_proof=read(wired/"ui.json")
    assert wired_proof==dict(read(HERE.parent/"v298/evidence/workbench.json")["wired"]["workbench"],
                             page=str(wired/"workbench/index.html"))
    terminal(NEW)
    checkpoints={}
    with (NEW/"history.csv").open() as stream:
        last=int(list(csv.reader(stream))[-1][2])
    for step in [*range(500,last,500),last]:
        for stem,ext in (("restart","csv"),("surface","vtk")):
            filename=f"{stem}_{step:06d}.{ext}"
            checkpoints[filename]=sha(NEW/filename)
            if step==last:
                assert checkpoints[filename]==sha(NEW/f"{stem}.{ext}")
    timing={}
    with zipfile.ZipFile(proof/"receipts.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for f in NEW.rglob("*.json"):
            if f.name!="body.json":
                z.write(f,"native/"+f.relative_to(NEW).as_posix())
        for name in ("case.cfg","history.csv","engine_stdout.log"):
            z.write(NEW/name,"native/"+name)
        for folder in (path,wired):
            rec=read(folder/"debugger_process.json")
            assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
            assert not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
            timing[folder.name]=rec["elapsedSeconds"]
            z.write(folder/"debugger_process.json",folder.name+"/debugger_process.json")
    save(proof/"checkpoints.json",dict(folder=NEW.relative_to(ROOT).as_posix(),
         initial=sha(NEW/"initial_flow.csv"),terminalIteration=last,files=checkpoints))
    save(proof/"workbench.json",dict(midbody=ui,wiredRegression=wired_proof,seconds=timing,
         allRowsCompared=len(delta),faceCSVSHA256=sha(path/"faces.csv"),cutCSVSHA256=sha(path/"section.csv"),
         dataSHA256=sha(path/"workbench/case-data.js")))


if __name__=="__main__":
    audit()
