"""按两份独立归档和新终帧核验云图、探针、截面导出及初场到终点的保存序列。"""
import argparse
import csv
import json
import zipfile
from pathlib import Path

import numpy as np

from mesh_audit import NEW, ROOT, REF_SMOOTH, WIRED, read, save, sha, shared_faces, terminal
from compare_meshes import table


def audit(ui_folder,out):
    ui=read(ui_folder/"ui.json")
    rec=read(ui_folder/"debugger_process.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
    assert ui["webglError"]==0 and ui["visiblePixelSamples"]>0
    assert ui["sourceFrame"]==sha(NEW/"restart.csv") and ui["smoothFrame"]==sha(REF_SMOOTH/"restart.csv")
    assert ui["coarseFrame"]==sha(WIRED/"restart.csv") and ui["meshSHA256"]==sha(NEW/"mesh.su2")
    cases=ROOT/"cases/thermistor_b"
    archive=cases/"leadsMidbody3.zip";assert sha(archive)==read(archive.with_suffix(".json"))["archiveSHA256"]
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None and z.read("case.json")== (NEW/"case.json").read_bytes()
        assert z.read("body.json")== (NEW/"body.json").read_bytes()
        for name in ("wall.csv","flow.csv","summary.json","flow.json"):
            assert z.read(name)==(NEW/"surface"/name).read_bytes()
        w=table(z.read("wall.csv"))
    with zipfile.ZipFile(cases/"leads3.zip") as z:base=table(z.read("wall.csv"))
    with zipfile.ZipFile(cases/"midbody3.zip") as z:smooth=table(z.read("wall.csv"))
    i,j=shared_faces(read(NEW/"body.json"),read(REF_SMOOTH/"body.json")).T
    body_delta=np.full(len(w),np.nan);body_delta[i]=w["heatFluxWm2"][i]-smooth["heatFluxWm2"][j]
    grid_delta=w["heatFluxWm2"]-base["heatFluxWm2"]
    raw=(ui_folder/"workbench/case-data.js").read_text("utf8").strip()
    assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
    data=json.loads(raw[len("window.ROCKET_CASE="):-1]);fields={f["key"]:f for f in data["fields"]}
    expected={"heatFlux":w["heatFluxWm2"]/1000,"pressure":w["pressurePa"]/1000,
              "nearWallCellTemperature":w["nearWallCellTemperatureK"],"leadGridDeltaHeat":grid_delta,
              "leadGridDeltaPressure":w["pressurePa"]-base["pressurePa"],"leadGridDeltaTemperature":
              w["nearWallCellTemperatureK"]-base["nearWallCellTemperatureK"],"pairedDeltaHeat":body_delta}
    for key,array in expected.items():
        f=fields[key]
        np.testing.assert_allclose(np.asarray(f["values"],float),array,rtol=0,atol=1e-9,equal_nan=True)
        np.testing.assert_array_equal(np.array(f["status"])>0,np.isfinite(array))
    assert all(v==2 for v in fields["wallInput"]["status"]) and not any(fields["wallSolved"]["status"])
    hashes={}
    for prefix,key,delta in (("grid","leadGridDeltaHeat",grid_delta),("geometry","pairedDeltaHeat",body_delta)):
        assert ui[prefix+"Peak"]==np.nanargmax(abs(delta))+1
        for kind in ("faces","cut"):
            path=ui_folder/f"{prefix}_{kind}.csv"
            hashes[path.name]=sha(path)
            with path.open(encoding="utf-8-sig",newline="") as stream:
                exported=list(csv.DictReader(stream))
            assert len(exported)==(len(w) if kind=="faces" else ui[prefix+"CutRows"])
            for ix,row in enumerate(exported):
                face=int(row["face_id"])-1
                assert 0<=face<len(w) and row["field_key"]==key and row["unit"]=="W/m²"
                if kind=="cut":assert float(row["section_x_m"])==ui["sectionX"]
                else:assert face==ix and abs(float(row["area_m2"])-w["areaM2"][face])<1e-18
                if np.isfinite(delta[face]):assert abs(float(row["value"])-delta[face])<1e-9
                else:assert row["value"]=="" and int(row["availability_code"])==0
    term=terminal(NEW)
    with (NEW/"history.csv").open() as stream:last=int(list(csv.reader(stream))[-1][2])
    frames={}
    for step in [*range(500,last,500),last]:
        for kind,ext in (("restart","csv"),("surface","vtk")):
            name=f"{kind}_{step:06d}.{ext}";frames[name]=sha(NEW/name)
            if step==last:assert frames[name]==sha(NEW/f"{kind}.{ext}")==term["outputSHA256"][f"{kind}.{ext}"]
    out.mkdir(exist_ok=False)
    save(out/"workbench.json",dict(ui=ui,seconds=rec["elapsedSeconds"],peakFamilyRSSBytes=rec["peakFamilyRSSBytes"],
         exports=hashes,offlineDataSHA256=sha(ui_folder/"workbench/case-data.js")))
    save(out/"frames.json",dict(directory=NEW.relative_to(ROOT).as_posix(),initial=sha(NEW/"initial_flow.csv"),
         iteration=last,checkpoints=frames))
    with zipfile.ZipFile(out/"receipts.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for f in NEW.rglob("*.json"):
            if f.name!="body.json":z.write(f,"native/"+f.relative_to(NEW).as_posix())
        for name in ("case.cfg","history.csv","engine_stdout.log"):z.write(NEW/name,"native/"+name)
        z.write(ui_folder/"debugger_process.json","ui/debugger_process.json")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ui",type=Path);parser.add_argument("output",type=Path)
    args=parser.parse_args();ui=args.ui.resolve();output=args.output.resolve()
    for folder in (ui,output):folder.relative_to(ROOT/"outputs")
    audit(ui,output)
