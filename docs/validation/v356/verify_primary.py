"""把工作台导出的每一面与本次实际终帧和上一版独立算出的细网场分别比较。"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

import numpy as np

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
sys.path.insert(0,str(DOC.parent/"v351"))
from verify_scopes import EXPORT


def read(path):
    return json.loads(path.read_bytes())


def verify(ui,web):
    r=read(ui/"ui.json");job=Path(r["runDirectory"])
    case=ROOT/"cases/simple_launcher/wake"
    with ZipFile(case/"wake_heat_grid.zip") as z:
        reference=np.loadtxt(z.open("fine_wall.csv"),delimiter=",",skiprows=1)
    with ZipFile(case/"wind_tunnel.zip") as z:
        coarse=np.loadtxt(z.open("wall.csv"),delimiter=",",skiprows=1)
    w=job/("heat/disturbed_half_2/wall.csv" if r["action"]=="solve" else "wall/fine_wall.csv")
    native=np.loadtxt(w,delimiter=",",skiprows=1)
    assert hashlib.sha256(w.read_bytes()).hexdigest()==r["wallSHA256"]
    np.testing.assert_array_equal(native[:,:3],reference[:,:3])
    np.testing.assert_allclose(native[:,3:],reference[:,3:],atol=1e-12,rtol=0)
    assert abs(-native[:,1]@native[:,3]-r["outwardW"])<1e-13
    assert abs(sum(r["regionW"].values())-r["outwardW"])<1e-13
    t=read(case/"condition.json")["elapsedPhysicalTimeS"]
    values={"heatFlux":(native[:,3],"W/m²"),"pressure":(native[:,2]/1000,"kPa"),
            "wallInput":(np.full(len(native),300.001),"K"),"nearWallCellTemperature":(native[:,4],"K"),
            "wakeCoarseHeat":(coarse[:,3],"W/m²"),"wakeHeatDifference":(native[:,3]-coarse[:,3],"W/m²")}
    counts={}
    for key in r["fields"]:
        field,unit=values[key]
        process=subprocess.run(["node","-e",EXPORT,str((ui/"workbench").resolve()),str(web.resolve()),key],
                               capture_output=True,check=True,timeout=45)
        table=list(csv.DictReader(io.StringIO(process.stdout.decode("utf-8-sig"))))
        assert len(table)==len(native)
        assert all(row["单位"]==unit and row["场量标识"]==key and row["算法"] and float(row["时间 / s"])==t for row in table)
        if key=="heatFlux":
            assert all("细热网" in row["算法"] and "细网" in row["物理量"] for row in table)
        np.testing.assert_array_equal([int(next(iter(row.values()))) for row in table],native[:,0])
        np.testing.assert_allclose([float(row["数值"]) for row in table],field,rtol=0,atol=1e-12)
        with (ui/(key+"_section.csv")).open(encoding="utf-8-sig",newline="") as stream:
            section=list(csv.DictReader(stream))
        assert section and all(float(row["elapsed_physical_time_s"])==t for row in section)
        np.testing.assert_allclose([float(row["value"]) for row in section],field[[int(row["face_id"])-1 for row in section]],rtol=0,atol=1e-12)
        counts[key]={"fullRows":len(table),"sectionRows":len(section)}
    with (ui/"wall_pair.csv").open(encoding="utf-8-sig",newline="") as stream:
        pair=list(csv.DictReader(stream))
    assert len(pair)==len(native) and all(row["source_traceable"]=="1" and row["joint_verified"]=="0" and
                                          "细热网" in row["heat_flux_method"] for row in pair)
    np.testing.assert_array_equal([int(row["face_id"]) for row in pair],native[:,0])
    np.testing.assert_allclose([float(row["pressure"])*1000 for row in pair],native[:,2],atol=1e-8,rtol=0)
    np.testing.assert_allclose([float(row["heat_flux"]) for row in pair],native[:,3],atol=1e-13,rtol=0)
    evidence=dict(ui=r,rows=counts,pairRows=len(pair),areaM2=float(native[:,1].sum()),
                  maxFluxChangeFromV355Wm2=float(max(abs(native[:,3]-reference[:,3]))),
                  maxFirstCellChangeFromV355K=float(max(abs(native[:,4]-reference[:,4]))))
    if r["action"]=="solve":
        done=read(job/"completed.json")
        assert done["summary"]["wallSHA256"]==r["wallSHA256"] and len(done["steps"])==2
        assert sorted(p.name for p in (job/"heat").iterdir())==["disturbed_half_1","disturbed_half_2"]
        evidence["newRun"]={key:done[key] for key in ("flowReusedAfterVerification","initialReusedAfterVerification")}
        for i,record in enumerate(done["steps"],1):
            folder=job/"heat"/f"disturbed_half_{i}"
            check=read(folder/"debugger_process.json")
            assert check["status"]=="terminal" and check["debuggerExitCode"]==0 and all(check["inputsStable"].values())
            assert not (check["supervisorError"] or check["stopReason"] or check["survivingProcessIDs"])
            h=record["heat"]
            assert h["accepted"] and hashlib.sha256((folder/"thermal.npz").read_bytes()).hexdigest()==h["fieldSHA256"]
            evidence["newRun"][folder.name]=dict(flowSHA256=h["sourceFrameSHA256"],initialSHA256=h["initialThetaSHA256"],
                thermalSHA256=h["fieldSHA256"],wallSHA256=h["wallSHA256"],elapsedPhysicalS=h["elapsedPhysicalS"],
                storageW=record["independentStorageW"],farW=h["farW"],outwardW=h["outwardW"],residual=h["originalResidual"],
                seconds=check["elapsedSeconds"],peakRSSBytes=check["peakFamilyRSSBytes"])
    return evidence


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("ui",type=Path);p.add_argument("web",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args()
    a.output.write_bytes((json.dumps(verify(a.ui,a.web),ensure_ascii=False,indent=2,allow_nan=False)+"\n").encode())
