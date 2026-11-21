"""核对两套工作台导出与独立原面表，归集外域对照的原生及界面证据。"""
import argparse
import csv
import io
import json
import shutil
import zipfile
from pathlib import Path

import numpy as np

from check_state import BASE,HERE,NEW,ROOT,read,save,sha
from compare_domains import load


def verify(ui,out):
    proof=read(ui/"ui.json");rec=read(ui/"debugger_process.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not (rec["survivingProcessIDs"] or rec["supervisorError"] or rec["stopReason"])
    comparison=read(NEW.parent/"comparison01/comparison.json")
    a,b,_=load(BASE);s,w,_=load(NEW)
    cases=ROOT/"cases/thermistor_b"
    assert sha(cases/"outer30.zip")==read(cases/"outer30.json")["archiveSHA256"]==sha(NEW/"surface/surface.zip")
    checks=read(NEW.parent/"checkpoints.json")
    assert checks["fieldChecksPassed"] and max(checks["residualWindow"].values())< -8
    for f in checks["frames"]:assert sha(NEW/f["file"])==f["sha256"]
    with np.load(NEW/"nodes.npz") as z:
        near=z["firstMap"];wall=np.unique(z["wallMap"]);xyz=z["xyz"]
    plus=np.loadtxt(NEW/"plus20/restart.csv",delimiter=",",skiprows=1)
    np.testing.assert_allclose(plus[:,1:4],xyz,rtol=0,atol=1e-13)
    assert len(near)==len(b) and np.all(plus[wall,9]==1) and plus[:,9].min()>=-1e-12 and plus[:,9].max()<=1+1e-12
    out.mkdir(exist_ok=False)
    modes={}
    with zipfile.ZipFile(out/"face_comparison.zip","x",compression=zipfile.ZIP_DEFLATED) as zipped:
        for mode in ("givenWall","givenPower"):
            p=proof[mode];summary=p["summary"];direct=ui/mode
            assert p["webglError"]==0 and p["visibleSamples"]>0 and p["faces"]==len(w)==10272
            assert summary["baseFrame"]==a["sourceFrameSHA256"] and summary["expandedFrame"]==s["sourceFrameSHA256"]
            assert summary["mode"]==mode and abs(summary["wall30K"]-comparison[mode]["wall30K"])<1e-10
            k0,k=(1.,1.) if mode=="givenWall" else (.004/a["totalOutwardW"],.004/s["totalOutwardW"])
            q0,q=k0*b["heatFluxWm2"],k*w["heatFluxWm2"]
            t0,t=[293.15+factor*(v["nearWallCellTemperatureK"]-293.15) for factor,v in ((k0,b),(k,w))]
            dq=q-q0;dp=w["pressurePa"]-b["pressurePa"]
            raw=(direct/"workbench/case-data.js").read_text("utf8").strip()
            assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
            data=json.loads(raw[len("window.ROCKET_CASE="):-1]);fields={f["key"]:f for f in data["fields"]}
            assert data["outerDomainComparison"]==summary
            assert read(direct/"faces.json")["outerDomainComparison"]==summary
            identity=data["evidenceIdentity"]["sourceFrame"]
            assert identity["sha256"]==a["sourceFrameSHA256"] and identity["outerDomainReference"]["frameSHA256"]==s["sourceFrameSHA256"]
            assert identity["outerDomainReference"]["mode"]==mode and identity["outerDomainReference"]["wallSHA256"]==s["wallSHA256"]
            for key,vals in (("heatFlux",q0/1000),("outerHeat30",q/1000),("outerDeltaPressure",dp),("outerDeltaHeat",dq),("outerDeltaFirst",t-t0)):
                np.testing.assert_allclose(fields[key]["values"],vals,rtol=0,atol=1e-9)
                assert all(code==1 for code in fields[key]["status"])
            assert p["probeFace"]==int(np.argmax(abs(dq)))+1
            for name in ("faces","section"):
                with (direct/(name+".csv")).open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
                assert len(rows)==(len(dq) if name=="faces" else p["cutRows"])
                for j,row in enumerate(rows):
                    i=int(row["face_id"])-1
                    assert 0<=i<len(w) and row["field_key"]=="outerDeltaHeat" and row["unit"]=="W/m²" and abs(float(row["value"])-dq[i])<1e-9
                    if name=="faces":assert i==j and abs(float(row["area_m2"])-w["areaM2"][i])<1e-18
                    else:assert float(row["section_x_m"])==p["sectionX"]
            stream=io.StringIO(newline="")
            writer=csv.writer(stream);writer.writerow(("face_id","region","area_m2","p20_Pa","p30_Pa","q20_into_W_m2","q30_into_W_m2",
                                                       "first20_K","first30_K","Tw20_K","Tw30_K"))
            for j,label in enumerate(data["region"]):
                writer.writerow((j+1,label,w["areaM2"][j],b["pressurePa"][j],w["pressurePa"][j],q0[j],q[j],t0[j],t[j],summary["wall20K"],summary["wall30K"]))
            zipped.writestr(mode+".csv",stream.getvalue().encode("utf-8-sig"))
            modes[mode]=dict(ui=p,faceCSVSHA256=sha(direct/"faces.csv"),cutCSVSHA256=sha(direct/"section.csv"),
                             pageSHA256=sha(direct/"workbench/case-data.js"))
        zipped.writestr("sources.json",json.dumps({m:proof[m]["summary"] for m in modes},ensure_ascii=False,indent=2).encode("utf8"))
    save(out/"workbench.json",dict(modes=modes,seconds=rec["elapsedSeconds"],maxRSSBytes=rec["peakFamilyRSSBytes"],
                                 onlyReferenceOutOfRange=proof.get("onlyReferenceOutOfRange"),faceArchiveSHA256=sha(out/"face_comparison.zip")))
    with zipfile.ZipFile(out/"receipts.zip","x",compression=zipfile.ZIP_DEFLATED) as z:
        for file in NEW.rglob("*.json"):
            if file.name!="body.json":z.write(file,"native/"+file.relative_to(NEW).as_posix())
        for n in ("history.csv","case.cfg","engine_stdout.log"):z.write(NEW/n,"native/"+n)
        z.write(ui/"debugger_process.json","ui/debugger_process.json")
    for src in (NEW.parent/"checkpoints.json",NEW/"inputs.json",NEW/"shell.json",NEW/"surface/summary.json",NEW.parent/"comparison01/comparison.json"):
        shutil.copyfile(src,out/src.name)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("ui",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();ui=a.ui.resolve();out=a.output.resolve()
    for path in (ui,out):path.relative_to(ROOT/"outputs")
    verify(ui,out)
