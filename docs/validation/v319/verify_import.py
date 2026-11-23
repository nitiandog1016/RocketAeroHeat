"""直接按二进制 STL 面记录和原生壁面表核对导出、热流积分、区域及截平面。"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np

from make_fixture import FACET, ROOT

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent / "v311"))
from verify_delivery import cut_values


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def archived(stem, member="wall.csv"):
    with ZipFile(ROOT/"cases/thermistor_b"/(stem+".zip")) as z:
        return np.genfromtxt(io.BytesIO(z.read(member)),delimiter=",",names=True)


def check_text(ui):
    proof=read(ui/"text-ui.json");source=ui/"translated-m.stl"
    receipt=read(ui/"debugger_process.json")
    assert receipt["debuggerExitCode"]==0 and all(receipt["inputsStable"].values()) and not receipt["survivingProcessIDs"]
    assert digest(source)==proof["sha256"] and proof["scale"]==1
    corners=np.array([s.split()[1:] for s in source.read_text().splitlines() if s.startswith("vertex ")],dtype=float)
    center=np.unique(corners,axis=0).mean(axis=0)
    tri=((corners-center)@np.array(proof["rotationMatrix"]).T+center).reshape(-1,3,3)
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    assert abs(area.sum()-proof["areaM2"])<1e-18
    order=np.loadtxt(ui.parent/"fixture/reference_ids.txt",dtype=int)-1
    p=archived("windPoint12")["pressurePa"][order]/1000
    segments=cut_values(ui/"text-section.csv",tri,proof["xM"],p,"pressure","kPa")
    return dict(faces=len(area),areaM2=float(area.sum()),sectionRecords=len(segments),xM=proof["xM"],
                maxVertexErrorM=proof["maxVertexErrorM"],maxRelativeAreaError=proof["maxAreaRelativeError"],
                sectionSHA256=digest(ui/"text-section.csv"),rejections=read(ui/"identity.json"))


def check(ui):
    if (ui/"text-ui.json").is_file():
        return check_text(ui)
    run = read(ui/"ui.json")
    fixture = ROOT/"outputs/v319_validation/fixture"
    path = fixture/"probe-mm.stl"
    with path.open("rb") as f:
        f.seek(80)
        n = int.from_bytes(f.read(4),"little")
        tri = np.fromfile(f,dtype=FACET)["xyz"].astype(float)*.001
    order = np.loadtxt(fixture/"reference_ids.txt",dtype=int)-1
    area = np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    assert n==10272 and np.array_equal(np.sort(order),np.arange(n))
    assert run["source"]["sha256"]==digest(path)
    prepared = Path(run["prepared"]); grid = read(prepared/"prepared.json")
    assert not grid["reusedGrid"] and grid["meshSHA256"]==run["meshSHA256"]==digest(prepared/"mesh.su2")
    body = read(prepared/"body.json")
    names = np.asarray(body["regionNames"])
    rid = np.asarray(body["regionID"])[order]
    old = archived("longWire12")
    factor = read(ROOT/"cases/thermistor_b/resistance12.json")["factor"]
    actions = ("native",) if "native" in run else ("solve","archive","uniform")
    if "native" in run:
        terminal=read(ui/"terminal.json")
        assert terminal["exitCode"]==0 and terminal["inputsStable"] and terminal["pid"]==read(ui/"invocation.json")["pid"]
    else:
        terminal=read(ui/"debugger_process.json")
        assert terminal["status"]=="terminal" and terminal["debuggerExitCode"]==0 and all(terminal["inputsStable"].values())
        assert not(terminal["stopReason"] or terminal["supervisorError"] or terminal["survivingProcessIDs"])
    output = {}
    for action in actions:
        proof = run[action]
        dest = ui/action
        raw = (dest/"workbench/case-data.js").read_text("utf8").strip()
        assert raw.startswith("window.ROCKET_CASE=") and raw.endswith(";")
        page = json.loads(raw[len("window.ROCKET_CASE="):-1])
        native = {}
        if action=="native":
            job = Path(proof["directory"]);work = job/"run"
            done = read(work/"completed.json");summary = read(work/"surface/summary.json")
            ack = read(job/"progress.json")
            assert read(job/("receipt-"+ack["invocation"]+".json"))==ack and ack["status"]=="complete"
            stage = work/done["stage"]
            assert proof["archiveSource"] is False and done["frame"]==proof["frame"]["sha256"]==digest(stage/"restart.csv")
            assert done["summarySHA256"]==digest(work/"surface/summary.json")
            assert done["wallSHA256"]==digest(work/"surface/wall.csv")
            assert summary["meshSHA256"]==run["meshSHA256"] and summary["bodySHA256"]==digest(prepared/"body.json")
            assert all(max(summary[k].values()) < -8 for k in ("window","plus20Window"))
            assert abs(summary["energy"]["relativeDifference"])<.005 and summary["relativeChange20"]<.002
            phases=[]
            for child in sorted(work.glob("steady[0-9][0-9]")):
                rec=read(child/"debugger_process.json");seal=read(child/"native.json")
                assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and not(rec["stopReason"] or rec["survivingProcessIDs"] or rec["supervisorError"])
                assert all(rec["inputsStable"].values()) and seal["outputs"]["restart.csv"]==digest(child/"restart.csv")
                assert seal["steps"]==320 and all((child/("restart_%06d.csv"%i)).is_file() for i in (100,200,300))
                phases.append(dict(stage=child.name,seconds=rec["elapsedSeconds"],frame=seal["outputs"]["restart.csv"],
                                   window=max(seal["last20"].values())))
            for name in ("reaction0","plus20","reaction20","step_control","step_observed"):
                child=stage/name;rec=read(child/"debugger_process.json");seal=read(child/"native.json")
                assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
                assert not(rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
                assert seal["outputs"]["restart.csv"]==digest(child/"restart.csv")
            assert digest(stage/"reaction0/restart.csv")==summary["reactionSHA256"]
            assert digest(stage/"plus20/restart.csv")==summary["plus20FrameSHA256"]
            assert digest(stage/"reaction20/restart.csv")==summary["reaction20SHA256"]
            assert digest(work/"prepared/mesh.su2")==run["meshSHA256"]
            native=dict(stages=phases,extra20Change=summary["relativeChange20"],energyError=summary["energy"]["relativeDifference"],
                        invocation=ack["invocation"],frame=done["frame"],reaction=summary["reactionSHA256"],wall=done["wallSHA256"])
            source=np.genfromtxt(work/"surface/wall.csv",delimiter=",",names=True);scale=1
        elif action=="uniform":
            source=archived("windPoint12","six_wall.csv");scale=2
        elif action=="archive":
            source=archived("resistance12");scale=1
        else:
            source=old;scale=factor
        source=source[order];ref_area=source["areaM2"]
        power=ref_area*source["heatFluxWm2"]*scale
        q=power/area
        tw=np.full(n,305.15) if action=="uniform" else 293.15+scale*(source["wallInputK"]-293.15)
        first=293.15+scale*(source["nearWallCellTemperatureK"]-293.15)
        fields={f["key"]:f for f in page["fields"]}
        def close(actual,expected,tol):
            np.testing.assert_allclose(actual,expected,rtol=0,atol=tol)
        close(np.array(page["vertices"])[page["faces"]],tri,1e-18)
        close(page["faceArea"],area,1e-22)
        close(fields["heatFlux"]["values"],q/1000,1e-10)
        close(fields["pressure"]["values"],source["pressurePa"]/1000,1e-11)
        close(fields["pressureGauge"]["values"],source["pressureGaugePa"],1e-11)
        close(fields["wallInput"]["values"],tw,1e-10)
        close(fields["nearWallCellTemperature"]["values"],first,1e-10)
        if action!="uniform":
            uniform_scale=factor if action in ("solve","archive") else 1
            dq=(scale*source["heatFluxWm2"]-uniform_scale*source["uniformHeatFluxWm2"])*ref_area/area
            close(fields["wireHeatFluxChange"]["values"],dq/1000,1e-10)
        if action in ("solve","archive"):
            dq=(scale*source["heatFluxWm2"]-old["heatFluxWm2"][order])*ref_area/area
            close(fields["thermometryHeatFluxChange"]["values"],dq/1000,1e-10)
        assert np.array_equal(page["region"],names[rid-1]) and proof["pixels"]>0 and proof["webGLError"]==0
        assert np.array_equal(page["thermistorImport"]["referenceFaceID"],order+1)
        with (dest/"wall.csv").open(encoding="utf-8-sig",newline="") as stream:
            rows=list(csv.DictReader(stream))
        assert len(rows)==n
        for i,row in enumerate(rows):
            assert int(row["face_id"])==i+1 and int(row["reference_face_id"])==order[i]+1 and row["region"]==names[rid[i]-1]
            assert row["stl_sha256"]==run["source"]["sha256"] and row["volume_mesh_sha256"]==run["meshSHA256"]
            assert abs(float(row["area_m2"])-area[i])<1e-22 and abs(float(row["reference_area_m2"])-ref_area[i])<1e-22
            assert abs(float(row["pressure"])*1000-source["pressurePa"][i])<1e-8
            assert abs(float(row["heat_flux"])*1000-q[i])<1e-7
            assert abs(float(row["wall_input_K"])-tw[i])<1e-10 and abs(float(row["first_fluid_K"])-first[i])<1e-10
            assert "STL" in row["heat_flux_method"] and row["pressure_method"] and row["joint_verified"] in ("0","false")
        for k,(region,probe) in enumerate(zip(proof["regions"],proof["probeIDs"]),1):
            take=rid==k
            assert rid[probe-1]==k and region["name"]==names[k-1] and region["faces"]==int(take.sum())
            close(region["areaM2"],area[take].sum(),1e-18);close(region["outwardW"],-power[take].sum(),1e-12)
        segments=cut_values(dest/"section.csv",tri,.000172,q/1000,"heatFlux","kW/m²")
        peak=int(np.argmin(q))
        assert proof["peakFaceID"]==peak+1 and proof["referencePeak"]==order[peak]+1
        close(proof["inwardW"],power.sum(),1e-12)
        residual=np.abs(np.array(fields["heatFlux"]["values"])*1000*area-power)
        output[action]=dict(faces=n,areaM2=float(area.sum()),wholeOutW=float(-power.sum()),beadOutW=float(-power[np.isin(rid,[3,4,5])].sum()),
                            peakFaceID=peak+1,referencePeak=int(order[peak])+1,sectionRecords=len(segments),regions=proof["regions"],
                            maxFacePowerDifferenceW=float(residual.max()),csvSHA256=digest(dest/"wall.csv"),
                            sectionSHA256=digest(dest/"section.csv"),fieldFrame=proof["frame"]["sha256"],native=native)
    return output


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("ui",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();out=check(a.ui);a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("x",encoding="utf8") as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write("\n")
