"""从参考面、实际 STL 面积、两份物理解和导出差图，逐一核对新时刻的两个工作台。"""
import argparse
import csv
import io
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CASE=ROOT/"cases/short_sphere"


def check(folder):
    ui=json.loads((folder/"ui.json").read_bytes());branch=ui["branch"]
    proc=json.loads((folder/"window_process.json").read_bytes())
    assert proc["exitCode"]==0 and not(proc["stopReason"] or proc["livePids"])
    assert ui["wrongTimeRejected"] and ui["faces"]==32480 and ui["pixels"]>0
    late=CASE/"later1";name="supportedTime65"+branch.title()+".zip"
    pair=json.loads((late/"supportedTime65.json").read_bytes());early=json.loads((CASE/"supportedTime65.json").read_bytes())
    assert ui["terminal"]==pair[branch+"Frame"] and ui["previous"]==pair["history"][branch][-2]["frame"]
    assert ui["earlierFrame"]==early[branch+"Frame"] and ui["timeS"]==pair["elapsedS"]
    dtype=[("normal","<f4",(3,)),("vertex","<f4",(3,3)),("tag","<u2")]
    verts=np.fromfile(folder/"input.stl",dtype=dtype,offset=84)["vertex"].astype(float)
    area=.5*np.linalg.norm(np.cross(verts[:,1]-verts[:,0],verts[:,2]-verts[:,0]),axis=1)
    with ZipFile(late/name) as z,ZipFile(CASE/name) as original:
        wall=np.genfromtxt(io.BytesIO(z.read("wall.csv")),delimiter=",",names=True)
        old=np.genfromtxt(io.BytesIO(original.read("wall.csv")),delimiter=",",names=True)
    order=np.arange(len(area)) if ui["mode"]=="archive" else np.arange(len(area)-1,-1,-1)
    ratio=wall["areaM2"][order]/area
    with (folder/"wall.csv").open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
    assert len(rows)==32480
    for i,r in enumerate(rows):
        assert int(r["face_id"])==i+1 and int(r["reference_face_id"])==order[i]+1
        assert r["source_frame_sha256"]==ui["terminal"] and r["reference_time_frame_sha256"]==ui["earlierFrame"]
        assert abs(float(r["reference_time_s"])-early["elapsedS"])<1e-14 and r["joint_verified"]=="0"
    def numbers(key):return np.array([float(r[key]) for r in rows])
    np.testing.assert_allclose(numbers("pressure")*1000,wall["absolutePa"][order],atol=1e-8,rtol=0)
    np.testing.assert_allclose(numbers("heat_flux")*1000,wall["fluxInWm2"][order]*ratio,atol=1e-8,rtol=1e-12)
    np.testing.assert_allclose(numbers("first_fluid_K"),wall["firstK"][order],atol=1e-10,rtol=0)
    np.testing.assert_allclose(numbers("wall_input_K"),wall["wallInputK"][order],atol=1e-10,rtol=0)
    page=json.loads((folder/"workbench/case-data.js").read_text("utf8").split("=",1)[1].strip().rstrip(";"))
    assert page["shortSupportTemporal"]["baseFrame"]==ui["earlierFrame"]
    assert ("4 个整步" if branch=="full" else "8 个半步") in page["resultCondition"]["timeCaption"]
    fields={f["key"]:f for f in page["fields"]}
    errors={}
    for field,csv_key,key,factor in (("timeFluxChange","delta_time_heat_flux_W_m2","fluxInWm2",ratio),
             ("timeGaugeChange","delta_time_gauge_Pa","gaugePa",1.),("timeFirstChange","delta_time_first_K","firstK",1.)):
        wanted=(wall[key]-old[key])[order]*factor
        errors[field]=max(abs(numbers(csv_key)-wanted))
        np.testing.assert_allclose(numbers(csv_key),wanted,rtol=1e-12,atol=1e-11)
        np.testing.assert_allclose(fields[field]["values"],wanted,rtol=1e-12,atol=1e-11)
        limit=max(pair["timeComparison"][b]["fields"][key]["maxAbs"] for b in ("full","half"))*np.max(factor)
        np.testing.assert_allclose(page["comparisonColorRanges"][field],[-limit,limit],rtol=1e-10,atol=1e-14)
        assert all("BDF1" in method for method in fields[field]["methodNames"])
    with (folder/"section.csv").open(encoding="utf-8-sig",newline="") as f:section=list(csv.DictReader(f))
    assert len(section)>0
    for row in section:
        assert row["field_key"]=="timeFluxChange" and row["unit"]=="W/m²"
        face=int(row["face_id"])-1
        assert 0<=face<len(area)
        assert abs(float(row["value"])-fields["timeFluxChange"]["values"][face])<1e-10
    assert abs(area @ numbers("delta_time_heat_flux_W_m2")+ui["powerDeltaW"])<1e-10
    for field in ("wallSolved","edgeTemperature","recoveryTemperature"):
        assert all(v is None for v in fields[field]["values"])
    with (folder/"stations.csv").open(encoding="utf-8-sig",newline="") as f:t=list(csv.DictReader(f))
    assert len(t)==7
    np.testing.assert_allclose([float(r["earlyNu"]) for r in t],early[branch+"Nu"],rtol=1e-12)
    np.testing.assert_allclose([float(r["computedNu"]) for r in t],pair[branch+"Nu"],rtol=1e-12)
    return dict(window=ui,rows=len(rows),sectionRows=len(section),maxDifferenceReadbackErrors=errors,seconds=proc["elapsedSeconds"],peakRSS=proc["peakSampledRSSBytes"])


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("archive",type=Path);p.add_argument("local",type=Path);p.add_argument("output",type=Path);a=p.parse_args()
    result=dict(archiveHalf=check(a.archive),localFull=check(a.local))
    with a.output.open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,indent=2,ensure_ascii=False);f.write("\n")
