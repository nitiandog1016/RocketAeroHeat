"""用实际 MATLAB 回调和独立归档逐行检查外域差、测站及周向截面。"""
import argparse
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

from outer_grid import HERE, ROOT, read, sha
sys.path.insert(0,str(HERE.parent/"v334"))
from run_import import fixture
from runSU2NativeDebugger import observe
from compare_domains import CASE,snapshot
import numpy as np
import pandas as pd


def check(out,branch,order,records):
    a,x,s,_=snapshot(CASE,branch,read(CASE/"supportedTime65.json"))
    p=CASE/"outer20";b,y,t,_=snapshot(p,branch,read(p/"supportedTime65.json"))
    file=pd.read_csv(out/"wall.csv",float_precision="round_trip");ui=read(out/"ui.json")
    np.testing.assert_array_equal(file.reference_face_id,order+1)
    v=records["v"].astype(float)*(.001 if branch=="half" else 1)
    area=np.linalg.norm(np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]),axis=1)/2
    np.testing.assert_allclose(file.area_m2,area,rtol=1e-11,atol=0)
    factor=x["areaM2"][order]/area
    delta=(y["fluxInWm2"]-x["fluxInWm2"])[order]*factor
    np.testing.assert_allclose(file.delta_heat_flux_W_m2,delta,rtol=0,atol=2e-11)
    np.testing.assert_allclose(file.delta_gauge_Pa,(y["gaugePa"]-x["gaugePa"])[order],rtol=0,atol=2e-14)
    np.testing.assert_allclose(file.heat_flux*1000,y["fluxInWm2"][order]*factor,rtol=2e-12)
    np.testing.assert_allclose(file.pressure*1000,y["absolutePa"][order],rtol=0,atol=3e-10)
    for col,key in (("wall_input_K","wallInputK"),("first_fluid_K","firstK")):
        np.testing.assert_allclose(file[col],y[key][order],rtol=0,atol=1e-10)
    assert file.outer_radius_D.eq(20).all() and file.reference_radius_D.eq(14).all()
    for column,expected in (("stl_sha256",sha(out/"specimen.stl")),("reference_domain_frame_sha256",a["frameSHA256"]),
                            ("source_frame_sha256",b["frameSHA256"]),("reference_domain_mesh_sha256",a["meshSHA256"]),
                            ("volume_mesh_sha256",b["meshSHA256"])):
        assert file[column].eq(expected).all()
    stations=np.array(ui["stations"])
    np.testing.assert_array_equal(order[stations[:,1].astype(int)-1]+1,t["originalFaceID"])
    np.testing.assert_allclose(stations[:,2:],np.column_stack((t["reportedNu"],s["computedNu"],t["computedNu"])),rtol=1e-13)
    assert abs(area@file.delta_heat_flux_W_m2+ui["powerW"]-ui["basePowerW"])<1e-10
    cut=pd.read_csv(out/"section.csv",float_precision="round_trip")
    assert len(cut)==72 and cut.field_key.eq("outerFluxChange").all() and cut.unit.eq("W/m²").all()
    np.testing.assert_allclose(cut.value,delta[cut.face_id-1],rtol=0,atol=2e-11)
    if branch=="half":
        js=(out/"workbench/case-data.js").read_text(encoding="utf8")
        model=json.loads(js.split("=",1)[1].strip().rstrip(";"));fields={f["key"]:f["values"] for f in model["fields"]}
        for key in ("wallSolved","edgeTemperature","recoveryTemperature"):assert all(v is None for v in fields[key])
        np.testing.assert_allclose(fields["outerFluxChange"],delta,rtol=0,atol=2e-11)
        assert model["shortSupportOuter"]["baseFrame"]==a["frameSHA256"]
    return dict(**ui,sectionRows=len(cut),maxDeltaFluxErrorWm2=float(np.max(abs(file.delta_heat_flux_W_m2-delta))),csv=sha(out/"wall.csv"))


def main(name,branch):
    assert name.isidentifier()
    out=ROOT/"outputs/v335_validation"/name/branch;out.mkdir(parents=True,exist_ok=False)
    order,recs=fixture(out,branch);front=out/"frontend"
    paths=["web","createAerothermalWebTab.m","createAerothermalEvidenceView.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v334",*paths],cwd=ROOT))) as z:z.extractall(front)
    for path in ("web/evidence.js","createAerothermalEvidenceView.m"):shutil.copyfile(ROOT/path,front/path)
    def q(f):return "'"+f.as_posix().replace("'","''")+"'"
    cmd=f"cd({q(out)});addpath({q(ROOT)},{q(HERE)});addpath({q(front)},'-begin');check_domain_ui({q(out)},'{branch}')"
    locked=[*front.rglob("*"),*ROOT.glob("*ShortSupport*.m"),ROOT/"shortSupportImportBoundary.m",ROOT/"buildAerothermalViewerData.m",
            ROOT/"exportWallPairEvidence.m",ROOT/"rocketAeroHeatWorkbench.m",ROOT/"aeroHeatWorkflowState.m",ROOT/"runAeroHeatSingleCase.m",
            HERE/"check_domain_ui.m",*list((CASE/"outer20").glob("*")),out/"specimen.stl"]
    r=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
              {str(f):sha(f) for f in locked if f.is_file()},240,300*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"])
    result=check(out,branch,order,recs);result.update(seconds=r["elapsedSeconds"],peakRSS=r["peakFamilyRSSBytes"])
    (out/"checked.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("name");p.add_argument("branch",choices=("full","half"))
    a=p.parse_args();main(a.name,a.branch)
