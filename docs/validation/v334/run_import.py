"""为真实五步回调制作已知面序的 STL，隔离前端并核对独立原表。"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASES = ROOT / "cases/short_sphere"
sys.path.insert(0, str(ROOT / "scripts"))
from runSU2NativeDebugger import observe, sha


def fixture(out, branch):
    with ZipFile(CASES / "supported65.zip") as z:
        raw = z.read("supported_m.stl")
    assert hashlib.sha256(raw).hexdigest() == "eb4e406088070a7b1b80141c2f650469e46477845107439dfa39842636dcf0ff"
    dtype = np.dtype([("normal", "<f4", 3), ("v", "<f4", (3, 3)), ("attribute", "<u2")])
    records = np.frombuffer(raw, dtype=dtype, offset=84).copy()
    order = np.arange(len(records))
    if branch == "half":
        order = np.random.default_rng(65).permutation(len(records))
        records = records[order]
        # 绕 Z 逆时针 90°，工作台需确认反向旋转；循环换位保留三角面法向。
        v = np.roll(records["v"].astype(float) * 1000, 1, axis=1)
        records["v"] = v[:, :, [1, 0, 2]] * [-1, 1, 1] + [.25, .125, -.5]
        normal = records["normal"].copy()
        records["normal"] = normal[:, [1, 0, 2]] * [-1, 1, 1]
    (out / "specimen.stl").write_bytes(raw[:84] + records.tobytes())
    if branch == "full":
        for name in ("reversed", "duplicate"):
            damaged = records.copy()
            if name == "reversed":
                damaged["v"][0] = damaged["v"][0, [0,2,1]]
            else:
                damaged[0] = records[1]
            (out / (name + ".stl")).write_bytes(raw[:84] + damaged.tobytes())
    np.savetxt(out / "expected_order.txt", order + 1, fmt="%d")
    return order, records


def verify(out, branch, order, records):
    import pandas as pd
    with ZipFile(CASES / f"supportedTime65{branch.title()}.zip") as z:
        with z.open("wall.csv") as stream:
            wall = pd.read_csv(stream, float_precision="round_trip").iloc[order]
        with z.open("stations.csv") as stream:
            stations = pd.read_csv(stream)
        with z.open("direction.csv") as stream:
            direction = pd.read_csv(stream).iloc[order]
        s = json.loads(z.read("summary.json"))
    t = pd.read_csv(out / "wall.csv", float_precision="round_trip")
    proof = json.loads((out / "ui.json").read_text(encoding="utf-8"))
    v = records["v"].astype(float) * (.001 if branch == "half" else 1.)
    area = np.linalg.norm(np.cross(v[:, 1]-v[:, 0], v[:, 2]-v[:, 0]), axis=1)/2
    np.testing.assert_array_equal(t.face_id, np.arange(len(order))+1)
    np.testing.assert_array_equal(t.reference_face_id, order+1)
    np.testing.assert_allclose(t.area_m2, area, rtol=1e-11, atol=0)
    np.testing.assert_allclose(t.pressure*1000, wall.absolutePa, rtol=0, atol=3e-10)
    np.testing.assert_allclose(t.wall_input_K, wall.wallInputK, rtol=0, atol=1e-10)
    np.testing.assert_allclose(t.first_fluid_K, wall.firstK, rtol=0, atol=1e-10)
    err = np.asarray(t.heat_flux*1000*area) - np.asarray(wall.fluxInWm2*wall.areaM2)
    assert abs(err).max()<1e-12 and abs(err.sum())<1e-10
    regional = [float(-np.sum(np.asarray(t.heat_flux*1000*area)[wall.regionID==k])) for k in range(1,6)]
    np.testing.assert_allclose(regional,[r["outW"] for r in s["regions"]],rtol=0,atol=1e-10)
    assert proof["peak"]==t.heat_flux.idxmin()+1 and proof["peakReference"]==order[proof["peak"]-1]+1
    assert proof["nativePeakReference"]==s["peakOutFace"]==order[proof["nativePeakFile"]-1]+1
    cut=pd.read_csv(out/"section.csv",float_precision="round_trip")
    assert len(cut)==72 and cut.field_key.eq("heatFlux").all() and cut.unit.eq("kW/m²").all()
    np.testing.assert_allclose(cut.value,t.heat_flux.to_numpy()[cut.face_id-1],rtol=0,atol=1e-13)
    np.testing.assert_array_equal(order[np.array(proof["stations"])-1]+1,stations.originalFaceID)
    np.testing.assert_array_equal(np.loadtxt(out/"reverse.txt",dtype=int),np.flatnonzero(direction.iloc[:,6]==-1)+1)
    assert t.source_frame_sha256.eq(s["frameSHA256"]).all() and t.stl_sha256.eq(sha(out/"specimen.stl")).all()
    assert np.allclose(t.physical_time_s,s["elapsedSeconds"],rtol=1e-12,atol=0) and np.allclose(t.physical_step_s,s["dtSeconds"],rtol=1e-12,atol=0)
    assert t.heat_flux_method.str.contains("功率守恒").all() and t.pressure_method.str.contains("SU2").all()
    if branch == "half":
        js=(out/"workbench/case-data.js").read_text(encoding="utf-8")
        model=json.loads(js[len("window.ROCKET_CASE="):].strip().rstrip(";"))
        fields={f["key"]:f["values"] for f in model["fields"]}
        np.testing.assert_array_equal(model["shortSupportImport"]["referenceFaceID"],order+1)
        np.testing.assert_allclose(fields["heatFlux"],t.heat_flux,rtol=2e-14)
        for key,col in (("ShearPa","shearPa"),("SpeedMps","speedMps"),("SpeedMagnitudeMps","speedMagnitudeMps"),("ShearMagnitudePa","shearMagnitudePa"),("Status","status")):
            expected=direction[col].to_numpy(dtype=float,copy=True)
            if col=="status":expected[expected==0]=np.nan
            np.testing.assert_allclose(np.array(fields["supportTimehalf"+key],dtype=float),expected,rtol=3e-13,atol=1e-15,equal_nan=True)
    return dict(rows=len(order), sectionRows=len(cut),maxFacePowerErrorW=float(abs(err).max()), regionalOutW=regional, **proof)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("name");p.add_argument("branch",choices=["half","full"]);a=p.parse_args()
    assert a.name.isidentifier()
    out=ROOT/"outputs/v334_validation"/a.name;out.mkdir(parents=True,exist_ok=False)
    order,records=fixture(out,a.branch)
    keep=["web","createAerothermalWebTab.m","createAerothermalEvidenceView.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    front=out/"frontend"
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v333",*keep],cwd=ROOT))) as z:
        z.extractall(front)
    shutil.copyfile(ROOT/"web/evidence.js",front/"web/evidence.js")
    def quote(path):
        return "'"+path.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_import({quote(out)},'{a.branch}')"
    files=[*front.rglob("*"),*ROOT.glob("*ShortSupport*.m"),ROOT/"shortSupportImportBoundary.m",ROOT/"runAeroHeatSingleCase.m",
           ROOT/"aeroHeatWorkflowState.m",ROOT/"rocketAeroHeatWorkbench.m",ROOT/"exportWallPairEvidence.m",ROOT/"buildAerothermalViewerData.m",
           HERE/"check_import.m",HERE/"run_import.py",HERE/"dispatch_probe.m.in",out/"specimen.stl"]
    receipt=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                    {str(f):sha(f) for f in files if f.is_file()},240,300*1024**2)
    assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and all(receipt["inputsStable"].values()), receipt["stopReason"]
    assert not(receipt["stopReason"] or receipt["supervisorError"] or receipt["survivingProcessIDs"])
    result=verify(out,a.branch,order,records)
    result.update(seconds=receipt["elapsedSeconds"],peakFamilyRSSBytes=receipt["peakFamilyRSSBytes"],csvSHA256=sha(out/"wall.csv"))
    (out/"checked.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
