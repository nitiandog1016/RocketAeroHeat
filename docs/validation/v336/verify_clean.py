"""核对独立安装及原生终点，逐原面比较相同步长的旧、新结果和两个实窗的导出。"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import tempfile
from zipfile import ZipFile

for var in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[var] = "1"

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "v333"))
from verify_time import verify, table
from build_time import read, sha


def compare_branch(job, branch):
    """用参考面面积比较同一步长，分别核对五区积分、七站、峰值和近壁反向面。"""
    baseline = ROOT / "cases/short_sphere" / ("supportedTime65" + branch.title() + ".zip")
    current = job / "export" / branch / "surface.zip"
    with ZipFile(baseline) as old, ZipFile(current) as new:
        a, b = (table(z, "wall.csv") for z in (old, new))
        ra, rb = (json.loads(z.read("summary.json")) for z in (old, new))
        sta, stb = (table(z, "stations.csv") for z in (old, new))
        da, db = (table(z, "direction.csv") for z in (old, new))
    assert ra["condition"] == rb["condition"] and ra["dtSeconds"] == rb["dtSeconds"]
    assert ra["elapsedSeconds"] == rb["elapsedSeconds"] and len(a) == len(b) == 32480
    for col in ("faceID", "regionID", "areaM2", "wallInputK"):
        np.testing.assert_array_equal(a[col], b[col])
    for col in ("angleDeg", "reportedNu", "originalFaceID"):
        np.testing.assert_array_equal(sta[col], stb[col])
    area = a["areaM2"]
    condition = ra["condition"]
    dynamic = .5 * condition["rho"] * condition["speedMps"] ** 2
    peak = max(abs(a["fluxInWm2"]))
    dp = b["gaugePa"] - a["gaugePa"]
    dq = b["fluxInWm2"] - a["fluxInWm2"]
    rms_p = np.sqrt(np.average(dp ** 2, weights=area))
    rms_q = np.sqrt(np.average(dq ** 2, weights=area))
    region_delta = []
    for region in range(1, 6):
        mask = a["regionID"] == region
        old_watts, new_watts = (-area[mask] @ v["fluxInWm2"][mask] for v in (a, b))
        region_delta.append(new_watts / old_watts - 1)
    nu_delta = stb["computedNu"] / sta["computedNu"] - 1
    power_delta = rb["energy"]["wallW"] / ra["energy"]["wallW"] - 1
    assert rms_p < .01 * dynamic and rms_q < .01 * peak and max(abs(dq)) < .05 * peak
    assert max(abs(np.array(region_delta))) < .01 and max(abs(nu_delta)) < .01 and abs(power_delta) < .01
    assert ra["peakOutFace"] == rb["peakOutFace"] and ra["jointExperimentalAreaM2"] == rb["jointExperimentalAreaM2"] == 0
    np.testing.assert_array_equal(da["faceID"], db["faceID"])
    np.testing.assert_array_equal(da["status"] == -1, db["status"] == -1)
    return dict(pressureRMSPa=rms_p, pressureMaxPa=max(abs(dp)), heatRMSWm2=rms_q, heatMaxWm2=max(abs(dq)),
                pressureRMSOverDynamic=rms_p/dynamic, heatRMSOverPeak=rms_q/peak, heatMaxOverPeak=max(abs(dq))/peak,
                regionalRelativeChange=region_delta, stationRelativeChange=nu_delta.tolist(), wallPowerRelativeChange=power_delta,
                firstTemperatureMaxK=max(abs(b["firstK"]-a["firstK"])), peakFace=rb["peakOutFace"],
                reverseFaces=int(np.count_nonzero(db["status"] == -1)), newFrame=rb["frameSHA256"],
                referenceFrame=ra["frameSHA256"], sameWallBytes=ra["files"]["wall.csv"] == rb["files"]["wall.csv"])


def compare_windows(job, solve, reopen, incomplete_supervisor=False):
    """将两个窗口实际导出的文件面与原生参考面对应，按 STL 实际面积核对热功率。"""
    pair = read(job / "export/pair.json")
    fixture=solve.parent/"supported_m.stl"
    triangles=np.fromfile(fixture, dtype=[("normal","<f4",(3,)), ("vertex","<f4",(3,3)), ("tag","<u2")], offset=84)
    vertices=triangles["vertex"].astype(float)
    file_area=.5*np.linalg.norm(np.cross(vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,0]),axis=1)
    with ZipFile(job / "export/half/surface.zip") as z:
        wall = table(z, "wall.csv")
        stations = table(z, "stations.csv")
        summary = json.loads(z.read("summary.json"))
    reports = []
    for folder, action in ((solve, "solve"), (reopen, "reopen")):
        ui = read(folder / "ui.json")
        if action=="solve" and incomplete_supervisor:
            failure=read(folder/"supervisor_failure.json")
            assert not (folder/"window_process.json").exists() and failure["pythonSupervisorExit"]==1
            assert failure["verifiedUIReceiptSHA256"]==sha(folder/"ui.json") and failure["matlabExitCode"] is None
            assert not failure["observedWindowPidsAfterExit"]
            resource=dict(seconds=None,peakBytes=None,matlabExitCode=None,supervisorFailure=failure["waitError"])
        else:
            proc = read(folder / "window_process.json")
            assert not (proc["exitCode"] or proc["livePids"] or proc["stopReason"])
            resource=dict(seconds=proc["elapsedSeconds"],peakBytes=proc["peakSampledRSSBytes"],matlabExitCode=proc["exitCode"])
        assert Path(ui["directory"]).resolve() == job and ui["action"] == action and ui["source"]["sha256"] == pair["halfFrame"]
        assert ui["source"]["previousSHA256"] == pair["firstHalfFrame"] and ui["source"]["initialSHA256"] == pair["seed"]
        assert ui["webGLError"] == 0 and ui["pixels"] > 0
        assert ui["faces"] == len(wall) and ui["peak"] == summary["peakOutFace"]
        assert ui["reverseFaces"] == pair["halfReverse"]["faces"]
        station_columns=("angleDeg","originalFaceID","reportedNu","computedNu")
        np.testing.assert_allclose(ui["stations"],np.column_stack([stations[k] for k in station_columns]),rtol=1e-12,atol=1e-11)
        np.testing.assert_allclose(ui["regionalOutW"], [r["outW"] for r in summary["regions"]], atol=1e-11, rtol=0)
        values = []
        with (folder / "wall.csv").open(encoding="utf-8-sig", newline="") as f:
            for i, row in enumerate(csv.DictReader(f), 1):
                assert int(row["face_id"]) == int(row["reference_face_id"]) == i
                assert row["source_frame_sha256"] == pair["halfFrame"] and row["volume_mesh_sha256"] == pair["meshSHA256"]
                assert row["pressure_unit"] == "kPa" and row["heat_flux_unit"] == "kW/m²"
                assert row["pressure_method"] and row["heat_flux_method"] and row["joint_verified"] == "0"
                values.append([float(row[k]) for k in ("area_m2","pressure","heat_flux","wall_input_K","first_fluid_K")])
        area, pressure, flux, given, first = np.array(values).T
        np.testing.assert_allclose(area,file_area,atol=0,rtol=1e-11)
        np.testing.assert_allclose(pressure*1000, wall["absolutePa"], rtol=0, atol=1e-8)
        np.testing.assert_allclose(flux*1000*area, wall["fluxInWm2"]*wall["areaM2"], rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose([given, first], [wall["wallInputK"], wall["firstK"]], rtol=0, atol=1e-10)
        np.testing.assert_allclose(-1000*(area@flux), ui["outW"], rtol=0, atol=1e-10)
        reports.append(dict(action=action, wallSHA256=sha(folder / "wall.csv"), pixels=ui["pixels"],**resource))
        if action == "reopen":
            file=folder/"workbench/case-data.js"
            page=json.loads(file.read_text("utf8").split("=",1)[1].strip().rstrip(";"))
            assert page["evidenceIdentity"]["sourceFrame"]["sha256"] == pair["halfFrame"]
            fields={field["key"]:field["values"] for field in page["fields"]}
            for key, expected in (("heatFlux", flux), ("pressure", pressure),
                                  ("wallInput", given), ("nearWallCellTemperature", first)):
                np.testing.assert_allclose(fields[key],expected,atol=1e-10,rtol=1e-12)
            for key in ("wallSolved","edgeTemperature","recoveryTemperature"):
                assert len(fields[key])==len(wall) and all(v is None for v in fields[key])
            reports[-1]["offlineDataSHA256"] = sha(file)
    assert reports[0]["wallSHA256"] == reports[1]["wallSHA256"]
    return reports


def audit(app, job, windows=None, reread=None, incomplete_supervisor=False):
    job.relative_to(app / "outputs")
    origin = read(job / "construction.json")
    steady = app / origin["steady"]
    steady.relative_to(app / "outputs")
    stage = read(steady / "converged.json")
    assert stage["frame"] == origin["steadyFrame"]
    with tempfile.TemporaryDirectory(dir=job.parent) as folder:
        data = Path(folder)
        for dest, source in (("supportedTime65.json", "pair.json"), ("supportedTime65Full.zip", "full/surface.zip"),
                             ("supportedTime65Half.zip", "half/surface.zip")):
            os.link(job / "export" / source, data / dest)
        raw = verify(data, job / "run", base_frame=stage["frame"], steady=steady / stage["stage"] / "restart_000199.csv")
    result = dict(raw=raw, regression={branch:compare_branch(job, branch) for branch in ("full", "half")})
    if windows:
        result["windows"] = compare_windows(job, *windows,incomplete_supervisor=incomplete_supervisor)
    if reread:
        receipt=read(reread/"completed.json")
        assert Path(receipt["source"]).resolve()==job/"run"
        hashes={}
        for name in ("construction.json","export/pair.json","export/full/surface.zip","export/half/surface.zip"):
            hashes[name]=sha(reread/name)
            assert hashes[name]==sha(job/name)
        assert hashes["export/pair.json"]==receipt["pairSHA256"]
        result["rereadSHA256"]=hashes
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("app", type=Path)
    parser.add_argument("job", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--windows", nargs=2, type=Path)
    parser.add_argument("--reread",type=Path)
    parser.add_argument("--record-solve-supervisor-failure",action="store_true",
                        help="保留 solve02 的独立窗口监督失败，逐行审核其实际导出")
    args = parser.parse_args()
    report = audit(args.app.resolve(), args.job.resolve(), args.windows, args.reread,args.record_solve_supervisor_failure)
    with args.output.open("x", encoding="utf8") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
