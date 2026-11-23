"""从原生节点反力和独立交叠积分复算 10,272 面，逐行核对侧栏、三维快照与导出。"""
import argparse
import csv
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT / "scripts"), str(HERE.parent / "v311"), str(HERE.parent / "v314")]
from audit_angular import independent_overlap
from verify_delivery import cut_values
from thermistorBead import read, sha, write_json
from thermistorReaction import terminal


def verify(job, controls, ui, evidence, reference_job=None):
    # 终态、已暂停的空任务、两份历史参照的侧栏入口分别核对。
    done = read(job / "completed.json")
    stage = job / done["stage"]
    for name in ("replay_pause.json","surface_pause.json"):
        pause = read(job.parent / name)
        assert pause["stage"] == stage.name and pause["pause"]["job"] == done["token"]
        assert all(sha(stage / file) == h for file,h in pause["recordedBeforeResume"].items())
    s = read(stage / "surface/summary.json")
    proof = read(ui / "ui.json")
    buttons = read(controls / "ui.json")
    paused = Path(buttons["pausedDirectory"])
    assert read(paused / "progress.json")["status"] == "paused" and not list(paused.glob("steady*"))
    assert sha(paused / "cold/initial_flow.csv") == done["initialSHA256"] == sha(job / "cold/initial_flow.csv")
    ref = read(ROOT / "cases/thermistor_b/windPoint12.json")
    assert buttons["solve"] == ref["sixKFrame"] and buttons["archive"] == ref["sourceFrameSHA256"]
    assert done["frameSHA256"] == proof["frame"] == terminal(stage)["outputSHA256"]["restart.csv"]
    assert done["frameSHA256"] != buttons["archive"] and proof["faces"] == 10272
    for folder in (controls, ui):
        r = read(folder / "debugger_process.json")
        assert r["status"] == "terminal" and r["debuggerExitCode"] == 0 and all(r["inputsStable"].values())
        assert not (r["survivingProcessIDs"] or r["supervisorError"] or r["stopReason"])
    assert proof["pixels"] > 0 and proof["webglError"] == 0
    spec = read(stage / "case.json")
    body = read(stage / "body.json")
    tri = np.array(body["V"])[np.array(body["F"]) - 1]
    def areas(faces):
        return np.linalg.norm(np.cross(faces[:, 1]-faces[:, 0], faces[:, 2]-faces[:, 0]), axis=1)/2
    area = areas(tri)
    with np.load(stage / "nodes.npz") as z:
        xyz, wet, first, parent, child = [z[k] for k in ("xyz","wallMap","originalFirstMap","overlapParent","overlapChild")]
    small = xyz[wet] * spec["diameterM"]
    native_area = areas(small)
    weights = independent_overlap(tri, small, parent, child) * native_area[child, None]
    def project(values):
        return np.bincount(parent, (weights * values[wet[child]]).sum(axis=1), minlength=len(area)) / area
    indices = np.unique(wet)
    lump = np.bincount(wet.ravel(), np.repeat(native_area/3, 3), minlength=len(xyz))
    with (stage / "reaction0/restart.csv").open(encoding="ascii") as stream:
        keys = next(csv.reader(stream))
        rhs = np.loadtxt(stream, delimiter=",", usecols=keys.index("Residual_Species_0"))
    density = np.zeros(len(xyz))
    density[indices] = rhs[indices] * spec["rho"]*3*spec["cp"]*24*spec["diameterM"]**2 / lump[indices]
    flow = np.loadtxt(stage / "restart.csv", delimiter=",", skiprows=1, usecols=(4,9))
    q = project(density)
    p = project(spec["rho"]*9*flow[:,0]) + 101325
    temp = 293.15 + 24*flow[first,1].mean(axis=1)
    def equal(a, b, tol):
        np.testing.assert_allclose(a, b, rtol=0, atol=tol)
    # 科研出口的单位以及给定壁温、共同首层流体温度、未求固体温度身份。
    wall = np.genfromtxt(stage / "surface/wall.csv", names=True, delimiter=",")
    np.testing.assert_array_equal(wall["faceID"], np.arange(len(q)))
    equal(q, wall["heatFluxWm2"], 1e-8);equal(p, wall["pressurePa"], 1e-9);equal(temp, wall["nearWallCellTemperatureK"], 1e-11)
    raw = (ui / "workbench/case-data.js").read_text("utf8").strip()
    head = "window.ROCKET_CASE="
    assert raw.startswith(head) and raw.endswith(";")
    data = json.loads(raw[len(head):-1])
    m = data["measuredPoint"]
    assert "独立冷启动" in data["resultCondition"]["caption"]
    assert m["coldSource"]["frame"] == proof["frame"] and m["coldSource"]["initialSHA256"] == done["initialSHA256"]
    actual = read(job / "comparison.json")
    assert m["coldComparison"]["frame"] == actual["frame"] and m["jointVerifiedM2"] == 0
    for name, refdata in actual["references"].items():
        shown = m["coldComparison"]["references"][name]
        assert shown["referenceFrame"] == refdata["referenceFrame"] and shown["withinTolerance"] == refdata["withinTolerance"]
        for field in ("pressurePa","fluxWm2","commonFirstK","relativeFlux"):
            equal([shown[field]["maxAbs"],shown[field]["areaRMS"]],[refdata[field]["maxAbs"],refdata[field]["areaRMS"]],1e-12)
            assert shown[field]["faceID"] == refdata[field]["faceID"]
    fields = {f["key"]:f["values"] for f in data["fields"]}
    for k, values in (("pressure",p/1000),("pressureGauge",p-101325),("heatFlux",q/1000),
                      ("nearWallCellTemperature",temp),("wallInput",305.15*np.ones(len(q)))):
        equal(fields[k], values, 1e-9)
    assert "wallSolved" not in fields or all(v is None for v in fields["wallSolved"])
    receipt = read(ui / "heat.json")
    assert receipt["measuredPoint"] == m and receipt["csvSHA256"] == sha(ui / "heat.csv")
    for filename in ("wall.csv","heat.csv"):
        with (ui / filename).open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == len(q)
        for i,row in enumerate(rows):
            assert int(row["face_id"]) == i+1
            if filename == "wall.csv":
                assert row["reference_frame_sha256"] == proof["frame"] and float(row["temperature_difference_ratio"]) == 1
                assert float(row["reference_wall_K"]) == 305.15 and row["experiment_doi"] == m["experiment"]["doi"]
                assert abs(1000*float(row["pressure"])-p[i]) < 1e-8 and abs(1000*float(row["heat_flux"])-q[i]) < 1e-7
            else:
                assert row["unit"] == "kW/m²" and row["qualification_code"] == "0" and row["availability_code"] == "2"
                assert abs(float(row["value"])-q[i]/1000) < 1e-9
    segments = cut_values(ui / "section.csv", tri, proof["sectionXM"], q/1000, "heatFlux", "kW/m²")
    group = np.array(body["regionID"])
    powers = [-float(area[group==i] @ q[group==i]) for i in range(1,8)]
    equal(powers, [r["outwardW"] for r in m["regions"]], 1e-13)
    equal(sum(powers), s["energy"]["wallReleaseW"], 1e-12)
    equal(sum(powers[2:5])/(12*3.094e-6), proof["h"], 1e-9)
    order = np.argsort(q)
    # 镜像面排序可能只相差浮点舍入，保留按原生原面表确定的实际探针面号。
    assert proof["peakFace"] == s["peakFaceID"] == int(np.argmin(wall["heatFluxWm2"]))+1
    assert abs(q[proof["peakFace"]-1]-min(q)) < 1e-8
    stages = []
    initial = done["initialSHA256"]
    with ZipFile(ROOT / "cases/thermistor_b/windPoint12.zip") as archive:
        native = np.genfromtxt(io.BytesIO(archive.read("wall.csv")), names=True, delimiter=",")
    # 用独立重建的量对照归档，按每个真实面只计一次超差面积。
    dp = p-native["pressurePa"]
    dq = q-native["heatFluxWm2"]
    fraction = dq/native["heatFluxWm2"]
    dt = temp-native["nearWallCellTemperatureK"]
    rejected = dict(pressure=abs(dp)>=.001, relativeHeatFlux=abs(fraction)>=.0005, firstTemperature=abs(dt)>=.001)
    missing = np.logical_or.reduce(list(rejected.values()))
    def coverage(mask):
        return dict(faces=int(mask.sum()),areaM2=float(area[mask].sum()),areaFraction=float(area[mask].sum()/area.sum()))
    agreement = dict(byVariable={k:coverage(v) for k,v in rejected.items()},union=coverage(missing),
                     byRegion=[dict(regionID=i,name=name,**coverage(missing & (group==i))) for i,name in enumerate(body["regionNames"],1)])
    # 两条分支各自已有的 +20 保存帧和对应壁反力，可直接检验再推进同样步数后局部差距是否缩小。
    # 仅作为证据复核的可选输入，实际冷态求解及工作台不读取参照体场。
    paired = None
    if reference_job is not None:
        def additional_fields(folder,source_frame):
            summary=read(folder/"surface/summary.json")
            plus, reaction=folder/"plus20",folder/"reaction20"
            step, released=terminal(plus),terminal(reaction)
            assert summary["temperatureSampling"]==s["temperatureSampling"] and summary["nodeMapSHA256"]==s["nodeMapSHA256"]
            assert step["inputSHA256"][str(plus/"initial_flow.csv")]==summary["sourceFrameSHA256"]==source_frame
            assert step["outputSHA256"]["restart.csv"]==summary["plus20FrameSHA256"]==released["inputSHA256"][str(reaction/"initial_flow.csv")]
            assert released["outputSHA256"]["restart.csv"]==summary["reaction20SHA256"]
            assert step["inputSHA256"][str(plus/"mesh.su2")]==released["inputSHA256"][str(reaction/"mesh.su2")]==s["meshSHA256"]
            assert summary["energy"]["temperatureScaleK"]==24
            with (reaction/"restart.csv").open(encoding="ascii") as stream:
                header=next(csv.reader(stream))
                residual=np.loadtxt(stream,delimiter=",",usecols=header.index("Residual_Species_0"))
            state=np.loadtxt(plus/"restart.csv",delimiter=",",skiprows=1,usecols=(4,9))
            face_density=np.zeros(len(xyz))
            face_density[indices]=residual[indices]*spec["rho"]*3*spec["cp"]*24*spec["diameterM"]**2/lump[indices]
            return project(spec["rho"]*9*state[:,0]),project(face_density),293.15+24*state[first,1].mean(axis=1)
        pc,qc,tc=additional_fields(stage,done["frameSHA256"])
        pr,qr,tr=additional_fields(reference_job,ref["sourceFrameSHA256"])
        mask=(abs(pc-pr)>=.001)|(abs(qc/qr-1)>=.0005)|(abs(tc-tr)>=.001)
        paired=dict(coldFrame=read(stage/"surface/summary.json")["plus20FrameSHA256"],
                    referenceFrame=read(reference_job/"surface/summary.json")["plus20FrameSHA256"],
                    maxPressureDifferencePa=float(max(abs(pc-pr))),maxRelativeFluxDifference=float(max(abs(qc/qr-1))),
                    maxFirstDifferenceK=float(max(abs(tc-tr))),relativePowerDifference=float(abs((area@qc)/(area@qr)-1)),
                    coldExtra20MaxFluxChangeWm2=float(max(abs(qc-q))),referenceExtra20MaxFluxChangeWm2=float(max(abs(qr-native["heatFluxWm2"]))),
                    union=coverage(mask))
    def iteration_difference(source):
        state = np.loadtxt(source, delimiter=",", skiprows=1, usecols=(4,9))
        gauge = project(spec["rho"]*9*state[:,0])
        first_k = 293.15+24*state[first,1].mean(axis=1)
        dp = gauge-native["pressureGaugePa"]
        return dict(maxGaugeDifferenceToV315Pa=float(max(abs(dp))), meanGaugeDifferenceToV315Pa=float(np.average(dp,weights=area)),
                    maxFirstDifferenceToV315K=float(max(abs(first_k-native["nearWallCellTemperatureK"]))),
                    maxGaugeDifferenceToColdFinalPa=float(max(abs(gauge-wall["pressureGaugePa"]))))
    for dest in sorted(job.glob("steady*")):
        a = read(dest / "field_audit.json")
        rec = terminal(dest)
        assert a["initialSHA256"] == initial == rec["inputSHA256"][str(dest / "initial_flow.csv")]
        initial = a["lastFrame"]
        stages.append(dict(stage=dest.name, steps=a["steps"], seconds=a["seconds"], end=initial, residuals=a["last20Residual"],
                           iterationComparison=iteration_difference(dest/"restart.csv"),
                           peakRSS=read(dest / "resources.json")["peakSampledFamilyRSSBytes"]))
    assert initial == done["frameSHA256"]
    evidence.mkdir(parents=True, exist_ok=False)
    with (evidence/"outside_tolerance.csv").open("w",encoding="utf8",newline="") as stream:
        out=csv.writer(stream)
        out.writerow(["faceID_one_based","regionID","areaM2","pressureDifferencePa","heatFluxDifferenceWm2","relativeHeatFluxDifference","commonFirstDifferenceK"])
        out.writerows((i+1,int(group[i]),area[i],dp[i],dq[i],fraction[i],dt[i]) for i in np.flatnonzero(missing))
    write_json(evidence / "verification.json", dict(stages=stages, energy=s["energy"], relativeChange20=s["relativeChange20"],
        residualWindow=s["residualWindow"], plus20Window=s["plus20Window"], peaks=[dict(faceID=int(i)+1,fluxWm2=q[i]) for i in order[:2]],
        wallRangePa=[min(p),max(p)], firstCellRangeK=[min(temp),max(temp)], commonSectionSegments=len(segments), powerByRegionW=powers,
        additional20Comparison=iteration_difference(stage/"plus20/restart.csv"), pageSHA256=sha(ui / "workbench/case-data.js"),
        sectionSHA256=sha(ui / "section.csv"), initialSHA256=done["initialSHA256"],v315NativeTolerance=agreement,
        pairedAdditional20=paired))
    with ZipFile(evidence / "receipts.zip","x",compression=ZIP_DEFLATED) as bundle:
        def save(folder, names, tag):
            for name in names:
                bundle.write(folder / name, tag + "/" + name)
        save(job, ["job.json","cold/cold.json","completed.json","comparison.json",*(p.name for p in job.glob("receipt-*.json"))], "job")
        bundle.write(job.parent / "replay_pause.json", "job/pause_at_replay.json")
        bundle.write(job.parent / "surface_pause.json", "job/pause_before_surface.json")
        for dest in sorted(job.glob("steady*")):
            save(dest, ["engine_request.json","engine_process.json","case.cfg","case.json","field_audit.json","history.csv","resources.json"], dest.name)
        for key in ("reaction0","plus20","reaction20"):
            save(stage/key, ["engine_request.json","engine_process.json","case.cfg","history.csv","resources.json"], key)
        for key in ("control","observed"):
            save(stage/"ledger"/key, ["step_process.json","checked.json"], "ledger/"+key)
        save(stage, ["ledger/energy.json","surface/summary.json","surface/wall.csv"], "result")
        save(ui, ["debugger_process.json","ui.json","section.csv"], "ui")
        save(controls, ["debugger_process.json","ui.json"], "controls")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    for arg in ("job","controls","ui","evidence"):ap.add_argument(arg, type=Path)
    ap.add_argument("--reference-job",type=Path,help="可选，v315 原生 12 K 任务根目录")
    args = ap.parse_args()
    verify(*(p.resolve() for p in (args.job,args.controls,args.ui,args.evidence)),
           reference_job=args.reference_job.resolve() if args.reference_job else None)
