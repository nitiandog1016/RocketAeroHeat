"""热敏珠独立冷态的终态审核、七区积分及与两份已发布风洞测点结果的同面比较。"""
import io
import json
from zipfile import ZipFile

import numpy as np

from thermistorBead import ROOT, read, sha, write_json
from thermistorCrossflowSurface import collect


def compare(stage):
    """比较相同物理面的数值，保留迭代位置、两个对称峰和所有区域的实际差值。"""
    case_dir = ROOT / "cases/thermistor_b"
    key = read(case_dir / "windPoint12.json")
    archive = case_dir / "windPoint12.zip"
    assert sha(archive) == key["archiveSHA256"]
    s = read(stage / "surface/summary.json")
    table = np.genfromtxt(stage / "surface/wall.csv", delimiter=",", names=True)
    area, group = table["areaM2"], table["regionID"]
    assert len(area) == 10272 and np.all(area > 0)
    np.testing.assert_array_equal(table["faceID"], np.arange(len(area)))
    results = {}
    with ZipFile(archive) as bundle:
        exp = json.loads(bundle.read("comparison.json"))["experiment"]
        for name, prefix, ratio in (("scaledSix", "six_", 2), ("nativeTwelve", "", 1)):
            ref = np.genfromtxt(io.BytesIO(bundle.read(prefix + "wall.csv")), delimiter=",", names=True)
            summary = json.loads(bundle.read(prefix + "summary.json"))
            for field in ("faceID", "areaM2", "regionID"):
                np.testing.assert_array_equal(ref[field], table[field])
            reference_q = ratio * ref["heatFluxWm2"]
            reference_t = 293.15 + ratio * (ref["nearWallCellTemperatureK"] - 293.15)
            dp = table["pressurePa"] - ref["pressurePa"]
            dt = table["nearWallCellTemperatureK"] - reference_t
            dq = table["heatFluxWm2"] - reference_q
            assert np.all(reference_q != 0)
            regions = []
            for i, region in enumerate(s["regions"], 1):
                mask = group == i
                watts = float(np.sum(-area[mask] * table["heatFluxWm2"][mask]))
                other = float(np.sum(-area[mask] * reference_q[mask]))
                assert abs(watts-region["outwardW"]) < 1e-12
                regions.append(dict(name=region["name"], watts=watts, referenceW=other, differenceW=watts-other))
            def metric(delta):
                idx = int(np.argmax(abs(delta)))
                return dict(maxAbs=float(abs(delta[idx])), faceID=idx+1,
                            areaRMS=float(np.sqrt(np.average(delta**2, weights=area))))
            pressure, heat, temp = metric(dp), metric(dq), metric(dt)
            relative = metric(dq/reference_q)
            power = abs(s["totalOutwardW"]/(ratio*summary["totalOutwardW"])-1)
            results[name] = dict(referenceFrame=summary["sourceFrameSHA256"], referenceSteps=summary["cumulativeSteps"],
                                 pressurePa=pressure, fluxWm2=heat, commonFirstK=temp, relativeFlux=relative,
                                 relativePower=power, regions=regions, referencePeakFaceID=int(np.argmin(reference_q))+1,
                                 referencePeakWm2=float(min(reference_q)),
                                 withinTolerance=bool(pressure["maxAbs"]<.001 and temp["maxAbs"]<.001 and
                                                      relative["maxAbs"]<.0005 and power<.0005))
    h = s["beadAndJointsW"] / (12 * exp["areaM2"])
    return dict(frame=s["sourceFrameSHA256"], coldIterations=s["cumulativeSteps"], references=results,
                interpretation="Independent cold start. The reference branches started from a developed field; iteration positions differ.",
                peakFaceID=s["peakFaceID"], peakWm2=float(min(table["heatFluxWm2"])), wholeW=s["totalOutwardW"],
                beadW=s["beadAndJointsW"], beadHWm2K=h, experimentHWm2K=exp["hWm2K"], differenceFromExperiment=h/exp["hWm2K"]-1,
                faces=len(area), areaM2=float(area.sum()), jointExperimentVerifiedM2=0)


def finish(run, stage, job, checkpoint):
    """独立重读每个回执后才发布完成标记；再次打开同一目录不会重复原生推进。"""
    collect(stage, checkpoint=checkpoint)
    comparison = compare(stage)
    path = run / "comparison.json"
    if path.exists():
        assert read(path) == comparison
    else:
        write_json(path, comparison)
    s = read(stage / "surface/summary.json")
    record = dict(token=job["token"], stage=stage.name, jobSHA256=sha(run / "job.json"),
                  initialSHA256=read(run / "cold/cold.json")["initial"]["sha256"], frameSHA256=s["sourceFrameSHA256"],
                  reactionSHA256=s["sourceReactionSHA256"], caseSHA256=s["caseSHA256"], wallSHA256=s["wallSHA256"],
                  meshSHA256=s["meshSHA256"], bodySHA256=s["bodySHA256"], nodeMapSHA256=s["nodeMapSHA256"],
                  summarySHA256=sha(stage / "surface/summary.json"), comparisonSHA256=sha(path))
    done = run / "completed.json"
    if done.exists():
        assert read(done) == record
    else:
        write_json(done, record)
