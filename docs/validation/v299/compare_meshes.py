"""核对两套混合网的同面结果，另按历史三角面守恒对照四个轴向带。"""
import io
import json
import zipfile

import numpy as np

from check_grid import HERE, NEW, OLD, ROOT, read, save, sha
from compare import project_to_original
from thermistorBead import A, AREA, DT


def table(raw):
    return np.genfromtxt(io.BytesIO(raw), names=True, delimiter=",")


def axial_position(body):
    return np.asarray(body["V"])[np.asarray(body["F"])-1, 0].mean(axis=1)/A


def band_powers(x, area, q_fields):
    bands = []
    for lo, hi in zip((-1, -.5, 0, .5), (-.5, 0, .5, 1)):
        mask = (x >= lo) & (x < hi)
        bands.append(dict(xOverA=[lo, hi], faces=int(mask.sum()), areaM2=float(area[mask].sum()),
                          outwardW={label:float(-area[mask] @ q[mask]) for label, q in q_fields.items()}))
    assert sum(b["faces"] for b in bands) == len(area)
    return bands


def audited_field(folder):
    surf = folder/"surface"
    s, flow = read(surf/"summary.json"), read(surf/"flow.json")
    assert sha(folder/"restart.csv") == s["sourceFrameSHA256"] == flow["sourceFrameSHA256"]
    assert sha(surf/"wall.csv") == s["wallSHA256"] and sha(surf/"flow.csv") == flow["tableSHA256"]
    w, f = [table((surf/name).read_bytes()) for name in ("wall.csv", "flow.csv")]
    assert np.array_equal(w["faceID"], np.arange(len(w))) and np.array_equal(w["faceID"]+1, f["face_id"])
    np.testing.assert_allclose(-w["areaM2"] @ w["heatFluxWm2"], s["totalOutwardW"], atol=1e-13, rtol=0)
    assert s["relativeChange20"] < .002 and abs(s["energy"]["relativeDifference"]) < .005
    return s, w, f


def main():
    proof = HERE/"evidence"
    grid = read(proof/"inputs.json")
    assert all(sha(NEW/k) == v for k, v in grid["inputs"].items())
    old, base, old_flow = audited_field(OLD)
    new, refined, flow = audited_field(NEW)
    assert old == read(HERE.parent/"v298/evidence/surface.json")
    assert old["bodySHA256"] == new["bodySHA256"] == sha(NEW/"body.json")
    np.testing.assert_array_equal(base["areaM2"], refined["areaM2"])
    area = base["areaM2"]
    body = read(NEW/"body.json")
    x = axial_position(body)
    native = band_powers(x, area, {"baseMixed":base["heatFluxWm2"], "refined":refined["heatFluxWm2"]})
    cases = ROOT/"cases/thermistor_b"
    with zipfile.ZipFile(cases/"surface.zip") as z:
        original = json.loads(z.read("body.json"))
    with zipfile.ZipFile(cases/"surface_7p5.zip") as z:
        structured = table(z.read("wall.csv"))
    with zipfile.ZipFile(cases/"paired3.zip") as z:
        mapped = table(z.read("historical.csv"))
    a = structured["areaM2"]
    np.testing.assert_array_equal(a, mapped["areaM2"])
    _, q, _, _ = project_to_original(NEW, NEW/"reaction0", original, a, True)
    assert abs(-a@q-new["totalOutwardW"]) < 1e-13
    historic = band_powers(axial_position(original), a, dict(structured=structured["heatFluxWm2"],
                         baseMixed=mapped["pairedHeatFluxWm2"], refined=q))
    prior = read(HERE.parent/"v298/evidence/spatial.json")["smoothPowerChangeByAxialBand"]
    for band, previous in zip(historic, prior):
        w = band["outwardW"]
        assert abs(w["baseMixed"]-w["structured"]-previous["extraOutwardW"]) < 1e-13
    mid = [b["outwardW"] for b in historic[1:3]]
    recovered = sum(w["baseMixed"]-w["refined"] for w in mid)
    original_gap = sum(w["baseMixed"]-w["structured"] for w in mid)
    reverse = flow["reverse_both"].astype(bool)
    changed = reverse != old_flow["reverse_both"].astype(bool)
    peaks = {}
    for label, w, f in (("baseMixed", base, old_flow), ("refined", refined, flow)):
        j = int(np.argmin(w["heatFluxWm2"]))
        peaks[label] = dict(faceID=j+1, xOverA=float(x[j]), meridianStrip=int(f["meridian_strip_0"][j]),
                            intoWallWm2=float(w["heatFluxWm2"][j]))
    result = dict(nativeBands=native, historicalBands=historic, middleOriginalGapW=original_gap,
                  middleReductionW=recovered, middleFractionRecovered=recovered/original_gap, peaks=peaks,
                  reverseChangedFaces=int(changed.sum()), reverseChangedAreaM2=float(area[changed].sum()),
                  oldReverseAreaM2=float(area[old_flow["reverse_both"]>0].sum()), newReverseAreaM2=float(area[reverse].sum()),
                  outwardChangeW=new["totalOutwardW"]-old["totalOutwardW"],
                  meanHWm2K=new["totalOutwardW"]/(AREA*DT), relativeExperimentError=new["beadMeanHWm2K"]/new["experimentMeanHWm2K"]-1,
                  baselineFrame=old["sourceFrameSHA256"], refinedFrame=new["sourceFrameSHA256"])
    save(proof/"comparison.json", result)
    save(proof/"surface.json", new)
    save(proof/"flow.json", read(NEW/"surface/flow.json"))
    # 随源码交付轻量原面归档，完整三维场保存在独立求解目录。
    archive = cases/"midbody3.zip"
    with zipfile.ZipFile(archive, "x") as z:
        for path in (NEW/"case.json", NEW/"body.json", *(NEW/"surface"/name for name in ("wall.csv","summary.json","flow.csv","flow.json"))):
            z.writestr(zipfile.ZipInfo(path.name, date_time=(2026, 1, 1, 0, 0, 0)), path.read_bytes(),
                       compress_type=zipfile.ZIP_DEFLATED)
    contract = {k:grid["geometry"][k] for k in ("nodes","cells","meshSHA256","bodySHA256","nodeMapSHA256")}
    contract.update(initialSHA256=sha(NEW/"initial_flow.csv"), configSHA256=sha(NEW/"case.cfg"),
                    caseSHA256=sha(NEW/"case.json"), archiveSHA256=sha(archive))
    save(cases/"midbody3.json", contract)


if __name__ == "__main__":
    main()
