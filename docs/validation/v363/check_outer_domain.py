"""由两份已封存的同形风洞面表，复算外域从 14D 至 22D 时的逐区差异。"""
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASE = ROOT / "cases/simple_launcher"


def compare():
    fine_path = ROOT/"outputs/v350_validation/adapt_fine02/wall.csv"
    geom = json.loads((CASE/"geometry.json").read_bytes())
    bc = json.loads((CASE/"tunnel.json").read_bytes())
    with zipfile.ZipFile(CASE/"wind_tunnel.zip") as bundle:
        source = json.loads(bundle.read("summary.json"))
        outer_bytes = bundle.read("wall.csv")
    inner_bytes = fine_path.read_bytes()
    assert source["accepted"] and source["condition"]==bc
    assert source["surfaceSHA256"]==geom["fileSHA256"]
    input_sha = {key:hashlib.sha256(raw).hexdigest() for key,raw in (("fine",inner_bytes),("outer",outer_bytes))}
    assert all(input_sha[k]==source["runs"][k]["wallSHA256"] for k in input_sha)
    inner,outer = (np.genfromtxt(io.BytesIO(raw),delimiter=",",names=True) for raw in (inner_bytes,outer_bytes))
    np.testing.assert_array_equal(inner["faceID"],np.arange(1,geom["faces"]+1))
    np.testing.assert_array_equal(inner["faceID"],outer["faceID"])
    np.testing.assert_array_equal(inner["areaM2"],outer["areaM2"])
    speed = float(np.linalg.norm(bc["velocityMps"]))
    dynamic_pressure = .5*bc["densityKgM3"]*speed**2
    delta_p = inner["pressurePa"]-outer["pressurePa"]
    regions = {}
    squared_error = 0.
    for r in geom["regions"]:
        key = r["name"]
        idx = slice(r["firstFace"]-1,r["lastFace"])
        area = outer["areaM2"][idx]
        assert abs(area.sum()-r["areaM2"])<1e-12
        pressure_rms = np.sqrt(np.average(delta_p[idx]**2,weights=area))/dynamic_pressure
        heat14,heat22 = (-float(area@data["heatFluxWm2"][idx]) for data in (inner,outer))
        power_diff = abs(heat14/heat22-1)
        np.testing.assert_allclose(power_diff,source["outerDomain"]["regionPowerDifferences"][key],atol=1e-13,rtol=0)
        regions[key] = dict(pressureRmsPercentDynamic=100*pressure_rms,
                            outwardPower14DW=heat14,outwardPower22DW=heat22,powerDifferencePercent=100*power_diff)
        squared_error += area.sum()*pressure_rms**2
    overall = np.sqrt(squared_error/geom["areaM2"])
    np.testing.assert_allclose(overall,source["outerDomain"]["pressureRmsOverDynamic"],atol=1e-13,rtol=0)
    return dict(surfaceSHA256=geom["fileSHA256"],speedMps=speed,wallGivenK=bc["wallTemperatureK"],
                dynamicPressurePa=dynamic_pressure,pressureRmsPercentDynamic=100*overall,
                wallTableSHA256=input_sha,regions=regions)


if __name__=="__main__":
    print(json.dumps(compare(),indent=2,allow_nan=False))
