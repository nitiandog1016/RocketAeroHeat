"""由已登记的 M0.30 来流和实测网格间距复算量纲参照，不生成面场。"""
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def read_json(path):
    return json.loads(path.read_bytes())


def assess():
    geometry = read_json(ROOT / "cases/simple_launcher/geometry.json")
    runs = read_json(HERE.parent / "v358/evidence/cases.json")
    case, = [r for r in runs["cases"] if r["mach"] == 0.30]
    if runs["surfaceSHA256"] != geometry["fileSHA256"] or case["altitudeM"] != 5000:
        raise ValueError("飞行工况或几何发生变化")

    # 与 standardAtmosphere.m 相同的干空气常数和 Sutherland 黏度。
    gas_constant = 287.05287
    gamma = 1.4
    temp = case["freeTemperatureK"]
    speed = case["speedMps"]
    density = case["freePressurePa"] / (gas_constant * temp)
    viscosity = 1.716e-5 * (temp / 273.15) ** 1.5 * (273.15 + 110.4) / (temp + 110.4)
    unit_re = density * speed / viscosity
    total_temp = temp * (1 + (gamma - 1) * case["mach"] ** 2 / 2)
    if not math.isclose(speed, case["mach"] * math.sqrt(gamma * gas_constant * temp)):
        raise ValueError("来流速度与登记的气体状态不一致")

    # 采用来流物性和平板剪力做尺度换算，实际壁面剪力及 y+ 须有本工况解。
    mesh = read_json(HERE.parent / "v357/flow_evidence.json")["mesh"]
    median_height = {grid: mesh[grid]["firstNormalMedianMm"] for grid in ("oldResolution", "refinedResolution")}
    reference = read_json(HERE / "plate02.json")["references"]
    estimates = {
        grid: {code: unit_re * height_mm * 1e-3 * math.sqrt(cf / 2) for code, cf in reference.items()}
        for grid, height_mm in median_height.items()
    }
    return {
        "mach": case["mach"], "speedMps": speed, "densityKgM3": density,
        "viscosityPaS": viscosity, "reDiameter": unit_re * geometry["maxDiameterM"],
        "reLength": unit_re * geometry["lengthM"], "totalTemperatureK": total_temp,
        "givenWallMinusTotalK": case["wallGivenK"] - total_temp,
        "tunnelMedianNormalHeightMm": median_height, "flatPlateScaleEstimateOnly": estimates,
        "actualFlightWallYPlus": None, "wholeBodyPowerW": None,
    }


if __name__ == "__main__":
    print(json.dumps(assess(), ensure_ascii=False, indent=2, allow_nan=False))
