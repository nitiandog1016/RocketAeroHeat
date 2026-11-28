"""复核 NASA 等温平板表 I(a)，比较工程公式的定物性、已转捩极限。"""
import argparse
import csv
import json
from math import isclose, isfinite, sqrt
from pathlib import Path
from statistics import mean


HERE = Path(__file__).resolve().parent
BTU_IT_J = 1055.05585262
FT_M = 0.3048


def compare(table):
    """保留每个测站的两列 St，按原作者给出的来流物性 Re_x 比较。"""
    if any(len(row) != len(table["columns"]) for row in table["rows"]):
        raise ValueError("原表列数与表头不符")
    rows = [dict(zip(table["columns"], row)) for row in table["rows"]]
    if [row["strip"] for row in rows] != list(range(2, 24)):
        raise ValueError("表 I(a) 最高速度档必须包含铜带 2 至 23")
    flux_factor = BTU_IT_J / (3600 * FT_M ** 2)
    free_k = (table["freeTemperatureF"] - 32) * 5 / 9 + 273.15
    stations = []
    for row in rows:
        re_x = row["reXMillions"] * 1e6
        st = row["stTimes1000"] / 1000
        corrected = row["propertyCorrectedStTimes1000"] / 1000
        delta_k = row["meanWallMinusFreeF"] * 5 / 9
        heat_out = flux_factor * row["outwardFluxBtuPerHrFt2"]
        h = flux_factor * 1.8 * row["hBtuPerHrFt2F"]
        # 此处只重算 calcSurfaceHeating.m 中局部湍流公式的定物性极限。
        # 试件在第一条铜带已转捩；不调用整机转捩判定，也不拟合常数。
        turbulent_st = 0.0592 / re_x ** 0.2 / (2 * 0.72 ** (2 / 3))
        stations.append({
            "strip": row["strip"], "reX": re_x, "stMeasured": st,
            "stPropertyCorrected": corrected, "stConstantPropertyFormula": turbulent_st,
            "errorVsCorrectedPercent": 100 * (turbulent_st / corrected - 1),
            "errorVsUncorrectedPercent": 100 * (turbulent_st / st - 1),
            "meanWallK": free_k + delta_k, "meanWallMinusFreeK": delta_k,
            "measuredQIntoWallWm2": -heat_out, "measuredHWm2K": h,
            "printedFluxDifferencePercent": 100 * (h * delta_k / heat_out - 1),
        })
    errors = [s["errorVsCorrectedPercent"] for s in stations]
    return {
        "speedMps": table["nominalSpeedFtPerS"] * FT_M, "freeTemperatureK": free_k,
        "prandtl": 0.72, "stationCount": len(stations),
        "meanSignedErrorVsCorrectedPercent": mean(errors),
        "rmsRelativeErrorVsCorrectedPercent": sqrt(mean(err ** 2 for err in errors)),
        "meanSignedErrorVsUncorrectedPercent": mean(s["errorVsUncorrectedPercent"] for s in stations),
        "maximumAbsoluteErrorVsCorrectedPercent": max(map(abs, errors)),
        "worstStrip": max(stations, key=lambda s: abs(s["errorVsCorrectedPercent"]))["strip"],
        "maxPrintedFluxDifferencePercent": max(abs(s["printedFluxDifferencePercent"]) for s in stations),
        "stations": stations,
    }


def check_heat_readback(result, csv_path):
    """按铜带号和原始量纲独立核对 MATLAB 的逐行有符号热流及相对偏差。"""
    with csv_path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != result["stationCount"]:
        raise ValueError("MATLAB 输出的铜带数与原表不一致")
    errors = []
    for source, raw in zip(result["stations"], rows):
        row = {k: float(v) for k, v in raw.items()}
        if not all(map(isfinite, row.values())) or source["strip"] != row["strip"]:
            raise ValueError("逐站读回的数值或铜带号不符")
        for measured, field in (("reX", "reX"), ("meanWallK", "givenWallK"),
                                ("measuredQIntoWallWm2", "measuredQIntoWallWm2")):
            if not isclose(source[measured], row[field], rel_tol=1e-11):
                raise ValueError(f"第 {source['strip']} 条铜带的 {field} 与原表不符")
        estimated = row["estimatedQIntoWallWm2"]
        error = 100 * (estimated / source["measuredQIntoWallWm2"] - 1)
        if estimated >= 0 or not isclose(error, row["relativeErrorPercent"], abs_tol=1e-9):
            raise ValueError("入壁正向或热流相对偏差不一致")
        source["matlabEstimatedQIntoWallWm2"] = estimated
        errors.append(error)
    return {
        "meanSignedRelativeErrorPercent": mean(errors),
        "rmsRelativeErrorPercent": sqrt(mean(err ** 2 for err in errors)),
        "minimumRelativeErrorPercent": min(errors),
        "maximumAbsoluteRelativeErrorPercent": max(map(abs, errors)),
        "worstStrip": result["stations"][max(range(len(errors)), key=lambda i: abs(errors[i]))]["strip"],
        "sameHeatFluxDirectionStations": sum(s["matlabEstimatedQIntoWallWm2"] * s["measuredQIntoWallWm2"] > 0
                                              for s in result["stations"]),
    }


def plot(result, path):
    """画出两种对照；热流保留入壁正向，横坐标沿用原表的来流物性定义。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stations = result["stations"]
    x = [s["reX"] / 1e6 for s in stations]
    fig, (heat, ax) = plt.subplots(2, 1, figsize=(7.2, 6.8), sharex=True, constrained_layout=True)
    heat.plot(x, [s["matlabEstimatedQIntoWallWm2"] for s in stations],
              color="#225883", label="Local formula with project air properties")
    heat.scatter(x, [s["measuredQIntoWallWm2"] for s in stations],
                 marker="o", s=24, facecolors="none", edgecolors="#ab492c", label="Experiment")
    heat.set(title="Heated flat plate at 38.71 m/s", ylabel=r"Heat flux into wall (W/m$^2$)")
    ax.plot(x, [1000*s["stConstantPropertyFormula"] for s in stations],
            color="#225883", label="Local turbulent formula (constant properties)")
    ax.scatter(x, [1000*s["stMeasured"] for s in stations], s=22, marker="x",
               color="#716b63", label="Measured St")
    ax.scatter(x, [1000*s["stPropertyCorrected"] for s in stations], s=24,
               facecolors="none", edgecolors="#ab492c", label="Authors' property-corrected St")
    ax.set(xlabel=r"Local Reynolds number $Re_x$ (millions)", ylabel=r"Local Stanton number $10^3 St$")
    for panel in (heat, ax):
        panel.grid(axis="y", alpha=0.2)
        panel.legend(frameon=False, fontsize=9)
    fig.savefig(path, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", type=Path)
    parser.add_argument("--heat-readback", type=Path, default=HERE / "plate_heat.csv")
    args = parser.parse_args()
    result = compare(json.loads((HERE / "table_ia_127.json").read_bytes()))
    result["heatFluxComparison"] = check_heat_readback(result, args.heat_readback)
    print(json.dumps(result, indent=2, allow_nan=False))
    if args.figure:
        plot(result, args.figure)
