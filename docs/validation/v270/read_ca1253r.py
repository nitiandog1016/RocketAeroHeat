"""从 NASA 原版扫描图复读三个可分辨的仪器散点。时间原点与仪器工况见 RESULTS.md。"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from PIL import Image
from scipy.ndimage import label


REPORT_SHA256 = "b0826991798736328c54549f9baf3107b86b28f207dc77186f2241709b36a799"
# 在同一张扫描图、紧邻测点的网格标定，左上角为像素原点。
TIME_AXIS = ((776.5, 4500), (834.0, 4540))
FLUX_AXIS = ((1859.5, 1.5), (1717.5, 2.0))
MARKS = ((795, 1815), (822, 1781), (822, 1791))
NOMINAL_TIME = 4520
SI_FACTOR = 1055.05585262 / 0.3048**2  # 国际蒸汽表 Btu；图中面积单位为 ft²。


def scale(pixel, axis, offset1=0, offset2=0):
    (p1, v1), (p2, v2) = axis
    return v1 + (pixel-p1-offset1)*(v2-v1)/(p2+offset2-p1-offset1)


def reading_bounds(pixels, axis):
    # 读图范围包括整个墨点以及两根扫描网格线各 ±3 像素的定位差。
    values = [scale(x, axis, d1, d2) for x, d1, d2 in itertools.product(
        (pixels.min()-.5, pixels.max()+.5), (-3, 3), (-3, 3))]
    return [min(values), max(values)]


def read_mark(gray, x, y):
    patch = gray[y-5:y+6, x-5:x+6] < 128
    ink, _ = label(patch, np.ones((3, 3)))
    mark = ink[5, 5]
    if not mark:
        raise ValueError("指定位置没有量热计散点。")
    row, col = np.where(ink == mark)
    if min(row.min(), col.min()) == 0 or max(row.max(), col.max()) == 10:
        raise ValueError("选定墨点与扫描网格或曲线相连，无法独立读数。")
    px, py = col+x-5, row+y-5
    time_s = scale(px.mean(), TIME_AXIS)
    flux = scale(py.mean(), FLUX_AXIS)
    return dict(centerPixel=[px.mean(), py.mean()], elapsedSeconds=time_s,
                timeReadingBoundsS=reading_bounds(px, TIME_AXIS),
                heatFluxBtuITPerFt2S=flux, heatFluxWm2=flux*SI_FACTOR,
                heatFluxReadingBoundsWm2=[q*SI_FACTOR for q in reading_bounds(py, FLUX_AXIS)])


def annotate(gray, points, destination):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    left, right, top, bottom = 728, 890, 1640, 1905
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    extent = (scale(left-.5, TIME_AXIS), scale(right-.5, TIME_AXIS),
              scale(bottom-.5, FLUX_AXIS), scale(top-.5, FLUX_AXIS))
    ax.imshow(gray[top:bottom, left:right], cmap="gray", vmin=0, vmax=255,
              origin="upper", extent=extent, aspect="auto", interpolation="nearest")
    ax.axvline(NOMINAL_TIME, color="#0072b2", ls="--", lw=1.2, label="Nominal comparison time")
    ax.scatter([p["elapsedSeconds"] for p in points], [p["heatFluxBtuITPerFt2S"] for p in points],
               s=115, facecolors="none", edgecolors="#bf4200", lw=1.5, label="Selected flight marks")
    ax.set(xlabel="Seconds from range zero", ylabel="Gauge heat flux [Btu/(ft² s)]",
           title="Apollo AS-202  ·  CA1253R")
    ax.legend(loc="upper right", framealpha=1, fontsize=8)
    fig.savefig(destination, dpi=160)
    plt.close(fig)


def main(pdf, out):
    if hashlib.sha256(pdf.read_bytes()).hexdigest() != REPORT_SHA256:
        raise ValueError("任务报告与审核版本不一致。")
    out.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory(dir=out) as temp:
        stem = Path(temp)/"page"
        subprocess.run(["pdftoppm", "-f", "172", "-l", "172", "-scale-to", "3200",
                        "-gray", "-singlefile", "-png", str(pdf), str(stem)], check=True, timeout=60)
        with Image.open(stem.with_suffix(".png")) as raster:
            if raster.size != (2465, 3200):
                raise ValueError("原扫描页的像素坐标不匹配。")
            gray = np.array(raster.convert("L"))
    points = [read_mark(gray, x, y) for x, y in MARKS]
    if not (points[0]["timeReadingBoundsS"][1] < NOMINAL_TIME <
            min(p["timeReadingBoundsS"][0] for p in points[1:])):
        raise ValueError("参考时刻与两侧读图范围不符。")
    record = dict(mission="Apollo AS-202", spacecraft="011", sensor="CA1253R", sensorLocationLabel="LOC 2",
                  surface="command module conical surface", instrument="asymptotic calorimeter",
                  observationSource=dict(report="MSC-A-R-66-5", printedPage="7-68", pdfPage=172,
                                         figure="7.3-5(a)", sha256=REPORT_SHA256),
                  clock="range zero (mission report notice, PDF p. 6)", nominalTimeS=NOMINAL_TIME,
                  exactTimeSampleWm2=None, graphSamples=points, conversionToWm2=SI_FACTOR,
                  uncertaintyScope="reading bounds include mark width and grid placement only",
                  experimentalUncertaintyWm2=None, sameStationPressurePa=None,
                  sensorFaceTemperatureK=None, adjacentSkinTemperatureK=None,
                  simulatedPressurePa=None, simulatedConvectiveHeatFluxWm2=None,
                  boundSurfaceFaceID=None, newComputedAreaM2=0, newQualifiedAreaM2=0)
    (out/"station.json").write_text(json.dumps(record, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    annotate(gray, points, out/"readings.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mission_pdf", type=Path)
    parser.add_argument("new_output_directory", type=Path)
    args = parser.parse_args()
    main(args.mission_pdf, args.new_output_directory)
