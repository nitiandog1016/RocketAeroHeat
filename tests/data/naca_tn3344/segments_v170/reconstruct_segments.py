"""按原图明确尺寸恢复名义电极段；保留旧图读数及未闭合的中心身份。"""

from pathlib import Path
import csv
import hashlib
import json
import math
import numpy as np


ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
DATA = OUT.parent
INCH = 0.0254
RADIUS = 2 * INCH
DIAMETER = 4 * INCH
SHOULDER = math.pi * RADIUS / 2
POSITION_TOLERANCE = 0.005 * INCH


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def area_primitive(s):
    """完整圆周的名义外表面积原函数，鼻尖为零。"""
    if s <= SHOULDER:
        return 2 * math.pi * RADIUS**2 * (1 - math.cos(s / RADIUS))
    return 2 * math.pi * RADIUS**2 + 2 * math.pi * RADIUS * (s - SHOULDER)


def independent_area(a, b):
    """独立高斯积分使用周长函数，并在球柱接点拆分。"""
    x, w = np.polynomial.legendre.leggauss(48)
    split = [a] + ([SHOULDER] if a < SHOULDER < b else []) + [b]
    result = 0.0
    for left, right in zip(split[:-1], split[1:]):
        s = (left + right) / 2 + (right - left) / 2 * x
        radius = np.where(s <= SHOULDER, RADIUS * np.sin(s / RADIUS), RADIUS)
        result += float(np.dot(w, 2 * math.pi * radius)) * (right - left) / 2
    return result


def write_csv(name, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


# 上下仅对应原图纸面；两套经线在圆周上相差180度。
stations = []
for j in range(16):
    theta = 6 * j
    stations.append({
        "station_id": j,
        "meridian": "upper" if j % 2 else "lower",
        "geometry": "hemisphere" if theta < 90 else "shoulder",
        "theta_deg": theta,
        "cylinder_offset_inch": 0.0,
        "s_m": RADIUS * math.radians(theta),
        "dimension_source": "Figure1(c): upper6+12n_deg; lower12n_deg",
    })
for offset, meridian in [(13, "lower"), (27, "upper"), (40, "lower"),
                         (54, "upper"), (67, "lower"), (81, "upper")]:
    stations.append({
        "station_id": len(stations),
        "meridian": meridian,
        "geometry": "cylinder",
        "theta_deg": "",
        "cylinder_offset_inch": offset / 64,
        "s_m": SHOULDER + offset / 64 * INCH,
        "dimension_source": "Figure1(c): upper27n/64_in; lower(13+27n)/64_in",
    })
for row in stations:
    row["s_over_D"] = row["s_m"] / DIAMETER
    row["station_position_tolerance_m"] = POSITION_TOLERANCE
    row["source_pdf_page"] = 43
    row["source_printed_page"] = 42

source_paths = [
    ROOT / "outputs/v168_validation/sources/NACA_TN_3344.pdf",
    ROOT / "outputs/v168_validation/sources/NACA_TN_3344.txt",
    DATA / "raw/figure1c.png", DATA / "raw/figure4b.png",
    DATA / "heat_parameter.csv", DATA / "digitization.json",
]
input_hashes = {str(p.relative_to(ROOT)).replace("\\", "/"): digest(p) for p in source_paths}
with (DATA / "heat_parameter.csv").open(encoding="utf-8-sig", newline="") as stream:
    heat = list(csv.DictReader(stream))
assert len(heat) == 20 and len(stations) == 22
segments = []
for j, original in enumerate(heat, 1):
    left, center, right = stations[j-1:j+2]
    assert left["meridian"] == right["meridian"] != center["meridian"]
    a, c, b = left["s_m"], center["s_m"], right["s_m"]
    assert a < c < b
    midpoint = (a+b)/2
    area = area_primitive(b) - area_primitive(a)
    area_low = area_primitive(b-POSITION_TOLERANCE) - area_primitive(a+POSITION_TOLERANCE)
    area_high = area_primitive(b+POSITION_TOLERANCE) - area_primitive(max(0, a-POSITION_TOLERANCE))
    sx = float(original["s_over_D"])
    tol = float(original["reading_uncertainty_s_over_D"])
    error_tc = sx - c/DIAMETER
    error_mid = sx - midpoint/DIAMETER
    matches = [k for k in range(1, 21)
               if abs(sx-stations[k]["s_over_D"]) < abs(error_tc)+1e-14]
    assert matches == [j]
    segments.append({
        "heat_row_1based": j,
        "left_voltage_station_id": j-1,
        "right_voltage_station_id": j+1,
        "voltage_meridian": left["meridian"],
        "opposite_temperature_station_id": j,
        "temperature_meridian": center["meridian"],
        "s_left_m": a, "s_right_m": b,
        "s_left_over_D": a/DIAMETER, "s_right_over_D": b/DIAMETER,
        "temperature_station_s_over_D": c/DIAMETER,
        "voltage_interval_arc_midpoint_s_over_D": midpoint/DIAMETER,
        "center_identity_same": abs(c-midpoint) < 1e-14,
        "area_full_circumference_m2": area,
        "area_quarter_domain_m2": area/4,
        "area_min_with_station_tolerance_m2": area_low,
        "area_max_with_station_tolerance_m2": area_high,
        "area_max_relative_station_sensitivity": max(area-area_low, area_high-area)/area,
        "independent_area_relative_error": abs(independent_area(a,b)-area)/area,
        "original_plot_s_over_D": sx,
        "original_Nu_s_over_sqrt_Re_s": float(original["Nu_s_over_sqrt_Re_s"]),
        "original_reading_uncertainty_s_over_D": tol,
        "original_reading_uncertainty_parameter": float(original["reading_uncertainty_parameter"]),
        "original_source_pixel_x": int(original["source_pixel_x"]),
        "original_source_pixel_y": int(original["source_pixel_y"]),
        "plot_minus_temperature_station_s_over_D": error_tc,
        "plot_minus_interval_midpoint_s_over_D": error_mid,
        "within_plot_reading_bound_at_temperature_station": abs(error_tc) <= tol,
        "within_plot_plus_station_bounds_at_temperature_station": abs(error_tc) <= tol+POSITION_TOLERANCE/DIAMETER,
        "within_plot_plus_station_bounds_at_interval_midpoint": abs(error_mid) <= tol+POSITION_TOLERANCE/DIAMETER,
        "nearest_temperature_station_ordinal_unique": True,
        "strict_local_abscissa_identity_closed": abs(c-midpoint) < 1e-14 and abs(error_tc) <= tol,
    })

assert all(input_hashes[str(p.relative_to(ROOT)).replace("\\", "/")] == digest(p) for p in source_paths)
write_csv("nominal_stations.csv", stations)
write_csv("nominal_segments.csv", segments)
summary = {
    "schema": "RocketAeroHeat.naca-nominal-electrode-segments.v170",
    "source_document": "NACA TN3344, printed pages10/11/13/42/48; PDF pages11/12/14/43/49",
    "selection": "M_inf1.97, Re_D2.28e6, Figure4(b) circle-with-center-dot",
    "nominal_geometric_station_count": len(stations),
    "nominal_voltage_segment_count": len(segments),
    "nominal_segment_geometry_identified": True,
    "nominal_segments_overlap": True,
    "strict_all20_local_abscissa_identity_closed": False,
    "temperature_station_and_interval_midpoint_identical_rows": [r["heat_row_1based"] for r in segments if r["center_identity_same"]],
    "plot_reading_bound_mismatch_rows_at_temperature_station": [r["heat_row_1based"] for r in segments if not r["within_plot_reading_bound_at_temperature_station"]],
    "plot_plus_station_bound_mismatch_rows_at_temperature_station": [r["heat_row_1based"] for r in segments if not r["within_plot_plus_station_bounds_at_temperature_station"]],
    "plot_plus_station_bound_mismatch_rows_at_interval_midpoint": [r["heat_row_1based"] for r in segments if not r["within_plot_plus_station_bounds_at_interval_midpoint"]],
    "sphere_radius_m": RADIUS,
    "diameter_m": DIAMETER,
    "shoulder_s_over_D": SHOULDER/DIAMETER,
    "covered_s_over_D": [0.0, stations[-1]["s_over_D"]],
    "union_area_full_circumference_m2": area_primitive(stations[-1]["s_m"]),
    "sum_overlapping_segment_areas_m2_not_union": sum(r["area_full_circumference_m2"] for r in segments),
    "maximum_independent_area_relative_error": max(r["independent_area_relative_error"] for r in segments),
    "maximum_station_position_area_relative_sensitivity": max(r["area_max_relative_station_sensitivity"] for r in segments),
    "station_position_tolerance_m": POSITION_TOLERANCE,
    "station_position_tolerance_over_D": POSITION_TOLERANCE/DIAMETER,
    "uncertainty_rule": "旧图读范围与原文站位公差分列；并列界是保守可行性检查，不是独立随机误差合成。原文Nu±15%、Re±1%另列，不能用于吞掉坐标差。",
    "heat_sign_rule": "原文电热向流体放热为正；CFD入壁为正时需取负；段内先积分有符号热流，再除以名义面积。",
    "stop_rule": "不改旧CSV，不用理论曲线填点；第15至20行须保留中心身份敏感性，第18行不得宣称在旧坐标范围内闭合。",
    "old_inputs_sha256_before_and_after_identical": input_hashes,
    "generator_sha256": digest(Path(__file__).resolve()),
    "generated_csv_sha256": {name: digest(OUT/name) for name in ["nominal_stations.csv", "nominal_segments.csv"]},
}
(OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
print(json.dumps({key: summary[key] for key in [
    "nominal_geometric_station_count", "nominal_voltage_segment_count",
    "strict_all20_local_abscissa_identity_closed", "maximum_independent_area_relative_error",
    "maximum_station_position_area_relative_sensitivity",
    "plot_plus_station_bound_mismatch_rows_at_temperature_station",
]}, ensure_ascii=False, indent=2))
