"""Independent byte/transcription/geometry checks; does not call the builder or CFD.

Usage: python tests/data/tmr_mach7/verify_contract.py [output_directory]
The verifier reads both acquired official distributions, their ZIP membership,
and all five original meridional grids. Passing is a data-contract result only.
"""
import gzip
import hashlib
import json
import math
import pathlib
import re
import sys
import zipfile

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent


def sha(data):
    return hashlib.sha256(data).hexdigest()


def numeric_rows(text, columns):
    rows = []
    for line in text.splitlines():
        if re.match(r"^\s*[-+\d.]", line):
            row = [float(value) for value in line.split()]
            assert len(row) == columns and all(math.isfinite(v) for v in row)
            rows.append(row)
    return np.asarray(rows)


def zone(text, name, columns):
    selected = []
    active = False
    for line in text.splitlines():
        match = re.search(r'zone\s+t="([^"]+)"', line, re.I)
        if match:
            if active:
                break
            active = match.group(1) == name
        elif active:
            selected.append(line)
    return numeric_rows("\n".join(selected), columns)


def close(actual, expected, tolerance=2e-12):
    actual, expected = np.asarray(actual), np.asarray(expected)
    assert actual.shape == expected.shape
    assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected))
    delta = np.abs(actual-expected)
    assert np.all(delta <= tolerance*np.maximum(1, np.abs(expected))), float(np.max(delta))


def read_grid(payload):
    values = np.fromstring(gzip.decompress(payload).decode("ascii"), sep=" ")
    assert values[0] == 2
    dims = values[1:5].astype(int).reshape(2, 2)
    offset, result = 5, []
    for ni, nj in dims:
        count = ni*nj
        x = values[offset:offset+count].reshape((ni, nj), order="F")*.01
        r = values[offset+count:offset+2*count].reshape((ni, nj), order="F")*.01
        assert np.all(np.isfinite(x)) and np.all(r > 0)
        result.append((x, r))
        offset += 2*count
    assert offset == len(values)
    return result


def main():
    contract = json.loads((ROOT/"contract.json").read_text(encoding="utf-8"))
    manifest = json.loads((ROOT/"source_manifest.json").read_text(encoding="utf-8"))
    close([contract["scope"]["flare_angle_deg"], contract["scope"]["angle_of_attack_deg"]], [20, 0])
    assert contract["scope"]["axisymmetric"] and not contract["scope"]["planar_2d_equivalent"]
    assert not contract["scope"]["solver_validation_complete"]
    assert contract["scope"]["new_physical_area_m2"] == 0
    assert contract["source_manifest_sha256"] == sha((ROOT/"source_manifest.json").read_bytes())
    assert contract["sources"]["original_report_local_sha256"] is None
    source_hashes, extracted_count = {}, 0
    for record in manifest["files"]:
        path = (ROOT/record["file"]).resolve()
        assert ROOT.resolve() in path.parents
        data = path.read_bytes()
        assert len(data) == record["bytes"] and sha(data) == record["sha256"]
        source_hashes[record["file"]] = sha(data)
        if record["archive_file"]:
            with zipfile.ZipFile(ROOT/record["archive_file"]) as archive:
                assert archive.read(record["archive_member"]) == data
            extracted_count += 1

    raw = ROOT/"raw"
    p = numeric_rows((raw/"pw_exp.dat").read_text(), 2)
    q = numeric_rows((raw/"qw_exp.dat").read_text(), 2)
    u = numeric_rows((raw/"upstream_exp.dat").read_text(), 3)
    archive_pq = (raw/"KussoyCompCorner.dat").read_text()
    author_p = zone(archive_pq, "KussoyP20", 4)
    author_q = zone(archive_pq, "KussoyQ20", 4)
    author_u = zone((raw/"KussoyBLProfile.dat").read_text(), "KussoyBL20_IncomingBL", 5)
    assert p.shape == (27, 2) and q.shape == (27, 2) and u.shape == (21, 3)
    unique = []
    seen = set()
    for row in p:
        key = tuple(row)
        unique.append(key not in seen)
        seen.add(key)
    assert sum(unique) == 24 and np.array_equal(np.flatnonzero(np.logical_not(unique))+1, [25, 26, 27])
    close(p[unique], author_p[:, :2], 1e-6)
    close(q, author_q[:, [0, 2]], 1e-6)
    close(u, author_u[:, [0, 3, 4]], 1e-6)
    assert np.all(author_p[:, 2:] == 0) and np.all(author_q[:, [1, 3]] == 0)

    conditions = (raw/"KussoyConditions.txt").read_text()
    def condition(name):
        return float(re.search(r"(?im)^\s*"+re.escape(name)+r"\s*=\s*([-+\d.Ee]+)", conditions).group(1))
    ref = contract["experimental_local_reference"]
    pairs = [("Mach", "Minf", 1), ("temperature_K", "Tinf", 1),
             ("pressure_Pa", "Pinf", 1), ("density_kg_m3", "RhoInf", 1),
             ("velocity_m_s", "Uinf", 1), ("wall_temperature_K", "TWall", 1),
             ("reported_wall_heating_magnitude_W_m2", "QWallInf", 1),
             ("wall_shear_Pa", "TauWallInf", 1), ("unit_Re_m_inv", "Re/m", 1),
             ("boundary_layer_thickness_m", "Delta0", .01),
             ("displacement_thickness_m", "Delta*0", .01),
             ("momentum_thickness_m", "Theta0", .01), ("Cf", "CfInf", 1)]
    for key, name, factor in pairs:
        close(ref[key], condition(name)*factor)
    assert ref["total_temperature_K"] is None
    inlet = contract["computational_truncated_inflow"]
    assert inlet["incoming_turbulence_profiles"] is None and inlet["not_a_measured_complete_profile"]
    wind = (raw/"wind_sa.inp").read_text()
    values = re.search(r"(?im)^freestream\s+total\s+([^\r\n]+)", wind).group(1).split()
    close([inlet[k] for k in ["Mach", "total_pressure_psi_WIND_input", "total_temperature_rankine_WIND_input", "angle_of_attack_deg", "sideslip_deg"]], [float(v) for v in values])
    close(inlet["total_pressure_Pa_WIND_input"], float(values[1])*6894.757293168)
    close(inlet["total_temperature_K_WIND_input"], float(values[2])*5/9)
    assert ref["Mach"] == 7.05 and inlet["Mach"] == 7.11
    close([ref["location_source_x_m"], inlet["source_x_m"]], [-.06, -.7839])
    for name, source, value_key, factor in [
        ("pressure_rows", p, "pressure_Pa", 576),
        ("heating_rows", q, "reported_heating_magnitude_W_m2", 9300)]:
        rows = contract["measurements"][name]
        close([r["source_row"] for r in rows], np.arange(1, 28))
        s = source[:, 0]*.01
        x = np.where(s <= 0, s, s*math.cos(math.radians(20)))
        radius = .1015+np.maximum(s, 0)*math.sin(math.radians(20))
        for key, actual in [("S_m", s), ("source_x_m", x), ("project_X_m", -x),
                            ("radius_m", radius), ("ratio", source[:, 1]),
                            (value_key, source[:, 1]*factor)]:
            close([r[key] for r in rows], actual)
    assert [r["independent_comparison_point"] for r in contract["measurements"]["pressure_rows"]] == unique
    profile = contract["measurements"]["upstream_profile"]
    for key, source in [("y_from_wall_m", u[:, 0]*.01), ("velocity_ratio", u[:, 1]),
                        ("total_temperature_ratio", u[:, 2]), ("measured_M", author_u[:, 1]),
                        ("pressure_ratio", author_u[:, 2])]:
        close([r[key] for r in profile], source)
    assert q[0, 0] < 0 and "12.06   0.98" in (raw/"KussoyAxiCompCorner.txt").read_text()

    grid_summary, previous = [], None
    with zipfile.ZipFile(raw/"Axiswblim7_grids.zip") as archive:
        for level in range(5):
            suffix = "kussoy_2zn_{}levdn.p2dfmt.gz".format(level)
            member = [n for n in archive.namelist() if n.endswith(suffix)]
            assert len(member) == 1
            grids = read_grid(archive.read(member[0]))
            close(grids[0][0][-1, :], grids[1][0][0, :])
            close(grids[0][1][-1, :], grids[1][1][0, :])
            blocks = []
            for block_index, (x, r) in enumerate(grids):
                ni, nj = x.shape
                assert (ni, nj) == ((321-1)//(2**level)+1, (401-1)//(2**level)+1)
                nominal = .10150094+np.maximum(x[:, 0], 0)*math.tan(math.radians(20))
                # Original grid is nominally 20 degrees: actual slope is about
                # 19.999867 degrees, differing by <=0.618 micrometre in radius.
                # Audit that source discrepancy; do not rescale source vertices.
                close(r[:, 0], nominal, 1e-6)
                if previous is not None:
                    close(x, previous[block_index][0][::2, ::2])
                    close(r, previous[block_index][1][::2, ::2])
                gap = np.hypot(x[:, 1]-x[:, 0], r[:, 1]-r[:, 0])
                flare = x[:, 0] > 0
                angle = None
                if np.count_nonzero(flare) >= 2:
                    slope = np.polyfit(x[flare, 0], r[flare, 0], 1)[0]
                    angle = float(np.rad2deg(np.arctan(slope)))
                blocks.append({"ni": ni, "nj": nj, "first_normal_node_distance_min_m": float(gap.min()),
                               "max_nominal_20deg_radius_difference_m": float(np.max(np.abs(r[:, 0]-nominal))),
                               "fitted_flare_angle_deg": angle})
            grid_summary.append({"level_down": level, "blocks": blocks})
            previous = grids

    result = {"status": "passed_independent_data_contract_only", "solver_run": False,
              "new_physical_area_m2": 0, "source_files_verified": len(manifest["files"]),
              "archive_members_byte_verified": extracted_count,
              "independent_pressure_points": 24, "raw_pressure_rows": 27,
              "heat_points": 27, "inlet_profile_points": 21,
              "grids_checked": grid_summary,
              "source_sha256": source_hashes,
              "contract_sha256": sha((ROOT/"contract.json").read_bytes()),
              "builder_sha256": sha((ROOT/"build_contract.py").read_bytes()),
              "verifier_sha256": sha(pathlib.Path(__file__).read_bytes()),
              "missing_physical_inputs": ["measured turbulence profiles", "complete transverse velocity profile", "experimental numerical total-temperature reference", "frozen transport and turbulence closure"],
              "limits": ["data/geometry checks are not NS/RANS or experimental validation", "original report local PDF byte hash is unavailable", "grid node spacing is not cell-center spacing or y-plus"]}
    if len(sys.argv) > 1:
        output = pathlib.Path(sys.argv[1])
        output.mkdir(parents=True, exist_ok=True)
        (output/"independent_contract_check.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print("Independent TMR contract verified: {} raw files, {} archive members, five grid levels; no solver run.".format(len(manifest["files"]), extracted_count))


if __name__ == "__main__":
    main()
