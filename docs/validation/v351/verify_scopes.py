"""回读网页实际导出的三个全表面 CSV，核对逐面取样范围和数值。"""
import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
from zipfile import ZipFile

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
CASE = ROOT / "cases/simple_launcher"
EXPORT = """
const fs = require('node:fs');
const path = require('node:path');
const [folder, lib, key] = process.argv.slice(1);
const text = fs.readFileSync(path.join(folder, 'case-data.js'), 'utf8').trim();
if (!text.startsWith('window.ROCKET_CASE=') || !text.endsWith(';')) throw Error('结果格式不符');
const model = JSON.parse(text.slice('window.ROCKET_CASE='.length,-1));
const field = model.fields.find(item => item.key === key);
process.stdout.write(require(path.join(lib, 'export.js')).csv(model, field));
"""


def check(folder, web):
    study = json.loads((CASE / "wind_tunnel_stations.json").read_bytes())
    with ZipFile(CASE / "wind_tunnel.zip") as bundle:
        with bundle.open("wall.csv") as f:
            wall = np.loadtxt(f, delimiter=",", skiprows=1)
    geometry = json.loads((CASE / "geometry.json").read_bytes())
    groups = {r["name"]: r for r in geometry["regions"]}
    units = {"pressure": "kPa", "pressureGauge": "Pa", "heatFlux": "W/m²"}
    ref_pressure = json.loads((CASE / "tunnel.json").read_bytes())["ambientPressurePa"]
    arrays = {"pressure": wall[:, 2]/1000, "pressureGauge": wall[:, 2]-ref_pressure, "heatFlux": wall[:, 3]}
    station_labels = [f"{r:.3f} R 处 {p:+.5f}" for r,p in
                      zip(study["tapRadiusOverBaseRadius"],study["grids"]["outer"]["tapGaugeMilliPa"])]
    counts = {}
    for key in units:
        run = subprocess.run(["node", "-e", EXPORT, str(folder), str(web), key],
                             check=True, capture_output=True, timeout=45)
        rows = csv.DictReader(io.StringIO(run.stdout.decode("utf-8-sig")))
        values = []
        for face, row in enumerate(rows, 1):
            assert int(row[rows.fieldnames[0]]) == face and row["单位"] == units[key]
            scope = row["逐面元适用范围"]
            on_base = groups["base"]["firstFace"] <= face <= groups["base"]["lastFace"]
            on_nose = groups["sphericalNose"]["firstFace"] <= face <= groups["sphericalNose"]["lastFace"]
            if on_base:
                assert "0 至 80.0 mm" in scope
            elif on_nose:
                assert "鼻端面平均" in scope
                if study["noseTipFaceIDs"][0] <= face <= study["noseTipFaceIDs"][1]:
                    assert "0.000 至 1.224 mm" in scope
            else:
                assert "原三角面面积平均" in scope
            assert ("径向点值" in scope) == (on_base and key != "heatFlux")
            if on_base and key != "heatFlux":
                assert all(label in scope for label in station_labels)
            values.append(float(row["数值"]))
        np.testing.assert_allclose(values, arrays[key], rtol=0, atol=1e-12)
        counts[key] = len(values)
    return counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbench", type=Path)
    parser.add_argument("--web", type=Path, help="独立的发布版网页脚本目录")
    args = parser.parse_args()
    print(json.dumps(check(args.workbench.resolve(), (args.web or args.workbench).resolve())))
