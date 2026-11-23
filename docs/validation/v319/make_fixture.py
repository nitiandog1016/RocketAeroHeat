"""按独立的面号真值生成可从工作台普通文件页打开的毫米制 STL。"""
import argparse
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
FACET = np.dtype([("normal", "<f4", 3), ("xyz", "<f4", (3, 3)), ("attribute", "<u2")])


def make(out):
    out.mkdir(parents=True, exist_ok=True)
    with ZipFile(ROOT / "cases/thermistor_b/leads3.zip") as z:
        body = json.loads(z.read("body.json"))
    triangles = np.asarray(body["V"])[np.asarray(body["F"]) - 1]
    rng = np.random.RandomState(319)
    order = rng.permutation(len(triangles))
    records = np.zeros(len(order), dtype=FACET)
    for i, j in enumerate(order):
        records["xyz"][i] = np.roll(triangles[j], rng.randint(3), axis=0) * 1000
    normals = np.cross(records["xyz"][:, 1] - records["xyz"][:, 0],
                       records["xyz"][:, 2] - records["xyz"][:, 0])
    records["normal"] = normals / np.linalg.norm(normals, axis=1)[:, None]
    with (out / "probe-mm.stl").open("xb") as f:
        f.write(b"thermistor B  |  mm".ljust(80))
        f.write(len(order).to_bytes(4, "little"))
        records.tofile(f)
    np.savetxt(out / "reference_ids.txt", order + 1, fmt="%d")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    make(parser.parse_args().out)
