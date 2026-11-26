"""从原 STL、发布的 CAD 区间及导出的 CSV 独立核对逐面壁温与面积。"""
import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CASE = ROOT / "cases/rfz_rec"


def read(path):
    return json.loads(path.read_bytes())


def wall_hash(mode, groups, temp, cavity):
    h = sha256(("fullCAD/" + mode).encode("utf8"))
    for data in (groups, temp, cavity):
        h.update(np.array([len(data), 1], dtype="<f8").tobytes())
        h.update(np.asarray(data, dtype="<f8").tobytes())
    return h.hexdigest()


def verify(directory, stl):
    model = read(CASE / "geometry.json")
    layers = read(CASE / "deg190_boundary.json")
    receipt = read(directory / "workbench.json")
    canonical = read(CASE / "full_cad_input.json")
    changed = read(directory / "RFZ-REC-input.json")
    assert sha256(stl.read_bytes()).hexdigest() == model["stlSHA256"]
    assert read(directory / "full_cad_input.json") == canonical
    with (CASE / "cad_faces.csv").open(encoding="utf8") as f:
        regions = list(csv.DictReader(f))
    groups = np.repeat([int(v["cadFace"]) for v in regions],
                       [int(v["lastFace"])-int(v["firstFace"])+1 for v in regions])
    cavity = np.isin(groups, [s for cap in layers["caps"] for s in cap["innerCADFaces"]])
    dtype = [("normal", "<f4", 3), ("vertices", "<f4", (3,3)), ("attr", "<u2")]
    tri = np.fromfile(stl, offset=84, dtype=dtype)["vertices"].astype(float)
    area = np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]), axis=1)/2
    csv_file = directory / "RFZ-REC-wall-input.csv"
    max_error = 0.
    temp = np.where(cavity, 315., 300.)
    with csv_file.open(encoding="utf-8-sig",newline="") as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames == ["faceID","cadFace","areaM2","givenWallK","temperatureSource"]
        count = 0
        for i, row in enumerate(reader):
            assert int(row["faceID"]) == i+1 and int(row["cadFace"]) == groups[i]
            assert float(row["givenWallK"]) == temp[i]
            assert row["temperatureSource"] == ("用户单独给定内壁温度" if cavity[i] else "论文给定外壁温度")
            max_error = max(max_error,abs(float(row["areaM2"])-area[i]))
            count += 1
    assert count == model["faces"] and max_error < 1e-12
    assert int(cavity.sum()) == receipt["innerFaces"] and abs(area[cavity].sum()-receipt["innerAreaM2"])<1e-10
    assert receipt["csvSHA256"] == sha256(csv_file.read_bytes()).hexdigest()
    assert canonical["wallSHA256"] == wall_hash("matchExterior",groups,np.full(count,300.),cavity)
    assert changed["wallSHA256"] == wall_hash("specified",groups,temp,cavity)
    assert canonical["meshSHA256"] != changed["meshSHA256"]
    assert receipt["computedFaces"] == 0
    return dict(faceRows=count,areaM2=float(area.sum()),cavityAreaM2=float(area[cavity].sum()),
                maximumAreaDifferenceM2=max_error,csvSHA256=receipt["csvSHA256"],
                canonicalWallSHA256=canonical["wallSHA256"],changedWallSHA256=changed["wallSHA256"])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    p.add_argument("--stl",type=Path,default=ROOT/"outputs/rfz_rec_geometry/RFZ-REC_m.stl")
    args=p.parse_args()
    print(json.dumps(verify(args.directory,args.stl),indent=2))
