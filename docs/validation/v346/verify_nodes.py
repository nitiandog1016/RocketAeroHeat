"""用原始 STL、边界分区和四面体逐行核对节点温界及工作台导出。"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[3]
CASE = ROOT / "cases/rfz_rec"
STATUS = ["内部对点待求温", "外盒对点", "四壁点同温", "四壁点异温", "节点温界冲突"]
LOCATION = ["空气内部", "外盒", "实体壁"]


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(case_file, asset, face_csv=None, node_csv=None):
    source = read(ROOT / "docs/validation/v345/evidence/runs.json")
    accepted = read(ROOT / "docs/validation/v345/evidence/air.json")
    mat = asset or (ROOT / source["assetDirectory"] / "air.mat")
    file_hash = digest(mat)
    if asset is None:
        assert file_hash == accepted["fileSHA256"]
    case = read(case_file)
    assert case["stlSHA256"] == accepted["stlSHA256"]
    stl = ROOT / "outputs/rfz_rec_geometry/RFZ-REC_m.stl"
    assert digest(stl) == case["stlSHA256"] and digest(CASE / "cad_faces.csv") == case["cadTableSHA256"]
    stl_type = [("normal", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2")]
    triangles = np.fromfile(stl, dtype=stl_type, offset=84)["v"].astype(float)
    area = np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]), axis=1)/2
    data = loadmat(mat)
    # MAT 文件头含写入时间。重建件按有序坐标、连接及相邻单元号核验。
    geometry = hashlib.sha256()
    for name in ("vertices", "tetrahedra", "wall", "outer", "wallOwner"):
        array = data[name]
        assert array.dtype == np.dtype("<f8" if name == "vertices" else "<i4")
        geometry.update(name.encode("ascii") + b"\x00")
        geometry.update(np.asarray(array.shape, dtype="<i8").tobytes())
        geometry.update(array.tobytes(order="C"))
    assert geometry.hexdigest() == "c55b1ae554fb2867579ad3a1b1779e2de24936ca65520dc2e888175c67af9c4f"
    faces = data["wall"]
    assert np.array_equal(data["vertices"][faces-1], triangles)
    owner = data["wallOwner"].ravel()
    adjacent = data["tetrahedra"][owner-1]
    extra = ~np.any(adjacent[:, :, None] == faces[:, None, :], axis=2)
    assert np.all(extra.sum(axis=1) == 1)
    opposite = adjacent[extra]
    with (CASE / "cad_faces.csv").open(encoding="utf8") as stream:
        regions = list(csv.DictReader(stream))
    groups = np.repeat([int(r["cadFace"]) for r in regions], [int(r["lastFace"])-int(r["firstFace"])+1 for r in regions])
    inside_ids = {i for cap in read(CASE / "deg190_boundary.json")["caps"] for i in cap["innerCADFaces"]}
    inside = np.isin(groups, list(inside_ids))
    given = np.where(inside, case["cavityWallK"], case["exteriorWallK"])
    outer_source = "论文给定外壁温度"
    inner_source = "项目假设：内壁随外壁同温" if case["cavityMode"] == "matchExterior" else "用户单独给定内壁温度"

    # 独立按温度分组建立节点集合，共享节点同时落入两个集合时保留多值。
    n = len(data["vertices"])
    values = np.unique(given)
    membership = np.zeros((len(values), n+1), dtype=bool)
    for i, temperature in enumerate(values):
        membership[i, faces[given == temperature]] = True
    wall = membership.any(axis=0)
    bounds = np.where(membership, values[:, None], np.nan)
    lo = np.full(n+1, np.nan)
    hi = lo.copy()
    lo[wall], hi[wall] = np.nanmin(bounds[:, wall], axis=0), np.nanmax(bounds[:, wall], axis=0)
    conflict = membership.sum(axis=0) > 1
    mixed = conflict[adjacent].sum(axis=1)
    far = np.isin(opposite, data["outer"])
    loc = np.where(wall[opposite], 2, np.where(far, 1, 0))
    status = np.where(mixed, 4, loc + ((loc == 2) & (lo[opposite] != given)))
    is13 = groups == 13
    stats = {}
    for label, mask in (("cad13", is13), ("all", np.ones(len(faces), dtype=bool))):
        assert abs(area[mask].sum() - (477.36662700004854 if label == "cad13" else 740.7290034031088)) < 1e-8
        stats[label] = [dict(status=name, faces=int(np.sum(mask & (status == i))),
                             areaM2=float(area[mask & (status == i)].sum())) for i, name in enumerate(STATUS)]
    vertex_ids = np.flatnonzero(wall)
    if face_csv:
        with face_csv.open(encoding="utf-8-sig", newline="") as stream:
            rows = csv.DictReader(stream)
            size = 0
            for i, row in enumerate(rows):
                assert [int(row[k]) for k in ("faceID", "cadFace", "ownerCellID", "oppositeNodeID", "mixedWallNodesInCell")] == [
                    i+1, groups[i], owner[i], opposite[i], mixed[i]]
                assert row["nodalConstraint"] == STATUS[status[i]] and row["oppositeLocation"] == LOCATION[loc[i]]
                assert abs(float(row["areaM2"])-area[i]) < 1e-12 and float(row["givenWallK"]) == given[i]
                assert row["temperatureSource"] == (inner_source if inside[i] else outer_source)
                for key, value in (("oppositeWallMinK", lo[opposite[i]]), ("oppositeWallMaxK", hi[opposite[i]])):
                    check_optional(row[key], value)
                size += 1
            assert size == len(faces)
    if node_csv:
        incident = np.bincount(faces.ravel(), minlength=n+1)
        in_counts = np.bincount(faces[inside].ravel(), minlength=n+1)
        with node_csv.open(encoding="utf-8-sig", newline="") as stream:
            size = 0
            for i, row in enumerate(csv.DictReader(stream)):
                vid = vertex_ids[i]
                assert int(row["vertexID"]) == vid and int(row["incidentFaces"]) == incident[vid]
                assert float(row["minimumGivenK"]) == lo[vid] and float(row["maximumGivenK"]) == hi[vid]
                check_optional(row["uniqueGivenK"], np.nan if conflict[vid] else lo[vid])
                expected_sources = set()
                if in_counts[vid]:
                    expected_sources.add(inner_source)
                if in_counts[vid] < incident[vid]:
                    expected_sources.add(outer_source)
                assert set(row["temperatureSource"].split(" / ")) == expected_sources
                size += 1
            assert size == len(vertex_ids)
    return dict(caseSHA256=digest(case_file), matSHA256=file_hash, orderedGeometrySHA256=geometry.hexdigest(), wallNodes=len(vertex_ids),
                ambiguousNodes=int(conflict.sum()), cases=stats, verifiedFaceCSV=bool(face_csv), verifiedNodeCSV=bool(node_csv))


def check_optional(token, value):
    actual = float(token or "nan")
    assert actual == value or (np.isnan(actual) and np.isnan(value))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", type=Path, help="工作台验收导出目录")
    parser.add_argument("--asset", type=Path, help="重新生成且有序几何一致的 air.mat")
    args = parser.parse_args()
    inputs = {"uniform": CASE / "full_cad_input.json",
              "split": ROOT / "docs/validation/v344/evidence/rotated_input.json"}
    result = {}
    for name, case_file in inputs.items():
        files = [args.export / (name+"-"+suffix+".csv") for suffix in ("faces", "nodes")] if args.export else []
        result[name] = check(case_file, args.asset, *files)
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
