"""从已保存 MAT、冻结 STL 和 CAD 分组独立复核空气域；不调用生成程序。"""
import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import sys

for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASE = ROOT / "cases/rfz_rec"
sys.path.insert(0, str(HERE.parent / "v344"))
from verify_surface_input import wall_hash


def read(path):
    return json.loads(path.read_bytes())


def face_keys(triangles, base):
    a, b, c = np.sort(triangles.astype(np.int64), axis=1).T
    return (a * base + b) * base + c


def verify(directory, first):
    summary = read(directory / "air.json")
    asset = directory / "air.mat"
    sha = sha256(asset.read_bytes()).hexdigest()
    assert sha == summary["fileSHA256"]
    data = loadmat(asset)
    points = data["vertices"]
    cells, wall, outer = [data[name] - 1 for name in ("tetrahedra", "wall", "outer")]
    owner = data["wallOwner"].ravel() - 1
    assert all(a.min() >= 0 and a.max() < len(points) for a in (cells, wall, outer))
    assert np.array_equal(np.unique(cells), np.arange(len(points)))
    assert len(np.unique(points, axis=0)) == len(points) and np.isfinite(points).all()
    assert len(np.unique(np.sort(cells, axis=1), axis=0)) == len(cells)

    # 从每一个四面体的四种三节点组合反查边界；内部同面两侧的对顶点须异侧。
    count = len(cells)
    ids = np.arange(4)
    faces = np.concatenate([cells[:, ids != i] for i in ids])
    opposite = cells.T.ravel()
    base = len(points)
    keys = face_keys(faces, base)
    order = np.argsort(keys)
    unique, starts, counts = np.unique(keys[order], return_index=True, return_counts=True)
    assert np.all((counts == 1) | (counts == 2))
    boundary = counts == 1
    wanted = face_keys(np.vstack((wall, outer)), base)
    assert np.array_equal(np.sort(wanted), unique[boundary])
    source = order[starts[np.searchsorted(unique, wanted)]]
    assert np.array_equal(owner, source[:len(wall)] % count)
    assert len(np.unique(wanted)) == len(wanted)
    pairs = starts[~boundary]
    left, right = order[pairs], order[pairs + 1]
    tri = points[faces[left]]
    normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    h1 = np.einsum("ij,ij->i", normal, points[opposite[left]] - tri[:, 0])
    h2 = np.einsum("ij,ij->i", normal, points[opposite[right]] - tri[:, 0])
    assert np.all(h1 * h2 < 0)
    graph = coo_matrix((np.ones(len(pairs), dtype=bool), (left % count, right % count)), shape=(count, count))
    assert connected_components(graph, directed=False, return_labels=False) == 1
    del tri, faces, graph

    # 原壁按 STL 顺序比较每个三角面的三个坐标；分别验算各 CAD 区和给定壁温指纹。
    case = read(CASE / "full_cad_input.json")
    stl = ROOT / "outputs/rfz_rec_geometry/RFZ-REC_m.stl"
    assert sha256(stl.read_bytes()).hexdigest() == case["stlSHA256"] == summary["stlSHA256"]
    stltype = [("normal", "<f4", 3), ("tri", "<f4", (3, 3)), ("attr", "<u2")]
    tri = np.fromfile(stl, offset=84, dtype=stltype)["tri"].astype(float)
    assert np.array_equal(points[wall], tri)
    normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = np.linalg.norm(normal, axis=1) / 2
    heights = np.einsum("ij,ij->i", normal, points[opposite[source[:len(wall)]]] - tri[:, 0]) / (2 * area)
    assert np.all(heights > 0)
    with (CASE / "cad_faces.csv").open(encoding="utf8") as stream:
        regions = list(csv.DictReader(stream))
    groups = np.zeros(len(wall), dtype=int)
    start = 0
    for index, row in enumerate(regions, 1):
        end = int(row["lastFace"])
        assert int(row["cadFace"]) == index and int(row["firstFace"]) == start + 1 and end > start
        groups[start:end] = index
        assert abs(area[start:end].sum() - float(row["stlAreaM2"])) < 1e-10
        start = end
    assert start == len(wall) and len(regions) == 204
    layers = read(CASE / "deg190_boundary.json")
    cavity = np.isin(groups, [v for cap in layers["caps"] for v in cap["innerCADFaces"]])
    temp = np.where(cavity, case["cavityWallK"], case["exteriorWallK"])
    assert wall_hash(case["cavityMode"], groups, temp, cavity) == case["wallSHA256"]
    assert int(cavity.sum()) == layers["groups"]["cavity"]["faces"]
    assert abs(area[cavity].sum() - layers["groups"]["cavity"]["areaM2"]) < 1e-9

    # 六个外盒面各自具有正确坐标、绕序和矩形总面积。
    extent = np.ptp(tri.reshape(-1, 3), axis=0).max()
    lo, hi = tri.min(axis=(0, 1)) - extent / 2, tri.max(axis=(0, 1)) + extent / 2
    assert np.array_equal(data["box"], np.array([lo, hi]))
    assert np.all((points >= lo - 1e-10) & (points <= hi + 1e-10))
    far = points[outer]
    vector = np.cross(far[:, 1] - far[:, 0], far[:, 2] - far[:, 0]) / 2
    box_faces = np.zeros(len(far), dtype=int)
    for axis in range(3):
        for value, sign in [(lo[axis], -1), (hi[axis], 1)]:
            select = np.all(abs(far[:, :, axis] - value) < 1e-10, axis=1)
            box_faces += select
            projected = sign * vector[select, axis]
            assert np.all(projected > 0)
            assert abs(projected.sum() / np.prod(np.delete(hi - lo, axis)) - 1) < 1e-12
    assert np.all(box_faces == 1)

    # 逐体行列式求体积。以等边四面体为参考，按 Frobenius 条件数独立求 SICN。
    corners = points[cells]
    edges = corners[:, 1:] - corners[:, :1]
    volumes = np.linalg.det(edges) / 6
    assert np.all(volumes > 0)
    solid = np.einsum("ij,ij->", tri.mean(axis=1), normal) / 6
    error = abs(volumes.sum() / (np.prod(hi - lo) - solid) - 1)
    assert error <= 2e-10 and abs(volumes.sum() - summary["airVolumeM3"]) < 1e-8
    reference = np.array([[1, 0, 0], [.5, np.sqrt(3)/2, 0], [.5, 1/np.sqrt(12), np.sqrt(2/3)]])
    jac = np.linalg.inv(reference) @ edges
    quality = 3 / (np.linalg.norm(jac, axis=(1, 2)) * np.linalg.norm(np.linalg.inv(jac), axis=(1, 2)))
    max_error = np.max(abs(quality - data["cellSICN"].ravel()))
    assert max_error < 1e-11
    assert abs(quality.min() - summary["quality"]["minSICN"]) < 1e-11
    assert (quality < .1).sum() == summary["quality"]["belowPointOne"]
    assert len(unique) - len(wanted) == summary["internalFaces"] and summary["fluidComponents"] == 1
    assert len(points) == summary["nodes"] and count == summary["tetrahedra"]
    assert (len(wall), len(outer)) == (summary["wallFaces"], summary["outerFaces"])
    assert abs(heights.min() - summary["minimumWallOppositeHeightM"]) < 1e-12
    result = dict(matSHA256=sha, vertices=len(points), cells=count, wall=len(wall), outer=len(outer),
                  internalFaces=len(pairs), components=1, regions=len(regions), cavityFaces=int(cavity.sum()),
                  cavityAreaM2=float(area[cavity].sum()), exteriorAreaM2=float(area[~cavity].sum()),
                  minVolumeM3=float(volumes.min()), airVolumeM3=float(volumes.sum()), relativeVolumeError=float(error),
                  minSICN=float(quality.min()), belowPointOne=int(np.sum(quality < .1)), maxSICNDifference=float(max_error),
                  minWallHeightM=float(heights.min()), computedFaces=0, computedAreaM2=0)
    low_owner = quality[owner] < .1
    result["wallWithOwnerBelowPointOne"] = {
        name: dict(faces=int(np.sum(mask & low_owner)), areaM2=float(area[mask & low_owner].sum()))
        for name, mask in (("exterior", ~cavity), ("cavity", cavity))
    }
    if first:
        original = loadmat(first)
        result["sameAsFirst"] = {key: np.array_equal(original[key], data[key])
                                 for key in ("wall", "vertices", "tetrahedra", "wallOwner", "outer", "box")}
        assert result["sameAsFirst"]["wall"]
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--first", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.directory, args.first), indent=2, allow_nan=False))
