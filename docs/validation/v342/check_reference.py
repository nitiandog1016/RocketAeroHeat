"""独立量取作者点云到冻结 STL 的距离。只读取几何，保留作者坐标与半域。"""
import argparse
import hashlib
import json
import os
from pathlib import Path

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[3]
COORDINATES = ("e511731392daafa3e97379a6a40c3dc307ce26705519fc1a843d5327eca4cc8d",
               "62f415a45076b7604d9feac76886ede6a3073678c606a3266139e4b219b59a02",
               "4e11f63efb2a4495a263a71ab049bce3a4f2415ef12966ca6df25d079cd01c38")


def dot(a, b):
    return np.einsum("ij,ij->i", a, b)


def squared_distance(points, triangles):
    # 投影位于三角形内时取到平面的距离，否则取到三条闭线段的最短距离。
    a, b, c = triangles.transpose(1, 0, 2)
    normal = np.cross(b-a, c-a)
    plane = dot(points-a, normal)**2 / dot(normal, normal)
    inside = np.ones(len(points), dtype=bool)
    segments = np.full(len(points), np.inf)
    for start, end in ((a,b), (b,c), (c,a)):
        edge = end-start
        delta = points-start
        inside &= dot(np.cross(edge, delta), normal) >= 0
        t = np.clip(dot(delta, edge)/dot(edge, edge), 0, 1)
        residual = delta - t[:, None]*edge
        segments = np.minimum(segments, dot(residual, residual))
    return np.where(inside, plane, segments)


def distances(points, tri):
    # 最近顶点给出上界；重心与三角形任一点相距不超过最长顶点半径。
    # 按半径分组，搜全可能改善当前上界的三角形，逐三角形计算真实距离。
    center = tri.mean(axis=1)
    radius = np.linalg.norm(tri-center[:, None, :], axis=2).max(axis=1)
    bins = np.floor(np.log2(radius)).astype(int)
    best = cKDTree(tri.reshape(-1, 3)).query(points)[0]**2
    for group in np.unique(bins):
        ids = np.flatnonzero(bins == group)
        tree = cKDTree(center[ids])
        for offset in range(0, len(points), 256):
            chunk = slice(offset, offset+256)
            bounds = np.sqrt(best[chunk]) + radius[ids].max() + 1e-12
            hits = tree.query_ball_point(points[chunk], bounds)
            counts = [len(h) for h in hits]
            if not any(counts):
                continue
            p = np.repeat(np.arange(len(hits)), counts)
            cells = ids[np.concatenate([h for h in hits if h])]
            dist2 = squared_distance(points[chunk][p], tri[cells])
            np.minimum.at(best[chunk], p, dist2)
    return np.sqrt(best)


def check(stl, source):
    model = json.loads((ROOT/"cases/rfz_rec/geometry.json").read_bytes())
    assert hashlib.sha256(stl.read_bytes()).hexdigest() == model["stlSHA256"]
    dtype = [("normal", "<f4", 3), ("vertex", "<f4", (3,3)), ("attr", "<u2")]
    tri = np.fromfile(stl, offset=84, dtype=dtype)["vertex"].astype(float)
    axes = []
    for name, expected in zip("XYZ", COORDINATES):
        file = source/(name+"_data.npy")
        assert hashlib.sha256(file.read_bytes()).hexdigest() == expected
        axes.append(np.load(file, allow_pickle=False).reshape(-1))
    points = np.column_stack(axes)
    assert points.shape == (72706, 3) and np.isfinite(points).all()
    d = distances(points, tri)
    assert np.isfinite(d).all()
    # 给小样本独立遍历全部三角面，检查空间搜索包含真正最近的单元。
    for i in np.linspace(0,len(points)-1,9,dtype=int):
        direct = np.sqrt(squared_distance(np.broadcast_to(points[i], (len(tri),3)), tri).min())
        assert abs(direct-d[i]) < 1e-11
    return dict(dataset="EUCASS_2023/deg190", coordinateSHA256=dict(zip("XYZ",COORDINATES)),
                referencePoints=len(d), stlSHA256=model["stlSHA256"],
                referenceBoundsM=[points.min(axis=0).tolist(),points.max(axis=0).tolist()],
                distanceQuantilesM={str(n):float(np.percentile(d,n)) for n in (0,50,90,95,99,100)},
                pointsWithinMm={str(mm):int(sum(d<=mm/1000)) for mm in (1,10,50)},
                comparedPhysicalFields=[], coordinateTransform="identity (m)")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("stl",type=Path);p.add_argument("source",type=Path);p.add_argument("result",type=Path)
    args=p.parse_args()
    result=check(args.stl,args.source)
    with args.result.open("x",encoding="utf8",newline="\n") as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write("\n")
    print(json.dumps(result))
