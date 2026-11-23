"""在 Short 65 的八条真实边及三个网格尺度上检验逐边 0.5 幅值限制。"""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import time
from zipfile import ZipFile

for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[variable] = "1"
import numpy as np
import psutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ARCHIVE = ROOT/"docs/validation/v326/evidence/receipts.zip"
REFERENCE = ROOT/"cases/short_sphere/supportedFlux65.json"
FIELDS = dict(constant=[0, 0, 0], x=[1, 0, 0], y=[0, 1, 0], z=[0, 0, 1], oblique=[1, -2, 3])
TOL = 1e-12


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sources():
    assert digest(ARCHIVE.read_bytes()) == "37b679a7181e5d1a45855242f7caef298f8387a1a47e4ebd6fa882562ca9dca5"
    assert digest(REFERENCE.read_bytes()) == "1c4fb8cff2b83f23ac4343b602b9b28a1c9d0f3026a4443928f366cc0b3c3a05"
    ref = json.loads(REFERENCE.read_bytes())
    with ZipFile(ARCHIVE) as archive:
        hashes = json.loads(archive.read("manifest.json"))["included"]
        def raw(name):
            data = archive.read(name)
            assert digest(data) == hashes[name]
            return data
        nodes = np.loadtxt(io.BytesIO(raw("first_nodes.csv")), delimiter=",", skiprows=1)
        edges = np.loadtxt(io.BytesIO(raw("outputs/v326_validation/detail01/scalar_edges.csv")), delimiter=",")
        prep = json.loads(raw("outputs/v325_validation/prep02/prepared/prepared.json"))
    assert prep["files"]["mesh.su2"] == ref["meshSHA256"] and prep["caseSHA256"] == ref["caseSHA256"]
    return ref, 1/prep["condition"]["pe"], nodes, edges[edges[:, 0].argsort()]


def dot(a, b):
    return np.sum(a*b, axis=-1)


def edge_flux(left, right, normal, jump, gradient, diffusivity):
    """复算 SU2 沿边修正；分别返回两个节点共享的原通量、受限通量和隐式两点系数。"""
    dx = right - left
    length2 = dot(dx, dx)
    if np.any(length2 <= 0):
        raise ValueError("节点重合")
    projection = dot(dx, normal)/length2
    diagonal = diffusivity*projection
    two = -diagonal*jump
    corrected_grad = gradient + ((jump-dot(dx, gradient))/length2)[:, None]*dx
    original = -diffusivity*dot(corrected_grad, normal)
    correction = original-two
    limited = two + np.sign(correction)*np.minimum(abs(correction), abs(two))
    return original, limited, two, diagonal


def compare(left, right, normal, gradient, diffusivity):
    # 以实际坐标生成两个端值，解析法向通量直接由面面积向量积分。
    jump = (0.25 + dot(right, gradient)) - (0.25 + dot(left, gradient))
    original, limited, two, diag = edge_flux(left, right, normal, jump, gradient, diffusivity)
    _, reverse, _, reverse_diag = edge_flux(right, left, -normal, -jump, gradient, diffusivity)
    exact = -diffusivity*dot(normal, gradient)
    scale = diffusivity*np.linalg.norm(normal, axis=1)*(np.linalg.norm(gradient) or 1.)
    assert np.all(scale > 0) and np.all(diag > 0)
    np.testing.assert_allclose(diag, reverse_diag, rtol=TOL, atol=0)
    assert np.max(abs(original-exact)/scale) < TOL and np.max(abs(reverse+limited)/scale) < TOL
    # 单独扰动两点沿边项，核对延迟修正下对角的符号和大小。
    test_increment = np.ones_like(jump)
    _, _, base, _ = edge_flux(left, right, normal, 0*jump, gradient, diffusivity)
    _, _, left_perturbed, _ = edge_flux(left, right, normal, -test_increment, gradient, diffusivity)
    np.testing.assert_allclose(left_perturbed-base, diag, rtol=TOL, atol=0)
    return exact, original, limited, two, diag, scale, reverse


def actual_star(ref, diffusivity, nodes, edges):
    xyz = {int(row[0]): row[1:4] for row in nodes}
    point = ref["first"]["worstID"]
    assert len(edges) == len(np.unique(edges[:, 0])) == ref["row"]["edges"]
    assert np.all((edges[:, 1] == point) ^ (edges[:, 2] == point))
    left, right = [np.array([xyz[int(n)] for n in edges[:, k]]) for k in (1, 2)]
    sign = np.where(edges[:, 1] == point, 1., -1.)
    normals = edges[:, 3:6]
    rows = []
    for name, values in FIELDS.items():
        g = np.array(values)
        exact, original, limited, two, diag, scale, reverse = compare(left, right, normals, g, diffusivity)
        for k, edge in enumerate(edges):
            rows.append(dict(field=name, edge=int(edge[0]), i=int(edge[1]), j=int(edge[2]), exact=exact[k],
                             original=original[k], twoPoint=two[k], limited=limited[k], reverse=reverse[k],
                             diagonal=diag[k], scale=scale[k], targetSign=int(sign[k]),
                             normalizedError=abs(limited[k]-exact[k])/scale[k]))
    worst = max(rows, key=lambda r: r["normalizedError"])
    return rows, dict(pointID=point, xyzD=xyz[point].tolist(), worst=worst,
                      segmentEndsD=np.stack((left, right), axis=1).tolist(),
                      edgeIDs=edges[:, 0].astype(int).tolist())


def median_cut(n, axis, shear):
    """在同一物理截面拼接节点中心对偶面。截面边界半片、角点四分片。"""
    assert n % 2 == 1
    a, b = (axis+1) % 3, (axis+2) % 3
    u, v = np.meshgrid(np.arange(n+1)/n, np.arange(n+1)/n, indexing="ij")
    cut = np.zeros((u.size, 3))
    cut[:, axis], cut[:, a], cut[:, b] = .5, u.ravel(), v.ravel()
    basis = np.eye(3)
    basis[1, 0] = shear  # (ξ,η,ζ) 映射为 (ξ+sη,η,ζ)
    half = .5/n
    left = (cut-half*np.eye(3)[axis]) @ basis
    right = (cut+half*np.eye(3)[axis]) @ basis
    quad = []
    for da, db in ((-half,-half), (half,-half), (half,half), (-half,half)):
        c = cut.copy()
        c[:, a] = np.clip(c[:, a]+da, 0, 1)
        c[:, b] = np.clip(c[:, b]+db, 0, 1)
        quad.append(c @ basis)
    normal = .5*np.cross(quad[2]-quad[0], quad[3]-quad[1])
    analytic_normal = np.cross(basis[a], basis[b])
    np.testing.assert_allclose(normal.sum(axis=0), analytic_normal, rtol=TOL, atol=TOL)
    return left, right, normal, analytic_normal


def refinement(diffusivity):
    rows, failures = [], []
    for shear in (0., .5):
        for axis in range(3):
            for field, values in FIELDS.items():
                g = np.array(values)
                series = []
                for n in (3, 9, 27):
                    l, r, normals, whole = median_cut(n, axis, shear)
                    exact, original, limited, two, _, scale, reverse = compare(l, r, normals, g, diffusivity)
                    analytic = -diffusivity*whole.dot(g)
                    assert abs(exact.sum()-analytic) < TOL*scale.sum()
                    error = abs(limited.sum()-analytic)/scale.sum()
                    if shear == 0 or field == "constant":
                        assert error < TOL
                    series.append(error)
                    rows.append(dict(shear=shear, axis=axis, field=field, divisions=n, faces=len(normals),
                                     analytic=analytic, original=original.sum(), twoPoint=two.sum(),
                                     limited=limited.sum(), naturalScale=scale.sum(), normalizedError=error,
                                     maxReversalError=float(np.max(abs(reverse+limited)/scale))))
                order_ok = all(fine <= max(TOL, coarse/3+TOL) for coarse, fine in zip(series, series[1:]))
                if not order_ok or series[-1] >= .01:
                    failures.append(dict(shear=shear, axis=axis, field=field, errors=series,
                                         orderPass=bool(order_ok), finestPass=bool(series[-1]<.01)))
    return rows, failures


def save_csv(path, rows):
    with path.open("x", encoding="utf8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def dump(path, value):
    with path.open("x", encoding="utf8", newline="\n") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def run(out):
    started = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    ref, diffusivity, nodes, edges = sources()
    actual, star = actual_star(ref, diffusivity, nodes, edges)
    cuts, failures = refinement(diffusivity)
    save_csv(out/"star.csv", actual)
    save_csv(out/"cuts.csv", cuts)
    result = dict(caseSHA256=ref["caseSHA256"], meshSHA256=ref["meshSHA256"], referenceSHA256=digest(REFERENCE.read_bytes()),
                  originalFaces=ref["originalFaces"], originalAreaM2=ref["originalAreaM2"],
                  diffusivity=diffusivity, limitCoefficient=.5, star=star, failedCuts=failures,
                  passed=not failures, cutRows=len(cuts), starRows=len(actual))
    dump(out/"check.json", result)
    process = psutil.Process()
    peak = getattr(process.memory_info(), "peak_wset", None)
    seconds = time.monotonic()-started
    hashes = {p.name: digest(p.read_bytes()) for p in sorted(out.iterdir())}
    assert seconds < 60 and (peak is None or peak < 1024**3)
    assert sum(p.stat().st_size for p in out.iterdir()) < 20*1024**2
    dump(out/"receipt.json", dict(methodGate=result["passed"], program= digest(Path(__file__).read_bytes()),
                                  elapsedSeconds=seconds, peakProcessBytes=peak, outputs=hashes))
    print("methodGate", result["passed"], "failed cuts", len(failures))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    run(parser.parse_args().output)
