"""为标准 RFZ 几何生成可复用的空气体网，或核对并打开已有网格。"""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from scipy.io import savemat

from auditAirVolume import audit
import installGmshRuntime
from meshOriginalSurfaceAir import mesh_air
from prepareRFZGeometry import CASE, DEFAULT_OUT, ROOT, STL_TYPE, digest, open_or_build


def mesh_identity(vertices, faces):
    """复算工作台的一基编号、MATLAB 列序、双精度原面指纹。"""
    fingerprint = sha256()
    for array in ([*vertices.shape, *faces.shape], vertices, faces + 1):
        fingerprint.update(np.asarray(array, dtype="<f8").tobytes(order="F"))
    return fingerprint.hexdigest()


def read_wall(stl, case):
    triangles = np.fromfile(stl, dtype=STL_TYPE, offset=84)["corners"].astype(float)
    points, first, inverse = np.unique(triangles.reshape(-1, 3), axis=0,
                                       return_index=True, return_inverse=True)
    order = np.argsort(first)
    vertices = points[order]
    faces = np.argsort(order)[inverse].reshape(-1, 3).astype(np.int32)
    assert len(faces) == case["faces"] and mesh_identity(vertices, faces) == case["meshSHA256"]
    return vertices, faces


def source_identity():
    scripts = ("prepareRFZAirVolume.py", "meshOriginalSurfaceAir.py", "auditAirVolume.py",
               "auditSU2HybridTopology.py", "prepareRFZGeometry.py", "installGmshRuntime.py")
    records = ("full_cad_input.json", "geometry.json", "cad_faces.csv",
               "deg190_boundary.json", "deg190_reference.json")
    return dict(inputs={name: digest(CASE / name) for name in records},
                scripts={name: digest(ROOT / "scripts" / name) for name in scripts},
                gmshVersion=installGmshRuntime.VERSION, wheelSHA256=installGmshRuntime.DIGEST)


def build(folder, case, geometry, identity):
    vertices, wall = read_wall(DEFAULT_OUT / "RFZ-REC_m.stl", case)
    points, cells, outer, box, quality = mesh_air(vertices, wall)
    owner, area, summary = audit(points, cells, wall, outer, box)
    assert abs(area.sum() - geometry["stlAreaM2"]) < 1e-8
    assert abs(summary["bodyVolumeM3"] - geometry["stlVolumeM3"]) < 1e-8
    asset = folder / "air.mat"
    savemat(asset, dict(vertices=points, tetrahedra=cells + 1, wall=wall + 1,
                        outer=outer + 1, box=box, wallOwner=owner + 1, cellSICN=quality), oned_as="column")
    summary.update(schema="rfz-cad-air-v1", source=identity,
                   quality=dict(minSICN=float(quality.min()), belowPointOne=int(np.sum(quality < .1))),
                   canonicalInputSHA256=identity["inputs"]["full_cad_input.json"],
                   canonicalMeshSHA256=case["meshSHA256"], fileSHA256=digest(asset),
                   stlSHA256=case["stlSHA256"], cadTableSHA256=case["cadTableSHA256"],
                   physicalFields=[], computedFaces=0, computedAreaM2=0.)
    (folder / "air.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf8")


def prepare(base):
    if not __debug__:
        raise RuntimeError("几何核对需要以常规 Python 模式运行。")
    case = json.loads((CASE / "full_cad_input.json").read_bytes())
    geometry = json.loads((CASE / "geometry.json").read_bytes())
    open_or_build(None, DEFAULT_OUT)
    source = source_identity()
    assert case["domain"] == "fullCAD" and case["stlSHA256"] == geometry["stlSHA256"]
    assert case["cadTableSHA256"] == source["inputs"]["cad_faces.csv"]
    assert case["referenceCaseSHA256"] == source["inputs"]["deg190_reference.json"]
    key = sha256(json.dumps(source, sort_keys=True).encode("utf8")).hexdigest()
    target = base / key
    base.mkdir(parents=True, exist_ok=True)
    if not target.is_dir():
        with installGmshRuntime.file_lock(base / (key + ".lock")):
            if not target.is_dir():
                work = Path(tempfile.mkdtemp(prefix="building-", dir=base))
                print("几何准备目录", work, flush=True)
                build(work, case, geometry, source)
                assert source_identity() == source
                work.rename(target)
    receipt = json.loads((target / "air.json").read_bytes())
    assert receipt["schema"] == "rfz-cad-air-v1" and receipt["source"] == source
    assert receipt["fileSHA256"] == digest(target / "air.mat")
    assert receipt["canonicalMeshSHA256"] == case["meshSHA256"]
    return target


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "outputs/rfz_air_volume")
    parser.add_argument("--receipt", type=Path, help="将此次已核对的资产路径写入新的 JSON，供工作台使用")
    args = parser.parse_args()
    base = args.out.resolve()
    base.relative_to(ROOT / "outputs")
    target = prepare(base)
    message = json.dumps(dict(directory=str(target), summarySHA256=digest(target / "air.json")))
    if args.receipt:
        with args.receipt.open("x", encoding="utf8") as stream:
            stream.write(message + "\n")
    print(message, flush=True)
