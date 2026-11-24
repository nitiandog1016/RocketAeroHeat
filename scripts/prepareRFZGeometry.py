"""从作者固定版本的 RFZ-REC STEP 生成米制工作台输入，并保留每个 CAD 面的三角面号。"""
import argparse
import csv
import hashlib
import json
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np

import installGmshRuntime

ROOT = Path(__file__).resolve().parents[1]
COMMIT = "3a174ee3d0b0626d6029e9608d862b674f3b266b"
STEP_SHA = "0442f31ef5cc0ffb0676879c17299da68245b35a3ac104d82d7d46419f2033db"
STEP_BYTES = 14081331
URL = f"https://raw.githubusercontent.com/rfz-model/RFZ_data/{COMMIT}/CAD/stp/RFZ-REC.step"
LOCAL_STEP = ROOT / "outputs/rfz_reference/author_data_20261124/files/CAD/stp/RFZ-REC.step"
DEFAULT_OUT = ROOT / "outputs/rfz_rec_geometry"
CASE = ROOT / "cases/rfz_rec"
STL_TYPE = np.dtype([("normal", "<f4", (3,)), ("corners", "<f4", (3, 3)), ("attribute", "<u2")])


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_file(requested):
    if requested is not None:
        source = requested
    elif LOCAL_STEP.is_file():
        source = LOCAL_STEP
    else:
        source = ROOT / "outputs/rfz_reference/RFZ-REC.step"
        source.parent.mkdir(parents=True, exist_ok=True)
        if not source.is_file():
            with installGmshRuntime.file_lock(source.with_suffix(".lock")):
                if not source.is_file():
                    installGmshRuntime.download_pinned(URL, source, source.with_suffix(".part"), STEP_BYTES, STEP_SHA)
    if source.stat().st_size != STEP_BYTES or digest(source) != STEP_SHA:
        raise ValueError("输入 STEP 与已核对的作者版本不符。")
    return source


def topology(corners):
    """按实际写入 STL 的单精度坐标合并共享顶点，检查所有边、绕序、连通性及重复面。"""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    vertices, ix = np.unique(corners.reshape(-1, 3), axis=0, return_inverse=True)
    faces = ix.reshape(-1, 3)
    assert len(np.unique(np.sort(faces, axis=1), axis=0)) == len(faces), "重复三角面"
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    undirected, links, counts = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)
    assert np.all(counts == 2), "开放边或非流形边"
    assert np.all(np.bincount(links, weights=np.sign(edges[:, 1] - edges[:, 0])) == 0), "相邻面绕序不一致"
    graph = coo_matrix((np.ones(len(undirected)), undirected.T), shape=(len(vertices), len(vertices)))
    assert connected_components(graph, directed=False, return_labels=False) == 1, "表面不连通"
    return len(vertices), len(undirected)


def prepare(source, out):
    # 建面仅使用 CAD。各曲面的三角面按 Gmsh 单线程生成的单元号顺序写入。
    installGmshRuntime.main()
    sys.path.insert(0, str(ROOT / "tools/runtime/gmsh" / installGmshRuntime.VERSION))
    import gmsh

    gmsh.initialize(readConfigFiles=False)
    try:
        for name, value in (("General.NumThreads", 1), ("Mesh.Algorithm", 6), ("Mesh.MeshSizeMin", .035),
                            ("Mesh.MeshSizeMax", .45), ("Mesh.MeshSizeFromCurvature", 30)):
            gmsh.option.setNumber(name, value)
        gmsh.option.setString("Geometry.OCCTargetUnit", "M")
        gmsh.model.occ.importShapes(str(source))
        gmsh.model.occ.synchronize()
        solids = gmsh.model.getEntities(3)
        assert len(solids) == 1
        surfaces = sorted(abs(tag) for dim, tag in gmsh.model.getBoundary(solids) if dim == 2)
        assert surfaces == [tag for _, tag in sorted(gmsh.model.getEntities(2))] and len(surfaces) == 204
        gmsh.model.mesh.generate(2)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        order = np.argsort(tags)
        node_ids = tags[order]
        points = xyz.reshape(-1, 3)[order]
        rows, blocks = [], []
        first = 1
        for surface in surfaces:
            kinds, elements, nodes = gmsh.model.mesh.getElements(2, surface)
            assert list(kinds) == [2] and len(elements[0]) > 0
            ids = nodes[0].reshape(-1, 3)[np.argsort(elements[0])]
            idx = np.searchsorted(node_ids, ids)
            assert np.array_equal(node_ids[idx], ids)
            blocks.append(points[idx].astype("<f4"))
            rows.append([surface, first, first + len(ids) - 1, gmsh.model.getType(2, surface),
                         gmsh.model.occ.getMass(2, surface)])
            first += len(ids)
        cad_volume = gmsh.model.occ.getMass(*solids[0])
    finally:
        gmsh.finalize()

    corners = np.concatenate(blocks).astype(float)
    nodes, edges = topology(corners)
    cross = np.cross(corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0])
    twice_area = np.linalg.norm(cross, axis=1)
    assert np.all(np.isfinite(corners)) and np.all(twice_area > 0)
    volume = np.einsum("ij,ij->", corners[:, 0], cross) / 6
    area = float(np.sum(twice_area) / 2)
    cad_area = sum(row[4] for row in rows)
    assert abs(area / cad_area - 1) < .01 and volume > 0 and abs(volume / cad_volume - 1) < .01
    stl = out / "RFZ-REC_m.stl"
    records = np.zeros(len(corners), dtype=STL_TYPE)
    records["corners"] = corners
    records["normal"] = cross / twice_area[:, None]
    with stl.open("xb") as file:
        file.write(b"DLR RFZ-REC | SI metres | CAD boundary | RocketAeroHeat".ljust(80, b" "))
        file.write(struct.pack("<I", len(records)))
        records.tofile(file)
    table = out / "cad_faces.csv"
    with table.open("x", encoding="utf8", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["cadFace", "firstFace", "lastFace", "surfaceType", "cadAreaM2", "stlAreaM2"])
        for row in rows:
            writer.writerow([*row, float(twice_area[row[1]-1:row[2]].sum() / 2)])
    summary = dict(schema="rfz-rec-geometry-v1", sourceURL=URL, stepSHA256=STEP_SHA,
                   gmshVersion=installGmshRuntime.VERSION, unit="m", cadFaces=len(rows),
                   faces=len(records), vertices=nodes, edges=edges, components=1,
                   boundsM=[corners.min(axis=(0, 1)).tolist(), corners.max(axis=(0, 1)).tolist()],
                   stlAreaM2=area, cadAreaM2=cad_area, stlVolumeM3=float(volume), cadVolumeM3=cad_volume,
                   stlSHA256=digest(stl), regionsSHA256=digest(table))
    with (out / "geometry.json").open("x", encoding="utf8", newline="\n") as f:
        json.dump(summary, f, indent=2, allow_nan=False)
        f.write("\n")
    print(json.dumps(summary), flush=True)


def verify_case(out):
    expected = json.loads((CASE / "geometry.json").read_bytes())
    actual = json.loads((out / "geometry.json").read_bytes())
    if actual != expected or digest(out / "RFZ-REC_m.stl") != expected["stlSHA256"] or digest(out / "cad_faces.csv") != expected["regionsSHA256"]:
        raise ValueError(f"几何或 CAD 面清单与冻结版本不符：{out}")


def open_or_build(step, out):
    if out.exists():
        if step is not None:
            source_file(step)
        verify_case(out)
        print(f"已核对 RFZ-REC 几何：{out}", flush=True)
        return
    # 新计算写在独立目录；只有全部文件与冻结版本一致时，才发布到工作台路径。
    out.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=out.name + "-", dir=out.parent))
    print(f"几何准备目录 {work}", flush=True)
    prepare(source_file(step), work)
    verify_case(work)
    try:
        work.rename(out)
    except FileExistsError:
        verify_case(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, help="作者原始 RFZ-REC STEP 文件")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="新的几何算例目录")
    args = parser.parse_args()
    open_or_build(args.step, args.out.resolve())
