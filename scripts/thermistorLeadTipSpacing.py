"""局部二分热敏珠前引线端面的第一层，保留其余原网和共同温度采样位置。"""
import argparse
from pathlib import Path
import shutil
import sys

import numpy as np

from thermistorBead import ROOT, read, sha, write_json
from thermistorLeadMesh import meridian, sweep
from thermistorOuterDomain import load_mesh, write_mesh
from thermistorReaction import terminal

HERE = ROOT / "docs/validation/v313"
SOURCE_FILES = {"mesh.su2": "meshSHA256", "nodes.npz": "nodeMapSHA256", "case.json": "caseSHA256",
                "body.json": "bodySHA256", "case.cfg": "configSHA256", "initial_flow.csv": "initialSHA256"}


def prepare(base, seed, dest):
    """从本机已核验的同工况解构建新网和初猜，旧节点号保持可直接对照。"""
    key = read(ROOT / "cases/thermistor_b/sidewind3.json")
    assert all(sha(base / name) == key[field] for name, field in SOURCE_FILES.items())
    rec = terminal(seed)
    summary = read(seed / "surface/summary.json")
    frame = sha(seed / "restart.csv")
    assert frame == rec["outputSHA256"]["restart.csv"] == summary["sourceFrameSHA256"]
    assert frame == summary["energy"]["inputFrameSHA256"] and frame == read(seed / "field_audit.json")["lastFrame"]
    assert summary["meshSHA256"] == key["meshSHA256"] and summary["caseSHA256"] == key["caseSHA256"]
    assert max(summary["residualWindow"].values()) < -8 and max(summary["plus20Window"].values()) < -8
    assert summary["relativeChange20"] < .002 and abs(summary["energy"]["relativeDifference"]) < .005

    old, coords, marks = load_mesh(base / "mesh.su2", ("far", "wall"))
    spec = read(base / "case.json")
    assert spec["outerRadiusOverD"] == 20 and spec["velocityMps"] == [0, -3, 0] and spec["wallK"] == 299.15
    with np.load(base / "nodes.npz") as src:
        rays, wall, common, area = (src[k] for k in ("ids", "wallMap", "firstMap", "area"))
        np.testing.assert_array_equal(coords, src["xyz"])
    line, groups = meridian()
    width = len(line)
    np.testing.assert_allclose(coords[rays[0,:width],:2], line, rtol=0, atol=1e-13)
    cap = np.flatnonzero(groups == 6)[0]
    start = cap - 2
    assert np.all(groups[start-1:cap] == 5) and np.all(groups[cap:] == 6)
    # 只在端面与两排引线的法向边上取中点；轴线和接缝引用同一编号。
    node_pairs = np.stack((rays[:,start:width], rays[:,width+start:2*width]), axis=-1)
    edges, inverse = np.unique(node_pairs.reshape(-1,2), axis=0, return_inverse=True)
    mid = len(coords) + inverse.reshape(node_pairs.shape[:2])
    nodes = np.vstack((coords, coords[edges].mean(axis=1)))
    split = np.column_stack((rays[:,:2*width], mid))
    replacements = {}
    for lo, hi in zip(split[:-1], split[1:]):
        for j in range(start-1, width-1):
            a, b, c, d = j, j+1, width+j+1, width+j
            right = 2*width + (b-start)
            if j == start-1:
                parts = ([a,b,right], [a,right,c,d])
            else:
                left = right-1
                parts = ([a,b,right,left], [left,right,c,d])
            parent = tuple(sweep([a,b,c,d], lo, hi))
            replacements[parent] = [tuple(sweep(part, lo, hi)) for part in parts]
    cells = []
    affected = []
    for i, cell in enumerate(old):
        if cell in replacements:
            affected.append(i)
            cells.extend(replacements.pop(cell))
        else:
            cells.append(cell)
    assert not replacements and len(affected) == (width-start)*48
    first_by_wall = np.full(len(coords), -1, dtype=int)
    first_by_wall[wall] = common
    first_by_wall[edges[:,0]] = np.arange(len(coords), len(nodes))
    first = first_by_wall[wall]
    dest.mkdir(parents=True, exist_ok=False)
    write_mesh(dest / "mesh.su2", cells, nodes, marks)
    np.savez_compressed(dest / "nodes.npz", xyz=nodes, ids=np.column_stack((rays,mid)), wallMap=wall,
                        area=area, firstMap=first, commonMap=common, newEdges=edges, parents=affected)
    for name in ("body.json", "case.cfg"):
        shutil.copyfile(base / name, dest / name)
    # 记录基础首层及局部比例，实际法向距离由读回程序逐面核验。
    spec["baseFirstLayerM"] = spec.pop("firstLayerM")
    spec.update(frontCapFirstLayerRatio=.5, transitionLeadRows=1, splitLeadRows=2)
    write_json(dest / "case.json", spec)
    orig = np.loadtxt(seed / "restart.csv", delimiter=",", skiprows=1)
    np.testing.assert_array_equal(orig[:,0], np.arange(len(coords)))
    np.testing.assert_allclose(orig[:,1:4], coords, atol=1e-13, rtol=0)
    state = np.column_stack((np.arange(len(nodes)), nodes, np.vstack((orig[:,4:], orig[edges,4:].mean(axis=1)))))
    with (base / "initial_flow.csv").open(encoding="ascii") as f:
        header = f.readline().strip()
    np.savetxt(dest / "initial_flow.csv", state, fmt=["%d"] + ["%.17g"] * 9, delimiter=",", header=header, comments="")
    write_json(dest / "inputs.json", dict(base=str(base), seed=str(seed), seedSHA256=frame,
               generatorSHA256=sha(Path(__file__)), baseMeshSHA256=key["meshSHA256"],
               files={k: sha(dest / k) for k in SOURCE_FILES}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("base", "seed", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    dest = args.output.resolve()
    dest.relative_to(ROOT / "outputs")
    prepare(args.base.resolve(), args.seed.resolve(), dest)
