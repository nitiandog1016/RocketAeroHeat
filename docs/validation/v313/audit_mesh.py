"""独立读回细网，按父单元、物理壁和真实距离检验局部二分。"""
import argparse
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(HERE.parent / "v311"))
from audit_mesh import face_keys, normals  # v311 的有向面编码
from thermistorBead import read, sha, write_json
from thermistorLeadMesh import SIDES
from thermistorOuterDomain import load_mesh


def topology(cells, xyz, marks):
    keys, signs = [], []
    volumes = np.zeros(len(cells))
    min_cone = np.inf
    widths = np.fromiter(map(len,cells), int)
    for width, sides in SIDES.items():
        indices = np.flatnonzero(widths == width)
        if not len(indices):
            continue
        group = np.asarray([cells[i] for i in indices])
        centers = xyz[group].mean(axis=1)
        for side in sides:
            face = group[:,side]
            points = xyz[face]
            dv = np.einsum("ij,ij->i", normals(points), points.mean(axis=1)-centers)/3
            assert dv.min() > 0
            volumes[indices] += dv
            min_cone = min(min_cone, dv.min())
            k, s = face_keys(face)
            keys.append(k)
            signs.append(s)
    for faces in marks.values():
        for width in (3,4):
            group = np.asarray([f for f in faces if len(f)==width])
            if len(group):
                k, s = face_keys(group)
                keys.append(k)
                signs.append(~s)
    k, s = np.vstack(keys), np.concatenate(signs)
    order = np.lexsort(k.T[::-1])
    left, right = order[::2], order[1::2]
    assert np.array_equal(k[left], k[right]) and np.all(s[left] != s[right])
    assert np.all(np.any(k[left[1:]] != k[left[:-1]], axis=1))
    return volumes, float(min_cone)


def audit(folder):
    locks = read(folder / "inputs.json")
    base, seed = Path(locks["base"]), Path(locks["seed"])
    key = read(ROOT / "cases/thermistor_b/sidewind3.json")
    assert all(sha(folder/name) == h for name,h in locks["files"].items())
    assert sha(base/"mesh.su2") == key["meshSHA256"] and sha(base/"nodes.npz") == key["nodeMapSHA256"]
    assert sha(folder/"body.json") == key["bodySHA256"] and sha(folder/"case.cfg") == key["configSHA256"]
    assert sha(seed/"restart.csv") == locks["seedSHA256"]
    old, original, marks = load_mesh(base/"mesh.su2", ("far","wall"))
    cells, xyz, fine_marks = load_mesh(folder/"mesh.su2", ("far","wall"))
    assert marks == fine_marks
    np.testing.assert_array_equal(xyz[:len(original)], original)
    with np.load(folder/"nodes.npz") as data, np.load(base/"nodes.npz") as ref:
        wall, first, common, area, edges, parents, rays = (data[k] for k in
                ("wallMap","firstMap","commonMap","area","newEdges","parents","ids"))
        for k in ("wallMap","area"):
            np.testing.assert_array_equal(data[k],ref[k])
        np.testing.assert_array_equal(common,ref["firstMap"])
        np.testing.assert_array_equal(rays[:,:ref["ids"].shape[1]],ref["ids"])
        np.testing.assert_array_equal(xyz,data["xyz"])
    np.testing.assert_array_equal(xyz[len(original):],original[edges].mean(axis=1))
    np.testing.assert_array_equal(rays[0],rays[-1])
    off_axis = np.linalg.norm(xyz[rays[0],1:],axis=1)>0
    assert len(np.unique(rays[:-1,off_axis])) == 48*off_axis.sum()
    assert np.all(rays[:,~off_axis]==rays[0,~off_axis]) and len(np.unique(rays))==len(xyz)
    volumes, min_cone = topology(cells, xyz, fine_marks)
    # 仅重算被替换父单元的体积，并逐个核对它两个子单元的顶点来源和体积和。
    changed = set(parents.tolist())
    j = 0
    max_error = 0.
    for i, cell in enumerate(old):
        if i not in changed:
            assert cells[j] == cell
            j += 1
            continue
        vertices = set(cells[j]) | set(cells[j+1])
        assert {v for v in vertices if v<len(original)} == set(cell)
        assert all(set(edges[v-len(original)]) <= set(cell) for v in vertices if v>=len(original))
        p = original[list(cell)]
        volume = 0.
        for side in SIDES[len(cell)]:
            f = p[list(side)]
            volume += (normals(f[None,:,:])[0] @ (f.mean(axis=0)-p.mean(axis=0)))/3
        error = abs(volumes[j:j+2].sum()/volume-1)
        assert error < 2e-12
        max_error = max(max_error,error)
        j += 2
    assert j == len(cells) and abs(volumes.sum()/read(base/"mesh_audit.json")["volumeD3"]-1)<2e-13
    body, case = read(folder/"body.json"), read(folder/"case.json")
    groups = np.asarray(body["regionID"])
    tri = np.asarray(body["V"])[np.asarray(body["F"])-1]
    scale = case["diameterM"]
    np.testing.assert_allclose(scale*xyz[wall],tri,rtol=0,atol=1e-16)
    normal = normals(tri)
    np.testing.assert_allclose(np.linalg.norm(normal,axis=1), area, rtol=1e-12,atol=0)
    normal /= area[:,None]
    def distance(ids):
        return scale * np.einsum("fvi,fi->fv",xyz[ids]-xyz[wall],normal).mean(axis=1)
    coarse, fine = distance(common), distance(first)
    mask = (first!=common).any(axis=1)
    assert np.all(coarse>0) and np.all(fine>0) and np.all(mask[groups==7]) and not np.any(mask[groups<6])
    np.testing.assert_allclose(fine[groups==7],coarse[groups==7]/2,rtol=1e-12,atol=1e-18)
    assert np.all(fine[mask] < coarse[mask])
    np.testing.assert_array_equal(coarse[~mask],fine[~mask])
    initial = np.loadtxt(folder/"initial_flow.csv",delimiter=",",skiprows=1)
    prev = np.loadtxt(seed/"restart.csv",delimiter=",",skiprows=1)
    np.testing.assert_array_equal(initial[:,0],np.arange(len(xyz)))
    np.testing.assert_array_equal(initial[:,1:4],xyz)
    np.testing.assert_array_equal(initial[:len(prev),4:],prev[:,4:])
    np.testing.assert_array_equal(initial[len(prev):,4:],prev[edges,4:].mean(axis=1))
    assert np.isfinite(initial).all() and initial[:,9].min()>=-1e-12 and initial[:,9].max()<=1+1e-12
    cap, lead = groups==7, (groups==6)&mask
    def span(vals):
        return [float(vals.min()),float(vals.max())]
    stats = dict(nodes=len(xyz),cells=len(cells),addedNodes=len(edges),splitParents=len(parents),
                 unchangedCells=len(old)-len(parents),changedFirstFaces=int(sum(mask)),capFaces=int(sum(cap)),
                 adjacentLeadFaces=int(sum(lead)),changedFaceXM=span(tri[mask,:,0]),
                 capBaseNormalM=span(coarse[cap]),capFineNormalM=span(fine[cap]),
                 adjacentLeadBaseNormalM=span(coarse[lead]),adjacentLeadFineNormalM=span(fine[lead]),
                 peakBaseNormalM=coarse[10266],peakFineNormalM=fine[10266],volumeD3=volumes.sum(),
                 minCellVolumeD3=volumes.min(),minFaceConeD3=min_cone,maxParentRelativeVolumeError=max_error,
                 axisAndSeamChecked=True,faces=len(area),areaM2=area.sum(),commonPositionNodesExact=True,
                 baseMeshSHA256=key["meshSHA256"],seedSHA256=locks["seedSHA256"],files=locks["files"])
    out = folder/"mesh_audit.json"
    assert not out.exists()
    write_json(out,stats)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared",type=Path)
    audit(p.parse_args().prepared.resolve())
