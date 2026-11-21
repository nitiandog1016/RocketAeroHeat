"""在带线热敏珠两个实际 SU2 体网上核对相同的近壁控制体、切口边和表面分带。"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from thermistorBead import CP, D, DT, RHO, U, read, sha
from thermistorLeadMesh import SIDES, meridian, original
from thermistorOuterDomain import load_mesh
from thermistorReaction import terminal

BASE = ROOT / "outputs/v305_validation"
GRIDS = (ROOT / "outputs/v300_validation/wired01", ROOT / "outputs/v304_validation/wired01")
LABELS = ["后引线端面","后引线侧面","后接头","珠面后极侧","珠面后中段","珠面前中段","珠面前极侧","前接头","前引线侧面","前引线端面"]
ENDS = np.array([4,16,22,38,54,70,86,92,104,108])
REGIONS = np.array([1,2,3,4,4,4,4,5,6,7])
SCALE = 12 * D**2 * RHO * U * CP * DT


def save(path, value):
    with path.open("x", encoding="utf8") as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)


def mapping(folder, ref, target):
    rec = terminal(folder)
    s = read(folder/"surface/summary.json")
    assert sha(folder/"restart.csv") == ref["frameSHA256"] == s["sourceFrameSHA256"]
    assert sha(folder/"mesh.su2") == ref["meshSHA256"] == s["meshSHA256"]
    assert sha(folder/"reaction0/restart.csv") == ref["reactionSHA256"] == s["sourceReactionSHA256"]
    assert sha(folder/"surface/wall.csv") == ref["wallSHA256"] == s["wallSHA256"]
    assert s["caseSHA256"] == sha(folder/"case.json") and s["bodySHA256"] == sha(folder/"body.json")
    assert rec["status"]=="terminal" and max(s["residualWindow"].values()) < -8 and abs(s["energy"]["relativeDifference"]) < .005
    cells, xyz, _ = load_mesh(folder/"mesh.su2")
    body, links, area = original()
    assert body == read(folder/"body.json")
    points = len(meridian()[0])
    assert points == 109 and ENDS[-1] == points-1
    with np.load(folder/"nodes.npz") as z:
        np.testing.assert_array_equal(xyz, z["xyz"])
        np.testing.assert_array_equal(area, z["area"])
        rings = z["ids"][:, :8*points].reshape(5, 8, points)
        wall = z["wallMap"]
    np.testing.assert_array_equal(wall, rings[links[:,:,0], 0, links[:,:,1]])
    strips = links[:,:,1].min(axis=1)
    face_group = np.searchsorted(ENDS, strips, side="right")
    np.testing.assert_array_equal(REGIONS[face_group], body["regionID"])
    # 由子午网节点和周向射线生成共同编号，轴上等价节点使用第一次出现的编号。
    shared, first = np.unique(rings.ravel(), return_index=True)
    canonical = np.full(len(xyz), -1, dtype=int)
    canonical[shared] = first
    group = np.full(len(xyz), -1, dtype=int)
    node_bands = np.minimum(len(ENDS)-1, np.searchsorted(ENDS, np.arange(points), side="right"))
    group[rings[:, :7]] = np.broadcast_to(node_bands, rings[:, :7].shape)
    assert np.array_equal(group[rings[:, :7]], np.broadcast_to(node_bands, rings[:, :7].shape))
    assert np.all(group[rings[:, 7]] == -1)
    touched = []
    edges = set()
    for cell in cells:
        ids = np.asarray(cell)
        if np.all(group[ids] == -1):
            continue
        numbered = canonical[ids]
        assert np.all(numbered >= 0)
        touched.append(tuple(numbered))
        for side in SIDES[len(cell)]:
            for a,b in zip(side, side[1:]+side[:1]):
                if group[ids[a]]>=0 or group[ids[b]]>=0:
                    edges.add(tuple(sorted((numbered[a],numbered[b]))))
    edge_list = np.array(sorted(edges))
    np.savez_compressed(target, canonical=canonical, group=group, expectedEdges=edge_list, faceBand=face_group)
    summary=dict(nodes=len(xyz), selectedNodes=int((group>=0).sum()), sharedTouchedCells=len(touched),
                 incidentEdges=len(edges), faces=len(area), meshSHA256=s["meshSHA256"], frameSHA256=s["sourceFrameSHA256"],
                 nodeMapSHA256=sha(folder/"nodes.npz"), mapSHA256=sha(target), faceTableSHA256=s["wallSHA256"])
    return xyz[rings], sorted(touched), edge_list, face_group, summary


def check(out):
    out.mkdir(parents=True,exist_ok=False)
    source=read(HERE.parent/"v304/evidence/comparison.json")["grids"][1:]
    pairs=[mapping(p, ref, out/f"grid{i}.npz") for i,(p,ref) in enumerate(zip(GRIDS,source),2)]
    a,b=pairs
    for j in (0,2,3):
        np.testing.assert_array_equal(a[j],b[j])
    assert a[1] == b[1] and a[-1]["selectedNodes"]==b[-1]["selectedNodes"]
    save(out/"geometry.json",dict(groups=LABELS,meridianStripEnds=ENDS.tolist(),grids=[a[-1],b[-1]],
         faceCount=np.bincount(a[3]).tolist(),levelOffsetsD=[0,.001,.002,.00325,.0045,.0075,.011,.015],
         scaleW=SCALE,sourceFiles={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__).resolve(),HERE/"PLAN.md",HERE.parent/"v304/evidence/comparison.json")}))
    print(a[-1]["selectedNodes"], len(a[1]), len(a[2]), np.bincount(a[3]))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("output",type=Path)
    out=p.parse_args().output.resolve();out.relative_to(BASE);check(out)
