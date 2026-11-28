"""使用已有等距热网生成器为同一固定边界建立较少倾斜的子午三角网。"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from simpleLauncherHeatMesh import heat_grid


def prepare(grid, output):
    runtime = ROOT/"tools/runtime/gmsh/4.15.2"
    for entry in json.loads((runtime/"runtime.json").read_bytes())["files"]:
        assert hashlib.sha256((runtime/entry["path"]).read_bytes()).hexdigest() == entry["sha256"]
    original = ROOT/f"outputs/v350_validation/{grid}01"
    with np.load(original/"initial.npz") as data:
        xy = data["xy"]  # 仅读取几何坐标。
    wall = np.array(json.loads((original/"condition.json").read_bytes())["wall"])
    ids = np.arange(len(xy)).reshape(-1,len(wall))
    np.testing.assert_array_equal(xy[ids[0]], wall)
    points, cells, wet, far, axis, station = heat_grid(xy,ids,.04)
    np.testing.assert_allclose(station, np.arange(len(wall)), atol=1e-10, rtol=0)
    np.testing.assert_allclose(points[wet], wall, atol=1e-12, rtol=0)
    assert (points[axis,1] == 0).all()
    # 原网与新网的远界共用同一条折线；核对每个顶点在原各段上单调前进。
    chains = xy[ids[-1]]
    seg = np.diff(chains, axis=0)
    delta = points[far,None,:] - chains[:-1]
    t = np.clip((delta*seg).sum(axis=2)/(seg**2).sum(axis=1),0,1)
    dist = np.linalg.norm(delta-t[:,:,None]*seg,axis=2)
    index = dist.argmin(axis=1)
    assert dist.min(axis=1).max() < 1e-10 and min(np.diff(index+t[np.arange(len(far)),index])) > 0
    assert max(np.min(np.linalg.norm(chains[:,None]-points[far],axis=2),axis=1)) < 1e-10
    assert np.all(np.cross(points[cells[:,1]]-points[cells[:,0]],points[cells[:,2]]-points[cells[:,0]])>0)
    left = axis[points[axis,0] <= 0]
    right = axis[points[axis,0] >= wall[-1,0]]
    assert len(left)+len(right)==len(axis)
    def pairs(nodes):
        return np.column_stack((nodes[:-1],nodes[1:]))
    axis_edges = np.vstack([pairs(side[np.argsort(points[side,0])]) for side in (left,right)])
    with output.open("xb") as stream:
        np.savez_compressed(stream,xy=points,cells=cells,wall=pairs(wet),far=pairs(far),axis=axis_edges)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("grid",choices=("fine","outer"))
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    prepare(args.grid,args.output)
