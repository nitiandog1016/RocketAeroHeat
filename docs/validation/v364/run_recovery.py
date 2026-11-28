"""保留首轮工况与原始字典，用独立目录执行更换子午离散的对照。"""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
from run_outer import BASE, BIN, GRIDS, HERE, ROOT, TOOL_SHA, flight, read_meridian, run, sha, write_wedge


def recovery_grid(grid):
    path=ROOT/f"outputs/v364_validation/{grid}.npz"
    with np.load(path) as data:
        xy, cells = data["xy"], data["cells"]
        edges = {key:data[key] for key in ("wall","far","axis")}
    source=ROOT/f"outputs/v350_validation/{grid}01/mesh.su2"
    assert sha(source)==GRIDS[grid]
    old, _, markers=read_meridian(source)
    np.testing.assert_allclose(xy[edges["wall"]],old[markers["wall"]],atol=1e-12,rtol=0)
    return path, xy, cells, edges


def execute(grid, original, job):
    report=json.loads((original/"potentialFoam/report.json").read_bytes())
    assert report["passed"] and sha(BIN/"potentialFoam")==TOOL_SHA["potentialFoam"]
    locks=[Path(p) for p in report["inputSHA256"]]
    assert all(sha(p)==report["inputSHA256"][str(p)] for p in locks)
    mesh_file, xy, cells, edges=recovery_grid(grid)
    job.mkdir(parents=True,exist_ok=False)
    case=job/"case"
    write_wedge(case/"constant/polyMesh",xy*flight()[1]["maxDiameterM"],cells,edges)
    shutil.copytree(original/"case/system",case/"system")
    shutil.copytree(original/"initial",job/"initial")
    shutil.copytree(job/"initial",case/"0")
    locks += [HERE/"RECOVERY.md",HERE/"remesh.py",ROOT/"scripts/simpleLauncherHeatMesh.py",Path(__file__),mesh_file,
              *(case/"system").iterdir(),*(job/"initial").iterdir(),*(case/"constant/polyMesh").iterdir()]
    for name,args in (("checkMesh",[]),("potentialFoam",["-initialiseUBCs","-writePhi"])):
        assert sha(BIN/name)==TOOL_SHA[name]
        receipt=run([BIN/name,*args],case,job/name,180,[*locks,BIN/name])
        log=(job/name/"stdout.txt").read_text()
        if not receipt["passed"] or not log.rstrip().endswith("End") or (name=="checkMesh" and "Mesh OK." not in log):
            raise RuntimeError(f"需要检查 {job/name}")
        print(f"{grid} {name} {receipt['elapsedSeconds']:.2f} s",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("grid",choices=GRIDS)
    parser.add_argument("original",type=Path)
    parser.add_argument("folder",type=Path)
    args=parser.parse_args()
    for p in (args.original,args.folder):
        p.resolve().relative_to(BASE)
    execute(args.grid,args.original,args.folder)
