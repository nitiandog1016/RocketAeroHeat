"""从独立写出的四张体网定位局部变化，确认带线新网可送入求解。"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/"v299"))
from check_grid import ROOT, A, D, cells, crossings, extent, meridian, read, save, sha, terminal
from check_mesh import shared_faces

WIRED = ROOT/"outputs/v297_validation/fine01"
SMOOTH = ROOT/"outputs/v298_validation/smooth01"
REF_SMOOTH = ROOT/"outputs/v299_validation/smooth01"
NEW = ROOT/"outputs/v300_validation/wired01"


def middle(cells_to_test):
    return {c for c in cells_to_test if all(abs(v[0])*D/A <= .5 for v in c[1])}


def check():
    base = read(HERE.parent/"v297/evidence/fine_inputs.json")
    assert all(sha(WIRED/k)==v for k,v in base["inputs"].items())
    for folder, version in ((SMOOTH,"v298"), (REF_SMOOTH,"v299")):
        recorded = read(HERE.parent/version/"evidence/inputs.json")["inputs"]
        assert all(sha(folder/k)==v for k,v in recorded.items())
    for folder in (WIRED, SMOOTH, REF_SMOOTH):
        terminal(folder)
    prep = read(NEW/"debugger_process.json")
    assert prep["status"]=="terminal" and prep["debuggerExitCode"]==0 and all(prep["inputsStable"].values())
    assert not (prep["stopReason"] or prep["supervisorError"] or prep["survivingProcessIDs"])
    inputs = {k:sha(NEW/k) for k in base["inputs"]}
    for k in ("case.cfg","body.json","case.json"):
        assert inputs[k]==base["inputs"][k]
    mapping = shared_faces(read(NEW/"body.json"),read(REF_SMOOTH/"body.json"))
    with np.load(NEW/"nodes.npz") as n, np.load(WIRED/"nodes.npz") as w, np.load(REF_SMOOTH/"nodes.npz") as s:
        for k in ("wallMap","firstMap"):
            np.testing.assert_array_equal(n["xyz"][n[k]],w["xyz"][w[k]])
        np.testing.assert_array_equal(n["area"],w["area"])
        np.testing.assert_array_equal(n["area"][mapping[:,0]],s["area"][mapping[:,1]])
        common_area=float(n["area"][mapping[:,0]].sum())
    old,new,sm,ref = (set.union(*cells(f)) for f in (WIRED,NEW,SMOOTH,REF_SMOOTH))
    gone,added=old-new,new-old
    core_old,core_new=middle(gone),middle(added)
    assert core_old and core_new and core_old==middle(sm-ref) and core_new==middle(ref-sm)
    wall,regions=meridian()
    # 接头最内侧的轴向位置取自真实母线；全部变化顶点须处于两接头之间。
    joints=np.flatnonzero(np.isin(regions,[2,4]))
    joint_start=float(np.min(abs(wall[np.r_[joints,joints+1],0]))*D/A)
    removed,inserted=extent(gone,wall),extent(added,wall)
    for changed in (removed,inserted):
        assert max(map(abs,changed["xOverA"])) < joint_start
        assert changed["vertexWallDistanceD"][0]>.014 and changed["vertexWallDistanceD"][1]<.25
    coarse,fine=crossings(WIRED),crossings(NEW)
    reference=read(HERE.parent/"v299/evidence/inputs.json")["refined"]
    for c,f,r in zip(coarse,fine,reference):
        if f["degFromUpstream"] in (60,85,120):
            assert f["maxCrossingIntervalOverD"] < c["maxCrossingIntervalOverD"]
            np.testing.assert_allclose(f["maxCrossingIntervalOverD"],r["maxCrossingIntervalOverD"],atol=.004,rtol=0)
    req=read(NEW/"engine_request.json")
    assert req["threads"]==1 and req["maximumSeconds"]==720
    save(HERE/"evidence/inputs.json",dict(geometry=read(NEW/"geometry.json"),inputs=inputs,coarse=coarse,refined=fine,
         smoothRefined=reference,removed=removed,added=inserted,unchangedCells=len(old&new),
         pairedCoreParents=len(core_old),pairedCoreChildren=len(core_new),matchedParentsAll=len(gone&(sm-ref)),
         matchedChildrenAll=len(added&(ref-sm)),jointInnerAbsXOverA=joint_start,
         commonFaces=len(mapping),commonAreaM2=common_area,sources={p:sha(ROOT/p) for p in req["sourceFiles"]}))


if __name__=="__main__":
    check()
