"""在通过 96 档映射核对的带线珠体上有限推进，按同域终帧继续。"""
import argparse
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read
from thermistorSidewindExtent import run

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared",type=Path);p.add_argument("out",type=Path);p.add_argument("--steps",type=int,required=True)
    p.add_argument("--parent",type=Path)
    a=p.parse_args();folder=a.prepared.resolve();out=a.out.resolve();out.relative_to(ROOT/"outputs")
    record=read(folder/"mesh_audit.json")
    assert record["reportFaces"]==10272 and record["commonPositionsExact"] and record["nativeWallTriangles"]==20544
    assert read(folder/"case.json")["angularCells"]==96 and record["maxNativeShapeIntegralFractionError"]<2e-11
    assert 1<=a.steps<=(500 if a.parent else 440)
    files=[HERE/k for k in ("PLAN.md","NATIVE.md","audit_angular.py","run_mesh.py","solve_angular.py")]
    files.extend((folder/"inputs.json",ROOT/"scripts/thermistorAngularSpacing.py",HERE.parent/"v313/audit_mesh.py",HERE.parent/"v311/audit_mesh.py"))
    run(folder,out,a.steps,a.parent.resolve() if a.parent else None,stage_files=files)
