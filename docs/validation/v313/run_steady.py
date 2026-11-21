"""复用侧吹原生求解和终帧核对，只接收已验收的 20D 前端面细网。"""
import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read, sha
from thermistorSidewindExtent import run


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--steps",type=int,required=True)
    parser.add_argument("--parent",type=Path)
    args=parser.parse_args()
    prepared, out = args.prepared.resolve(), args.output.resolve()
    out.relative_to(ROOT/"outputs")
    audit=read(prepared/"mesh_audit.json")
    assert read(prepared/"case.json")["outerRadiusOverD"]==20 and audit["commonPositionNodesExact"]
    assert audit["capFaces"]==336 and audit["changedFirstFaces"]==624 and audit["splitParents"]==336
    source=ROOT/"scripts/thermistorLeadTipSpacing.py"
    assert sha(source)==read(prepared/"inputs.json")["generatorSHA256"]
    assert 1<=args.steps<=(800 if args.parent else 1000)
    files=[HERE/"PLAN.md",HERE/"NATIVE.md",HERE/"audit_mesh.py",HERE/"run_mesh.py",
           Path(__file__).resolve(),source,prepared/"inputs.json",HERE.parent/"v311/audit_mesh.py"]
    run(prepared,out,args.steps,args.parent.resolve() if args.parent else None,stage_files=files)
