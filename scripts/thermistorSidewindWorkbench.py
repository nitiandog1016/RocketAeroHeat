"""为统一侧栏准备 B 珠全周侧吹体网，并从已确认的冷态运行同一求解与热量核对。"""
import argparse
from pathlib import Path
import sys

from thermistorBead import ROOT, read, sha, write_json
from thermistorLeadCrossflow import expand
from thermistorLeadWorkbench import REFINED_LEADS, check_grid, geometry
from thermistorCrossflowSurface import collect, observer_binary

sys.path.insert(0,str(ROOT/"docs/validation/v311"))
from audit_mesh import audit
from run_steady import run


def prepare(out):
    sector=ROOT/"outputs/v300_validation/wired01"
    if not sector.is_dir():
        sector=out.parent/"axial_grid"
        geometry(sector,REFINED_LEADS)
    check_grid(sector,REFINED_LEADS)
    expand(sector,out)
    info=audit(out);write_json(out/"mesh_audit.json",info)
    key=read(ROOT/"cases/thermistor_b/sidewind3.json");body=read(out/"body.json")
    for name,field in (("mesh.su2","meshSHA256"),("nodes.npz","nodeMapSHA256"),("body.json","bodySHA256"),
                       ("case.json","caseSHA256"),("initial_flow.csv","initialSHA256"),("case.cfg","configSHA256")):
        assert info["files"][name]==key[field]
    write_json(out/"prepared.json",dict(directory=str(out),nodes=info["nodes"],cells=info["cells"],angularCells=48,
               wallFaces=info["physicalFaces"]["wall"],meshSHA256=key["meshSHA256"],bodySHA256=key["bodySHA256"],
               nodeMapSHA256=key["nodeMapSHA256"],caseSHA256=key["caseSHA256"],meshAuditSHA256=sha(out/"mesh_audit.json"),
               preview=dict(V=body["V"],F=body["F"])))


def calculate(prepared,out):
    spec=read(prepared/"prepared.json")
    assert sha(prepared/"mesh_audit.json")==spec["meshAuditSHA256"]
    out.mkdir(parents=True,exist_ok=False)
    observer=observer_binary(out/"observer_build")
    first=out/"steady01";run(prepared,first,400,None,4)
    accepted=first
    if not read(first/"field_audit.json")["allResidualsBelowMinus8"]:
        accepted=out/"steady02";run(prepared,accepted,1000,first,4)
    if not read(accepted/"field_audit.json")["allResidualsBelowMinus8"]:
        raise RuntimeError("六项残差尚未全部达标，保存帧已留在本次计算目录。")
    collect(accepted,observer)
    write_json(out/"completed.json",dict(stage=accepted.name,frameSHA256=sha(accepted/"restart.csv"),
               wallSHA256=sha(accepted/"surface/wall.csv"),meshSHA256=sha(accepted/"mesh.su2")))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("action",choices=("prepare","calculate"))
    parser.add_argument("folder",type=Path);parser.add_argument("--prepared",type=Path);args=parser.parse_args()
    out=args.folder.resolve();out.relative_to(ROOT/"outputs")
    if args.action=="prepare":prepare(out)
    else:
        if args.prepared is None:parser.error("需要本次已确认的全周体网目录")
        folder=args.prepared.resolve();folder.relative_to(ROOT/"outputs");calculate(folder,out)
