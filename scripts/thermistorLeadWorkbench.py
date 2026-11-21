"""五步工作台使用的带引线热敏珠入口。建网只准备输入，显式计算才推进并审核流场。"""
import argparse
import os
from pathlib import Path
import uuid
import numpy as np

from thermistorBead import D, ROOT, read, sha, solve, write_json
from thermistorLeadMesh import WEDGE, SECTORS
from runThermistorLeadCase import configuration, prepare, probes
from thermistorLeadSurface import collect

CONTRACT=ROOT/"cases/thermistor_b/leads3.json"
MIDBODY=CONTRACT.with_name("midbody3.json")
REFINED_LEADS=CONTRACT.with_name("leadsMidbody3.json")
THIRD_LEADS=CONTRACT.with_name("leadsThird3.json")


def check_grid(folder,contract=CONTRACT):
    fixed=read(contract)
    for name,key in (("mesh.su2","meshSHA256"),("nodes.npz","nodeMapSHA256"),("body.json","bodySHA256"),
                     ("initial_flow.csv","initialSHA256"),("case.cfg","configSHA256"),("case.json","caseSHA256")):
        if sha(folder/name)!=fixed[key]:
            raise ValueError("热敏珠混合网输入与登记来源不符："+name)
    g=read(folder/"geometry.json")
    assert all(g[key]==fixed[key] for key in ("meshSHA256","nodeMapSHA256","bodySHA256","nodes","cells"))
    return fixed


def geometry(folder,contract=CONTRACT):
    prepare(folder,.001,smooth=contract==MIDBODY,midbody=contract in (MIDBODY,REFINED_LEADS),third_grid=contract==THIRD_LEADS)
    fixed=check_grid(folder,contract)
    with np.load(folder/"nodes.npz") as arr:
        wall=arr["wallMap"]
        nodes,f=np.unique(wall[:len(wall)*WEDGE//SECTORS],return_inverse=True)
        v=arr["xyz"][nodes]*D
    g=read(folder/"geometry.json")
    write_json(folder/"prepared.json",dict(directory=str(folder),nodes=g["nodes"],cells=g["cells"],
               wallFaces=g["boundaryFaces"]["wall"],angularCells=WEDGE,meshSHA256=fixed["meshSHA256"],
               nodeMapSHA256=fixed["nodeMapSHA256"],bodySHA256=fixed["bodySHA256"],
               preview=dict(V=v.tolist(),F=(f.reshape(-1,3)+1).tolist())))


def calculate(prepared,folder,contract=CONTRACT):
    check_grid(prepared,contract)
    assert (prepared/"case.cfg").read_text("ascii")==configuration()
    kit=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",ROOT/"tools/runtime/su2-thermal-build"))
    required=("su2.zip","build-v117/compile_commands.json","build-v117/build.ninja","toolchain/w64devkit/bin/g++.exe")
    if not all((kit/name).is_file() for name in required):
        raise FileNotFoundError("本工况需要固定 SU2 构建工具来核对壁面热流。请将 RAH_THERMISTOR_SU2_BUILD 指向已安装的工具目录。")
    folder.mkdir(parents=True,exist_ok=False)
    names=("mesh.su2","nodes.npz","body.json","initial_flow.csv","case.cfg","case.json","geometry.json")
    for name in names:
        os.link(prepared/name,folder/name)
    template=read(prepared/"engine_request.json")
    sources=[ROOT/s for s in template["sourceFiles"] if (ROOT/s).parent!=prepared]
    sources.extend(folder/n for n in names)
    sources.extend([contract,Path(__file__).resolve(),prepared/"prepared.json"])
    template.update(workDirectory=str(folder),jobToken=str(uuid.uuid4()),
                    sourceFiles=[p.relative_to(ROOT).as_posix() for p in sources])
    write_json(folder/"engine_request.json",template)
    check_grid(folder,contract)
    solve(folder)
    probes(folder)
    collect(folder,folder/"surface")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("action",choices=("prepare","calculate"))
    cli.add_argument("folder",type=Path)
    cli.add_argument("--prepared",type=Path)
    meshes=cli.add_mutually_exclusive_group()
    meshes.add_argument("--midbody",action="store_true",help="使用光滑热敏珠中段加密网")
    meshes.add_argument("--refine-leads",action="store_true",help="使用带线热敏珠中段加密网")
    meshes.add_argument("--third-midbody",action="store_true",help="使用带线热敏珠第三档中段网")
    args=cli.parse_args()
    folder=args.folder.resolve();folder.relative_to(ROOT/"outputs")
    contract=MIDBODY if args.midbody else CONTRACT
    if args.refine_leads:contract=REFINED_LEADS
    if args.third_midbody:contract=THIRD_LEADS
    if args.action=="prepare":geometry(folder,contract)
    else:
        if args.prepared is None:cli.error("请指定 --prepared 及已确认的体网目录")
        source=args.prepared.resolve();source.relative_to(ROOT/"outputs")
        calculate(source,folder,contract)
