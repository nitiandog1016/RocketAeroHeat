"""建立固定几何的滑移球势流对照，在 WSL 内监督原生网格、求解与导出。"""
import argparse
import json
import math
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = Path.home() / ".cache/RocketAeroHeat/v362_validation"
BIN = Path("/opt/openfoam13/platforms/linux64GccDPInt32Opt/bin")
SOURCE = Path("/opt/openfoam13/applications/solvers/potentialFoam")
TOOL_SHA = {
    "potentialFoam": "a3dac7d96f7f7e25ec772aa4621f4a891de63348476c321930804c79a897df41",
    "blockMesh": "567bd6275c8b8e96f1d09d812ec0cb2c69808ad1a4e2b65a4c1be0e4423f05c6",
    "checkMesh": "6352014a4e465a46a801722e244570e96ae6523bfe3d48962327ff93b5e1bdac",
    "foamToVTK": "fd444f20d1178e4d99d4021d6fd59a18a34903ee7bab0f27ea7029788b36c53d",
}


def foam(name, body, kind="dictionary"):
    return f"FoamFile {{ version 2.0; format ascii; class {kind}; object {name}; }}\n{body}\n"


def mesh_dict():
    """每个极角面带一个径向分级块，前后轴点由两切面共用。"""
    count, radial = 128, 72
    lo, hi = 1., 2.
    for _ in range(60):
        g = (lo + hi) / 2
        if .002 * sum(g**j for j in range(radial)) < 13.5:
            lo = g
        else:
            hi = g
    grading = ((lo + hi) / 2)**(radial-1)
    points, blocks, ids = [], [], {}
    delta = math.radians(1)
    for k, radius in enumerate((.5, 14.)):
        for j in range(count+1):
            theta = math.pi*j/count
            for side in (-1, 1):
                if j in (0, count) and side == 1:
                    ids[k, j, side] = ids[k, j, -1]
                    continue
                r = 0. if j in (0, count) else radius*math.sin(theta)
                ids[k, j, side] = len(points)
                xyz = (radius*math.cos(theta), r*math.cos(delta), side*r*math.sin(delta))
                points.append("("+" ".join(f"{v:.17g}" for v in xyz)+")")
    patches = {key: [] for key in ("body", "far", "minus", "plus")}
    def face(nodes):
        return "("+" ".join(map(str, nodes))+")"
    for j in range(count):
        low = [ids[k, j+n, -1] for k, n in ((0,0),(1,0),(1,1),(0,1))]
        high = [ids[k, j+n, 1] for k, n in ((0,0),(1,0),(1,1),(0,1))]
        a,b,c,d = low
        e,f,g,h = high
        blocks.append(f"hex {face(low+high)} ({radial} 1 1) simpleGrading ({grading:.16g} 1 1)")
        patches["body"].append(face([a,e,h,d]))
        patches["far"].append(face([b,c,g,f]))
        patches["minus"].append(face(low[::-1]))
        patches["plus"].append(face(high))
    boundary = []
    for name, faces in patches.items():
        kind = {"body":"wall","far":"patch"}.get(name, "wedge")
        boundary.append(f"{name} {{ type {kind}; faces (\n"+"\n".join(faces)+"\n); }")
    return foam("blockMeshDict", "convertToMeters 1;\nvertices (\n"+"\n".join(points)+
                "\n);\nblocks (\n"+"\n".join(blocks)+"\n);\nboundary (\n"+"\n".join(boundary)+"\n);")


def initial_field(name, kind, dimensions, initial, body, far):
    return foam(name, f"dimensions {dimensions};\ninternalField uniform {initial};\nboundaryField\n{{\n"
                f"body {{ {body} }}\nfar {{ {far} }}\nminus {{ type wedge; }}\nplus {{ type wedge; }}\n}}", kind)


def prepare(job):
    case = job / "case"
    (case / "system").mkdir(parents=True)
    (case / "constant").mkdir()
    initial = job / "initial"
    initial.mkdir()
    defs = {
        "blockMeshDict": mesh_dict(),
        "controlDict": foam("controlDict", "application potentialFoam; startFrom startTime; startTime 0; "
                            "stopAt endTime; endTime 1; deltaT 1; writeControl timeStep; writeInterval 1; "
                            "writeFormat ascii; writePrecision 16; writeCompression off; runTimeModifiable false;"),
        "fvSchemes": foam("fvSchemes", """ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss linear; div(div(phi,U)) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }"""),
        "fvSolution": foam("fvSolution", """solvers
{ "(Phi|p)" { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; maxIter 10000; } }
potentialFlow { nNonOrthogonalCorrectors 4; PhiRefCell 0; PhiRefValue 0; }"""),
    }
    for name, contents in defs.items():
        (case / "system" / name).write_text(contents, encoding="ascii")
    for name, contents in {
        "U": initial_field("U","volVectorField","[0 1 -1 0 0 0 0]","(-1 0 0)","type slip;","type fixedValue; value uniform (-1 0 0);"),
        "Phi": initial_field("Phi","volScalarField","[0 2 -1 0 0 0 0]","0","type zeroGradient;","type zeroGradient;"),
        "p": initial_field("p","volScalarField","[0 2 -2 0 0 0 0]","0","type zeroGradient;","type fixedValue; value uniform 0;"),
    }.items():
        (initial / name).write_text(contents, encoding="ascii")
    shutil.copytree(initial, case / "0")


def execute(job):
    """每阶段继承只读输入锁和固定资源上限；失败回执与完整运行目录原样保留。"""
    sys.path.insert(0, str(ROOT/"scripts"))
    from runOpenFoamProcess import run
    from runEilmerProcess import sha
    for name, expected in TOOL_SHA.items():
        if sha(BIN/name) != expected:
            raise ValueError(f"{name} 程序指纹不符")
    case = job / "case"
    locks = [HERE/"PLAN.md", HERE/"RECOVERY.md", Path(__file__).resolve(), SOURCE/"potentialFoam.C", SOURCE/"createFields.H",
             *sorted((case/"system").iterdir()), *sorted((job/"initial").iterdir())]
    commands = [("blockMesh", []), ("checkMesh", []),
                ("potentialFoam", ["-initialiseUBCs","-writePhi","-writep"]),
                ("foamToVTK", ["-ascii","-noInternal","-noPointValues","-fields","(p U)","-time","0"])]
    for name, args in commands:
        record = job / name
        if name == "checkMesh":
            locks.extend(p for p in (case/"constant/polyMesh").iterdir() if p.is_file())
        if name == "potentialFoam":
            if (case/"0/phi").exists():
                raise ValueError("初场只含速度势 Phi；请选择能分别保存面流量 phi 的目录")
            for p in (job/"initial").iterdir():
                assert sha(p) == sha(case/"0"/p.name)
        if name == "foamToVTK":
            if (case/"0/phi").samefile(case/"0/Phi"):
                raise ValueError("面通量与势文件发生名称冲突")
            locks.extend(case/"0"/field for field in ("U","p","phi","Phi"))
        receipt = run([BIN/name, *args],case,record,180,locks=[*locks,BIN/name],output_mib=96,memory_mib=2048)
        log = (record/"stdout.txt").read_text()
        if not receipt["passed"] or not log.rstrip().endswith("End"):
            raise RuntimeError(f"{name} 原生阶段未正常结束，见 {record}")
        if name == "checkMesh" and "Mesh OK." not in log:
            raise RuntimeError(f"网格未通过，见 {record}")
        print(json.dumps({"stage":name,"seconds":receipt["elapsedSeconds"]}), flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder",type=Path)
    folder=parser.parse_args().folder.resolve()
    folder.relative_to(BASE)
    folder.mkdir(parents=True,exist_ok=False)
    prepare(folder)
    execute(folder)
