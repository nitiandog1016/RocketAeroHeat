"""在单一低速基态上复用 SU2 的原生组装，读取两种重建的实际面值。"""
import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import sys

DOC = Path(__file__).resolve().parent
ROOT = DOC.parents[2]
PREVIOUS = ROOT / "outputs/v367_validation"
OUT = ROOT / "outputs/v368_validation"
sys.path.insert(0, str(DOC.parent / "v367"))
from build_observer import BIN, BUILD, invoke, replace, sha
from match_state import observe, read, save, terminal

MATCHED = PREVIOUS / "matched01"
LOCKS = {
    "build01/CNumericsSIMD.cpp": "e02810ed36c35b618a47c725010495a354a22ff30a922c3679e62cd29297b764",
    "build01/observed_upwind.hpp": "90a30cf8ba6d4c517ba9e3f2e08b282c8c0bcb4345370e7f17017ce2f0a157d9",
    "build_match01/link_process.json": "9d74a2e0f169a380eaa540b717718913669bfdc4e0c98e5b59836b6bede5e0e3",
    "matched01/case.cfg": "4fa82adc166c2ddc51c7a5d9766b041f068252d91080f026c58c849aa99ef078",
    "matched01/seed_00000.csv": "998600a5f54616298783e53f78370ecd51b203f64e639835f598d4d8a7871346",
    "matched01/mesh.su2": "f5656636b80e5a47a7d0cdb04b52f7a97e24fefb274f298ae1c31f3a42e4da68",
    "matched01/assembly_state.bin": "c5991167e43217a612451a06eeea92fa8788c541b4256db2be903be5963f2643",
    "matched01/bound_state.bin": "27819d594fbc46587fdf12fe69392a9a209039dff956cca2a38dc7c46b907d9a",
    "official_prefix01/restart_00001.csv": "7092b9a754666270f8a5918bcb2afbccdd41bbadddef0f9c6fbffc7cd362fa49",
}


def build(folder):
    """只重编读数所在的 SIMD 单元，其他对象沿用已经对齐的版本。"""
    base = PREVIOUS / "build01"
    matched = PREVIOUS / "build_match01"
    cached = json.loads(gzip.decompress((DOC.parent / "v272/evidence/build.json.gz").read_bytes()))
    assert all(sha(Path(p)) == h for p, h in cached["cachedBuildInputs"].items())
    assert sha(BIN / "g++.exe") == cached["compilerSHA256"]
    objects = read(matched / "build.json")
    assert sha(matched / "CNewtonIntegration.obj") == objects["newtonObjectSHA256"]
    assert sha(base / "CIntegration.obj") == objects["reusedObjects"]["CIntegration.obj"]
    assert sha(base / "assembly.hpp") == read(base / "build.json")["observerSHA256"]
    folder.mkdir(parents=True, exist_ok=False)
    for name in ("CNumericsSIMD.cpp", "assembly.hpp"):
        shutil.copyfile(base / name, folder / name)
    header = (base / "observed_upwind.hpp").read_text(encoding="utf8")
    marker = "/*--- Finalize in derived class (static polymorphism). ---*/"
    header = replace(header, marker, "static_assert(Double::Size <= 16);\n    " + marker)
    marker = "template<class Derived, class Base>"
    # 两个不同成员函数共享当前线程、当前向量的中心质量通量。
    header = replace(header, marker, "inline thread_local std::array<double,16> centralMassSample{};\n" + marker)
    marker = "TailBudget::central[k]=flux(3)[k];"
    header = replace(header, marker, "{ " + marker + " centralMassSample[k]=flux(0)[k]; }")
    marker = "const auto energyBeforeViscosity=flux(nDim+1);"
    header = replace(header, marker, (DOC / "face_probe.inc").read_text(encoding="utf8") + "\n    " + marker)
    (folder / "observed_upwind.hpp").write_text(header, encoding="utf8")
    cmd = read(base / "CNumericsSIMD_process.json")["args"]
    obj = folder / "CNumericsSIMD.obj"
    for flag, dest in (("-MQ", obj), ("-MF", str(obj) + ".d"), ("-o", obj)):
        cmd[cmd.index(flag) + 1] = str(dest)
    cmd[-1] = str(folder / "CNumericsSIMD.cpp")
    invoke("compile", cmd, BUILD, folder, 180, 3)
    link = read(matched / "link_process.json")["args"]
    link[link.index(str(base / obj.name))] = str(obj)
    exe = folder / "probe.exe"
    link[2] = str(exe)
    invoke("link", link, BUILD, folder, 180, 3)
    assert sum(p.stat().st_size for p in folder.iterdir()) < 300*1024**2
    save(folder / "build.json", {"executableSHA256": sha(exe), "probeSHA256": sha(DOC / "face_probe.inc"),
                                 "recipeSHA256": sha(Path(__file__))})


def run(folder, executable_folder, method):
    """绑定官方工作状态和物理上一步，仅读取一次原生基态组装。"""
    exe = executable_folder / "probe.exe"
    assert sha(exe) == read(executable_folder / "build.json")["executableSHA256"]
    terminal(MATCHED)
    folder.mkdir(parents=True, exist_ok=False)
    files = [exe, executable_folder / "build.json", DOC / "PLAN.md", DOC / "face_probe.inc", Path(__file__)]
    for name in ("seed_00000.csv", "mesh.su2", "assembly_state.bin", "case.cfg"):
        shutil.copyfile(MATCHED / name, folder / name)
        files.append(folder / name)
    if method == "green_gauss":
        with (folder / "case.cfg").open("a", encoding="ascii") as f:
            f.write("NUM_METHOD_GRAD_RECON= GREEN_GAUSS\n")
    os.environ.update(PATH=str(BIN) + os.pathsep + os.environ["PATH"],
                      RAH_M030_TAIL_BUDGET="1", RAH_M030_ASSEMBLY_STATE="1")
    observe([str(exe), "case.cfg", "--threads", "2"], folder, {str(p): sha(p) for p in files}, 150, 64*1024**2)
    terminal(folder)
    assert sha(folder / "bound_state.bin") == LOCKS["matched01/bound_state.bin"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["build", "weighted", "green_gauss"])
    parser.add_argument("folder", type=Path)
    parser.add_argument("--build-folder", type=Path)
    args = parser.parse_args()
    dest = args.folder.resolve()
    dest.relative_to(OUT)
    assert all(sha(PREVIOUS / name) == h for name, h in LOCKS.items())
    if args.action == "build":
        build(dest)
    else:
        run(dest, args.build_folder.resolve(), args.action)
