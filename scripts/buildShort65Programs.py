"""将银球已验收的热壁与端点取样模块链接到本安装从上游源码生成的库。"""
import json
from pathlib import Path
import re
import sys
import tempfile
from zipfile import ZipFile

from installShort65Kernel import ROOT, HOME_DIR, LIBS, DIGESTS, digest
from buildThermistorWireDriver import MAIN, APPLY
from runEmbeddedSU2 import atomic_write_json
sys.path.insert(0,str(ROOT/"docs/validation/v333"))
from build_time import integration, scalar, header
from build_shared import SPECIES, SCALAR, FLOW, write_sources, arguments, invoke, once

INTEGRATION="SU2_CFD/src/integration/CIntegration.cpp"
WRITER="SU2_CFD/src/output/filewriter/CSU2FileWriter.cpp"


def check_kernel(prepared):
    sdk=HOME_DIR;src=sdk/"source";wc=sdk/"build";cc=sdk/"toolchain/w64devkit/bin/g++.exe"
    kernel=json.loads((sdk/"kernel.json").read_bytes());p=json.loads((prepared/"prepared.json").read_bytes())
    assert kernel["schema"]=="short65-kernel-v1" and kernel["upstreamCommit"]=="7048fd35f56d4bf52e586429d6746a36a10255f8"
    assert kernel["archives"] == {name: value[0] for name, value in DIGESTS.items()}
    assert digest(cc)==kernel["compiler"]=="315c23bc298aba9eaf2aece174fa7fcb0e404247176663e6a4c7707a60f8e444"
    assert all(digest(wc/k)==kernel["libraries"][k] for k in LIBS)
    assert p["files"]["mesh.su2"]==digest(prepared/"mesh.su2")=="a56a4c3dee8fd8bb74bc23ae4e85b71eaf7d5552dd99fccf8a6c8cb6a40789cd"
    assert (p["geometry"]["nodes"],p["geometry"]["farNodes"])==(261708,2690)
    # 用同一上游归档复查参与替换的文件；共用历史版本的生成函数与物理实现。
    archive=sdk/"downloads/su2.zip";assert digest(archive)==DIGESTS["su2"][0]
    with ZipFile(archive) as z:
        for rel in (MAIN,INTEGRATION,FLOW,SPECIES,SCALAR,WRITER):
            member,=[k for k in z.namelist() if k.endswith("/"+rel)]
            assert z.read(member)==(src/rel).read_bytes()
    return p


def generate_sources(out,src,p):
    """从固定上游生成适配源码，保存文件按双精度往返所需的位数写出。"""
    out.mkdir(parents=True,exist_ok=False)
    steady=out/"steady";time=out/"time";steady.mkdir();time.mkdir()
    write_sources(steady,src,p["geometry"]["nodes"],p["geometry"]["farNodes"])
    bc=once(APPLY,'  std::ifstream boundary(', '  if (count != 11906) throw std::runtime_error("Short 65 wall size");\n  std::ifstream boundary(')
    (steady/"SU2_CFD.cpp").write_text('#include <fstream>\n'+once((src/MAIN).read_text("utf8"),"  driver->StartSolver();",bc+"  driver->StartSolver();"),encoding="utf8")
    writer = once((src/WRITER).read_text("utf8"), "  data.precision(15);",
                  "  data.precision(std::numeric_limits<double>::max_digits10 - 1);")
    (steady/"CSU2FileWriter.cpp").write_text("#include <limits>\n"+writer, encoding="utf8")
    (time/"CIntegration.cpp").write_text(integration((src/INTEGRATION).read_text("utf8")),encoding="utf8")
    (time/"CScalarSolver.inl").write_text(scalar((steady/"CScalarSolver.inl").read_text("utf8")),encoding="utf8")
    (time/"silver.hpp").write_text(header((steady/"silver.hpp").read_text("utf8")),encoding="utf8")
    for name in ("CSpeciesSolver.cpp","far_mass.hpp","capacity.hpp"):(time/name).write_bytes((steady/name).read_bytes())
    return steady,time


def check_programs(native, prepared):
    """复用缓存前按现行算法重新生成小型源码，同时核对库、目标文件和四个可执行文件。"""
    spec = check_kernel(prepared)
    steady = json.loads((native / "steady/build.json").read_bytes())
    time = json.loads((native / "time/build.json").read_bytes())
    assert steady["kernel"] == digest(HOME_DIR / "kernel.json")
    assert steady["preparedSHA256"] == digest(prepared / "prepared.json")
    assert time["baseline"] == digest(native / "steady/build.json")
    with tempfile.TemporaryDirectory(dir=native.parent) as folder:
        fresh = Path(folder) / "generated"
        generate_sources(fresh, HOME_DIR / "source", spec)
        for name, description in (("steady", steady), ("time", time)):
            files = description["files"]
            sources = {p.name for p in (fresh / name).iterdir()}
            assert sources == {f for f in files if Path(f).suffix in (".cpp", ".hpp", ".inl")}
            for f in sources:
                if digest(fresh / name / f) != files[f]:
                    raise ValueError("银球程序的源码已经改变，请为当前版本准备新的计算缓存：" + f)
            for f, expected in files.items():
                if digest(native / name / f) != expected:
                    raise ValueError("银球程序缓存与编译记录不符：" + f)


def build(out,prepared):
    sdk=HOME_DIR;src=sdk/"source";wc=sdk/"build";bin_dir=sdk/"toolchain/w64devkit/bin";cc=bin_dir/"g++.exe"
    spec=check_kernel(prepared)
    steady,time=generate_sources(out,src,spec)
    entries=json.loads((wc/"compile_commands.json").read_bytes())
    def compile(folder,rel,extra=()):
        entry,=[e for e in entries if e["file"].replace("\\","/").endswith("/"+rel)]
        argv=arguments(entry["command"]);argv[0]=cc;obj=folder/(Path(rel).stem+".obj")
        for flag,v in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):argv[argv.index(flag)+1]=v
        argv[-1:]=["-iquote",(src/rel).parent,*extra,folder/Path(rel).name]
        invoke(obj.stem,argv,wc,folder,180,3,bin_directory=bin_dir);return obj
    # 编译器在实际 .cpp 所在目录优先查找本次生成的 .inl 和头文件。
    main=compile(steady,MAIN);flow=compile(steady,FLOW);writer=compile(steady,WRITER)
    include=("-iquote",(src/SCALAR).parent)
    csteady=compile(steady,SPECIES,include);ci=compile(time,INTEGRATION);ctime=compile(time,SPECIES,include)
    line=re.search(r"^build SU2_CFD/src/SU2_CFD\.exe: cpp_LINKER [^\n]+\n LINK_ARGS = ([^\n]+)",(wc/"build.ninja").read_text("utf8"),re.M)
    assert line;link_args=arguments(line[1])
    for folder,name,objs in ((steady,"wire",[main]),(steady,"shared",[main,csteady,flow]),
                             (time,"physical",[main,ci,ctime,flow]),(time,"release",[main,ci])):
        invoke(name,[cc,"-o",folder/(name+".exe"),*objs,writer,*link_args],wc,folder,180,3,bin_directory=bin_dir)
    def files(folder):
        return {q.name:digest(q) for q in folder.iterdir() if q.suffix in (".hpp",".cpp",".inl",".exe",".obj")}
    atomic_write_json(steady/"build.json",dict(schema="short65-source-build-v1",kernel=digest(sdk/"kernel.json"),
                      preparedSHA256=digest(prepared/"prepared.json"),exe=digest(steady/"shared.exe"),files=files(steady)))
    atomic_write_json(time/"build.json",dict(baseline=digest(steady/"build.json"),files=files(time)))
