"""从固定的上游包在本项目内构建 Short Test 65 的 Windows SU2 静态库。"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import urllib.request
from zipfile import ZipFile

from prepareSU2ThermalBuild import SOURCES
from runEmbeddedSU2 import atomic_write_json
from runThermistorMeasuredCold import exclusive

ROOT=Path(__file__).resolve().parents[1]
HOME_DIR=ROOT/"tools/runtime/short65-14d"
sys.path.insert(0,str(ROOT/"docs/validation/v272"))
from observe_flux import invoke

DIGESTS={
    "compiler":("54c2bb9962e8a30c8270e34c5efc3d5fee579f53087fa7d435b95fbba2af7b2b",68801240),
    "su2":("c45f77e16af4ce326e0fca3adf21009d61e4ecc0d6cac40e8dae22654e1c83dc",21404534),
    "meson":("4e2b6d43a766cc7f8ea2922339b358f9f325cc8e3dd406debb16a83382e9b532",7237778),
    "mel":("03c2b721d8d2d0aff26d02bee6a92619f2aaf60708c5fafd0268bb6fde1dd100",18560),
    "eigen":("d42e8453beab5c8d7b8bdc32ba7452604f711f64b53b0a11991ee13ea08dc918",3857918)}
LIBS=("SU2_CFD/src/libSU2core.a","Common/src/libSU2Common.a")
WINDOWS_HEADERS={
    "SU2_CFD/src/interfaces/fsi/CFlowTractionInterface.cpp":"interface_interpolation/CInterpolator.hpp",
    "SU2_CFD/include/interfaces/cfd/CMixingPlaneInterface.hpp":"containers/C2DContainer.hpp",
    "SU2_CFD/src/interfaces/cfd/CMixingPlaneInterface.cpp":"interface_interpolation/CInterpolator.hpp",
}


def digest(file):
    h=hashlib.sha256()
    with file.open("rb") as stream:
        for part in iter(lambda:stream.read(1024*1024),b""):h.update(part)
    return h.hexdigest()


def acquire(name,url,cache,dest):
    file=dest/(name+(".7z.exe" if name=="compiler" else ".zip"))
    expected,size=DIGESTS[name]
    if not file.exists():
        source=cache/file.name if cache else None
        with tempfile.NamedTemporaryFile(dir=str(dest),delete=False) as stream:
            tmp=Path(stream.name)
            try:
                reader=source.open("rb") if source and source.is_file() else urllib.request.urlopen(url,timeout=60)
                with reader:
                    while True:
                        block=reader.read(1024*1024)
                        if not block:break
                        stream.write(block)
                        if stream.tell()>size:raise ValueError("上游包超过登记大小："+name)
                stream.close()
                if tmp.stat().st_size!=size or digest(tmp)!=expected:raise ValueError("上游包内容与固定版本不符："+name)
                tmp.replace(file)
            finally:
                stream.close()
                if tmp.exists():tmp.unlink()
    if file.stat().st_size!=size or digest(file)!=expected:raise ValueError("缓存校验失败："+name)
    return file


def extract(archive,dest,skip_tests=False,fix_headers=False):
    mark=dest/".short65-source.json";key=digest(archive)
    if mark.is_file() and json.loads(mark.read_bytes())["archive"]!=key:
        raise ValueError("暂存源码来自其他归档："+str(dest))
    dest.mkdir(parents=True,exist_ok=True)
    with ZipFile(archive) as z:
        top=z.namelist()[0].split("/")[0]
        for item in z.infolist():
            rel=Path(item.filename).relative_to(top)
            if item.is_dir() or not rel.parts or (skip_tests and rel.parts[0] in ("tests","test cases")):continue
            target=(dest/rel).resolve();target.relative_to(dest.resolve())
            assert (item.external_attr>>16)&0o170000!=0o120000
            target.parent.mkdir(parents=True,exist_ok=True)
            original=z.read(item);data=original
            header=WINDOWS_HEADERS.get(rel.as_posix()) if fix_headers else None
            if header:
                old=('"../../../Common/include/'+header+'"').encode()
                assert original.count(old)==1
                data=original.replace(old,('"../../../../Common/include/'+header+'"').encode())
            if target.exists():
                actual=target.read_bytes()
                if actual==data:continue
                if actual!=original:raise ValueError("源码暂存目录存在不同内容："+str(rel))
            target.write_bytes(data)
    atomic_write_json(mark,dict(archive=key))


def install(cache=None):
    home=HOME_DIR;home.mkdir(parents=True,exist_ok=True)
    with exclusive(home):
        receipts=None
        packages=home/"downloads";packages.mkdir(exist_ok=True)
        archives={name:acquire(name,url,cache,packages) for name,url,_,_ in SOURCES}
        src=home/"source";build=home/"build";bin_dir=home/"toolchain/w64devkit/bin"
        if not (bin_dir/"g++.exe").is_file():
            receipts=Path(tempfile.mkdtemp(prefix="build_",dir=str(home)))
            tool=home/"toolchain";tool.mkdir(exist_ok=True)
            unpack=shutil.which("tar")
            if not unpack:raise FileNotFoundError("需要 Windows 自带的 tar 来解包编译器。")
            invoke("extract_compiler",[unpack,"-xf",archives["compiler"],"-C",tool],tool,receipts,180,2,bin_directory=tool)
        compiler=bin_dir/"g++.exe";ninja=bin_dir/"ninja.exe"
        for f in (compiler,ninja):
            if not f.is_file():raise FileNotFoundError("编译工具未完整解包："+str(f))
        assert digest(compiler)=="315c23bc298aba9eaf2aece174fa7fcb0e404247176663e6a4c7707a60f8e444"
        assert digest(ninja)=="c00843865feeb97d51aed8c8dadd4c3d9beeb1e3f045001d73ee0ec42cff0146"
        extract(archives["su2"],src,fix_headers=True);extract(archives["meson"],home/"meson",True)
        for name in ("eigen","mel"):extract(archives[name],src/"externals"/name)
        record_file=home/"kernel.json"
        if record_file.is_file():
            record=json.loads(record_file.read_bytes())
            assert record["schema"]=="short65-kernel-v1" and record["upstreamCommit"]=="7048fd35f56d4bf52e586429d6746a36a10255f8"
            assert record["archives"]=={name:v[0] for name,v in DIGESTS.items()}
            assert record["compiler"]==digest(compiler) and record["ninja"]==digest(ninja)
            assert all(record["libraries"][name]==digest(build/name) for name in LIBS)
            # 已通过的库保留原编译回执，重验安装不会使已有基流和专用程序失效。
            return record
        if receipts is None:receipts=Path(tempfile.mkdtemp(prefix="build_",dir=str(home)))
        if not (src/"su2preconfig.timestamp").exists():
            args=["--no-"+k for k in ("codi","medi","opdi","mpp","coolprop","mel","fado","mlpcpp","eigen")]
            invoke("preconfigure",[sys.executable,"preconfigure.py",*args],src,receipts,180,2,bin_directory=bin_dir)
        os.environ.update(CC=str(bin_dir/"gcc.exe"),CXX=str(compiler))
        if not (build/"build.ninja").is_file():
            opts=["-Dwith-mpi=disabled","-Dwith-omp=true","-Denable-cgns=false","-Denable-tecio=false",
                  "-Denable-pywrapper=false","-Denable-tests=false","-Dcpu-arch=x86-64","--buildtype=release"]
            invoke("configure",[sys.executable,home/"meson/meson.py","setup",build,src,*opts],home,receipts,180,2,bin_directory=bin_dir)
        invoke("libraries",[ninja,"-C",build,"-j","2",*LIBS],home,receipts,4200,3,bin_directory=bin_dir)
        record=dict(schema="short65-kernel-v1",upstreamCommit="7048fd35f56d4bf52e586429d6746a36a10255f8",
                    archives={name:v[0] for name,v in DIGESTS.items()},compiler=digest(compiler),ninja=digest(ninja),
                    libraries={name:digest(build/name) for name in LIBS},compilation=receipts.name,license="LGPL-2.1-or-later")
        atomic_write_json(record_file,record)
        return record


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--downloads",type=Path,help="可选，已下载的固定上游原始包")
    a=p.parse_args();print(json.dumps(install(a.downloads),indent=2))
