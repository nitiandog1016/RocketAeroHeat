"""摘录独立安装的真实构建、计算、热量账与窗口回执；大型保存帧只登记位置和指纹。"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
from zipfile import ZipFile, ZIP_DEFLATED

HERE = Path(__file__).resolve().parent


def read(file):
    return json.loads(file.read_bytes())


def pack(app, job, windows, dest):
    origin = read(job / "construction.json")
    cache = app / "outputs/short65_inputs/14D"
    native = app / origin["programs"]
    steady = app / origin["steady"]
    kernel = app / "tools/runtime/short65-14d"
    paths = [app.parent / (app.name+"_installed.json"), app.parent / (app.name+"_source.json"),
             app.parent / "public_downloads/receipt.json", app / "tools/runtime/gmsh/4.15.2/runtime.json",
             kernel / "kernel.json", app / origin["prepared"] / "prepared.json"]
    paths.extend(cache.glob("*.json"))
    post_sources=app.parent/"reopen_source.json"
    if post_sources.is_file():paths.append(post_sources)
    paths.extend((app.parent/"audits").glob("*_process.json"))
    # 归档本次验收的早期失败及更正。普通用户第一次成功时没有这些历史文件。
    initial_window=app.parent/"windows/solve01/window_process.json"
    if initial_window.is_file() and read(initial_window)["exitCode"]!=0:
        paths.extend([initial_window,app.parent/"solve02_source.json"])
        for failure in (app/"outputs").glob("short65_time_*.log"):
            if b"AssertionError" in failure.read_bytes():
                paths.append(failure)
        for previous in cache.glob("base_*/steady"):
            if previous!=steady and not (previous/"thermal.json").is_file():
                paths.extend(previous.glob("*/record.json"))
                paths.extend(previous.glob("*/debugger_process.json"))
                paths.extend(previous.glob("converged.json"))
    paths.extend(kernel.glob("build_*/*_process.json"))
    paths.extend(native.glob("*/*.json"))
    paths.extend(steady.glob("*.json"))
    paths.extend(job.glob("*.json"))
    paths.extend([job / "export/pair.json", job / "run/job.json"])
    for branch in ("full","half"):
        paths.append(job / "export" / branch / "summary.json")
    for window in windows:
        paths.append(window / "ui.json")
        path=window/"window_process.json"
        paths.append(path if path.exists() else window/"supervisor_failure.json")
    frames = {}
    for folder in (steady, job / "run"):
        for proc in folder.rglob("debugger_process.json"):
            info = read(proc)
            exe = Path(info["command"][0])
            exe.relative_to(native)
            for name in info["inputSHA256"]:
                Path(name).relative_to(app)
            assert info["status"] == "terminal" and info["debuggerExitCode"] == 0 and all(info["inputsStable"].values())
            assert not (info["stopReason"] or info["supervisorError"] or info["survivingProcessIDs"])
            paths.append(proc)
            for name in ("record.json","case.cfg","native.json","energy.json"):
                file = proc.parent / name
                if file.is_file():paths.append(file)
            if (proc.parent/"record.json").exists():
                record = read(proc.parent/"record.json")
                for name, digest in record.get("frames", {"restart.csv": record["end"]}).items():
                    frames[(proc.parent/name).relative_to(app).as_posix()] = digest
            if (proc.parent/"native.json").exists():
                native_step = read(proc.parent/"native.json")
                name = f"restart_{native_step['step']:05d}.csv"
                frames[(proc.parent/name).relative_to(app).as_posix()] = native_step["frame"]
    for receipt in native.glob("*/*_process.json"):
        data=read(receipt)
        assert not (data["exitCode"] or data["stopReason"] or data["livePids"])
        assert Path(data["cwd"]) == kernel / "build"
        assert Path(data["args"][0]) == kernel / "toolchain/w64devkit/bin/g++.exe"
        for arg in data["args"]:
            path=Path(arg)
            if path.is_absolute():path.relative_to(app)
    imports={}
    for exe in sorted(native.glob("*/*.exe")):
        output=subprocess.check_output([str(kernel/"toolchain/w64devkit/bin/objdump.exe"),"-p",str(exe)])
        names=re.findall(r"DLL Name:\s*(\S+)",output.decode("ascii",errors="replace"))
        assert {name.lower() for name in names}=={"kernel32.dll","msvcrt.dll"}
        imports[exe.name]=names
    assert set(imports)=={"shared.exe","wire.exe","physical.exe","release.exe"}
    with ZipFile(dest, "x", compression=ZIP_DEFLATED) as z:
        digests = {}
        for file in sorted(set(paths)):
            key=file.relative_to(app.parent).as_posix()
            data=file.read_bytes()
            z.writestr(key,data)
            digests[key]=hashlib.sha256(data).hexdigest()
        z.writestr("manifest.json", json.dumps(dict(files=digests, largeLocalFrames=frames,windowsDLLs=imports),indent=2))


def check(archive):
    with ZipFile(archive) as z:
        manifest=json.loads(z.read("manifest.json"))
        for name, expected in manifest["files"].items():
            assert hashlib.sha256(z.read(name)).hexdigest()==expected
        for name in (n for n in manifest["files"] if n.endswith("/energy.json")):
            e=json.loads(z.read(name))
            error=e["wallW"]-e["farW"]-e["storageW"]
            assert abs(error/e["wallW"])<.005 and abs(error-e["wallW"]*e["relativeBalance"])<1e-10
        assert len([n for n in manifest["files"] if n.endswith("/energy.json")])==3


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive",type=Path)
    parser.add_argument("--app",type=Path)
    parser.add_argument("--job",type=Path)
    parser.add_argument("--windows",nargs=2,type=Path)
    args=parser.parse_args()
    if args.app:
        pack(args.app.resolve(),args.job.resolve(),[p.resolve() for p in args.windows],args.archive)
    check(args.archive)
    print("已核对构建及原生回执指纹和三个端点能量账。")
