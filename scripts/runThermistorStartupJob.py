"""工作台调用的 B 型热敏珠 6 m/s 全周计算。复用两半步的原生推进和独立审核。"""
import argparse
from contextlib import contextmanager
import gzip
import msvcrt
import os
from pathlib import Path
import shutil
import subprocess
import sys

for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import psutil

ROOT = Path(__file__).resolve().parents[1]
RECIPE = ROOT / "docs/validation/v294"
sys.path.insert(0, str(RECIPE))
from advance import BUILD, ENGINE_SHA, GRID, OFFICIAL, OLD, PREVIOUS, RUN, grid_info, read, settings, sha
from prepare import bound_source
from compare import BASE, BASE_SHA
from runEmbeddedSU2 import atomic_write_json

JOBS = RUN / "workbench"
EVIDENCE = RECIPE / "evidence"
CONDITION = dict(speedMps=6, airK=293.15, airPa=101325, wallK=303.15,
                 elapsedPhysicalTimeS=7.80936086544745e-6)
GRID_SPEC = read(PREVIOUS / "evidence/grid.json")
GEOMETRY = "5e66bb2f70f0774b8ebb87a77b76061fc84b89b9f379d54e90879e81628bf62f"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def prepare(out):
    """从实际全周体网取壁面预览并绑定米制几何、初场和节点编号。"""
    info = grid_info()
    geometry = read(bound_source("geometry.json"))
    with np.load(GRID / "nodes.npz") as data:
        wall = data["ids"][:, :65]
        quads = np.stack([wall[:-1, :-1], wall[1:, :-1], wall[1:, 1:], wall[:-1, 1:]], axis=-1).reshape(-1, 4)
        triangles = np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
        triangles = triangles[np.all(np.diff(np.sort(triangles, axis=1), axis=1) != 0, axis=1)]
        nodes, indices = np.unique(triangles, return_inverse=True)
        vertices = data["xyz"][nodes] * geometry["transverseDiameterM"]
    require(info == GRID_SPEC and len(nodes) == 3026 and len(triangles) == 6048, "全周壁面预览与体网不符")
    atomic_write_json(out, dict(directory=str(GRID), meshSHA256=info["gridSHA256"],
        bodySHA256=geometry["originalBodySHA256"], seedSHA256=info["seedSHA256"], nodeMapSHA256=info["nodeMapSHA256"],
        angularCells=48, radialCells=56, cells=info["cells"], nodes=info["nodes"],
        wallFaces=len(quads), preview=dict(V=vertices.tolist(), F=(indices.reshape(-1, 3)+1).tolist())))


def preflight():
    """在创建新计算目录前核对该具名配方真正读取的本机资产。"""
    grid_info()
    build = read(PREVIOUS / "evidence/readout_build.json")
    require(sha(OFFICIAL) == ENGINE_SHA and sha(BUILD / "observe_bead.exe") == build["observerSHA256"]
            and read(BUILD / "build.json") == build and sha(BUILD / "CIntegration.cpp") == build["sourceSHA256"],
            "SU2 或固定端点观察器指纹有变化")
    require(sha(BASE) == BASE_SHA, "48 层参照场指纹有变化")
    for name in ("body.json", "nodes.npz", "geometry.json"):
        bound_source(name)
    assets = read(PREVIOUS / "evidence/local_assets.json")
    for name in ("case.cfg", "history_00001.csv"):
        path = OLD / "first01" / name
        require(sha(path) == assets[path.relative_to(ROOT).as_posix()]["sha256"], "v293 首步参照输入已变化")
    settings(1)
    settings(2)
    require(shutil.disk_usage(RUN).free >= 3 * 1024**3, "需要至少 3 GiB 可用磁盘空间")


def stages():
    for i, label in enumerate(("first", "second"), 1):
        yield label + "01", ["advance.py", str(i)]
        for stage in ("state", "observe", "release", "energy"):
            yield label + "_" + stage, ["verify.py", stage, str(i)]
    yield "compare", ["compare.py"]


def terminal(area, name):
    path = area / name / "debugger_process.json"
    require(path.is_file(), name + " 缺少终态回执。计算目录 " + str(area))
    receipt = read(path)
    require(receipt.get("status") == "terminal" and receipt.get("debuggerExitCode") == 0
            and not receipt["stopReason"] and not receipt["supervisorError"] and not receipt["survivingProcessIDs"]
            and receipt["inputsStable"] and all(receipt["inputsStable"].values()),
            "%s 尚无成功终态回执（状态 %s，停止原因 %s），现场 %s" %
            (name, receipt.get("status"), receipt.get("stopReason") or receipt.get("supervisorError") or "见原生日志", path.parent))
    return receipt


def check_step(area, label, i, old_sha):
    """回读实际旧新场和逐步审核链，防止借用其他计算目录的结果。"""
    src = area / (label + "01")
    native = terminal(area, src.name)
    require(native["command"] == [str(OFFICIAL), "case.cfg", "--threads", "1"]
            and native["inputSHA256"].get(str(OFFICIAL)) == ENGINE_SHA, "原生程序或线程设置不符")
    def bound(receipt, path):
        digest = sha(path)
        require(receipt["inputSHA256"].get(str(path)) == digest, "阶段输入或目录不符：" + str(path))
        return digest
    state_rec = terminal(area, label + "_state")
    state_path = area / (label + "_state/verified.json")
    state = read(state_path)
    frame = src / ("restart_%05d.csv" % i)
    history = src / ("history_%05d.csv" % i)
    config_sha = read(EVIDENCE / "local_assets.json")["outputs/v294_validation/" + src.name + "/case.cfg"]["sha256"]
    require(bound(native, src / "case.cfg") == config_sha and bound(native, src / "mesh.su2") == GRID_SPEC["gridSHA256"]
            and bound(native, src / ("seed_%05d.csv" % (i-1))) == old_sha, "正式推进的配置、体网或旧帧不符")
    bound(state_rec, src / "debugger_process.json")
    require(state["passed"] and state["normalNativeTermination"] and state["step"] == i
            and state["gridSHA256"] == GRID_SPEC["gridSHA256"] and state["oldFrameSHA256"] == old_sha
            and state["frameSHA256"] == bound(state_rec, frame) and state["historySHA256"] == bound(state_rec, history)
            and state["trailingPassingRows"] >= 20 and max(state["last20MaxLogRMS"].values()) < -8
            and state["nodeRHSMaxLogError"] < 1e-7 and state["elapsedPhysicalTimeS"] == i*CONDITION["elapsedPhysicalTimeS"]/2,
            "保存帧与六分量连续窗、物理时间或节点复算不符")
    if i == 1:
        prefix = gzip.decompress((PREVIOUS / "evidence/first_history.csv.gz").read_bytes())
        require(state["verifiedPrefixRows"] == 560 and history.read_bytes().startswith(prefix), "首半步前缀不符")
    energy_rec = terminal(area, label + "_energy")
    bound(energy_rec, state_path)
    for stage in ("observe", "release"):
        folder = area / (label + "_" + stage)
        terminal(area, folder.name)
        bound(energy_rec, folder / "debugger_process.json")
    result = read(area / (label + "_energy/accepted.json"))
    e = result["energy"]
    require(result["accepted"] and result["step"] == i and Path(result["acceptedFrame"]) == frame
            and result["elapsedPhysicalTimeS"] == state["elapsedPhysicalTimeS"] and e["oldFrameSHA256"] == old_sha
            and e["acceptedFrameSHA256"] == state["frameSHA256"] and e["dtSeconds"] == state["dtSeconds"]
            and result["gridSHA256"] == state["gridSHA256"] and abs(e["relativeImbalance"]) < .005
            and abs(e["airStorageW"]-e["independentStorageW"]) < 1e-9 and e["releasePeakRelativeError"] < 1e-5
            and e["relativeWallChange"] < .002, "壁反力、储热或接受帧不符")
    return result


def inspect(area):
    first = check_step(area, "first", 1, GRID_SPEC["seedSHA256"])
    second = check_step(area, "second", 2, first["energy"]["acceptedFrameSHA256"])
    comparison = terminal(area, "compare")
    summary = read(area / "compare/summary.json")
    require(summary["firstHalfEnergy"] == first["energy"] and summary["fineSurface"] == second["surface"]
            and summary["gridTransfer"] == GRID_SPEC and summary["originalFaces"] == 5216
            and summary["bodySHA256"] == GEOMETRY and summary["elapsedPhysicalTimeS"] == CONDITION["elapsedPhysicalTimeS"],
            "近壁对照和实际端点不符")
    files = [area / "first_energy/accepted.json", area / "second_energy/accepted.json",
             area / "second_energy/wall.csv", area / "second_energy/matched_temperature.csv"]
    require(all(comparison["inputSHA256"].get(str(p)) == sha(p) for p in files)
            and sha(files[2]) == second["surface"]["wallTableSHA256"]
            and sha(files[3]) == second["matchedTableSHA256"] == summary["matchedTableSHA256"], "同刻原面表绑定失败")
    return dict(directory=str(area), condition=CONDITION, frameSHA256=second["energy"]["acceptedFrameSHA256"],
                summarySHA256=sha(area / "compare/summary.json"), wallSHA256=sha(files[2]))


@contextmanager
def exclusive(area):
    with (area / "workbench.lock").open("a+b") as lock:
        lock.seek(0)
        try:
            msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
        except OSError as exc:
            raise ValueError("该计算目录已有工作台任务运行") from exc
        try:
            yield
        finally:
            lock.seek(0)
            msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)


def preserve_interrupted_release(area, stage):
    """接续时保留已停止、尚无保存帧的单次独立壁反力读数，原终态仍记未知。"""
    folder = area / stage
    path = folder / "debugger_process.json"
    if not stage.endswith("_release") or not path.is_file():
        return
    record = read(path)
    if record.get("status") != "running":
        return
    require(folder.resolve() == folder and record["command"] == [str(OFFICIAL), "case.cfg", "--threads", "1"]
            and record["inputSHA256"].get(str(OFFICIAL)) == ENGINE_SHA
            and not any(folder.glob("restart_*.csv")), "中断读数的目录、命令或终场需要人工核对")
    for key in ("processID", "supervisorProcessID"):
        try:
            proc = psutil.Process(record[key])
            require(proc.create_time() > path.stat().st_mtime + 2, "原读数进程仍可能活动，保留现场")
        except psutil.NoSuchProcess:
            pass
    require(all(sha(Path(k)) == v for k, v in record["inputSHA256"].items()), "中断读数的输入指纹发生变化")
    digest = sha(path)
    archive = folder.with_name(stage + ".interrupted-" + digest[:16])
    require(archive.parent.resolve() == area and not archive.exists(), "中断现场归档目录已有内容")
    folder.rename(archive)


def run(area, resume):
    require(area.parent == JOBS, "新任务须位于本工况专用输出目录")
    preflight()
    request = dict(directory=str(area), condition=CONDITION, meshSHA256=GRID_SPEC["gridSHA256"],
                   seedSHA256=GRID_SPEC["seedSHA256"])
    if resume:
        require(read(area / "request.json") == request, "所选目录与本工况不符")
    else:
        area.mkdir(parents=True, exist_ok=False)
        atomic_write_json(area / "request.json", request)
    with exclusive(area):
        progress = area / "progress.json"
        stage = "preflight"
        try:
            if resume:
                existing_names = [name for name, _ in stages() if (area / name).exists()]
                existing = [name in existing_names for name, _ in stages()]
                require(existing == sorted(existing, reverse=True), "保存阶段中间有缺项。已保留目录供核查")
            for stage, args in stages():
                directory = area / stage
                if resume and directory.exists():
                    if stage == existing_names[-1]:
                        preserve_interrupted_release(area, stage)
                if resume and directory.exists():
                    old = terminal(area, stage)
                    require(all(sha(Path(k)) == v for k, v in old["inputSHA256"].items()),
                            "接续时已保存阶段的输入已变化：" + stage)
                    if stage.endswith("_state"):
                        require(read(directory / "verified.json")["passed"], "保存的状态审核未通过")
                    if stage.endswith("_energy"):
                        first = stage == "first_energy"
                        previous = GRID_SPEC["seedSHA256"] if first else read(area/"first_energy/accepted.json")["energy"]["acceptedFrameSHA256"]
                        check_step(area, "first" if first else "second", 1 if first else 2, previous)
                    continue
                atomic_write_json(progress, dict(status="running", stage=stage, pid=os.getpid()))
                process = subprocess.run([sys.executable, str(RECIPE / args[0]), *args[1:], "--root", str(area)])
                if process.returncode:
                    log = directory / "debugger.log"
                    raise RuntimeError("%s 审核或推进失败（退出 %d）。现场 %s" % (stage, process.returncode, log))
            result = inspect(area)
            atomic_write_json(progress, dict(status="accepted", stage="compare", result=result))
            return result
        except Exception as exc:
            atomic_write_json(progress, dict(status="incomplete", stage=stage, reason=str(exc)))
            raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "resume", "read"))
    parser.add_argument("path", type=Path)
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    path = args.path.resolve()
    require(__debug__, "物理审核需要启用 Python 断言")
    path.relative_to(ROOT / "outputs")
    if args.action == "prepare":
        prepare(path)
    elif args.action == "read":
        args.receipt.resolve().relative_to(ROOT / "outputs")
        atomic_write_json(args.receipt, inspect(path))
    else:
        result = run(path, args.action == "resume")
        if args.receipt:
            args.receipt.resolve().relative_to(ROOT / "outputs")
            atomic_write_json(args.receipt, result)
