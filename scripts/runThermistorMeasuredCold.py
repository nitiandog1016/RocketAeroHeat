"""统一工作台的热敏珠冷启动任务。按同一任务的正常阶段接续，支持在保存后暂停。"""
import argparse
from contextlib import contextmanager
import msvcrt
from pathlib import Path
import shutil
import uuid

from thermistorMeasuredCold import GRID_KEY, MESH_FILES, STEPS, advance, audit, frame, grid, seed
from thermistorBead import ROOT, read, sha, write_json
from runEmbeddedSU2 import atomic_write_json
import numpy as np


@contextmanager
def exclusive(run):
    """由操作系统在进程退出时释放锁，同目录同时只能有一个推进者。"""
    with (run / "job.lock").open("a+b") as stream:
        stream.seek(0)
        try:
            msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
        except OSError as error:
            raise RuntimeError("该目录仍有本任务运行，请查看原任务进度。") from error
        try:
            yield
        finally:
            stream.seek(0)
            msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)


def check_seed(run, spec):
    folder = run / "cold"
    rec = read(folder / "cold.json")
    key = read(GRID_KEY)
    for name, field in MESH_FILES.items():
        assert sha(folder / name) == key[field]
    assert read(folder / "case.json") == spec and sha(folder / "case.json") == rec["caseSHA256"]
    with np.load(folder / "nodes.npz") as mapping:
        actual = frame(folder / "initial_flow.csv", mapping["xyz"], np.unique(mapping["wallMap"]), cold=True)
    assert rec["initial"] == actual and rec["initialSource"] == "geometry-and-boundary"
    return rec


class SavedPause(Exception):
    """当前操作已保存，可以从同一任务目录继续。"""


def execute(action, run, prepared, stage_limit=None, invocation=None, parent_pause=None):
    spec = grid(prepared)
    if action == "cold":
        run.mkdir(parents=True, exist_ok=False)
    elif not (run / "job.json").is_file():
        raise ValueError("所选目录没有本工况的冷启动任务记录。")
    with exclusive(run):
        identity = dict(kind="thermistor-b-measured-12K-cold", directory=str(run), geometry=sha(GRID_KEY),
                        pressurePa=spec["airPa"], velocityMps=spec["velocityMps"], wallK=spec["wallK"], airK=spec["airK"])
        if action == "cold":
            if shutil.disk_usage(run).free < 16 * 1024**3:
                raise OSError("冷启动需要 16 GiB 可用磁盘空间。")
            seed(prepared, run / "cold", spec)
            write_json(run / "job.json", dict(identity, token=str(uuid.uuid4()), seedSHA256=sha(run / "cold/cold.json")))
        job = read(run / "job.json")
        assert all(job[k] == v for k, v in identity.items()) and job["seedSHA256"] == sha(run / "cold/cold.json")
        check_seed(run, spec)
        token = str(uuid.UUID(invocation)) if invocation else str(uuid.uuid4())
        pause = run.with_name(run.name + "." + token + ".pause")
        def report(stage, status):
            progress = dict(stage=stage, status=status, invocation=token, job=job["token"])
            atomic_write_json(run / "progress.json", progress)
            if status in ("complete", "paused"):
                write_json(run / ("receipt-" + token + ".json"), progress)
        def checkpoint(stage):
            if pause.exists() or (parent_pause is not None and parent_pause.exists()):
                report(stage, "paused")
                raise SavedPause
            report(stage, "running")
        initial = run / "cold/initial_flow.csv"
        previous = 0
        maxima = []
        for number in range(1, 9):
            dest = run / f"steady{number:02d}"
            checkpoint(dest.name)
            if not dest.exists():
                advance(run, dest, initial, STEPS)
            result = audit(run, dest, initial, previous, STEPS)
            previous, initial = result["cumulativeSteps"], dest / "restart.csv"
            maxima.append(max(result["last20Residual"].values()))
            print(f"{dest.name}  steps={previous}  max(log10 RMS)={maxima[-1]:.6f}", flush=True)
            if result["allResidualsBelowMinus8"]:
                if stage_limit == number:
                    report(dest.name, "paused")
                    return
                from thermistorMeasuredCollect import finish
                finish(run, dest, job, checkpoint)
                report(dest.name, "complete")
                return
            if len(maxima) >= 3 and maxima[-1] >= maxima[-2] >= maxima[-3]:
                raise RuntimeError("连续两个阶段末窗未改善，已保留各段终态。")
            if stage_limit == number:
                report(dest.name, "paused")
                return
        raise RuntimeError("已达登记的总步数上限，六分量末窗尚未全部达标。")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("cold", "resume"))
    parser.add_argument("directory", type=Path)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--through-stage", type=int, choices=range(1, 9))
    parser.add_argument("--invocation", help="本次侧栏操作的 UUID，用于隔离暂停信号")
    args = parser.parse_args()
    run, prepared = args.directory.resolve(), args.prepared.resolve()
    for path in (run, prepared):
        path.relative_to(ROOT / "outputs")
    try:
        execute(args.action, run, prepared, args.through_stage, args.invocation)
    except SavedPause:
        print("当前阶段已保存，可从同一目录接续。")
