"""对实际热敏珠工作台任务检查并发拒绝、阶段缺失、终态恢复和原生数据身份。"""
import argparse
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import uuid
import zipfile

import psutil
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from runThermistorStartupJob import ENGINE_SHA, JOBS, OFFICIAL, RECIPE, inspect, preserve_interrupted_release, read, sha, stages


def invoke(action, area, proof):
    cmd = [sys.executable, str(ROOT / "scripts/runThermistorStartupJob.py"), action,
           str(area), "--receipt", str(proof)]
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def release_guards(out):
    """只在测试目录构造中断读数，核对活动进程、已有终场和变化输入都能保住现场。"""
    child = subprocess.Popen([sys.executable, "-c", "pass"])
    assert child.wait(timeout=15) == 0 and not psutil.pid_exists(child.pid)
    rejected = {}
    for case, message in (("live", "进程仍可能活动"), ("saved", "终场需要人工核对"), ("changed", "输入指纹发生变化")):
        area = out / case
        folder = area / "first_release"
        folder.mkdir(parents=True)
        seed = folder / "seed.csv"
        seed.write_bytes(b"first")
        owner = os.getpid() if case == "live" else child.pid
        receipt = dict(status="running", command=[str(OFFICIAL),"case.cfg","--threads","1"],
                       processID=owner, supervisorProcessID=owner,
                       inputSHA256={str(OFFICIAL): ENGINE_SHA, str(seed): sha(seed)})
        path = folder / "debugger_process.json"
        path.write_text(json.dumps(receipt), encoding="utf8")
        if case == "saved":
            (folder / "restart_00001.csv").touch()
        elif case == "changed":
            seed.write_bytes(b"changed")
        before = {p.name: sha(p) for p in folder.iterdir()}
        try:
            preserve_interrupted_release(area, "first_release")
        except ValueError as exc:
            assert message in str(exc), str(exc)
        else:
            raise AssertionError("应当保留现场：" + case)
        assert before == {p.name: sha(p) for p in folder.iterdir()} and len(list(area.iterdir())) == 1
        rejected[case] = True
    return rejected


def live_checks(area, out):
    progress = read(area / "progress.json")
    worker = psutil.Process(progress["pid"])
    assert any(Path(arg).name == "runThermistorStartupJob.py" for arg in worker.cmdline())
    first = read(area / "first01/debugger_process.json")
    cfd = psutil.Process(first["processID"])
    assert Path(cfd.cwd()) == area / "first01" and cfd in worker.children(recursive=True)
    old_request = sha(area / "request.json")
    result = {}
    for action in ("read", "resume"):
        target = out / (action + ".json")
        trial = invoke(action, area, target)
        (out / (action + ".log")).write_bytes(trial.stdout)
        assert trial.returncode and not target.exists()
        result[action + "Rejected"] = True
    assert cfd.is_running() and worker.is_running() and read(area / "progress.json") == progress
    assert sha(area / "request.json") == old_request
    result.update(activeChild=cfd.pid, activeSupervisor=worker.pid, activeStage=progress["stage"],
                  duplicateNativeStarts=False, directory=str(area), stageGapRejected=False)
    # 用独立临时任务检查缺项拒绝，保留正在运行的真实目录。
    gap = JOBS / ("recovery-gap-" + str(uuid.uuid4()))
    gap.mkdir()
    request = read(area / "request.json")
    request["directory"] = str(gap)
    (gap / "request.json").write_text(json.dumps(request), encoding="utf8")
    (gap / "first_state").mkdir()
    trial = invoke("resume", gap, out / "gap_result.json")
    (out / "gap.log").write_bytes(trial.stdout)
    assert trial.returncode and not (gap / "first01").exists()
    result["stageGapRejected"] = True
    return result


def recovered(area, out, interruption=None):
    before = {name: sha(area / name / "debugger_process.json") for name, _ in stages()}
    actual = inspect(area)
    receipt = out / "resumed_result.json"
    trial = invoke("resume", area, receipt)
    (out / "resume.log").write_bytes(trial.stdout)
    assert trial.returncode == 0 and read(receipt) == actual
    assert before == {name: sha(area / name / "debugger_process.json") for name, _ in stages()}
    reference = read(RECIPE / "evidence/comparison.json")["fineSurface"]
    current = read(area / "second_energy/accepted.json")["surface"]
    assert reference["bodySHA256"] == current["bodySHA256"]
    package = ROOT / "cases/thermistor_b/startup_6ms_nearwall56.zip"
    with zipfile.ZipFile(package) as archive:
        baseline = np.genfromtxt(io.BytesIO(archive.read("wall.csv")), names=True, delimiter=",")
    observed = np.genfromtxt(str(area / "second_energy/wall.csv"), names=True, delimiter=",")
    assert np.array_equal(baseline["faceID"], observed["faceID"]) and np.array_equal(baseline["areaM2"], observed["areaM2"])
    max_diff = {key: float(np.max(abs(observed[key]-baseline[key])))
                for key in ("pressurePa", "heatFluxWm2", "nearWallCellTemperatureK")}
    result = dict(result=actual, stagesPreserved=before, originalFaces=current["originalFaces"],
                maximumAbsoluteDifference=max_diff, referencePackageSHA256=sha(package),
                sameWallBytes=reference["wallTableSHA256"] == current["wallTableSHA256"],
                relativePowerChange=current["energy"]["wallReleaseW"]/reference["energy"]["wallReleaseW"] - 1)
    if interruption is not None:
        result["releaseGuardsPreservedOriginals"] = release_guards(out)
        prior = read(interruption)
        assert Path(prior["directory"]) == area
        digest = prior["originalReleaseReceiptSHA256"]
        archive = area / ("first_release.interrupted-" + digest[:16])
        assert sha(archive / "debugger_process.json") == digest != before["first_release"]
        assert sorted(p.name for p in archive.iterdir()) == prior["originalReleaseFileNames"]
        assert all(before[k] == v for k, v in prior["savedFirstStages"].items())
        result["reusedStagesAfterInterruption"] = prior["savedFirstStages"]
        result["interruptedReleaseReceiptSHA256"] = digest
    return result


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("mode", choices=("live", "recovered"))
    cli.add_argument("area", type=Path)
    cli.add_argument("out", type=Path)
    cli.add_argument("--interruption", type=Path)
    a = cli.parse_args()
    area, out = a.area.resolve(), a.out.resolve()
    assert area.parent == JOBS
    out.relative_to(ROOT / "outputs/v295_validation")
    out.mkdir(exist_ok=False)
    result = live_checks(area, out) if a.mode == "live" else recovered(area, out, a.interruption)
    (out / "checks.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf8")
