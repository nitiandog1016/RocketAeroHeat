"""接续验收：在同一任务的官方一步回放期间请求暂停，保存各已核验子阶段的身份。"""
import argparse
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from thermistorBead import read, sha, write_json


def test(run, invocation, result, after):
    job = read(run / "job.json")["token"]
    flag = run.with_name(run.name+"."+invocation+".pause")
    started = time.monotonic()
    triggered = False
    while time.monotonic() - started < 900:
        p = read(run / "progress.json")
        assert p["job"] == job
        if p["invocation"] != invocation:
            time.sleep(.5)
            continue
        if p["stage"] == after and not triggered:
            # 入口先核对调用身份，再观察官方程序创建的原生日志。
            logs = list(run.glob(f"steady*/ledger/{after}/step.log"))
            if logs:
                assert len(logs) == 1 and not (logs[0].parent / "step_process.json").exists()
                flag.touch(exist_ok=False)
                triggered = True
        receipt = run / ("receipt-"+invocation+".json")
        if triggered and receipt.is_file():
            pause = read(receipt)
            expected = "observed" if after == "control" else "surface"
            assert pause["status"] == "paused" and pause["stage"] == expected
            assert not (run / "completed.json").exists()
            missing = "ledger/observed" if after == "control" else "surface"
            assert not list(run.glob("steady*/"+missing))
            folders = list(run.glob(f"steady*/ledger/{after}/checked.json"))
            assert len(folders) == 1
            stage = folders[0].parents[2]
            names = ["engine_process.json","reaction0/engine_process.json","plus20/engine_process.json",
                     "reaction20/engine_process.json","ledger/control/checked.json","ledger/control/step_process.json"]
            if after == "observed":names += ["ledger/observed/checked.json","ledger/observed/step_process.json","ledger/energy.json"]
            write_json(result, dict(pause=pause, stage=stage.name, recordedBeforeResume={n:sha(stage/n) for n in names}))
            return
        if p["status"] == "complete":raise RuntimeError("未在完成前触发暂停")
        time.sleep(.5)
    raise TimeoutError("未在登记观察时间内得到本次暂停回执")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path);p.add_argument("invocation");p.add_argument("result",type=Path)
    p.add_argument("--after",choices=("control","observed"),default="control")
    args = p.parse_args()
    test(args.run.resolve(),args.invocation,args.result.resolve(),args.after)
