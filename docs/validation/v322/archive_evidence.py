"""归档银球试验的小型原始回执，绘制实际时间层与测热时刻的比较图。"""
from pathlib import Path
import json
from zipfile import ZIP_DEFLATED, ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
WORK = ROOT / "outputs/v322_validation/run01"
RESULT = WORK / "result_064"
OUT = HERE / "evidence"


def read(path):
    return json.loads(path.read_text("utf-8-sig"))


def archive():
    costs = {}
    with ZipFile(OUT / "receipts.zip", "x", ZIP_DEFLATED) as z:
        def save(path):
            z.write(path, path.relative_to(ROOT / "outputs").as_posix())

        def count(kind, seconds, memory):
            v = costs.setdefault(kind, dict(count=0, seconds=0, maxRSSBytes=0))
            v["count"] += 1
            v["seconds"] += seconds
            v["maxRSSBytes"] = max(v["maxRSSBytes"], memory)

        for path in sorted(WORK.rglob("debugger_process.json")):
            d = path.parent
            r = read(path)
            assert r["status"] == "terminal" and r["debuggerExitCode"] == 0
            assert all(r["inputsStable"].values()) and not (r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"])
            assert "Exit Success (SU2_CFD)" in (d / "debugger.log").read_text("ascii", errors="replace")
            kind = d.name if d.parent.name == "readout" else d.parent.name
            if (d / "prefix_rejected.json").exists():
                kind += "_rejected"
            count(kind, r["elapsedSeconds"], r["peakFamilyRSSBytes"])
            save(path)
            for file in ["case.cfg", "native.json", "stages.csv", "readout_updates.csv", "prefix_rejected.json"]:
                if (d / file).is_file():
                    save(d / file)
            save(next(d.glob("history_*.csv")))
        for d in [WORK / "observer", WORK / "observer4"]:
            for path in sorted(d.glob("*.json")):
                save(path)
                if path.name.endswith("_process.json"):
                    r = read(path)
                    assert r["exitCode"] == 0 and not (r["stopReason"] or r["livePids"])
                    count("build", r["elapsedSeconds"], r["peakSampledRSSBytes"])
        for path in sorted(WORK.rglob("accepted.json")):
            save(path)
        for path in [WORK/"job.json", WORK/"result_064.json", RESULT/"comparison.json"]:
            save(path)
        for name in ["compute01", "read01"]:
            d = WORK.parent / name
            save(d/"ui.json")
            save(d/"debugger_process.json")
            r = read(d/"debugger_process.json")
            count(name, r["elapsedSeconds"], r["peakFamilyRSSBytes"])
    with (OUT/"resources.json").open("x", encoding="utf8") as f:
        json.dump(costs, f, indent=2)
        f.write("\n")


def plot_response():
    match = read(RESULT / "comparison.json")
    check = read(OUT / "terminal.json")
    base = read(RESULT / "fine/summary.json")["baselineNu"]
    plt.rcParams.update({"font.family": "Microsoft YaHei", "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
    for name, style, color in [("coarse", "--", "#b46024"), ("fine", "-", "#254f75")]:
        trace = match["trace"][name]
        axes[0].plot([r["timeDU"] for r in trace], [r["nonAxisymmetricRmsOverU"]*1000 for r in trace],
                     style, color=color, label="Δt U/D = "+("0.25" if name=="coarse" else "0.125"))
    axes[0].axhline(check["wakeFine"]["steadyWakeRmsOverU"]*1000, color="#72767a", ls=":", label="扰动前基流")
    axes[0].set_ylabel(r"尾迹非轴对称速度 RMS / U  ($\times 10^{-3}$)")
    stations = match["samples"]["fine"]
    for station, label, marker, color in [(5, "119.8°", "o", "#254f75"), (6, "166.7°", "s", "#b46024")]:
        axes[1].plot([r["timeDU"] for r in stations], [1000*(r["localNu"][station]-base[station]) for r in stations],
                     color=color, marker=marker, ms=4, lw=1, label=label)
    axes[1].axhline(0, lw=.7, color="#72767a")
    axes[1].set_ylabel(r"瞬时 Nu − 基流 Nu  ($\times 10^{-3}$)")
    axes[0].set_title("尾迹速度  ·  每一物理步")
    axes[1].set_title("背风换热  ·  细步的六个采样时刻")
    for ax in axes:
        ax.set_xlabel("距扰动初场的时间  t U/D")
        ax.set_xticks(np.arange(0, 9, 2))
        ax.grid(axis="y", alpha=.18)
        ax.legend(frameon=False)
    fig.savefig(OUT / "wake.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    archive()
    plot_response()
