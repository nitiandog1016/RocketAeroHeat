"""复核低速试件终场、错误观测的拒绝和原面工作台回读。仅读取已有场。"""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/"scripts"))
from readCM011Tunnel import evaluate, read_energy_observation
from runEmbeddedSU2 import sha


def read_page(path):
    text = path.read_text(encoding="utf8").strip()
    prefix = "window.ROCKET_CASE="
    assert text.startswith(prefix) and text.endswith(";")
    return json.loads(text[len(prefix):-1])


def verify(base, prior, output):
    output.mkdir(parents=True, exist_ok=False)
    obs = base/"check01/fine02.json"
    source = prior/"fine02"
    stats, _ = evaluate(source)
    accepted = read_energy_observation(source, stats, obs)
    outer_stats, _ = evaluate(prior/"outer02")
    other = base/"check01/outer02.json"
    read_energy_observation(prior/"outer02", outer_stats, other)
    records = {}

    def expect_rejection(label, path):
        try:
            read_energy_observation(source, stats, path)
        except ValueError as error:
            records[label] = str(error)
        else:
            raise AssertionError("错误观测被采用：" + label)

    expect_rejection("otherGrid", other)
    data = json.loads(obs.read_bytes())
    for name in ("changedRestart", "compensatingTerms", "changedStage", "nonFiniteHeat"):
        bad = copy.deepcopy(data)
        if name == "changedRestart":
            bad["sourceSHA256"]["restart.csv"] = "0"*64
        elif name == "compensatingTerms":
            bad["gapTerms"]["farAdvectionChoiceW"] += 1e-6
            bad["gapTerms"]["farDiffusionChoiceW"] -= 1e-6
        elif name == "changedStage":
            bad["stageIncrements"]["source"]["excessEnergyW"] += 1e-6
        else:
            bad["nativeFarExcessOutW"] = float("nan")
        path = output/(name + ".json")
        path.write_text(json.dumps(bad), encoding="utf8")
        expect_rejection(name, path)

    # 实际走不带观测参数的原生读表入口，逐字节核对两种读法返回的壁面表。
    legacy = output/"no_observation"
    subprocess.run([sys.executable, str(ROOT/"scripts/readCM011Tunnel.py"), str(source), str(legacy)],
                   check=True, timeout=30)
    work = base/"workbench01"
    assert "nativeEnergyBalance" not in json.loads((legacy/"summary.json").read_bytes())
    assert sha(legacy/"wall.csv") == sha(work/"audit/wall.csv")
    summary = json.loads((work/"summary.json").read_bytes())
    assert summary["reusedNativeRun"] and summary["native"]["acceptedNumericalWall"]
    assert summary["native"]["nativeEnergyBalance"] == accepted

    # 同一批原面的几何、物理量、状态、方法、覆盖及热功率均应与上一版完全相同。
    current_page = work/"workbench/case-data.js"
    previous_page = prior/"fine_view01/workbench_final/case-data.js"
    new, old = read_page(current_page), read_page(previous_page)
    assert new.pop("nativeEnergyBalance") == accepted
    heat = next(f for f in new["fields"] if f["key"] == "heatFlux")
    old_note = next(f["stateNote"] for f in old["fields"] if f["key"] == "heatFlux")
    note = heat["stateNote"]
    assert note.startswith(old_note)
    for key in ("nativeWallReleaseW", "nativeFarExcessOutW", "nativeSurfaceImbalanceW"):
        assert f'{1e6*accepted[key]:.3f} µW' in note
    assert f'{100*accepted["nativeSurfaceImbalanceW"]/accepted["nativeWallReleaseW"]:.2f}%' in note
    heat["stateNote"] = old_note
    assert new == old
    coverage = new["wallCoverage"]
    assert coverage["fullSurfaceComputed"] and coverage["jointComputedFaceCount"] == len(new["faces"]) == 14140
    return dict(fieldKeys=[f["key"] for f in new["fields"]], originalFaces=len(new["faces"]),
                areaM2=coverage["jointComputedAreaM2"], powerIntoWallW=coverage["wholeSurfaceHeatPowerW"],
                priorPageSHA256=sha(previous_page), pageSHA256=sha(current_page), heatFluxStateNote=note,
                nativeEnergyBalance=accepted, fieldsAndCoverageMatch=True, invalidObservationsRejected=records,
                legacyWallCsvSHA256=sha(legacy/"wall.csv"), legacyObservationAbsent=True, browserTested=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("new_output", type=Path)
    args = parser.parse_args()
    result = verify(ROOT/"outputs/v272_validation", ROOT/"outputs/v271_validation", args.new_output.resolve())
    (args.new_output/"report.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf8")
