"""模块用途。只读复核小型保存场的有量纲尺度与独立热状态，保留原逐列失败。"""
import json
from pathlib import Path

import numpy as np

from auditSU2IdealAirState import audit, load_state
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save


def main():
    folder = ROOT/'outputs/v164_validation/small01'
    output = ROOT/'outputs/v164_validation/small_report.json'
    assert not output.exists()
    manifest = read(folder/'case_manifest.json')
    receipt = read(folder/'engine_process.json')
    assert receipt['exitCode'] == 0 and all(receipt['inputsStable'].values())
    assert not (folder/'admissibility.csv').exists()
    speed = float(np.linalg.norm(manifest['expectedVelocity']))
    assert speed > 0
    direct, scaled, details, locks = {}, {}, {}, {}
    for step in (0, 2, 3):
        name = 'restart_flow_%06d.csv' % step
        old_path = ROOT/'outputs/v163_validation/small02'/name
        new_path = folder/name
        old, a, nd = load_state(old_path, True, captured=True)
        new, b, new_nd = load_state(new_path, True, captured=True)
        assert list(a) == list(b) and nd == new_nd == 3 and old.shape == new.shape
        absolute = abs(new-old)
        denominator = np.maximum(1, abs(old))
        raw_error = absolute/denominator
        direct[str(step)] = float(raw_error.max())
        for index, field in enumerate(a):
            if field in ('Velocity_x', 'Velocity_y', 'Velocity_z'):
                denominator[:, index] = np.maximum(denominator[:, index], speed)
        scaled[str(step)] = float((absolute/denominator).max())
        assert scaled[str(step)] < 2e-13
        details[str(step)] = {field: dict(absoluteMaximum=float(absolute[:, i].max()),
            originalRelativeMaximum=float(raw_error[:, i].max())) for i, field in enumerate(a)}
        for path in (old_path, new_path):
            locks[str(path)] = sha(path)
    assert direct['3'] > 2e-13, '原首次失败必须保持'
    thermal = audit(folder)
    save(output, dict(passed=True, originalAllColumnGatePassed=False,
        originalAllColumnGate=2e-13, originalAllColumnErrors=direct,
        dimensionalVelocityScale=speed, scaledFieldErrors=scaled, perFieldDifferences=details,
        thermalErrors=thermal['errors'], sourceSHA256=locks, newFallbackEvents=0,
        actualNativeRuns=1, numericalBitwiseEquality=False,
        contractRejections=read(folder/'contract_rejections.json')['rejected'],
        supplementPlanSHA256=sha(ROOT/'docs/validation/v164/PLAN.md'),
        reviewSourceSHA256=sha(Path(__file__).resolve()), physicsQualified=False, newQualifiedArea=0))
    print(json.dumps(dict(passed=True, originalAllColumnGatePassed=False,
        scaledFieldErrors=scaled, thermalErrors=thermal['errors'])))


if __name__ == '__main__':
    main()
