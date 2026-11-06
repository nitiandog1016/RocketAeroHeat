"""核对原生线性历史，分开记录内层误差、非线性残差及物理资格。"""
import argparse
import csv
from pathlib import Path

import numpy as np

from runEmbeddedSU2 import atomic_write_json
from su2IdealAirInterface import settings
from su2ThermalRuntime import ROOT, read, sha


def review(folder):
    """模块一。只接受终态锁定的完整历史，保留原生相对残差的log10定义。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    receipt = read(folder/'engine_process.json')
    assert receipt['status'] == 'terminal'
    assert all(receipt['inputsStable'].values())
    normal = receipt['exitCode'] == 0 and not any(receipt.get(k, False)
        for k in ('timedOut', 'cancelled', 'supervisorError', 'outputBudgetExceeded'))
    assert sha(folder/'history.csv') == receipt['outputSHA256']['history.csv']
    assert sha(folder/'case.cfg') == receipt['inputSHA256'][str(folder/'case.cfg')]
    config = settings((folder/'case.cfg').read_text('utf8'))
    assert 'LINSOL' in config['HISTORY_OUTPUT']
    with (folder/'history.csv').open() as stream:
        rows = list(csv.reader(stream, skipinitialspace=True))
    keys = [s.strip() for s in rows[0]]
    data = np.asarray(rows[1:], dtype=float)
    assert data.ndim == 2 and data.shape[1] == len(keys) and np.isfinite(data).all()
    assert np.array_equal(data[:, keys.index('Inner_Iter')], np.arange(len(data)))
    limit = int(config['LINEAR_SOLVER_ITER'])
    target_log = float(np.log10(float(config['LINEAR_SOLVER_ERROR'])))
    groups = {}
    for name, iteration, error in [('flow', 'Linear_Solver_Iterations', 'LinSolRes'),
                                   ('turbulence', 'LinSolIterTurb', 'LinSolResTurb')]:
        counts, residuals = data[:, keys.index(iteration)], data[:, keys.index(error)]
        assert np.all(counts == counts.astype(int)) and min(counts) >= 0 and max(counts) <= limit
        failed_at_limit = (counts == limit) & (residuals > target_log)
        groups[name] = dict(iterationRange=[int(min(counts)), int(max(counts))],
            log10RelativeResidualRange=[float(min(residuals)), float(max(residuals))],
            targetReachedUpdates=int(sum(residuals <= target_log)),
            atLimitAboveTargetUpdates=int(sum(failed_at_limit)),
            atLimitAboveTargetFraction=float(np.mean(failed_at_limit)))
    fields = ['rms[Rho]', 'rms[RhoE]', 'rms[k]', 'rms[w]']
    residual = data[:, [keys.index(k) for k in fields]]
    full = normal and len(data) == int(config['ITER']) == 50
    # 模块二。候选触发遵循事前50更新合同；内层改善仍需独立物理验收。
    eligible = full and any(g['atLimitAboveTargetUpdates'] > 25 for g in groups.values())
    return dict(passedHistoryIdentity=True, updates=len(data), full50UpdateWindow=full,
        nativeNormalTerminal=normal, nativeExitCode=receipt['exitCode'],
        countDefinition='native BCGSTAB returns zero-based loop index i; zero may mean one executed iteration',
        linearSolver=config['LINEAR_SOLVER'], preconditioner=config['LINEAR_SOLVER_PREC'],
        linearIterationLimit=limit, targetLog10RelativeResidual=target_log, groups=groups,
        residualFields=fields, finalResiduals=residual[-1].tolist(),
        averagingWindowUpdates=min(20, len(data)), last20MeanResiduals=np.mean(residual[-20:], axis=0).tolist(),
        originalResidualThresholdPassed=bool(np.all(residual[-1] < -8)),
        conditionalCandidateEligible=eligible, nativeElapsedSeconds=receipt['elapsedSeconds'],
        inputSHA256={str(folder/name):sha(folder/name) for name in ('engine_process.json', 'history.csv', 'case.cfg')},
        reviewerSHA256=sha(Path(__file__)), physicsQualified=False, newQualifiedArea=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.relative_to(ROOT/'outputs')
    assert not output.exists()
    result = review(args.folder)
    atomic_write_json(output, result)
    print(result)
