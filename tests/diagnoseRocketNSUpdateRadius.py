# -*- coding: utf-8 -*-
"""定位真实守恒更新的状态可行半径及局部残差变化；只读已存NS证据。"""
from pathlib import Path
import hashlib
import json
import os
import sys
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def first_root(c0, c1, c2):
    """从严格合法初态出发，计算每个二次不等式首次离开可行域的正根。"""
    assert np.all(c0 > 0)
    roots = np.full(c0.shape, np.inf)
    linear = c2 == 0
    use = linear & (c1 < 0)
    roots[use] = -c0[use] / c1[use]
    disc = c1 * c1 - 4 * c2 * c0
    use = (~linear) & (disc > 0)
    indices = np.flatnonzero(use)
    b = c1[use]
    q = -.5 * (b + np.where(b >= 0, 1, -1) * np.sqrt(disc[use]))
    r1, r2 = q / c2[use], c0[use] / q
    # c(0)>0，正根处导数<0才是首次进入非法一侧，不将切触点当越界。
    for candidate in (r1, r2):
        valid = (candidate > 0) & (c1[use] + 2 * c2[use] * candidate < 0)
        ii = indices[valid]
        roots[ii] = np.minimum(roots[ii], candidate[valid])
    return roots


def run(case):
    assert case.replace('_', '').isalnum()
    target = VAL / (case + '_update_radius.json')
    assert not target.exists()
    source = VAL / case / 'wall_velocity_krylov_evidence.mat'
    domain_file = ROOT / 'outputs/v80_validation/domain02/shell_domain_evidence.mat'
    bindings = {p.relative_to(ROOT).as_posix(): sha(p) for p in (source, domain_file)}
    s = loadmat(source, simplify_cells=True, variable_names=['result', 'gas'])
    r, gas = s['result'], s['gas']
    domain = loadmat(domain_file, simplify_cells=True, variable_names=['domain'])['domain']
    nfaces = int(domain['quality']['sourceFaces'])
    nlayer = int(domain['quality']['layerCount'])
    rows = lambda x: [x] if isinstance(x, dict) else list(x)
    assert len(rows(r['trace'])) == 1, '本诊断要求显式选择单次接受的真实研究段。'
    tr = rows(r['trace'])[0]
    q, scale = tr['initialQ'], r['residualScale']
    dq = tr['linearUpdateU'] * scale
    rho, drho = q[:, 0], dq[:, 0]
    m, dm, energy, de = q[:, 1:4], dq[:, 1:4], q[:, 4], dq[:, 4]
    cv = gas['R'] / (gas['gamma'] - 1)
    low_t, high_t = gas['temperatureRange']
    u = m / rho[:, None]
    temp = (energy - .5 * np.sum(m * m, axis=1) / rho) / (rho * cv)
    assert np.all((rho > 0) & (temp > low_t) & (temp < high_t))

    def location(index):
        return dict(cellID=int(index + 1), layer=int(index // nfaces + 1),
                    sourceFacet=int(index % nfaces + 1), centroidM=domain['cells']['centroid'][index].tolist())

    density_radius = np.full(rho.shape, np.inf)
    density_radius[drho < 0] = -rho[drho < 0] / drho[drho < 0]
    radii = {'positiveDensity': density_radius}
    for name, limit, sign in [('minimumTemperature', low_t, 1), ('maximumTemperature', high_t, -1)]:
        c0 = rho * energy - .5 * np.sum(m * m, axis=1) - cv * limit * rho * rho
        c1 = rho * de + drho * energy - np.sum(m * dm, axis=1) - 2 * cv * limit * rho * drho
        c2 = drho * de - .5 * np.sum(dm * dm, axis=1) - cv * limit * drho * drho
        radii[name] = first_root(sign * c0, sign * c1, sign * c2)
    limits = []
    for name, values in radii.items():
        index = int(np.argmin(values))
        limits.append(dict(constraint=name, firstExitAlpha=float(values[index]), location=location(index)))
    limiting = min(limits, key=lambda x: x['firstExitAlpha'])
    radius = limiting['firstExitAlpha']

    def legal(alpha):
        candidate = q + alpha * dq
        rr = candidate[:, 0]
        if np.any(rr <= 0):
            return False
        tt = (candidate[:, 4] - .5 * np.sum(candidate[:, 1:4]**2, axis=1) / rr) / (rr * cv)
        return bool(np.all((tt >= low_t) & (tt <= high_t)))

    # 独立直接五量回代二分检查，未重算面重构/边界或完整NS。
    left, right = radius * .5, radius * 1.5
    assert legal(left) and not legal(right)
    for _ in range(65):
        mid = (left + right) * .5
        if legal(mid):
            left = mid
        else:
            right = mid
    relative_radius_difference = abs(radius - (left + right) * .5) / radius
    assert relative_radius_difference < 2e-11
    np.testing.assert_array_equal(q + tr['alpha'] * dq, r['Q'])

    sound = np.sqrt(gas['gamma'] * gas['R'] * temp)
    du = (dm - u * drho[:, None]) / rho[:, None]
    dt = (de - np.sum(u * dm, axis=1) + (.5 * np.sum(u * u, axis=1) - cv * temp) * drho) / (rho * cv)
    metrics = {}
    for name, values in [('densityFraction', np.abs(drho) / rho),
                         ('velocityOverLocalSound', np.linalg.norm(du, axis=1) / sound),
                         ('temperatureFraction', np.abs(dt) / temp)]:
        index = int(np.argmax(values))
        metrics[name] = dict(maximumFullDirection=float(values[index]), location=location(index),
                             maximumAtAcceptedAlpha=float(tr['alpha'] * values[index]),
                             quantilesFullDirection=np.quantile(values, [.5, .9, .99, .999]).tolist())

    initial = tr['initialRate'] / scale
    final = r['rate'] / scale
    ev = rows(r['linearLog'])[0]['steadyDirectionEvidence']
    np.testing.assert_array_equal(ev['directionU'], tr['linearUpdateU'])
    predicted = initial + tr['alpha'] * ev['value']
    nonlinear = final - predicted
    denom = np.linalg.norm(initial)
    expected_defect = np.linalg.norm(initial + ev['value']) / denom
    np.testing.assert_allclose(expected_defect, rows(r['linearLog'])[0]['fullSystemRelativeDefect'], rtol=1e-10)
    layer_energy = np.sum(nonlinear * nonlinear, axis=1).reshape(nlayer, nfaces).sum(axis=1)
    fraction = layer_energy / layer_energy.sum()
    hotspots = []
    for name, values in [('initial', initial), ('accepted', final), ('acceptedNonlinearRemainder', nonlinear)]:
        flat = np.argpartition(np.abs(values).ravel(), -12)[-12:]
        flat = flat[np.argsort(-np.abs(values).ravel()[flat])]
        cells, equations = np.unravel_index(flat, values.shape)
        hotspots.append(dict(kind=name, largestEntries=[dict(location=location(int(i)),
                            equation=['mass', 'xMomentum', 'yMomentum', 'zMomentum', 'energy'][j],
                            signedScaledResidualPerSecond=float(values[i, j])) for i, j in zip(cells, equations)]))
    trials = rows(r['lineLog'])
    report = dict(status='verified-real-conservative-state-radius-and-accepted-local-residual-diagnostic',
                  case=case, processID=os.getpid(), constraints=limits, limitingConstraint=limiting,
                  independentBisectionRelativeDifference=relative_radius_difference,
                  acceptedAlpha=float(tr['alpha']), acceptedAlphaFractionOfCellStateRadius=float(tr['alpha'] / radius),
                  fullDirectionRelativeStateMetrics=metrics,
                  originalNSDirectionDefect=float(expected_defect),
                  acceptedNonlinearRemainderRelativeToBaseResidual=float(np.linalg.norm(nonlinear) / denom),
                  acceptedNonlinearRemainderRelativeToLinearChange=float(np.linalg.norm(nonlinear) / np.linalg.norm(tr['alpha'] * ev['value'])),
                  nonlinearRemainderLayerSquareFractions=fraction.tolist(), residualHotspots=hotspots,
                  candidateTrials=[dict(alpha=float(x['alpha']), accepted=bool(x['accepted']),
                       candidateRMS=None if np.isnan(x['candidateRmsMerit']) else float(x['candidateRmsMerit']),
                       candidateMaximum=None if np.isnan(x['candidateMaximumResidual']) else float(x['candidateMaximumResidual']),
                       errorIdentifier=x['errorIdentifier'] if isinstance(x['errorIdentifier'], str) else '') for x in trials],
                  newActualNSEvaluations=0, flowStateChanged=False, interfaceOrBoundaryRecomputed=False,
                  physicalSolution=False, newProductionSolvedAreaM2=0, inputSHA256=bindings, scriptSHA256=sha(__file__))
    assert all(sha(ROOT / p) == h for p, h in bindings.items())
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf8')
    print(json.dumps({k: report[k] for k in ['status', 'limitingConstraint', 'acceptedAlphaFractionOfCellStateRadius',
                      'fullDirectionRelativeStateMetrics', 'acceptedNonlinearRemainderRelativeToLinearChange',
                      'nonlinearRemainderLayerSquareFractions']}))


if __name__ == '__main__':
    run(sys.argv[1])
