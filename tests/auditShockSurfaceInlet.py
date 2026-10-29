# -*- coding: utf-8 -*-
"""用闭式完全气体与独立Brent热力根，审核共同三维激波跃迁和物理入口梯度。"""
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq

from auditThermalPressureSurfaceFlow import air, sha
from auditPressureDrivenSurfaceFlow import rows

ROOT = Path(__file__).resolve().parents[1]


def thermal(T, p, gas):
    if gas['thermalPropertyMode'] == 'thermallyPerfect':
        cp, h, s = air(float(T))
        reference_s = air(float(gas['referenceTemperatureK']))[2]
    else:
        cp = gas['gamma']*gas['R']/(gas['gamma']-1)
        h, s, reference_s = cp*T, cp*np.log(T), cp*np.log(gas['referenceTemperatureK'])
    return cp, h, s-reference_s-gas['R']*np.log(p/gas['referencePressurePa'])


def jump(upstream, gas, A, orientation):
    """模块一：给定几何法向，用压缩体积分数求根；不读MATLAB密度比/激波角。"""
    p, T, rho = [float(upstream[k]) for k in ['pressure', 'temperature', 'density']]
    u = np.asarray(upstream['velocityVector']).reshape(3)
    n = orientation*np.cross(A[:, 0], A[:, 1]); n /= np.linalg.norm(n)
    w = float(u@n); cp, h, S = thermal(T, p, gas); R = gas['R']
    Mn = w/np.sqrt(cp/(cp-R)*R*T)
    assert Mn > 1
    if gas['thermalPropertyMode'] == 'caloricallyPerfect':
        g = gas['gamma']; r = (g+1)*Mn**2/((g-1)*Mn**2+2)
    else:
        H = h+w*w/2
        assert thermal(gas['temperatureRangeK'][1], p, gas)[1] >= H

        def energy_fraction(fraction):
            v = 1-fraction
            t2 = (p+rho*w*w*fraction)*v/(rho*R)
            if t2 < gas['temperatureRangeK'][0]:
                return -1.  # 焓和法向动能均低于上游；不外推物性。
            if t2 > gas['temperatureRangeK'][1]:
                return 1.  # 单调焓大于已核对的总焓；不外推物性。
            return (thermal(t2, p, gas)[1]+.5*(w*v)**2-H)/fraction

        fraction = brentq(energy_fraction, 1e-6, .9999, xtol=2e-14, rtol=2e-14)
        r = 1/(1-fraction)
    p2 = p+rho*w*w*(1-1/r); rho2 = rho*r; T2 = p2/(rho2*R)
    u2 = u+(w/r-w)*n; speed = np.linalg.norm(u2)
    cp2, h2, S2 = thermal(T2, p2, gas)
    mass = rho2*(u2@n); momentum = mass*u2+p2*n
    np.testing.assert_allclose(momentum, rho*w*u+p*n, rtol=3e-13, atol=2e-9)
    np.testing.assert_allclose(h2+speed**2/2, h+u@u/2, rtol=4e-13, atol=2e-7)
    assert S2 > S
    return np.array([p2, T2, rho2, speed, S2]), u2, n


def surface(q, K, tip):
    return tip+np.r_[-.5*q@K@q, q], np.vstack([-K@q, np.eye(2)])


def audit(evidence, process, target):
    """模块二：真实源身份、全部共同曲面点及两向导数独立回算。"""
    assert not target.exists()
    p = json.loads(process.read_text(encoding='utf-8-sig'))
    assert p['status'] == 'terminal' and p['exitCode'] == 0 and all(p['sourceStableAfter'].values())
    s = loadmat(evidence, simplify_cells=True)
    geometry = s['geometry']; assert geometry['F'].shape == (19752, 3)
    assert s['fs']['alpha'] == 20 and s['fs']['Mach'] == 5
    assert not s['report']['oldExternalStateUsedAsInlet'] and not s['report']['bodyMatched']
    errors, derivative_errors = [], []
    for inlet in rows(s['cases']):
        gas = inlet['gas']; up = inlet['uniformUpstream']; states = rows(inlet['states'])
        assert not inlet['shockShapeSolved'] and not inlet['bodyMatched']
        assert not inlet['boundaryLayerEdgeAvailable'] and not inlet['productionQualified']
        for q, saved, stored_geometry in zip(inlet['coordinates'], states, rows(inlet['geometrySamples'])):
            position, A = surface(q, s['K'], s['tip'])
            np.testing.assert_allclose(position, stored_geometry['position'], rtol=0, atol=2e-17)
            np.testing.assert_allclose(A, stored_geometry['firstDerivatives'], rtol=0, atol=2e-16)
            expected, velocity, n = jump(up, gas, A, inlet['normalOrientation'])
            actual = np.array([saved[k] for k in ['pressure', 'temperature', 'density', 'speed', 'entropy']])
            error = float(np.max(abs(actual-expected)/np.maximum(1, abs(expected))))
            assert error < 2e-10; errors.append(error)
            np.testing.assert_allclose(saved['velocityVector'], velocity, rtol=3e-12, atol=3e-9)
            np.testing.assert_allclose(saved['normal'], n, rtol=0, atol=2e-15)
            # 独立五点差分来自重新求解的热力根，不对保存的输出作插值。
            numerical = np.zeros((5, 2)); h = 2e-5
            for axis in range(2):
                offset = np.eye(2)[axis]*h; values = []
                for factor in [2, 1, -1, -2]:
                    _, local_A = surface(q+factor*offset, s['K'], s['tip'])
                    values.append(jump(up, gas, local_A, inlet['normalOrientation'])[0])
                numerical[:, axis] = (-values[0]+8*values[1]-8*values[2]+values[3])/(12*h)
            partial = np.vstack([saved['partialDerivatives'][k] for k in ['pressure', 'temperature', 'density', 'speed', 'entropy']])
            de = np.linalg.norm(partial-numerical, axis=1)/np.maximum(1, np.linalg.norm(partial, axis=1))
            assert np.max(de) < 2e-7; derivative_errors.append(float(np.max(de)))
            # QR最小范数解独立恢复曲面物理梯度，不能把标签导数直接输出为梯度。
            for k, name in enumerate(['pressure', 'temperature', 'density', 'speed', 'entropy']):
                reference = np.linalg.lstsq(A.T, partial[k], rcond=None)[0]
                np.testing.assert_allclose(saved['physicalGradients'][name], reference, rtol=3e-12, atol=3e-9)
                assert abs(np.dot(reference, n)) < 1e-9*max(1, np.linalg.norm(reference))
            np.testing.assert_allclose(saved['partialDerivatives']['totalEnthalpy'], 0, rtol=0, atol=2e-7)
            total = saved['inletTotalState']; T0, p0 = total['totalTemperature'], total['totalPressure']
            _, H, entropy0 = thermal(T0, p0, gas)
            np.testing.assert_allclose(H, saved['totalEnthalpy'], rtol=2e-12, atol=2e-7)
            np.testing.assert_allclose(entropy0, saved['entropy'], rtol=2e-12, atol=2e-8)
            assert not saved['bodyMatched'] and not saved['boundaryLayerEdgeAvailable']
    for plane in rows(s['planeCases']):
        expected, v, n = jump(plane['upstream'], plane['gas'], plane['geometry']['firstDerivatives'], 1)
        actual = np.array([plane['state'][k] for k in ['pressure', 'temperature', 'density', 'speed', 'entropy']])
        np.testing.assert_allclose(actual, expected, rtol=3e-11, atol=3e-8)
    # 当前源与实际任务绑定；冻结数据重放时仅验证明确的非输出来源。
    source_verified = 0
    for rel, expected_sha in p['sourceSHA256'].items():
        file = ROOT/rel
        if rel.startswith('outputs/'):
            continue
        assert file.is_file() and sha(file) == expected_sha, rel
        source_verified += 1
    report = dict(status='passed-independent-common-3D-shock-state-and-gradient-audit',
                  evidenceSHA256=sha(evidence), processSHA256=sha(process),
                  commonShockPoints=len(errors), planarAnalyticCases=len(rows(s['planeCases'])),
                  maximumRelativeStateDifference=max(errors), maximumRelativeDerivativeDifference=max(derivative_errors),
                  sourceSHARecordsVerified=source_verified, bodyMatched=False,
                  boundaryLayerEdgeAvailable=False, experimentalValidation=False, newProductionSolvedAreaM2=0)
    target.write_bytes((json.dumps(report, indent=2)+'\n').encode('utf-8'))
    print(json.dumps(report))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('evidence', type=Path); parser.add_argument('process', type=Path)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args(); audit(args.evidence, args.process, args.report)
