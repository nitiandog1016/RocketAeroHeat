# -*- coding: utf-8 -*-
"""独立复核原始变量度量、总方向增广问题和完整父矩阵实际剩余，不晋升物理解。"""
from pathlib import Path
import os
import json
import sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT, VAL, sha, local
from screenRocketParentImageScale import setup

mode, case = sys.argv[1:3]
assert mode in ('unit', 'real') and case.replace('_', '').isalnum()
folder = VAL / case
target = VAL / (case + '_independent_regularized_direction.json')
assert not target.exists()
bindings = {}


def read(path, variables=None):
    path = Path(path)
    path = path if path.is_absolute() else ROOT / path
    bindings[path.relative_to(ROOT).as_posix()] = sha(path)
    return loadmat(path, simplify_cells=True, variable_names=variables)


def primitive(q, gas):
    rho = q[:, 0]
    u = q[:, 1:4] / rho[:, None]
    t = (q[:, 4] / rho - .5 * np.sum(u * u, axis=1)) * (gas['gamma'] - 1) / gas['R']
    return np.column_stack([rho, u, t])


def weight(q, scale, gas, value):
    p = primitive(q, gas)
    dq = value * scale
    rho, u, t = p[:, 0], p[:, 1:4], p[:, 4]
    cv = gas['R'] / (gas['gamma'] - 1)
    dr = dq[:, 0]
    du = (dq[:, 1:4] - dr[:, None] * u) / rho[:, None]
    dt = (dq[:, 4] - np.sum(u * dq[:, 1:4], axis=1) + (.5*np.sum(u*u, axis=1)-cv*t)*dr) / (rho*cv)
    return np.column_stack([dr/rho, du/np.sqrt(gas['gamma']*gas['R']*t)[:, None], dt/t])


def rel(a, b):
    return float(np.linalg.norm(a-b)/max(1, np.linalg.norm(b)))


results = []
if mode == 'unit':
    s = read(folder/'regularized_parent_unit_evidence.mat')
    for a in s['cases']:
        q, gas, scale = a['Q'], a['gas'], a['scale']
        n, N = len(q), q.size
        p = primitive(q, gas)
        normalization = np.column_stack([p[:, 0], np.tile(np.sqrt(gas['gamma']*gas['R']*p[:, 4])[:, None], (1, 3)), p[:, 4]])
        Wfd, W = np.empty((N, N)), np.empty((N, N))
        for col in range(N):
            v = np.zeros((n, 5)); v.ravel(order='F')  # 显式保留MATLAB按列的未知量顺序。
            v[col % n, col // n] = 1
            h = 1e-5
            Wfd[:, col] = (((primitive(q+h*v*scale, gas)-primitive(q-h*v*scale, gas))/(2*h))/normalization).ravel(order='F')
            W[:, col] = weight(q, scale, gas, v).ravel(order='F')
        fd_error, matrix_error = rel(Wfd, W), rel(a['W'], W)
        assert fd_error < 1e-8 and matrix_error < 1e-12
        detail = a['detail']; penalty = detail['penaltyInOriginalResidualUnits']
        A = np.row_stack([a['C'], penalty*W]); b = np.concatenate([a['rhs'].ravel(order='F'), np.zeros(N)])
        dense = np.linalg.lstsq(A, b, rcond=None)[0]
        dense_error = rel(a['value'].ravel(order='F'), dense)
        assert dense_error < 1e-8
        lift = a['initial'].ravel(order='F') + a['M'] @ (a['columnScale']*a['coordinates']).ravel(order='F')
        lift_error = rel(a['value'].ravel(order='F'), lift); assert lift_error < 1e-12
        initial_metric = weight(q, scale, gas, a['initial'])
        denominator = np.linalg.norm(a['rhs'])
        expected_penalty = detail['options']['regularizationStrength']*denominator/max(1, np.linalg.norm(initial_metric))
        np.testing.assert_allclose(penalty, expected_penalty, rtol=1e-13)
        relative_c = np.linalg.norm(a['C']@a['value'].ravel(order='F')-a['rhs'].ravel(order='F'))/denominator
        np.testing.assert_allclose(relative_c, detail['relativeCResidual'], rtol=1e-10, atol=1e-14)
        results.append(dict(allPrimitiveColumnsFDDifference=fd_error, metricMatrixDifference=matrix_error,
                            denseAugmentedSolutionDifference=dense_error, totalCoordinateLiftDifference=lift_error,
                            originalCRelativeResidual=relative_c))
else:
    s = read(folder/'regularized_parent_candidate_evidence.mat')
    base = Path(s['base'])
    packed = read(base/'current_packed_evidence.mat', ['packed'])['packed']
    factor = read(base/'current_factor_evidence.mat', ['factor'])['factor']
    A, plan, order = setup(packed, factor); del packed, factor
    cscale = read(base/'current_image_scale_evidence.mat', ['columnScale'])['columnScale']
    correction = local(plan, (s['coordinates'][order]*cscale[order]).ravel()).reshape(len(order), 5)
    lift = s['initial'].copy(); lift[order] += correction
    lift_error = rel(lift, s['value']); assert lift_error < 1e-8
    image = np.empty_like(lift); image[order] = (A@s['value'][order].ravel()).reshape(len(order), 5)
    relative_c = float(np.linalg.norm(s['rhs']-image)/np.linalg.norm(s['rhs']))
    np.testing.assert_allclose(relative_c, s['detail']['relativeCResidual'], rtol=1e-12)
    source = read(s['cacheFile'], ['Q', 'scale', 'gas'])
    metric_value = weight(source['Q'], source['scale'], source['gas'], s['value'])
    initial_metric = weight(source['Q'], source['scale'], source['gas'], s['initial'])
    detail = s['detail']; penalty = detail['penaltyInOriginalResidualUnits']; denom = np.linalg.norm(s['rhs'])
    normed_objective = relative_c**2 + (penalty*np.linalg.norm(metric_value)/denom)**2
    for key, expected in [('initialPrimitiveUpdateNorm', np.linalg.norm(initial_metric)),
                          ('primitiveUpdateNorm', np.linalg.norm(metric_value)),
                          ('maximumDensityFraction', np.max(abs(metric_value[:, 0]))),
                          ('maximumVelocityOverSound', np.max(np.linalg.norm(metric_value[:, 1:4], axis=1))),
                          ('maximumTemperatureFraction', np.max(abs(metric_value[:, 4])))]:
        np.testing.assert_allclose(detail[key], expected, rtol=1e-12)
    # 两个已独立限制为1e-12的分量范数经平方相加，需要传播误差，不能复用分量相对门槛。
    # 保持上方完整C和原始变量范数的1e-12检查；不改MATLAB原C/NS目标。
    dc = 1e-12*abs(detail['relativeCResidual'])
    dw = 1e-12*abs(detail['primitiveUpdateNorm'])
    objective_bound = (2*abs(detail['relativeCResidual'])*dc + dc*dc +
                       (penalty/denom)**2*(2*abs(detail['primitiveUpdateNorm'])*dw + dw*dw) +
                       16*np.finfo(float).eps*abs(detail['normalizedAugmentedObjective']))
    assert abs(normed_objective-detail['normalizedAugmentedObjective']) <= objective_bound
    model = read(folder/'regularized_parent_model_evidence.mat')
    np.testing.assert_array_equal(model['value'], s['value'])
    full = float(np.linalg.norm(s['rhs']+model['fullModel'])/denom)
    np.testing.assert_allclose(full, model['report']['completeSelectedModelRelativeDefect'], rtol=1e-12)
    assert bool(model['report']['modelPointOneMet']) == (full <= .1)
    assert not model['report']['physicalStateChanged'] and model['report']['newActualNSEvaluations'] == 0
    results.append(dict(originalCVs=len(order), totalCoordinateLiftDifference=lift_error,
                        originalCRelativeResidual=relative_c, primitiveUpdateNorm=float(np.linalg.norm(metric_value)),
                        normalizedAugmentedObjective=normed_objective, fullSelectedModelRelativeDefect=full,
                        objectiveReferenceAbsoluteDifference=abs(normed_objective-detail['normalizedAugmentedObjective']),
                        propagatedComponentObjectiveAbsoluteBound=objective_bound,
                        fullModelPointOneMet=full <= .1, allRealBidiagonalizationColumnsIndependentlyVerified=False))

assert all(sha(ROOT/p) == h for p, h in bindings.items())
out = dict(status='passed-independent-primitive-metric-and-total-augmented-parent-direction', case=case, mode=mode,
           processID=os.getpid(), checks=results, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
           physicalSolution=False, inputSHA256=bindings, scriptSHA256=sha(__file__),
           helperSourceSHA256={name:sha(ROOT/'tests'/name) for name in ('screenRocketFullFLeastSquares.py','screenRocketParentImageScale.py')})
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256','helperSourceSHA256')}))
