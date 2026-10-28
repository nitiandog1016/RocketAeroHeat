# -*- coding: utf-8 -*-
"""独立核对原CV局部像、整列预算分配、保守/EOS提升及完整模型/原NS证据。"""
from pathlib import Path
import hashlib, json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]; VAL = ROOT / 'outputs/v83_validation'
case = sys.argv[1]; assert case.replace('_', '').isalnum()
target = VAL / (case + '_independent_local_coarse.json'); assert not target.exists()

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

file = VAL / case / 'local_coarse_evidence.mat'
s = loadmat(file, simplify_cells=True); c = s['coarse']; report = s['report']
parent_file = Path(s['parentFile']); parent = loadmat(parent_file, simplify_cells=True, variable_names=['coarse'])['coarse']
bindings = {str(f.relative_to(ROOT)): sha(f) for f in (file, parent_file)}
T = c['image'].tocsc(); L = c['candidateImage'].tocsc(); n = len(c['baseQ'])
retained = np.atleast_1d(c['retainedMixedColumns']).astype(int) - 1
removed = np.atleast_1d(c['removedMixedColumns']).astype(int) - 1
selected = np.atleast_1d(c['selectedCandidateColumns']).astype(int) - 1
cells = np.atleast_1d(c['localCellIDs']).astype(int) - 1
equations = np.atleast_1d(c['localEquationIDs']).astype(int) - 1
candidate_cells = np.atleast_1d(c['candidateCellIDs']).astype(int) - 1
candidate_equations = np.atleast_1d(c['candidateEquationIDs']).astype(int) - 1
nr = len(retained); nl = len(cells); nc = int(c['baseConservativeColumns'])
assert T.shape == (5 * n, nr + nl) and len(set(cells + n * equations)) == nl
assert set(retained) | set(removed) == set(range(parent['image'].shape[1]))
assert not set(retained) & set(removed) and set(range(nc)) <= set(retained)
assert (removed >= nc).all() and (cells >= 0).all() and (cells < n).all()
assert (equations >= 0).all() and (equations < 5).all()
assert int(c['imageBytes']) == 16 * T.nnz + 8 * (T.shape[1] + 1) <= c['maximumImageBytes'] == 268435456
assert c['constructionModelApplications'] == len(candidate_cells) == report['constructionModelApplications']
assert report['candidateCells'] * 5 == len(candidate_cells)
np.testing.assert_array_equal(cells, candidate_cells[selected]); np.testing.assert_array_equal(equations, candidate_equations[selected])
for a, b in ((T[:, :nr], parent['image'][:, retained]), (T[:, nr:], L[:, selected])):
    difference = a - b
    assert difference.nnz == 0 or np.max(np.abs(difference.data)) == 0
np.testing.assert_array_equal(c['baseQ'], parent['baseQ'])
np.testing.assert_array_equal(c['space']['geometry']['cells']['volume'], parent['space']['geometry']['cells']['volume'])
np.testing.assert_array_equal(c['columnLiftScale'][:nr], parent['columnLiftScale'][retained])
np.testing.assert_array_equal(c['localColumnLiftScale'], c['candidateColumnLiftScale'][selected])
np.testing.assert_array_equal(c['columnLiftScale'][nr:], c['localColumnLiftScale'])
G = (T.T @ T).toarray(); np.testing.assert_allclose(G, c['normalMatrix'], rtol=5e-10, atol=5e-13)
np.testing.assert_allclose(c['upper'].T @ c['upper'], G, rtol=5e-10, atol=5e-13)
np.testing.assert_allclose(np.diag(G), 1, rtol=1e-10, atol=1e-10)
assert np.linalg.eigvalsh(G)[0] > 0
P = c['space']['prolongation']; PP = c['mixedLift']['primitiveSpace']['prolongation']
ng = P.shape[1]; npg = PP.shape[1]
parent_selected = np.atleast_1d(c['mixedLift']['selectedPrimitiveColumns']).astype(int) - 1
Q = c['baseQ']; rho = Q[:, 0]; u = Q[:, 1:4] / rho[:, None]
cv = c['gas']['R'] / (c['gas']['gamma'] - 1)

def lift(weights):
    full = np.zeros(parent['image'].shape[1]); full[retained] = weights[:nr]
    full *= c['mixedLift']['columnLiftScale']
    base = P @ full[:nc].reshape((ng, 5), order='F')
    primitive = np.zeros(5 * npg); primitive[parent_selected] = full[nc:]
    fine = PP @ primitive.reshape((npg, 5), order='F')
    dr, du, dt = fine[:, 0], fine[:, 1:4], fine[:, 4]
    dq = np.column_stack((dr, u * dr[:, None] + rho[:, None] * du,
        Q[:, 4] / rho * dr + rho * (np.sum(u * du, axis=1) + cv * dt)))
    result = base + dq / c['scale']
    result[cells, equations] += weights[nr:] * c['localColumnLiftScale']
    return result

remaining = s['remaining'].ravel(order='F')
independent_weights = np.linalg.solve(G, T.T @ remaining)
# 保存因子重放用来核对提升实现；独立解/SVD另列。近相关小域不能声称系数逐值精确。
weights = np.linalg.solve(c['upper'], np.linalg.solve(c['upper'].T, T.T @ remaining))
expected = remaining - T @ weights
lift_error = float(np.linalg.norm(lift(weights) - s['correction']) / max(1, np.linalg.norm(s['correction'])))
assert lift_error < 1e-7
prediction_error = float(np.linalg.norm(expected - s['predicted'].ravel(order='F')) / max(1, np.linalg.norm(remaining)))
assert prediction_error < 5e-8
np.testing.assert_allclose(lift(s['probeCoefficients']), s['probeDirection'], rtol=5e-11, atol=1e-12)
probe_error = float(np.linalg.norm(T @ s['probeCoefficients'] - s['probeActual'].ravel(order='F')) / max(1, np.linalg.norm(s['probeActual'])))
assert probe_error < 1e-9
linearity = float(np.linalg.norm(s['combined'] - (.7 * s['correction'] - 1.3 * s['other'])) / max(1, np.linalg.norm(s['combined'])))
assert linearity < 1e-8

# 按每个释放步骤复算小Gram投影损失/释放字节比，确认只移除完整低收益原始变量列。
active = list(range(parent['image'].shape[1])); source_nnz = np.diff(parent['image'].tocsc().indptr)
prune_checks = []
for row in np.asarray(c['pruneHistory']).reshape((-1, 4)):
    pg = parent['normalMatrix'][np.ix_(active, active)]
    w = np.linalg.solve(pg, parent['image'][:, active].T @ remaining)
    inverse = np.linalg.inv(pg); losses = w * w / np.diag(inverse)
    eligible = np.array([j for j, identity in enumerate(active) if identity >= nc])
    scores = losses[eligible] / (16 * source_nnz[np.array(active)[eligible]] + 8)
    best = eligible[np.argmin(scores)]; removed_id = active[best]
    assert int(row[0]) - 1 == removed_id
    np.testing.assert_allclose(row[1], losses[best], rtol=2e-5, atol=1e-8)
    assert int(row[2]) == 16 * source_nnz[removed_id] + 8
    active.pop(best); assert len(active) == int(row[3]); prune_checks.append(removed_id + 1)
assert active == retained.tolist()
nf = len(c['space']['geometry']['inputGeometry']['F'])
np.testing.assert_array_equal(c['localOriginalFaceIDs'], cells % nf + 1)
np.testing.assert_array_equal(c['localLayerIDs'], cells // nf + 1)
for name, digest in c['liftSourceIdentity'].items():
    assert sha(ROOT / 'internalFlow/cutCell' / (name + '.m')) == digest
for name, digest in c['sourceIdentity'].items():
    path = ROOT / 'internalFlow/cutCell' / (name + '.m')
    if name == 'resolveNSCharacteristicBoundary3D':
        path = ROOT / 'internalFlow' / (name + '.m')
    assert sha(path) == digest

details = dict(retainedMixedColumns=nr, releasedPrimitiveParentColumns=prune_checks,
    selectedLocalColumns=nl, totalColumns=T.shape[1], imageBytes=int(c['imageBytes']),
    independentLiftRelativeDifference=lift_error, independentPredictionRelativeDifference=prediction_error,
    probeRelativeDifference=probe_error, fixedLinearityRelativeDifference=linearity,
    independentNormalSolveImageDifference=float(np.linalg.norm(T @ (independent_weights - weights)) / max(1, np.linalg.norm(remaining))),
    localPredictedRelativeDefect=float(np.linalg.norm(expected) / np.linalg.norm(s['rhs'])))
if report['kind'] == 'unit':
    reference = np.linalg.lstsq(T.toarray(), remaining, rcond=None)[0]
    error = float(np.linalg.norm(T @ (reference - weights)) / max(1, np.linalg.norm(remaining)))
    assert error < 2e-6; details['independentSVDImageDifference'] = error
    details['unitIndependentImageTolerance'] = 2e-6
    details['unitSmallGramNearlyDependent'] = bool(c['reciprocalConditionEstimate'] < 1e-10)
    details['coefficientWiseHighPrecisionClaimed'] = False
    errors = []
    for plus, minus, model in zip(s['validationPlus'], s['validationMinus'], s['validationModel']):
        measured = (plus - minus) / c['scale'] / (2 * s['h'])
        error = float(np.linalg.norm(measured - model) / max(1, np.linalg.norm(measured)))
        assert error < 2e-6; errors.append(error)
    assert len(errors) == 2 and report['actualOriginalNSEvaluations'] == 5
    details['originalNSDirectionFDErrors'] = errors
else:
    source = Path(s['sourceFile']); bindings[str(source.relative_to(ROOT))] = sha(source)
    r = loadmat(source, simplify_cells=True, variable_names=['result'])['result']
    logs = r['linearLog']; log = logs if isinstance(logs, dict) else list(logs)[-1]
    np.testing.assert_array_equal(s['remaining'], r['initialRate'] / r['residualScale'] + log['steadyDirectionEvidence']['value'])
    measured = (s['validationPlus'] - s['validationMinus']) / c['scale'] / (2 * s['h'])
    error = float(np.linalg.norm(measured - s['validationModel']) / max(1, np.linalg.norm(measured)))
    assert error < 2e-6 and report['actualOriginalNSEvaluations'] == 2
    details['originalNSDirectionFDErrors'] = [error]
if case in ('localunit02', 'localrocket02'):
    assert np.asarray(s['cellIDs']).dtype.kind == 'i'
    # v7会把整数值double紧缩为uint存储；需用mat_dtype=True检查MATLAB类，不能误判紧缩类型。
    if case == 'localunit02':
        typed = loadmat(file, simplify_cells=True, mat_dtype=True, variable_names=['coarse'])['coarse']
        assert np.asarray(typed['localCellIDs']).dtype.kind == 'f'
        details['integerInputCanonicalizationRegressionPassed'] = True
assert not report['physicalMeshChanged'] and not report['physicalSolution'] and report['newProductionSolvedAreaM2'] == 0
assert all(sha(ROOT / path) == digest for path, digest in bindings.items())
out = dict(status='passed-independent-original-CV-local-enrichment-lift-budget-projection-and-model-evidence',
    case=case, processID=os.getpid(), details=details, actualFlowUpdates=0,
    physicalSolution=False, newProductionSolvedAreaM2=0, inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2) + '\n', encoding='utf8'); print(json.dumps(out))
