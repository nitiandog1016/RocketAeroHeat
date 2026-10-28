"""独立公开Riemann参考核对单元/界面前瞻、阻挡分区、集合及真实温度导数。"""
from pathlib import Path
import json
import os
import sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix
from riemannReservoirReference import ROOT, sha, exact, PUBLIC_BINDINGS
from auditRocketTemperatureProjection import row_reference


def vector(value, dtype=None):
    return np.asarray(value, dtype=dtype).ravel()


def primitive(q, gas):
    rho = q[0]; velocity = q[1:4] / rho
    temperature = (q[4] - .5*rho*(velocity @ velocity)) / (rho*gas['R']/(gas['gamma']-1))
    return np.r_[rho, velocity, temperature]


def preview_reference(candidate, cut, gas, report):
    """独立逐面求解全部可检查输入；非法单元只进入阻挡集，不能当作合法面。"""
    f = cut['openFaces']; minus = vector(f['cellMinus'], int); plus = vector(f['cellPlus'], int)
    ids = np.flatnonzero((vector(f['area']) > 0) & (minus > 0) & (plus > 0))
    normals = np.asarray(f['normal']).reshape(-1, 3)
    states = np.row_stack([primitive(q, gas) for q in candidate]); low, high = gas['temperatureRange']
    valid = (states[:, 4] >= low) & (states[:, 4] <= high)
    eligible = valid[minus[ids]-1] & valid[plus[ids]-1]
    assessed = ids[eligible]; blocked = ids[~eligible]; failed = []; temperatures = []
    for face in assessed:
        sample, _, _ = exact(states[minus[face]-1], states[plus[face]-1], normals[face], gas)
        if sample[4] < low or sample[4] > high:
            failed.append(face+1); temperatures.append(float(sample[4]))
    np.testing.assert_array_equal(vector(report['assessedFaceIDs'], int), assessed+1)
    np.testing.assert_array_equal(vector(report['blockedFaceIDs'], int), blocked+1)
    np.testing.assert_array_equal(vector(report['blockedByOutOfRangeCellIDs'], int), np.flatnonzero(~valid)+1)
    np.testing.assert_array_equal(vector(report['predictedViolatingInterfaceFaceIDs'], int), failed)
    np.testing.assert_allclose(vector(report['predictedInterfaceTemperatureK']), temperatures, rtol=1e-11, atol=1e-9)
    assert report['completeInteriorFaceCount'] == len(ids) and report['allInputFacesAccountedFor']
    assert bool(report['allInteriorFacesAssessed']) == (not len(blocked))
    assert not report['temperatureDomainRelaxed'] and not report['nonTemperatureFailuresSuppressed']
    assert report['newActualNSEvaluations'] == report['newFlowAcceptances'] == 0
    assert report['originalKernelCalls'] <= 2*max(1, len(ids))+1
    assert report['kernelFaceInputsIncludingRechecks'] <= 16*max(1, len(ids))
    assert report['additionalArrayBoundBytes'] <= 256*1024**2
    return dict(interiorFaces=len(ids), independentlyAssessedFaces=len(assessed), blockedFaces=(blocked+1).tolist(),
                violatingFaces=failed, publicViolatingTemperatureK=temperatures,
                allAssessableFacesIndependentlyRecomputed=True)


def selected_reference(Q, candidate, info, cut, gas, rows, detail):
    """从原采样/候选单元与独立界面参考重建并集，并核对原CV稀疏导数映射。"""
    preview = preview_reference(candidate, cut, gas, detail['interfacePreview'])
    audit = info['exactEulerFluxAudit']; ids = vector(audit['faceID'], int)
    sample = np.asarray(audit['samplePrimitive']).reshape(-1, 5); low, high = gas['temperatureRange']
    margin = np.minimum(sample[:, 4]-low, high-sample[:, 4]); original = ids[margin <= 10]
    states = np.row_stack([primitive(q, gas) for q in candidate])
    bad = np.flatnonzero((states[:, 4] < low) | (states[:, 4] > high))+1
    required = vector(detail['requiredOriginalStateFaceIDs'], int)
    np.testing.assert_array_equal(vector(detail['originalBandFaceIDs'], int), original)
    np.testing.assert_array_equal(vector(detail['predictedViolatingCellIDs'], int), bad)
    f = cut['openFaces']; minus = vector(f['cellMinus'], int); plus = vector(f['cellPlus'], int)
    slots = np.searchsorted(ids, required); np.testing.assert_array_equal(ids[slots], required)
    regions = vector(audit['sampleRegion'], int)[slots]; assert set(regions) <= {1, 4}
    np.testing.assert_array_equal(np.where(regions == 1, minus[required-1], plus[required-1]), bad)
    expected = np.union1d(np.union1d(original, required), preview['violatingFaces']).astype(int)
    actual = vector(rows['faceIDs'], int); np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(vector(detail['selectedBandFaceIDs'], int), expected)
    assert not detail['constraintRankUsedToSelectFaces'] and detail['originalRowsExactlyRetained']
    assert detail['allPredictedViolatingInterfacesRepresented'] and detail['allPredictedViolatingCellsRepresented']
    np.testing.assert_array_equal(vector(rows['ownerCV'], int), minus[expected-1])
    np.testing.assert_array_equal(vector(rows['neighbourCV'], int), plus[expected-1])
    np.testing.assert_array_equal(np.asarray(rows['QL']).reshape(-1, 5), Q[minus[expected-1]-1])
    np.testing.assert_array_equal(np.asarray(rows['QR']).reshape(-1, 5), Q[plus[expected-1]-1])
    np.testing.assert_array_equal(np.asarray(rows['normal']).reshape(-1, 3), np.asarray(f['normal']).reshape(-1, 3)[expected-1])
    signs = {int(i): (-1 if sample[np.searchsorted(ids, i), 4]-low <= 10 else 1) for i in original}
    for i, cell in zip(required, bad):
        sign = -1 if states[cell-1, 4] < low else 1
        assert i not in signs or signs[i] == sign; signs[int(i)] = sign
    for i, temperature in zip(preview['violatingFaces'], preview['publicViolatingTemperatureK']):
        sign = -1 if temperature < low else 1
        assert i not in signs or signs[i] == sign; signs[int(i)] = sign
    np.testing.assert_array_equal(vector(rows['temperatureBoundSigns']), [signs[int(i)] for i in expected])
    n = len(Q); count = len(actual); column = np.arange(5)*n
    left = np.asarray(rows['leftRowsK']).reshape(-1, 5); right = np.asarray(rows['rightRowsK']).reshape(-1, 5)
    r = np.repeat(np.arange(count), 5)
    sparse = coo_matrix((np.r_[left.ravel(), right.ravel()],
                       (np.r_[r, r], np.r_[(minus[actual-1, None]-1+column).ravel(),
                                          (plus[actual-1, None]-1+column).ravel()])), shape=(count, 5*n)).tocsr()
    assert (rows['H']-sparse).nnz == 0
    return dict(preview=preview, originalRows=len(original), predictedViolatingCells=len(bad),
                selectedRows=count, currentOriginalCVSparseMappingExactlyVerified=True)


def contact_boundary(pair):
    """零接触速度是分支边界，显式展示分支切换，不伪称通过同波系中心差分。"""
    Q = np.asarray(pair['Q']); rows = pair['rows']; gas = pair['gas']; normal = vector(rows['normal'])
    scale = vector(rows['scale']); changes = []
    for side in range(2):
        for field in range(5):
            h = 1e-7*max(1., abs(Q[side, field]/scale[field])); plus = Q.copy(); minus = Q.copy()
            plus[side, field] += h*scale[field]; minus[side, field] -= h*scale[field]
            _, _, a = exact(primitive(plus[0], gas), primitive(plus[1], gas), normal, gas)
            _, _, b = exact(primitive(minus[0], gas), primitive(minus[1], gas), normal, gas)
            if a != b: changes.append(dict(side=side+1, field=field+1, plusRegion=a, minusRegion=b))
    assert changes, '必须实际观察到接触分支切换，不能用猜测豁免导数核对。'
    return dict(status='explicit-zero-contact-speed-nonsmooth-branch-boundary', switchingColumns=changes,
                sameSelectedRegionCenteredDerivativeClaimed=False)


def main(folder):
    folder = Path(folder).resolve(); source = folder/'cell_interface_lookahead_unit.mat'; before = sha(source)
    data = loadmat(source, simplify_cells=True); pairs = list(data['pairs']); checks = []
    assert len(pairs) == 4 and data['identifier'] == 'CutExactEuler:Vacuum'
    for i, pair in enumerate(pairs):
        result = selected_reference(np.asarray(pair['Q']), np.asarray(pair['candidate']), pair['info'], pair['cut'],
                                    pair['gas'], pair['rows'], pair['detail'])
        if i == 2:
            result['derivative'] = contact_boundary(pair)
        else:
            result['derivative'] = [row_reference(pair['rows'], h) for h in (3e-7, 1e-7)]
        assert not np.size(pair['oldRows']['faceIDs'])
        checks.append(result)
    domains = []
    for case in data['domains']:
        c = case['source']; rows = case['rows']; old = case['oldRows']
        result = selected_reference(c['Q'], case['candidate'], c['info'], c['cut'], c['gas'], rows, case['detail'])
        positions = np.searchsorted(vector(rows['faceIDs'], int), vector(old['faceIDs'], int))
        assert (old['H']-rows['H'][positions]).nnz == 0
        np.testing.assert_array_equal(vector(old['temperatureBoundSigns']), vector(rows['temperatureBoundSigns'])[positions])
        result['publicTemperatureDerivatives'] = [row_reference(rows, h) for h in (3e-7, 1e-7)]
        domains.append(result)
    batch = data['batchCase']
    batched = [preview_reference(batch['candidate'], batch['cut'], batch['gas'], report) for report in batch['reports']]
    assert batched[0] == batched[1] == batched[2]
    assert sha(source) == before
    inputs = {**PUBLIC_BINDINGS, source.relative_to(ROOT).as_posix(): before}
    for path in ('tests/auditRocketTemperatureProjection.py', 'tests/riemannReservoirReference.py'):
        inputs[path] = sha(ROOT/path)
    result = dict(status='passed-independent-public-cell-interface-preview-union-blockage-H-and-same-branch-derivatives',
                  pairs=checks, manufacturedDomains=domains, batchSizes=[1, 3, 4096], batchReference=batched[0],
                  zeroContactBoundaryExplicitlyClassified=True, derivativeThresholdUnchanged=1e-7,
                  inputSHA256=inputs, scriptSHA256=sha(__file__), newActualNSEvaluations=0,
                  newFlowAcceptances=0, physicalSolution=False, newProductionSolvedAreaM2=0)
    target = folder/'independent_cell_interface_lookahead.json'; assert not target.exists()
    target.write_bytes((json.dumps(result, indent=2)+'\n').encode('utf8'))
    print(json.dumps(dict(status=result['status'], pairs=len(checks), domains=len(domains), batchFaces=batched[0]['interiorFaces'])))


if __name__ == '__main__':
    main(sys.argv[1])
