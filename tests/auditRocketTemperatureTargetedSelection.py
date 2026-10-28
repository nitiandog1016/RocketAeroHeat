"""独立核对目标驱动集合、原条件保持、元数据身份和公开Riemann导数。"""
from pathlib import Path
import json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT, sha
from auditRocketTemperatureProjection import row_reference


def vector(value):
    return np.atleast_1d(value).ravel()


def check(Q, info, cut, gas, scale, direction, rows, detail):
    """从当前候选守恒量和原采样独立构造预期集合，不使用投影秩或NS结果。"""
    alpha = detail['targetAlpha']; candidate = Q + alpha * direction * scale
    rho = candidate[:, 0]; u = candidate[:, 1:4] / rho[:, None]
    T = (candidate[:, 4] - .5 * rho * np.sum(u*u, axis=1)) / (rho * gas['R'] / (gas['gamma'] - 1))
    low, high = gas['temperatureRange']; bad = np.flatnonzero((T < low) | (T > high))
    np.testing.assert_array_equal(bad+1, vector(detail['predictedViolatingCellIDs']))
    audit = info['exactEulerFluxAudit']; ids = vector(audit['faceID']).astype(int)
    samples = np.asarray(audit['samplePrimitive']); margin = np.minimum(samples[:, 4]-low, high-samples[:, 4])
    original = ids[margin <= 10]; required = vector(detail['requiredOriginalStateFaceIDs']).astype(int)
    slots = np.searchsorted(ids, required); np.testing.assert_array_equal(ids[slots], required)
    region = vector(audit['sampleRegion'])[slots]; assert set(region.tolist()) <= {1, 4}
    f = cut['openFaces']; minus = vector(f['cellMinus']).astype(int); plus = vector(f['cellPlus']).astype(int)
    cv = np.where(region == 1, minus[required-1], plus[required-1])
    np.testing.assert_array_equal(cv, bad+1)
    expected = np.union1d(original, required); actual = vector(rows['faceIDs']).astype(int)
    np.testing.assert_array_equal(actual, expected)
    K = max(10., float(np.max(margin[slots])))
    np.testing.assert_array_equal(vector(detail['candidateExtendedBandFaceIDs']), ids[margin <= K])
    assert K == detail['selectedTemperatureBandK']
    np.testing.assert_array_equal(vector(detail['originalBandFaceIDs']), original)
    np.testing.assert_array_equal(vector(rows['ownerCV']), minus[actual-1])
    np.testing.assert_array_equal(vector(rows['neighbourCV']), plus[actual-1])
    np.testing.assert_array_equal(np.atleast_2d(rows['QL']), Q[minus[actual-1]-1])
    np.testing.assert_array_equal(np.atleast_2d(rows['QR']), Q[plus[actual-1]-1])
    np.testing.assert_array_equal(np.atleast_2d(rows['normal']), np.asarray(f['normal'])[actual-1])
    selected_slots = np.searchsorted(ids, actual)
    np.testing.assert_allclose(vector(rows['baseTemperatureK']), samples[selected_slots, 4], rtol=0, atol=1e-10)
    expected_sign = np.where(samples[selected_slots, 4]-low <= K, -1, 1)
    np.testing.assert_array_equal(vector(rows['temperatureBoundSigns']), expected_sign)
    assert detail['originalRowsExactlyRetained'] and detail['allPredictedViolatingCellsRepresented']
    assert not detail['constraintRankUsedToSelectFaces'] and not detail['allExtendedBandFacesSelected']
    assert not detail['nonlinearTemperatureBoundGuaranteed'] and not detail['physicalBoundaryConditionAdded']
    public = [row_reference(rows, h) for h in (3e-7, 1e-7)]
    return dict(predictedViolatingCells=len(bad), originalRows=len(original), requiredOriginalStateFaces=len(required),
                targetedRows=len(actual), candidateBandRows=int(np.sum(margin <= K)),
                allRowMetadataMappedToCurrentCVAndFace=True, publicTemperatureDerivatives=public)


def main(folder):
    folder = Path(folder).resolve(); file = folder/'temperature_targeted_unit.mat'; digest = sha(file)
    data = loadmat(file, simplify_cells=True); checks = []
    for case in data['cases']:
        s = case['source']
        checks.append(check(s['Q'], s['info'], s['cut'], s['gas'], s['scale'], s['value'], case['rows'], case['detail']))
        old = s['base']; new = case['rows']; positions = np.searchsorted(vector(new['faceIDs']), vector(old['faceIDs']))
        assert (old['H'] - new['H'][positions]).nnz == 0
        np.testing.assert_array_equal(vector(old['temperatureBoundSigns']), vector(new['temperatureBoundSigns'])[positions])
        assert case['legacyDefaultExactlyReplayed']
    state = ROOT/'outputs/v83_validation/krylov37/wall_velocity_krylov_evidence.mat'
    prefix = ROOT/'outputs/v90_validation/continuous01/step01/linear_checkpoints/lookahead_free_prefix.mat'
    geometry = ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat'
    bindings = {str(p.relative_to(ROOT)): sha(p) for p in (state, prefix, geometry, file)}
    s = loadmat(state, simplify_cells=True, variable_names=['result', 'gas'])
    f = loadmat(prefix, simplify_cells=True, variable_names=['free', 'originalRows'])
    cut = loadmat(geometry, simplify_cells=True, variable_names=['domain'])['domain']
    real = check(s['result']['Q'], s['result']['info'], cut, s['gas'], s['result']['residualScale'], f['free'], data['realRows'], data['realReport'])
    old = f['originalRows']; new = data['realRows']; positions = np.searchsorted(vector(new['faceIDs']), vector(old['faceIDs']))
    assert (old['H'] - new['H'][positions]).nnz == 0
    np.testing.assert_array_equal(vector(old['temperatureBoundSigns']), vector(new['temperatureBoundSigns'])[positions])
    assert (real['predictedViolatingCells'], real['originalRows'], real['targetedRows'], real['candidateBandRows']) == (25, 124, 131, 193)
    assert digest == sha(file) and all(sha(ROOT/p) == h for p, h in bindings.items())
    out = dict(status='passed-independent-targeted-union-all-predicted-CVs-original-H-signs-current-metadata-and-public-derivatives',
               manufacturedCases=checks, currentQ37=real, inputSHA256=bindings, scriptSHA256=sha(__file__),
               derivativeReferenceSHA256=sha(ROOT/'tests/auditRocketTemperatureProjection.py'),
               newActualNSEvaluations=0, newFlowAcceptances=0, physicalSolution=False, newProductionSolvedAreaM2=0)
    target = folder/'independent_targeted_selection.json'; assert not target.exists()
    target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8'); print(json.dumps(dict(status=out['status'], currentQ37=real)))


if __name__ == '__main__':
    main(sys.argv[1])
