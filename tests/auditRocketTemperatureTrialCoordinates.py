"""独立核对候选坐标、全部小域界面与逐轮稠密约束方向；不授予物理求解资格。"""
from pathlib import Path
import json
import os
import sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from scipy.linalg import null_space
from riemannReservoirReference import ROOT, sha, PUBLIC_BINDINGS
from rocketTrialCoordinatesReference import trial_candidate
from auditRocketCellInterfaceLookahead import selected_reference, vector
from auditRocketTemperatureProjection import row_reference
from auditRocketTemperatureRepreview import cells, subset, repreview_arithmetic


def main(folder):
    folder = Path(folder).resolve()
    file = folder/'temperature_trial_coordinates_unit.mat'
    original = sha(file)
    data = loadmat(file, simplify_cells=True)
    pairs = []
    for pair in cells(data['pairs']):
        candidate = trial_candidate(pair['Q'], pair['value'], pair['scale'], .125,
                                    pair['gas'], pair['coordinates'])
        np.testing.assert_allclose(candidate, pair['candidate'], rtol=3e-15, atol=0)
        detail = pair['detail']
        assert detail.get('trialStateCoordinates', 'conservative') == pair['coordinates']
        check = selected_reference(pair['Q'], candidate, pair['info'], pair['cut'],
                                   pair['gas'], pair['rows'], detail)
        check['coordinates'] = pair['coordinates']
        check['publicTemperatureDerivatives'] = [row_reference(pair['rows'], h) for h in (3e-7, 1e-7)]
        pairs.append(check)
    results = []
    for case in cells(data['cases']):
        evidence = case['evidence']; report = case['detail']; n = len(case['Q'])
        assert report['repreview']['trialStateCoordinates'] == 'primitive-tangent'
        arithmetic = repreview_arithmetic(evidence['rows'], evidence, report)
        A = np.row_stack([case['C'], case['factor']['penaltyInOriginalResidualUnits']*case['W']])
        rhs = np.r_[case['rhs'].ravel(order='F'), np.zeros(5*n)]
        free = np.linalg.lstsq(A, rhs, rcond=None)[0]
        free_error = float(np.linalg.norm(evidence['free'].ravel(order='F')-free)/max(1., np.linalg.norm(free)))
        assert free_error < 1e-8
        cycles = cells(evidence['cycles'])
        directions = [evidence['free']] + [c['value'] for c in cycles]
        previews = []
        for value, selection in zip(directions, cells(evidence['previewSelections'])):
            assert selection['trialStateCoordinates'] == 'primitive-tangent'
            candidate = trial_candidate(case['Q'], value, case['scale'],
                                        report['repreview']['targetAlpha'], case['gas'], 'primitive-tangent')
            selected = subset(evidence['rows'], vector(selection['selectedBandFaceIDs'], int))
            previews.append(selected_reference(case['Q'], candidate, case['info'], case['cut'],
                                               case['gas'], selected, selection))
        errors = []
        for cycle in cycles:
            active = vector(cycle['activeRowIndices'], int)-1
            Z = null_space(evidence['rows']['H'][active].toarray())
            truth = Z @ np.linalg.lstsq(A @ Z, rhs, rcond=None)[0]
            error = float(np.linalg.norm(cycle['value'].ravel(order='F')-truth)/max(1., np.linalg.norm(truth)))
            assert error < 1e-8
            objective = float(np.linalg.norm(A @ cycle['value'].ravel(order='F')-rhs)**2/np.linalg.norm(case['rhs'])**2)
            np.testing.assert_allclose(objective, cycle['detail']['normalizedAugmentedObjective'], rtol=1e-9)
            errors.append(error)
        results.append(dict(arithmetic=arithmetic, freeDenseRelativeDifference=free_error,
                            projectedDenseRelativeDifferences=errors, previews=previews,
                            publicTemperatureDerivatives=[row_reference(evidence['rows'], h) for h in (3e-7, 1e-7)]))
    empty = data['actualEmpty']; wider = data['wider']
    assert not np.size(empty['faceIDs']) and empty['H'].shape[0] == 0
    np.testing.assert_array_equal(empty['scale'], wider['scale'])
    assert set(empty['gas']) == set(wider['gas']) and np.all(data['emptyRejects'])
    for key in empty['gas']:
        np.testing.assert_array_equal(empty['gas'][key], wider['gas'][key])
    assert data['report']['unknownCoordinatesRejected'] and data['report']['actualEmptySelectorMergeAndIdentityPassed']
    assert sha(file) == original
    result = dict(status='passed-independent-coordinate-path-all-small-interfaces-and-dense-repreview-audit',
                  pairs=pairs, cases=results, actualEmptySelectorIdentityVerified=True,
                  inputSHA256={file.relative_to(ROOT).as_posix(): original}, scriptSHA256=sha(__file__),
                  publicReferenceSHA256=PUBLIC_BINDINGS, newActualNSEvaluations=0, newFlowAcceptances=0,
                  manufacturedFacesNotClosedNSGeometry=True, manufacturedCNotPhysicalNSJacobian=True,
                  physicalSolution=False, newProductionSolvedAreaM2=0)
    target = folder/'independent_trial_coordinates.json'
    assert not target.exists()
    target.write_bytes((json.dumps(result, indent=2)+'\n').encode('utf8'))
    print(json.dumps(dict(status=result['status'], cases=len(results),
                         maximumDenseDifference=max(max(c['projectedDenseRelativeDifferences'], default=0) for c in results))))


if __name__ == '__main__':
    main(sys.argv[1])
