"""独立核对逐轮集合、强制活跃、全部稠密约束方向和公开温度导数，不授予物理资格。"""
from pathlib import Path
import json, os, sys
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.io import loadmat
from scipy.linalg import null_space
from riemannReservoirReference import ROOT, sha, PUBLIC_BINDINGS
from auditRocketCellInterfaceLookahead import selected_reference, vector
from auditRocketTemperatureProjection import row_reference


def cells(x):
    return [x] if isinstance(x, dict) else list(np.asarray(x, dtype=object).ravel())


def subset(rows, ids):
    positions = np.searchsorted(vector(rows['faceIDs'], int), ids); out = dict(rows)
    keys = ('faceIDs', 'ownerCV', 'neighbourCV', 'normal', 'QL', 'QR', 'leftRowsK', 'rightRowsK',
            'baseSamplePrimitive', 'baseTemperatureK', 'remainingUpperTemperatureK', 'sampleRegion', 'rowTwoNormK', 'temperatureBoundSigns')
    for key in keys:
        a = np.asarray(rows[key]); out[key] = np.atleast_1d(a)[positions] if a.ndim <= 1 else a[positions]
    out['H'] = rows['H'][positions]
    return out


def repreview_arithmetic(rows, evidence, report):
    """每轮只使用当时可见的行，独立重建强制面与朝外行并集，不把未来行倒灌到旧循环。"""
    ids = vector(rows['faceIDs'], int); H = rows['H']; signs = vector(rows['temperatureBoundSigns'])
    history = cells(report['repreview']['history']); selections = cells(evidence['previewSelections']); cycles = cells(evidence['cycles'])
    assert len(history) == len(selections) == len(cycles)+1
    available = vector(report['lookaheadSelection']['originalBandFaceIDs'], int)
    forced = np.array([], int); active = np.array([], int); value = evidence['free']; checks = []
    iterations = int(evidence['freeDetail']['nativeIterations'])
    for i, (selection, summary) in enumerate(zip(selections, history)):
        previous = available.copy(); available = np.union1d(available, vector(selection['selectedBandFaceIDs'], int))
        required = vector(selection['requiredOriginalStateFaceIDs'], int)
        interfaces = vector(selection.get('predictedViolatingInterfaceFaceIDs', []), int)
        forced = np.union1d(forced, np.union1d(required, interfaces))
        positions = np.searchsorted(ids, available); t = np.atleast_1d(H[positions] @ value.ravel(order='F'))
        b = 128*np.finfo(float).eps*np.atleast_1d(abs(H[positions]) @ abs(value.ravel(order='F')))
        outward = available[signs[positions]*t > b]
        needed = np.union1d(np.union1d(active, forced), outward); new = np.setdiff1d(needed, active)
        for key, expected in [('selectedFaceIDs', available), ('addedFaceIDs', np.setdiff1d(available, previous)),
                              ('forcedActiveFaceIDs', forced), ('outwardFaceIDs', outward),
                              ('activeFacesBeforeSolve', active), ('newActiveFaceIDs', new)]:
            np.testing.assert_array_equal(vector(summary[key], int), expected)
        ip = selection.get('interfacePreview'); fully = isinstance(ip, dict) and bool(ip['allInteriorFacesAssessed'])
        feasible = bool(selection['previewCellStateEvaluable']) and not np.size(selection['predictedViolatingCellIDs']) and fully and not len(interfaces)
        assert bool(summary['targetCandidateThermallyFeasible']) == feasible
        if i < len(cycles):
            cycle = cycles[i]; np.testing.assert_array_equal(vector(cycle['activeFaceIDs'], int), needed)
            np.testing.assert_array_equal(vector(cycle['availableFaceIDs'], int), available)
            np.testing.assert_array_equal(ids[vector(cycle['activeRowIndices'], int)-1], needed)
            np.testing.assert_array_equal(ids[vector(cycle['availableRowIndices'], int)-1], available)
            np.testing.assert_array_equal(ids[vector(cycle['newRowIndices'], int)-1], new)
            value = cycle['value']; iterations += int(cycle['detail']['nativeIterations']); active = needed
            after = np.atleast_1d(H[positions] @ value.ravel(order='F'))
            np.testing.assert_allclose(after, vector(cycle['temperatureDerivativeK']), atol=1e-8, rtol=1e-8)
            after_bound = 128*np.finfo(float).eps*np.atleast_1d(abs(H[positions]) @ abs(value.ravel(order='F')))
            np.testing.assert_allclose(after_bound, vector(cycle['localDotRoundoffBoundK']), rtol=1e-12, atol=0)
            np.testing.assert_array_equal(vector(cycle['remainingViolatingFaceIDs'], int), available[signs[positions]*after > after_bound])
            projector = cycle['projector']
            if 'originalToUniqueRow' in projector:
                representatives = vector(projector['exactRepresentativeRowIndices'], int)-1
                mapping = vector(projector['originalToUniqueRow'], int)-1
                selected_H = H[np.searchsorted(ids, needed)]
                assert (selected_H-selected_H[representatives[mapping]]).nnz == 0
                assert projector['allOriginalRowsExactlyRepresented'] and not projector['rankTruncationUsed']
        checks.append(dict(round=i, selectedRows=len(available), newActiveFaceIDs=new.tolist(),
                           targetCandidateThermallyFeasible=feasible))
    np.testing.assert_array_equal(available, ids)
    last = np.atleast_1d(H @ value.ravel(order='F')); bound = 128*np.finfo(float).eps*np.atleast_1d(abs(H) @ abs(value.ravel(order='F')))
    gate = bool(np.all(signs*last <= bound) and np.all(np.isin(forced, active)))
    assert bool(report['bandDirectionFeasible']) == gate and int(report['totalNativeIterations']) == iterations
    np.testing.assert_array_equal(vector(report['activeFaceIDs'], int), active)
    np.testing.assert_allclose(vector(report['temperatureDerivativeK']), last, rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(vector(report['localDotRoundoffBoundK']), bound, rtol=1e-12, atol=0)
    np.testing.assert_array_equal(vector(report['remainingViolatingRowIndices'], int), np.flatnonzero(signs*last > bound)+1)
    assert not report['oldActiveSetInherited'] and not report['constrainedOptimumClaimed'] and not report['physicalBoundaryConditionAdded']
    assert bool(report['repreview']['forcedRequiredRowsActuallySolved']) == bool(np.all(np.isin(forced, active)))
    assert bool(report['repreview']['targetCandidateThermallyFeasible']) == checks[-1]['targetCandidateThermallyFeasible']
    assert report['activeSolveCount'] <= report['options']['maximumActiveSolves'] <= 3 and len(history) <= 4
    return dict(rounds=checks, actualProjectedSolves=len(cycles), totalNativeIterations=iterations,
                original128EpsDotGatePassed=gate, targetThermalCheckPassed=checks[-1]['targetCandidateThermallyFeasible'])


def unit(folder):
    file = folder/'temperature_repreview_unit.mat'; before = sha(file); data = loadmat(file, simplify_cells=True); results = []
    for case in cells(data['cases']):
        variants = [('closing', case['Q'], case['info'], case['rhs'], case['evidence'], case['detail']),
                    ('one-solve-budget', case['Q'], case['info'], case['rhs'], case['budgetEvidence'], case['budget']),
                    ('curvature-stall', case['highQ'], case['highInfo'], case['stallRhs'], case['stallEvidence'], case['stall']),
                    ('high-kinetic-dot-gate-stress', case['highQ'], case['highInfo'], case['stressRhs'], case['stressEvidence'], case['stress'])]
        A = np.row_stack([case['C'], case['factor']['penaltyInOriginalResidualUnits']*case['W']]); n = len(case['Q'])
        for name, Q, info, rhs, evidence, report in variants:
            arithmetic = repreview_arithmetic(evidence['rows'], evidence, report); dense = []; previews = []
            directions = [evidence['free']] + [c['value'] for c in cells(evidence['cycles'])]
            for value, selection in zip(directions, cells(evidence['previewSelections'])):
                candidate = Q + report['repreview']['targetAlpha']*value*case['scale']
                selected = subset(evidence['rows'], vector(selection['selectedBandFaceIDs'], int))
                previews.append(selected_reference(Q, candidate, info, case['cut'], case['gas'], selected, selection))
            for cycle in cells(evidence['cycles']):
                active = vector(cycle['activeRowIndices'], int)-1; Z = null_space(evidence['rows']['H'][active].toarray())
                reference = Z @ np.linalg.lstsq(A @ Z, np.r_[rhs.ravel(order='F'), np.zeros(5*n)], rcond=None)[0]
                difference = float(np.linalg.norm(cycle['value'].ravel(order='F')-reference)/max(1., np.linalg.norm(reference)))
                assert difference < 1e-8
                objective = float(np.linalg.norm(A @ cycle['value'].ravel(order='F')-np.r_[rhs.ravel(order='F'), np.zeros(5*n)])**2/np.linalg.norm(rhs)**2)
                np.testing.assert_allclose(objective, cycle['detail']['normalizedAugmentedObjective'], rtol=1e-9)
                dense.append(difference)
            public = [row_reference(evidence['rows'], h) for h in (3e-7, 1e-7)]
            results.append(dict(variant=name, arithmetic=arithmetic, denseDirectionDifferences=dense, previews=previews,
                                publicTemperatureDerivatives=public, physicalSolution=False))
        assert case['sideConflictRejected'] and case['zeroReport']['activeSolveCount'] == 0
    assert sha(file) == before
    return dict(cases=results, inputSHA256={file.relative_to(ROOT).as_posix(): before}, manufacturedLinearSystemsNotPhysicalNSJacobian=True)


def main(folder):
    folder = Path(folder).resolve(); result = unit(folder)
    result.update(status='passed-independent-repreview-all-round-unions-forced-sets-dense-directions-public-temperatures-and-gates',
                  scriptSHA256=sha(__file__), publicReferenceSHA256=PUBLIC_BINDINGS,
                  newActualNSEvaluations=0, newFlowAcceptances=0, physicalSolution=False, newProductionSolvedAreaM2=0)
    target = folder/'independent_temperature_repreview.json'; assert not target.exists()
    target.write_bytes((json.dumps(result, indent=2)+'\n').encode())
    print(json.dumps(dict(status=result['status'], variants=len(result['cases']), maximumDenseDifference=max(max(c['denseDirectionDifferences'], default=0) for c in result['cases']))))


if __name__ == '__main__':
    main(sys.argv[1])
