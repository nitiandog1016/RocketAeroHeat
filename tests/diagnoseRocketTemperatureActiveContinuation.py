"""核对连续各态的全部内部采样温域，并用公开Riemann参考复算最热和最冷面。"""
from pathlib import Path
import json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from riemannReservoirReference import ROOT, sha, exact, PUBLIC_BINDINGS

OUT = ROOT / 'outputs/v87_validation'
if len(sys.argv) > 1:
    supplied = Path(sys.argv[1])
    OUT = (supplied if supplied.is_absolute() else ROOT / supplied).resolve()
    OUT.relative_to(ROOT.resolve())


def main():
    study_file = OUT / 'continuous01/active_continuation.json'
    study = json.loads(study_file.read_text('utf8'))
    geometry = ROOT / 'outputs/v80_validation/domain02/shell_domain_evidence.mat'
    cut = loadmat(geometry, simplify_cells=True, variable_names=['domain'])['domain']
    faces = cut['openFaces']; source_faces = len(cut['inputGeometry']['F'])
    bindings = {str(p.relative_to(ROOT)): sha(p) for p in (study_file, geometry)}
    bindings.update(PUBLIC_BINDINGS)
    steps = study['steps'] if isinstance(study['steps'], list) else [study['steps']]
    files = [Path(study['initialStateFile'])] + [Path(s['resultFolder']) / 'wall_velocity_krylov_evidence.mat'
        for s in steps if s['main'] and s['main']['acceptedNonlinearIterations']]
    results = []
    expected_ids = np.flatnonzero((faces['area'] > 0) & (faces['cellMinus'] > 0) & (faces['cellPlus'] > 0))
    for state_number, file in enumerate(files):
        bindings[str(file.relative_to(ROOT))] = sha(file)
        saved = loadmat(file, simplify_cells=True, variable_names=['result', 'gas'])
        result = saved['result']; gas = saved['gas']; audit = result['info']['exactEulerFluxAudit']
        ids = np.asarray(audit['faceID'], int) - 1
        np.testing.assert_array_equal(ids, expected_ids)
        T = np.asarray(audit['samplePrimitive'])[:, 4]
        low, high = map(float, gas['temperatureRange'])
        assert np.all(np.isfinite(T)) and np.all(T >= low) and np.all(T <= high)
        lower, upper = T - low, high - T
        sample_slots = np.unique(np.r_[np.argsort(upper)[:8], np.argsort(lower)[:8]])
        Q = result['Q']; cv = gas['R'] / (gas['gamma'] - 1)
        def primitive(q):
            rho = q[0]; u = q[1:4] / rho
            return np.r_[rho, u, (q[4] - .5 * rho * (u @ u)) / (rho * cv)]
        selected = []; maximum_error = 0.
        for slot in sample_slots:
            face = ids[slot]; owner = int(faces['cellMinus'][face]) - 1
            neighbour = int(faces['cellPlus'][face]) - 1
            sample, _, _ = exact(primitive(Q[owner]), primitive(Q[neighbour]), faces['normal'][face], gas)
            error = abs(sample[4] - T[slot]); maximum_error = max(error, maximum_error)
            assert error <= 1e-8 * max(1., abs(T[slot]))
            selected.append(dict(faceID=int(face + 1), temperatureK=float(T[slot]),
                lowerMarginK=float(lower[slot]), upperMarginK=float(upper[slot]),
                ownerSourceFacet=int(owner % source_faces + 1), ownerLayer=int(owner // source_faces + 1),
                neighbourSourceFacet=int(neighbour % source_faces + 1), neighbourLayer=int(neighbour // source_faces + 1),
                centroidM=faces['centroid'][face].tolist()))
        if state_number < len(steps):
            direction = steps[state_number]['directionReport']
            band_ids = np.asarray(direction['bandFaceIDs'], int).ravel()
            # 原十K集合仍须逐面一致；显式前瞻时再核对实际集合与全部所需面的并集。
            original_ids = ids[(upper <= 10.) | (lower <= 10.)] + 1
            selection = direction.get('lookaheadSelection')
            if selection:
                np.testing.assert_array_equal(np.atleast_1d(selection['originalBandFaceIDs']), original_ids)
                if selection.get('selectionPolicy') == 'original-and-predicted-required-faces':
                    expected = np.union1d(original_ids, np.atleast_1d(selection['requiredOriginalStateFaceIDs']))
                elif selection.get('selectionPolicy') in ('original-and-predicted-cell-and-interface-faces',
                                                          'original-and-iteratively-predicted-cell-and-interface-faces'):
                    expected = np.union1d(np.union1d(original_ids, np.atleast_1d(selection['requiredOriginalStateFaceIDs'])),
                                          np.atleast_1d(selection['predictedViolatingInterfaceFaceIDs']))
                    previews = [selection['interfacePreview']]
                    if 'repreview' in direction:
                        direction_path = Path(steps[state_number]['directionFile'])
                        bindings[str(direction_path.relative_to(ROOT))] = sha(direction_path)
                        saved_direction = loadmat(direction_path, simplify_cells=True,
                                                  variable_names=['evidence'])['evidence']
                        sequence = saved_direction['previewSelections']
                        sequence = [sequence] if isinstance(sequence, dict) else list(sequence)
                        previews = [p['interfacePreview'] for p in sequence if p['previewCellStateEvaluable']]
                        assert len(sequence) == direction['repreview']['previewCount']
                    for preview in previews:
                        assessed = np.atleast_1d(preview['assessedFaceIDs']).astype(int)
                        blocked = np.atleast_1d(preview['blockedFaceIDs']).astype(int)
                        assert not np.intersect1d(assessed, blocked).size
                        np.testing.assert_array_equal(np.sort(np.r_[assessed, blocked]), ids+1)
                        assert preview['allInputFacesAccountedFor'] and not preview['temperatureDomainRelaxed']
                        assert bool(preview['allInteriorFacesAssessed']) == (not len(blocked))
                else:
                    assert selection.get('selectionPolicy', 'all-extended-band') == 'all-extended-band'
                    width = selection['selectedTemperatureBandK']
                    expected = ids[(upper <= width) | (lower <= width)]+1
                np.testing.assert_array_equal(band_ids, expected)
                np.testing.assert_array_equal(band_ids, np.atleast_1d(selection['selectedBandFaceIDs']))
            else:
                np.testing.assert_array_equal(band_ids, original_ids)
        result_row = dict(state=state_number, sourceFile=str(file), interiorFaceCount=len(ids),
            interfaceTemperatureRangeK=[float(T.min()), float(T.max())],
            upperMarginCounts={str(b): int(np.sum(upper <= b)) for b in (.1, 1., 5., 10., 50.)},
            lowerMarginCounts={str(b): int(np.sum(lower <= b)) for b in (.1, 1., 5., 10., 50.)},
            actualCellTemperatureRangeK=[float(result['info']['primitive'][:, 4].min()), float(result['info']['primitive'][:, 4].max())],
            publicReferenceSampleFaces=len(sample_slots), publicReferenceMaximumTemperatureDifferenceK=maximum_error,
            extremeFaces=selected, prospectiveDirectionOrStepRadiusClaimed=False)
        results.append(result_row)
    assert all(sha(ROOT / name) == digest for name, digest in bindings.items())
    out = dict(status='passed-all-saved-current-temperature-margins-and-public-extreme-face-samples',
        states=results, allInteriorSamplesIndependentlyRecomputed=False,
        newActualNSEvaluations=0, newFlowAcceptances=0, physicalSolution=False,
        newProductionSolvedAreaM2=0, inputSHA256=bindings, scriptSHA256=sha(__file__))
    target = OUT / 'current_temperature_margins.json'; assert not target.exists()
    target.write_text(json.dumps(out, indent=2) + '\n', encoding='utf8')
    print(json.dumps(dict(status=out['status'], states=len(results), finalRangeK=results[-1]['interfaceTemperatureRangeK'])))


if __name__ == '__main__':
    main()
