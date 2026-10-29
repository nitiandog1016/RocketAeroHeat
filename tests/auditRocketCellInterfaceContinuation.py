"""逐轮复算真实候选全部CV/界面输入分区和必要面并集；公开解只声明实际抽样范围。"""
from pathlib import Path
import json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from riemannReservoirReference import ROOT, sha, exact, PUBLIC_BINDINGS
from auditRocketCellInterfaceLookahead import primitive, vector
from auditRocketTemperatureRepreview import cells, repreview_arithmetic


def audit_preview(Q, scale, gas, cut, value, selection, rows, original):
    """不借用未来集合证明本轮选面；全部CV温度独立复算，界面原解只抽样。"""
    low, high = gas['temperatureRange']; candidate = Q+selection['targetAlpha']*value*scale
    rho = candidate[:, 0]; positive_density = bool(np.all(rho > 0))
    internal = candidate[:, 4]-.5*np.sum(candidate[:, 1:4]**2, axis=1)/rho if positive_density else np.array([-1.])
    evaluable = positive_density and bool(np.all(internal > 0))
    assert bool(selection['previewCellStateEvaluable']) == evaluable
    if not evaluable:
        assert not isinstance(selection['interfacePreview'], dict)
        return dict(previewCellStateEvaluable=False, internalThermalCheckPassed=False,
                    publicSamples=[], allAssessableRiemannStatesIndependentlyRecomputed=False)
    T = internal/(rho*gas['R']/(gas['gamma']-1)); bad = np.flatnonzero((T < low) | (T > high))+1
    np.testing.assert_array_equal(bad, vector(selection['predictedViolatingCellIDs'], int))
    f = cut['openFaces']; minus = vector(f['cellMinus'], int); plus = vector(f['cellPlus'], int)
    ids = np.flatnonzero((f['area'] > 0) & (minus > 0) & (plus > 0))
    valid = (T >= low) & (T <= high); eligible = valid[minus[ids]-1] & valid[plus[ids]-1]
    assessed = ids[eligible]+1; blocked = ids[~eligible]+1; preview = selection['interfacePreview']
    np.testing.assert_array_equal(assessed, vector(preview['assessedFaceIDs'], int))
    np.testing.assert_array_equal(blocked, vector(preview['blockedFaceIDs'], int))
    np.testing.assert_array_equal(bad, vector(preview['blockedByOutOfRangeCellIDs'], int))
    assert bool(preview['allInteriorFacesAssessed']) == (not len(blocked)) and preview['allInputFacesAccountedFor']
    assert not preview['temperatureDomainRelaxed'] and not preview['nonTemperatureFailuresSuppressed']
    required = vector(selection['requiredOriginalStateFaceIDs'], int)
    interfaces = vector(selection.get('predictedViolatingInterfaceFaceIDs', []), int)
    expected = np.union1d(np.union1d(original, required), interfaces)
    np.testing.assert_array_equal(vector(selection['originalBandFaceIDs'], int), original)
    np.testing.assert_array_equal(vector(selection['selectedBandFaceIDs'], int), expected)
    selected = vector(rows['faceIDs'], int); positions = np.searchsorted(selected, required)
    np.testing.assert_array_equal(selected[positions], required)
    regions = vector(rows['sampleRegion'], int)[positions]; assert set(regions) <= {1, 4}
    np.testing.assert_array_equal(np.where(regions == 1, minus[required-1], plus[required-1]), bad)
    if len(bad):
        np.testing.assert_array_equal(vector(rows['temperatureBoundSigns'])[positions], np.where(T[bad-1] < low, -1, 1))
    # 每轮覆盖报告越温面、风险前128面、均匀64面及本轮选择面；不称独立重算全部大域Riemann解。
    margin = np.minimum(T-low, high-T)
    risk = np.minimum(margin[minus[assessed-1]-1], margin[plus[assessed-1]-1])
    uniform = assessed[np.linspace(0, len(assessed)-1, 64, dtype=int)] if len(assessed) else np.array([], int)
    samples = np.unique(np.r_[interfaces, assessed[np.argsort(risk, kind='stable')[:128]], uniform,
                              np.intersect1d(assessed, expected)]).astype(int)
    checked = []
    for face in samples:
        sample, _, region = exact(primitive(candidate[minus[face-1]-1], gas), primitive(candidate[plus[face-1]-1], gas), f['normal'][face-1], gas)
        violating = sample[4] < low or sample[4] > high
        assert violating == (face in interfaces)
        if violating:
            slot = np.flatnonzero(interfaces == face)[0]
            np.testing.assert_allclose(sample[4], vector(preview['predictedInterfaceTemperatureK'])[slot], rtol=1e-9)
        checked.append(dict(faceID=int(face), publicTemperatureK=float(sample[4]), region=region, violating=bool(violating)))
    return dict(previewCellStateEvaluable=True, cellCount=len(Q), cellTemperatureRangeK=[float(T.min()), float(T.max())],
                interiorFaceCount=len(ids), predictedViolatingCells=len(bad), assessedInterfaceCount=len(assessed),
                blockedInterfaceCount=len(blocked), localSelectedRows=len(expected), predictedViolatingInterfaces=len(interfaces),
                allCVTemperaturesAndInterfaceInputPartitionIndependentlyRecomputed=True,
                allAssessableRiemannStatesIndependentlyRecomputed=False, publicSampleCount=len(checked), publicSamples=checked,
                internalThermalCheckPassed=not len(bad) and not len(blocked) and not len(interfaces))


def main(folder):
    folder = Path(folder).resolve(); step = folder/'continuous01/step01'
    direction = step/'active_direction.mat'; prefix = step/'linear_checkpoints/lookahead_free_prefix.mat'
    context = loadmat(direction, simplify_cells=True, variable_names=['rows', 'sourceBinding', 'preparationFile', 'preparationSHA256', 'geometryBinding', 'directionReport', 'evidence'])
    p = json.loads(Path(context['preparationFile']).read_text('utf8')); assert sha(context['preparationFile']) == context['preparationSHA256']
    assert sha(context['sourceBinding']['sourceFile']) == context['sourceBinding']['sourceFileSHA256'] == p['sourceBinding']['sourceFileSHA256']
    assert sha(context['geometryBinding']['file']) == context['geometryBinding']['sha256']
    s = loadmat(p['artifacts']['cache']['file'], simplify_cells=True, variable_names=['Q', 'gas', 'scale'])
    cut = loadmat(context['geometryBinding']['file'], simplify_cells=True, variable_names=['domain'])['domain']
    saved = loadmat(prefix, simplify_cells=True); report = context['directionReport']; selection = report['lookaheadSelection']
    iterative = selection['selectionPolicy'] == 'original-and-iteratively-predicted-cell-and-interface-faces'
    assert iterative or selection['selectionPolicy'] == 'original-and-predicted-cell-and-interface-faces'
    original = vector(saved['originalRows']['faceIDs'], int); rows = context['rows']; selected = vector(rows['faceIDs'], int)
    positions = np.searchsorted(selected, original); np.testing.assert_array_equal(selected[positions], original)
    assert (saved['originalRows']['H']-rows['H'][positions]).nnz == 0
    np.testing.assert_array_equal(vector(saved['originalRows']['temperatureBoundSigns']), vector(rows['temperatureBoundSigns'])[positions])
    Q = s['Q']; f = cut['openFaces']; minus = vector(f['cellMinus'], int); plus = vector(f['cellPlus'], int)
    np.testing.assert_array_equal(vector(rows['ownerCV'], int), minus[selected-1])
    np.testing.assert_array_equal(vector(rows['neighbourCV'], int), plus[selected-1])
    np.testing.assert_array_equal(np.asarray(rows['QL']).reshape(-1, 5), Q[minus[selected-1]-1])
    np.testing.assert_array_equal(np.asarray(rows['QR']).reshape(-1, 5), Q[plus[selected-1]-1])
    np.testing.assert_array_equal(np.asarray(rows['normal']).reshape(-1, 3), f['normal'][selected-1])
    if iterative:
        evidence = context['evidence']; arithmetic = repreview_arithmetic(rows, evidence, report)
        selections = cells(evidence['previewSelections']); cycles = cells(evidence['cycles'])
        values = [evidence['free']]+[c['value'] for c in cycles]
        np.testing.assert_array_equal(saved['free'], evidence['free'])
    else:
        arithmetic = None; selections = [selection]; values = [saved['free']]
    checks = []; union = original.copy()
    for i, (value, preview) in enumerate(zip(values, selections)):
        result = audit_preview(Q, s['scale'], s['gas'], cut, value, preview, rows, original)
        union = np.union1d(union, vector(preview['selectedBandFaceIDs'], int))
        if iterative:
            h = cells(report['repreview']['history'])[i]
            np.testing.assert_array_equal(union, vector(h['selectedFaceIDs'], int))
            assert result['internalThermalCheckPassed'] == bool(h['targetCandidateThermallyFeasible'])
        result['round'] = i; checks.append(result)
    np.testing.assert_array_equal(union, selected)
    inputs = [direction, prefix, Path(context['preparationFile']), Path(p['artifacts']['cache']['file']), Path(context['sourceBinding']['sourceFile']), Path(context['geometryBinding']['file'])]
    bindings = {**PUBLIC_BINDINGS, **{path.relative_to(ROOT).as_posix(): sha(path) for path in inputs}}
    out = dict(status='passed-real-every-preview-all-CV-partitions-required-unions-original-H-and-public-interface-samples',
               iterativeRepreview=iterative, previews=checks, arithmetic=arithmetic,
               originalRows=len(original), finalSelectedRows=len(selected),
               allAssessableRiemannStatesIndependentlyRecomputed=False, originalRowsExactlyRetained=True,
               inputSHA256=bindings, scriptSHA256=sha(__file__),
               newActualNSEvaluations=0, newFlowAcceptances=0, physicalSolution=False, newProductionSolvedAreaM2=0)
    target = folder/'independent_cell_interface_continuation.json'; assert not target.exists()
    target.write_bytes((json.dumps(out, indent=2)+'\n').encode('utf8'))
    print(json.dumps(dict(status=out['status'], previews=len(checks), finalSelectedRows=len(selected),
                          publicSamplesPerRound=[c.get('publicSampleCount', 0) for c in checks])))


if __name__ == '__main__':
    main(sys.argv[1])
