"""核对v90真实同态恢复、目标驱动选面和自由前缀复用，不重算NS。"""
from pathlib import Path
import json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT, sha
from auditRocketTemperatureTargetedSelection import check, vector


def main(folder):
    folder = Path(folder).resolve(); direction_file = folder/'active_direction.mat'
    data = loadmat(direction_file, simplify_cells=True)
    r = data['directionReport']; replay = r['freeCheckpointReplay']; prefix = Path(replay['file'])
    assert sha(prefix) == replay['sha256'] and sha(replay['contextFile']) == replay['contextSHA256']
    saved = loadmat(prefix, simplify_cells=True)
    for key in ('free', 'freeCoordinates'):
        np.testing.assert_array_equal(saved[key], data['evidence'][key])
    assert replay['reusedNativeIterations'] == r['reusedFreeNativeIterations'] == saved['freeDetail']['nativeIterations'] == 128
    assert replay['newFreeNativeIterations'] == replay['newActualNSEvaluations'] == 0
    assert r['totalNativeIterations'] == 128+r['newNativeIterations']
    initial = saved['originalRows']; rows = data['rows']; positions = np.searchsorted(vector(rows['faceIDs']), vector(initial['faceIDs']))
    assert (initial['H'] - rows['H'][positions]).nnz == 0
    np.testing.assert_array_equal(vector(initial['temperatureBoundSigns']), vector(rows['temperatureBoundSigns'])[positions])
    s = loadmat(data['sourceBinding']['sourceFile'], simplify_cells=True, variable_names=['result', 'gas'])
    g = loadmat(data['geometryBinding']['file'], simplify_cells=True, variable_names=['domain'])['domain']
    selected = check(s['result']['Q'], s['result']['info'], g, s['gas'], s['result']['residualScale'], saved['free'], rows, r['lookaheadSelection'])
    failure_file = ROOT/'outputs/v83_validation/lookaheadfull90_process.json'
    failure = json.loads(failure_file.read_text('utf8'))
    assert failure['status'] == 'terminal' and failure['exitCode'] == 1 and all(failure['sourceStableAfter'].values())
    failed_sources = {k.replace('\\', '/'): v for k, v in failure['sourcesBefore'].items()}
    epochs = []
    for name in ('selectRocketTemperatureLookaheadRows', 'solveRocketTemperatureActiveDirection'):
        old = ROOT/('outputs/v90_validation/source_epochs/'+name+'_before_targeted_selection.m')
        assert sha(old) == failed_sources['internalFlow/cutCell/'+name+'.m']; epochs.append(old)
    audit_file = ROOT/'outputs/v90_validation/continuous01/independent_temperature_active.json'
    audit = json.loads(audit_file.read_text('utf8')); verified = audit['checks']['steps'][0]
    free = next(x for x in verified['directions'] if x['label'] == 'free')
    np.testing.assert_allclose(free['originalCResidual'], replay['currentOriginalCResidual'], rtol=1e-10)
    np.testing.assert_allclose(free['normalizedAugmentedObjective'], replay['currentAugmentedObjective'], rtol=1e-10)
    assert verified['activeArithmetic']['allCurrentBandRowsFeasible'] == r['bandDirectionFeasible']
    out = dict(status='passed-independent-current-Q37-free-prefix-targeted-all-required-CVs-and-original-failure-binding',
               reusedFreeIterations=128, newFreeIterations=0, newProjectedIterations=int(r['newNativeIterations']),
               completeCurrentObjectiveIndependentlyChecked=True, targetedSelection=selected,
               originalFailurePreserved=True, original124RowsAndSignsPreserved=True, rankTruncationUsed=False,
               newActualNSEvaluations=0, newFlowAcceptances=0, physicalSolution=False, newProductionSolvedAreaM2=0,
               inputSHA256={str(p.relative_to(ROOT)): sha(p) for p in [direction_file, prefix, failure_file, audit_file]+epochs},
               scriptSHA256=sha(__file__))
    target = folder/'independent_lookahead_recovery.json'; assert not target.exists()
    target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
    print(json.dumps({k: v for k, v in out.items() if k not in ('targetedSelection', 'inputSHA256')}))


if __name__ == '__main__':
    main(sys.argv[1])
