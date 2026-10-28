"""核对失败任务、同态自由检查点复用及原目标独立证据，不重复计算NS。"""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT, sha

out_folder = ROOT / 'outputs/v88_validation'
direction_file = out_folder / 'recovery01/active_direction.mat'
s = loadmat(direction_file, simplify_cells=True,
            variable_names=['evidence', 'directionReport', 'preparationSHA256'])
replay = s['directionReport']['freeCheckpointReplay']
checkpoint = Path(replay['file']); assert sha(checkpoint) == replay['sha256']
original = loadmat(checkpoint, simplify_cells=True, variable_names=['free', 'freeCoordinates', 'freeDetail'])
for name in ('free', 'freeCoordinates'):
    np.testing.assert_array_equal(s['evidence'][name], original[name])
assert s['preparationSHA256'] == replay['currentPreparationSHA256']
assert original['freeDetail']['nativeIterations'] == s['directionReport']['reusedFreeNativeIterations'] == 128
assert s['directionReport']['newNativeIterations'] == s['directionReport']['linearReport']['nativeIterations'] == 128
assert s['directionReport']['totalNativeIterations'] == 256 and replay['newActualNSEvaluations'] == 0
audit_file = out_folder / 'continuous01/independent_temperature_active.json'
audit = json.loads(audit_file.read_text('utf8'))
last = audit['checks']['steps'][-1]
free = next(x for x in last['directions'] if x['label'] == 'free')
np.testing.assert_allclose(free['originalCResidual'], replay['currentOriginalCResidual'], rtol=1e-10)
np.testing.assert_allclose(free['normalizedAugmentedObjective'], replay['currentAugmentedObjective'], rtol=1e-10)
representation = last['activeArithmetic']['cycles'][0]['exactDuplicateRepresentation']
assert representation['originalRows'] == 7 and representation['independentConditions'] == 3
assert representation['allOriginalRowsExactlyReconstructed'] and not representation['rankTruncationUsed']
failure_file = ROOT / 'outputs/v83_validation/activefull88_process.json'
failure = json.loads(failure_file.read_text('utf8'))
assert failure['status'] == 'terminal' and failure['exitCode'] == 1 and all(failure['sourceStableAfter'].values())
old_source = out_folder / 'source_epochs/buildRocketTemperatureProjector_before_exact_duplicates.m'
bound = next(v for k, v in failure['sourcesBefore'].items()
             if k.replace('\\', '/') == 'internalFlow/cutCell/buildRocketTemperatureProjector.m')
assert sha(old_source) == bound
proof = dict(status='passed-failed-source-current-free-replay-and-seven-original-row-identity',
    reusedFreeIterations=128, newlyExecutedProjectedIterations=128,
    currentOriginalCAndObjectiveIndependentlyChecked=True,
    allSevenOriginalRowsRetained=True, independentConditions=3,
    originalRankFailurePreserved=True, newActualNSEvaluations=0,
    physicalSolution=False, newProductionSolvedAreaM2=0,
    inputSHA256={str(f.relative_to(ROOT)): sha(f) for f in
                 (direction_file, checkpoint, audit_file, failure_file, old_source)}, scriptSHA256=sha(__file__))
target = out_folder / 'independent_recovery_identity.json'; assert not target.exists()
target.write_text(json.dumps(proof, indent=2)+'\n', encoding='utf8')
print(json.dumps({k: v for k, v in proof.items() if k != 'inputSHA256'}))
