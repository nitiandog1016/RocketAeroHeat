"""制造一个未接受的控制记录，验证工作台保留真实前驱，不能把候选方向冒充已更新结果。"""
from pathlib import Path
import json, sys
import buildRocketTemperatureActiveWorkbench as workbench
from screenRocketFullFLeastSquares import ROOT, sha


def main(folder):
    folder = Path(folder).resolve(); folder.relative_to(ROOT)
    original = folder/'continuous01/active_continuation.json'
    study = json.loads(original.read_text('utf8')); steps = study['steps']
    steps = [steps] if isinstance(steps, dict) else steps
    assert len(steps) == 1 and steps[0]['main']['acceptedNonlinearIterations'] == 1
    actual_endpoint = Path(steps[0]['resultFolder'])/'wall_velocity_krylov_evidence.mat'
    # 只制造控制元数据；真实MAT、实际接受结果和历史不变。该文件不是一次物理试算。
    step = steps[0]; step['main'] = []; step['resultFolder'] = ''; step['newMainActualNS'] = 0
    step['stopReason'] = 'manufactured-unaccepted-control-record-only'
    study['steps'] = steps; study['manufacturedControlRecordNotPhysicalRun'] = True
    target = folder/'noaccept_binding_unit'; assert not target.exists()
    (target/'continuous01').mkdir(parents=True)
    synthetic = target/'continuous01/active_continuation.json'
    synthetic.write_bytes((json.dumps(study)+'\n').encode())
    workbench.configure(target); workbench.aliases()
    selection = json.loads((target/'view_source_binding.json').read_text('utf8'))['selected']
    assert selection['displayedStateUnchanged'] and selection['displayCase'] == Path(study['initialStateFile']).parent.name
    proposal = workbench.bind(synthetic, step['directionFile'], study['initialStateFile'], step['step'])
    assert proposal['displayedStateUnchanged'] and proposal['proposedTotalDirectionBoundToCurrentState']
    assert not proposal['usedTotalDirectionExactlyVerified']
    rejected = False
    try:
        workbench.bind(synthetic, step['directionFile'], actual_endpoint, step['step'])
    except AssertionError:
        rejected = True
    assert rejected, '未接受的记录不能展示其他已更新端点。'
    result = dict(status='passed-manufactured-noaccept-control-current-state-identity-and-wrong-endpoint-rejection',
                  manufacturedControlRecordNotPhysicalRun=True, displayedStateUnchanged=True,
                  originalStudySHA256=sha(original), syntheticStudySHA256=sha(synthetic), scriptSHA256=sha(__file__),
                  physicalSolution=False, newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0)
    (target/'source_binding_test.json').write_bytes((json.dumps(result, indent=2)+'\n').encode())
    print(json.dumps(result))


if __name__ == '__main__': main(sys.argv[1])
