# -*- coding: utf-8 -*-
"""归档真实末态来源、Riemann导数、当前C、krylov21及有界状态方向研究。"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import subprocess

ROOT=Path(__file__).resolve().parents[1]; VAL=ROOT/'outputs/v83_validation'; DOC=ROOT/'docs/validation/v83'
target=DOC/'progress_current_state_direction.json'; assert not target.exists()


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''): h.update(b)
    return h.hexdigest()


def read(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def receipt(path): return dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),bytes=path.stat().st_size)


prior=DOC/'progress_riemann_reservoir.json'
assert sha(prior)=='ca2c3a8aa24535c2dde4c10a21f738e01c0767b998856b9e6608e404265916fb'
processes={}
for name in ['nsproberesultunit01','riemannreservoircontinue01','reservoirjacunit01','reservoirjacreal01',
             'reservoircurrentparent01','updateradiusunit01','updateradiusreal01','updateradiusns01','updateradiusui02']:
    path=VAL/(name+'_process.json'); p=read(path)
    assert p['status']=='terminal' and p['exitCode']==0 and all(p['sourceStableAfter'].values())
    processes[name]=dict(exitCode=0,sourceStable=True,engineProcessIDs=p['engineProcessIDs'],receipt=receipt(path))
checks=[]
for name in ['nsproberesultunit01_independent_direct_ns_probe.json','reservoirjacunit01_independent_reservoir_selected_direction.json',
             'reservoirjacreal01_independent_real_reservoir_direction.json','currentparent01_independent_current_matrix.json',
             'currentparent01_independent_current_candidate.json','krylov21_independent_riemann_reservoir.json',
             'krylov21_independent_direct_ns_probe.json','krylov21_independent_newton_directions.json',
             'krylov22_independent_direct_ns_probe.json','krylov22_independent_newton_directions.json','krylov22_independent_riemann_reservoir.json',
             'updateradiusunit01_independent_regularized_direction.json','updateradiusreal01_independent_regularized_direction.json',
             'workbench34/node_wall_velocity_summary.json','workbench34/independent_wall_velocity_view.json',
             'workbench35/node_wall_velocity_summary.json','workbench35/independent_wall_velocity_view.json']:
    path=VAL/name; r=read(path); assert r['status'].startswith('passed');epoch_bindings=[]
    for key in ['inputSHA256','sourceSHA256']:
        for input_name,digest in r.get(key,{}).items():
            if sha(ROOT/input_name)!=digest:
                # 只允许明确已归档的旧result来源单测；其旧资格不能改绑成新16项源码。
                assert path.name=='nsproberesultunit01_independent_direct_ns_probe.json'
                assert Path(input_name).name in ['buildCutCellSelectedDirectionCache.m','cutCellBoundaryDirectionalState.m','cutCellExactEulerBranchJacobian.m']
                old=VAL/'source_epochs'/('pre_reservoir_jac_'+Path(input_name).name);assert sha(old)==digest
                epoch_bindings.append(dict(originalBoundPath=input_name,originalSource=receipt(old),currentSourceQualificationClaimed=False))
    checks.append(dict(receipt=receipt(path),status=r['status'],originalSourceEpochBindings=epoch_bindings))
sources=['internalFlow/cutCell/'+name+'.m' for name in ['readRocketNSProbeBase','loadRocketDirectNSDirectionProbe',
    'cutCellExactEulerBranchJacobian','buildCutCellRiemannReservoirJacobian','buildCutCellSelectedDirectionCache',
    'cutCellBoundaryDirectionalState','rocketSelectedDirectionSourceIdentity','buildRocketPrimitiveUpdateMetric',
    'applyRocketPrimitiveUpdateMetric','solveRocketRegularizedParentDirection']]
sources+=['tests/'+name for name in ['createRocketNSDirectionProbe.m','continueRocketWallVelocityKrylov.m',
    'testRocketNSProbeResultSource.m','testRiemannReservoirSelectedDirection.m','studyRiemannReservoirSelectedDirection.m',
    'studyRocketCurrentReservoirParent.m','testRocketRegularizedParentDirection.m','studyRocketRegularizedParentDirection.m','studyRocketRegularizedDirectionNS.m',
    'diagnoseRocketCurrentNSDirection.py','diagnoseRocketNSUpdateRadius.py','auditRiemannReservoirSelectedDirection.py',
    'auditRiemannReservoirRealDirection.py','auditRocketCurrentParentMatrix.py','auditRocketCurrentParentCandidate.py',
    'auditRocketDirectNSDirectionProbe.py','auditNSRiemannReservoir3D.py','testRocketWallVelocityWorkbench.js',
    'auditRocketWallVelocityWorkbench.py','auditRocketRegularizedParentDirection.py','recordRocketCurrentDirectionProgress.py','packageRocketCurrentDirectionResearch.py']]
sources+=['exportRocketWallVelocityStudy.m','web/rocket-wall-velocity.js','docs/CURRENT_STATE_DIRECTION_V83.md',
          'docs/PROJECT_GOALS.md','docs/GITHUB_STAGE_READINESS.md']
hashes={name:sha(ROOT/name) for name in sources}
epochs=[receipt(p) for pattern in ['pre_result_source_*','pre_reservoir_jac_*','pre_reservoir_continue_*','pre_reservoir_norm_roundoff_*',
         'pre_residual_trend_20261026/*','pre_regularized_objective_roundoff_*','pre_current_progress_epoch_*'] for p in sorted((VAL/'source_epochs').glob(pattern)) if p.is_file()]
endpoint=read(VAL/'krylov22/wall_velocity_krylov_summary.json')
assert endpoint['totalNewEquationAcceptedIterations']==2 and endpoint['totalNewEquationResidualEvaluations']==9
assert endpoint['acceptedNonlinearIterations']==1 and endpoint['candidateTrials']==13
assert endpoint['maximumResidual']>endpoint['initialMaximumResidual'] and endpoint['rmsMerit']<endpoint['initialRmsMerit']
assert not endpoint['steadyConverged'] and not endpoint['physicalSolution'] and endpoint['newProductionSolvedAreaM2']==0
regularized=read(VAL/'updateradiusreal01/regularized_parent_summary.json')
assert not regularized['physicalStateChanged'] and regularized['newActualNSEvaluations']==regularized['newFlowAcceptances']==0
original_ns_study=read(VAL/'updateradiusns01/regularized_direction_ns_summary.json')
assert original_ns_study['alternatePredecessorCase']=='krylov20' and original_ns_study['referenceAcceptedCase']=='krylov21'
assert not original_ns_study['alternateHistoriesMerged'] and original_ns_study['newProductionSolvedAreaM2']==0
assert original_ns_study['main']==endpoint
reference=read(VAL/'krylov21/wall_velocity_krylov_summary.json')
assert original_ns_study['referenceRmsMerit']==reference['rmsMerit'] and original_ns_study['referenceMaximumResidual']==reference['maximumResidual']
assert original_ns_study['candidateRMSBetterThanReference'] and original_ns_study['candidateMaximumBetterThanReference']
assert sha(original_ns_study['candidateDirectionSource'])==original_ns_study['candidateDirectionSourceSHA256']
diagnostic=read(VAL/'krylov21_update_radius.json')
assert diagnostic['newActualNSEvaluations']==0 and not diagnostic['flowStateChanged']
assert all(sha(ROOT/p)==h for p,h in diagnostic['inputSHA256'].items())
current_diagnostic=read(VAL/'krylov22_update_radius.json')
assert current_diagnostic['newActualNSEvaluations']==0 and not current_diagnostic['flowStateChanged']
assert all(sha(ROOT/p)==h for p,h in current_diagnostic['inputSHA256'].items())
view=read(VAL/'workbench35/wall-velocity-study-data.json'); assert view['cases'][2]['summary']['resultEvidenceFolder']=='krylov22'
trend=view['cases'][2]['solverDiagnostics']['segmentResidualTrend']; assert trend['rmsDecreased'] and trend['maximumIncreased'] and not trend['steadyNumericalTargetsMet']
head=subprocess.check_output(['git','rev-parse','--short','HEAD'],cwd=ROOT,text=True).strip(); assert head=='42601c9'
paths=['krylov21/wall_velocity_krylov_evidence.mat','directnsprobereal03/direct_ns_direction_probe.mat',
       'krylov22/wall_velocity_krylov_evidence.mat','directnsprobereal04/direct_ns_direction_probe.mat',
       'reservoirjacreal01/riemann_selected_direction_real_evidence.mat','currentparent01/current_parent_candidate_evidence.mat',
       'updateradiusunit01/regularized_parent_unit_evidence.mat','updateradiusreal01/regularized_parent_candidate_evidence.mat',
       'updateradiusreal01/regularized_parent_model_evidence.mat','krylov21_update_radius.json','krylov22_update_radius.json',
       'workbench34/wall-velocity-study-data.json','workbench34/static_wall_velocity_comparison.png',
       'workbench34/static_wall_velocity_residuals.png','workbench34/static_wall_velocity_balance.png',
       'workbench35/wall-velocity-study-data.json','workbench35/static_wall_velocity_comparison.png',
       'workbench35/static_wall_velocity_residuals.png','workbench35/static_wall_velocity_balance.png']
out=dict(version=83,status='tested-current-result-state-Riemann-Jacobian-fresh-parent-and-regularized-krylov22-direction',
    recordedUTC=datetime.now(timezone.utc).isoformat(),processID=os.getpid(),priorProgress=receipt(prior),
    previousGoalTurnClassification='no-progress-status-answer-to-user-GitHub-version-question',
    currentGoalTurnClassification='progress-verified-previous-terminals-local-state-radius-new-regularized-direction-and-original-facet-UI',
    currentSourcesSHA256=hashes,sourceEpochs=epochs,matlabProcesses=processes,independentChecks=checks,
    evidence=[receipt(VAL/p) for p in paths],latestAcceptedEndpoint=endpoint,
    frozenQ21StateRadiusDiagnostic=diagnostic,localStateRadiusDiagnostic=current_diagnostic,
    regularizedDirectionResearch=regularized,regularizedDirectionOriginalNSStudy=original_ns_study,
    retainedReferenceFailures=[dict(receipt=receipt(VAL/name),report=read(VAL/name)) for name in
        ['updateradiusreal01_objective_reference_failure.json','current_direction_recorder_source_epoch_failure.json']],
    workAccounting=dict(oldCharacteristicProblemAccepted=13,oldCharacteristicProblemMainNS=258,
        frozenQ21BranchAccepted=2,frozenQ21BranchMainNS=9,selectedQ22BranchAccepted=2,selectedQ22BranchMainNS=9,
        regularizedLinearResearchNewNS=0,parallelBranchesMerged=False,
        alternateQ20BranchProbeNS=original_ns_study['newProbeActualNSEvaluations'],
        alternateQ20BranchMainNS=original_ns_study['newMainActualNSEvaluations'],
        oldAndNewProblemHistoriesMerged=False,independentDerivativeAndProbeNSCallsKeptSeparate=True),
    userInterface=dict(workbench='outputs/v83_validation/workbench35/index.html',latestCase='krylov22',
        originalFaces=19752,segmentRMSDownMaximumUpExplicit=True,scientificPNGsViewed=True,browserVisualOrWebGLRuntimeVerified=False),
    production=dict(formalQualifiedHeatAreaPercent=66.5626327424,newProductionSolvedAreaM2=0,physicalSolution=False,
        latestFormalVersion=82,version83IsolatedCompleteEngineeringDelivery=False,wholeGoalComplete=False,goalStatus='active'),
    delivery=dict(gitHEAD=head,gitCommitCreated=False,githubPushed=False,smallVersionLocalOnly=True,
        researchSourceSnapshotIsFormalDelivery=False),
    nextAction='Do not repeat the same weak 0.01 L2-regularized 128-step configuration: Q22 improves the same-base residuals but still takes 1/4096 and its maximum residual rises 6.05 percent from Q20. The limiting CV75951 is unchanged; the state radius increases only to 0.000728819 and nonlinear remainder still concentrates in layers 4/5. Investigate explicit control of rare large local primitive updates, with an accurately matched current-state metric, total-direction consistency and unchanged original model/NS 0.1 gates before any nonlinear update. Q22 is the selected endpoint, but Q20/Q21 numeric Jacobians or probe qualifications cannot be inherited to Q22; refresh current sources and state. Keep parallel Q21 and Q22 histories separate. Verify regional steady or justified statistical convergence, near-wall/grid/exterior/boundary/physics and public benchmarks; full-area single-condition results and scientific UI remain highest priority.',
    recorderSHA256=sha(__file__))
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status=out['status'],manifest=str(target),sha256=sha(target))))
