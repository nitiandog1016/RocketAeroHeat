"""记录v90真实研究/独立审计/浏览器/隔离与交付身份，保留物理资格限制。"""
from pathlib import Path
from datetime import datetime,timezone
import json
from packageRocketTemperatureLookahead import ROOT,OUT,DEL,CHANGED
from screenRocketFullFLeastSquares import sha

def read(file):return json.loads(Path(file).read_text(encoding='utf8'))
processes={}
for name in ('lookaheadunit90','lookaheadintegration90','lookaheadfull90','lookaheadrank90',
             'targetedunit90','targetedrecovery90','activeview90','activeisolated90'):
    file=ROOT/'outputs/v83_validation'/(name+'_process.json');p=read(file)
    assert p['status']=='terminal' and p['exitCode']==(1 if name=='lookaheadfull90' else 0)
    assert all(p['sourceStableAfter'].values())
    processes[name]=dict(file=str(file.relative_to(ROOT)),sha256=sha(file),exitCode=p['exitCode'],sourceStableDuringRun=True,
                         engineProcessIDs=p['engineProcessIDs'])
files={
 'unitSelection':OUT/'targeted01/independent_targeted_selection.json',
 'realDirection':OUT/'continuous01/independent_temperature_active.json',
 'sameStateRecovery':OUT/'recovery01/independent_lookahead_recovery.json',
 'fullNSAndWall':ROOT/'outputs/v83_validation/krylov38_independent_audit.json',
 'freshNSProbe':ROOT/'outputs/v83_validation/krylov38_independent_direct_ns_probe.json',
 'allOpenReservoirFaces':ROOT/'outputs/v83_validation/krylov38_independent_riemann_reservoir.json',
 'currentTemperatureMargins':OUT/'current_temperature_margins.json',
 'rankDiagnosis':OUT/'rank01/temperature_lookahead_rank.json',
 'distinctStarNormalIdentity':OUT/'rank01/star_pair_input_identity.json',
 'rejectedInterface':OUT/'interface01/interface_rejection.json',
 'publicRejectedInterface':OUT/'interface01/independent_interface_rejection.json',
 'originalFacetWorkbench':OUT/'workbench01/independent_wall_velocity_view.json',
 'nodeCSVSectionsAndWorkflow':OUT/'workbench01/node_wall_velocity_summary.json',
 'viewBinding':OUT/'workbench01/temperature_active_view_binding.json',
 'wrongViewSourceRejections':OUT/'view_binding_rejections.json',
 'browserReview':OUT/'browser_review.json','visualReview':OUT/'visual_review.json',
 'correctedScientificPlot':OUT/'scientific_plot.json','preservedFirstUIFailure':OUT/'ui_build_failure01/failure.json',
 'isolatedEntry':DEL/'isolated_comparison.json'}
evidence={k:dict(file=str(p.relative_to(ROOT)),sha256=sha(p),result=read(p)) for k,p in files.items()}
study=read(OUT/'continuous01/active_continuation.json');step=study['steps'][0];actual=step['main']
assert actual['totalNewEquationAcceptedIterations']==16 and actual['totalNewEquationResidualEvaluations']==70
assert actual['acceptedNonlinearIterations']==1 and not actual['steadyConverged']
assert step['newPreparationActualNS']==1 and step['newPreparationActualNSInRecovery']==0
assert step['newProbeActualNS']==2 and step['newMainActualNS']==5
assert not actual['physicalSolution'] and actual['newProductionSolvedAreaM2']==0
delivery=DEL/'package_verification.json'
out=dict(version=90,status='tested-current-Q37-lookahead-recovery-Q38-and-original-facet-UI-with-isolated-default-entry',
    recordedUTC=datetime.now(timezone.utc).isoformat(),previousGoalTurnClassification='progress',
    selectedState=dict(case='krylov38',file='outputs/v83_validation/krylov38/wall_velocity_krylov_evidence.mat',
        sha256=sha(ROOT/'outputs/v83_validation/krylov38/wall_velocity_krylov_evidence.mat'),
        acceptedAlpha=.0625,rmsResidualPerSecond=actual['rmsMerit'],maximumResidualPerSecond=actual['maximumResidual'],
        steadyBoundaryFluxBalance=actual['steadyBoundaryFluxBalance'],steadyConverged=False),
    sameEquationTotals=dict(acceptedIterations=16,actualNS=70),
    currentStudyWork=dict(preparationNSInFailedAttempt=1,newPreparationNSInRecovery=0,freshProbeNS=2,newMainNS=5,
        reusedFreeIterations=128,newProjectedIterations=128,newFreeIterationsInRecovery=0),
    targetedSelection=dict(originalRows=124,predictedViolatingCells=25,selectedFaces=131,outwardActiveRows=130,distinctConditions=48,
        allOriginalHAndSignsRetained=True,allPredictedRequiredCVsRepresented=True,rankUsedForSelection=False,rankTruncationUsed=False),
    limitation=dict(largerAlphaNotImproved=True,rejectedAlpha=.125,rejectedInterfaceFaceID=97923,
        rejectedStarTemperatureK=2007.5494927613981,cellOnlyPreviewDoesNotProtectAllRiemannSamples=True,
        steadyPhysicalRegionNotQualified=True,freshFullEntryBranchFromQ38NotSeparatelyExecuted=True,
        fullDomainRecomputedInIsolation=False),
    processes=processes,evidence=evidence,
    completeSourceFiles={n:sha(ROOT/n) for n in CHANGED if not n.endswith('progress_temperature_lookahead.json')},
    priorFrozenDelivery=dict(file='outputs/v89_validation/delivery/package_verification.json',sha256=sha(ROOT/'outputs/v89_validation/delivery/package_verification.json')),
    delivery=read(delivery) if delivery.exists() else dict(status='isolated-passed-archive-awaiting-actual-receipt'),
    formalHeatAreaFraction=.665626327424,newProductionSolvedAreaM2=0,physicalSolution=False,
    objectiveComplete=False,githubPushed=False)
target=ROOT/'docs/validation/v90/progress_temperature_lookahead.json';target.parent.mkdir(exist_ok=True)
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status=out['status'],delivery=out['delivery']['status'],selectedState=out['selectedState'],physicalSolution=False,newProductionSolvedAreaM2=0)))
