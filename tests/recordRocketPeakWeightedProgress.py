# -*- coding: utf-8 -*-
"""归档同方向坐标对照与当前状态峰值约束，保留未解资格及并列历史。"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,subprocess

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
target=ROOT/'docs/validation/v83/progress_peak_weighted_direction.json';assert not target.exists()


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def receipt(path):return dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),bytes=path.stat().st_size)


prior=ROOT/'docs/validation/v83/progress_current_state_direction.json'
assert sha(prior)=='a3d727a9cb3d5fcbc39b6ca3230f001f98f9338d1ac69e5ebe016e14d595cbb3'
process_names=['primitivereservoir01','peakweightunit01','peakweightreal01','peakzerounit01','peakzeroreal01','peakzeroreal02']
ns_file=VAL/'peakweightns01/peak_direction_ns_summary.json'
if ns_file.exists():process_names.append('peakweightns01')
if (VAL/'peakweightui01_process.json').exists():process_names.append('peakweightui01')
if (VAL/'peakweightui02_process.json').exists():process_names.append('peakweightui02')
processes={}
for name in process_names:
    path=VAL/(name+'_process.json');p=read(path)
    assert p['status']=='terminal' and p['exitCode']==0 and all(p['sourceStableAfter'].values())
    processes[name]=dict(exitCode=0,sourceStable=True,engineProcessIDs=p['engineProcessIDs'],receipt=receipt(path))
names=['krylov23_independent_riemann_reservoir.json','krylov23_independent_direct_ns_probe.json',
       'krylov23_independent_newton_directions.json','primitivereservoir01_independent_same_direction_trial.json',
       'peakweightunit01_independent_peak_direction.json','peakweightreal01_independent_current_matrix.json',
       'peakweightreal01_independent_peak_direction.json','workbench36/node_wall_velocity_summary.json',
       'workbench36/independent_wall_velocity_view.json','peakzerounit01_independent_zero_start.json','peakzeroreal01_independent_zero_start.json',
       'peakzeroreal02_independent_zero_start.json']
ns=read(ns_file) if ns_file.exists() else None
if ns and ns['main']:
    names+=['krylov24_independent_riemann_reservoir.json','krylov24_independent_direct_ns_probe.json','krylov24_independent_newton_directions.json']
if (VAL/'workbench38').exists():names+=['workbench38/node_wall_velocity_summary.json','workbench38/independent_wall_velocity_view.json']
checks=[]
for name in names:
    p=VAL/name;r=read(p);assert r['status'].startswith('passed');epoch_bindings=[]
    for key in ['inputSHA256','sourceSHA256']:
        for path,digest in r.get(key,{}).items():
            if sha(ROOT/path)!=digest:
                assert name=='workbench36/node_wall_velocity_summary.json' and Path(path).name=='testRocketWallVelocityWorkbench.js'
                old=VAL/'source_epochs/pre_peak_direction_provenance_20261026/testRocketWallVelocityWorkbench.js';assert sha(old)==digest
                epoch_bindings.append(dict(originalPath=path,originalSource=receipt(old),currentCodeQualificationClaimed=False))
    checks.append(dict(receipt=receipt(p),status=r['status'],originalSourceEpochBindings=epoch_bindings))
peak=read(VAL/'peakweightreal01/current_peak_summary.json')
assert peak['sourceCase']=='krylov22' and peak['currentStateCompleteNSExactlyReplayed']
assert not peak['oldNumericMatricesFactorsOrDerivativesReused'] and peak['newActualNSEvaluations']==1
assert peak['newFlowAcceptances']==peak['newProductionSolvedAreaM2']==0 and not peak['physicalStateChanged']
zero=read(VAL/'peakzeroreal01/peak_zero_summary.json')
assert zero['samePenaltyAsWarmRun'] and zero['newActualNSEvaluations']==zero['newFlowAcceptances']==0 and not zero['physicalStateChanged']
continued=read(VAL/'peakzeroreal02/peak_zero_summary.json')
assert continued['samePenaltyAsWarmRun'] and not continued['unrestartedLSQRClaimed']
assert continued['newActualNSEvaluations']==continued['newFlowAcceptances']==0
q22=read(VAL/'krylov22/wall_velocity_krylov_summary.json');q23=read(VAL/'krylov23/wall_velocity_krylov_summary.json')
assert q22['totalNewEquationAcceptedIterations']==q23['totalNewEquationAcceptedIterations']==2
assert q22['totalNewEquationResidualEvaluations']==9 and q23['totalNewEquationResidualEvaluations']==10
assert q23['rmsMerit']>q22['rmsMerit'] and q23['maximumResidual']<q22['maximumResidual']
selected='krylov24' if ns and ns['main'] else 'krylov22'
selected_summary=read(VAL/selected/'wall_velocity_krylov_summary.json')
if ns and ns['main']:
    assert ns['main']==selected_summary and ns['predecessorCase']=='krylov22' and not ns['parallelQ23HistoryAppended']
    assert sha(ns['candidateFile'])==ns['candidateFileSHA256'] and ns['probe']['originalNSDirectionTargetMet']
source_names=['internalFlow/cutCell/'+n+'.m' for n in ['reweightRocketPrimitiveUpdateMetric','solveRocketRegularizedParentDirection',
    'buildRocketPrimitiveUpdateMetric','applyRocketPrimitiveUpdateMetric','rocketShellTrialState']]
source_names+=['tests/'+n for n in ['studyRocketPrimitiveReservoirTrial.m','diagnoseRocketPrimitiveTrialRadius.py',
    'testRocketPeakWeightedDirection.m','studyRocketCurrentPeakDirection.m','auditRocketPeakWeightedDirection.py',
    'auditRocketCurrentParentMatrix.py','auditRocketPrimitiveReservoirTrial.py','auditRocketSteadyNewtonDirections.py',
    'plotRocketTrialCoordinateComparison.py','testRocketWallVelocityWorkbench.js','auditRocketWallVelocityWorkbench.py',
    'testRocketRegularizedZeroStart.m','studyRocketPeakZeroStart.m','continueRocketPeakZeroStart.m','auditRocketPeakZeroStart.py',
    'recordRocketPeakWeightedProgress.py','packageRocketPeakResearch.py']]
if ns_file.exists():source_names.append('tests/studyRocketPeakDirectionNS.m')
source_names.append('tests/testRocketDirectionPreparationBinding.m')
source_names+=['exportRocketWallVelocityStudy.m','web/rocket-wall-velocity.js','docs/PEAK_WEIGHTED_DIRECTION_V83.md',
    'docs/PROJECT_GOALS.md','docs/GITHUB_STAGE_READINESS.md']
hashes={name:sha(ROOT/name) for name in source_names}
epochs=[receipt(p) for pattern in ['pre_primitive_reservoir_trial_20261026/*','pre_peak_raw_metric_labels_*','pre_peak_matrix_layout_*',
                                 'pre_peak_norm_reference_*','pre_peak_fixed_normalization_*','pre_peak_zero_continue_*','pre_peak_option_type_*',
                                 'pre_peak_direction_provenance_20261026/*']
        for p in sorted((VAL/'source_epochs').glob(pattern)) if p.is_file()]
figure=VAL/'primitivereservoir01/trial_coordinate_figure_receipt.json';fr=read(figure)
assert all(sha(ROOT/path)==digest for path,digest in fr['inputSHA256'].items())
assert sha(VAL/'primitivereservoir01/trial_coordinate_comparison.png')==fr['figureSHA256']
diagnostic=read(VAL/'krylov24_update_radius.json')
assert diagnostic['newActualNSEvaluations']==0 and not diagnostic['flowStateChanged']
assert all(sha(ROOT/path)==digest for path,digest in diagnostic['inputSHA256'].items())
preserved=read(VAL/'peak_provenance_original_fields_preservation.json')
assert preserved['allDatasetBytesIdentical']
for row in preserved['cases']:assert sha(ROOT/row['path37'])==sha(ROOT/row['path38'])==row['sha256']
assert subprocess.check_output(['git','rev-parse','--short','HEAD'],cwd=ROOT,text=True).strip()=='42601c9'
out=dict(version=83,status='tested-same-direction-primitive-coordinate-comparison-and-fresh-current-peak-weighted-research',
    recordedUTC=datetime.now(timezone.utc).isoformat(),processID=os.getpid(),priorProgress=receipt(prior),
    previousGoalTurnClassification='progress-current-state-result-source-Jacobian-fresh-direction-Q22-and-original-facet-workbench',
    currentGoalTurnClassification='progress-same-direction-nonlinear-path-test-fresh-Q22-peak-weighted-operator-and-direction',
    currentSourcesSHA256=hashes,sourceEpochs=epochs,matlabProcesses=processes,independentChecks=checks,
    primitiveCandidateStudy=read(VAL/'primitivereservoir01/primitive_reservoir_trial_summary.json'),
    currentPeakDirectionStudy=peak,fixedObjectiveZeroStartStudy=zero,fixedObjectiveRestartedContinuation=continued,
    peakOriginalNSStudy=ns,selectedEndpointCase=selected,selectedEndpoint=selected_summary,
    selectedEndpointEvidence=receipt(VAL/selected/'wall_velocity_krylov_evidence.mat'),
    selectedConservativeUpdateRadiusDiagnostic=diagnostic,sourceMetadataPreservesFacetFields=preserved,
    selectedStateChangedFromQ22=bool(ns and ns['main'] and ns['main']['acceptedNonlinearIterations']>0),
    retainedReferenceFailure=dict(receipt=receipt(VAL/'peakweightreal01_norm_reference_failure.json'),report=read(VAL/'peakweightreal01_norm_reference_failure.json')),
    parallelQ20Branches=dict(Q22Accepted=2,Q22MainNS=9,Q23Accepted=2,Q23MainNS=10,historiesMerged=False),
    workAccounting=dict(newPrimitiveCoordinateMainNS=q23['residualEvaluations'],newCurrentStateReplayNS=peak['newActualNSEvaluations'],
        newZeroStartAndRestartNS=zero['newActualNSEvaluations']+continued['newActualNSEvaluations'],
        newQualifiedDirectionProbeNS=ns['newProbeActualNS'] if ns else 0,newSelectedMainNS=ns['newMainActualNS'] if ns else 0,
        selectedSameProblemMainNS=selected_summary['totalNewEquationResidualEvaluations'],
        selectedSameProblemAcceptances=selected_summary['totalNewEquationAcceptedIterations'],parallelAndProbeCallsMergedIntoMainHistory=False),
    userInterface=dict(selectedWorkbench='outputs/v83_validation/workbench38/index.html' if selected=='krylov24' else 'outputs/v83_validation/workbench35/index.html',
        actualDirectionPreparationSourceBound=True,
        coordinateStudyWorkbench='outputs/v83_validation/workbench36/index.html',primitiveCandidatePathExplicit=True,
        sameDirectionScientificFigure=receipt(VAL/'primitivereservoir01/trial_coordinate_comparison.png'),scientificPNGsViewed=True,browserVisualOrWebGLRuntimeVerified=False),
    production=dict(formalQualifiedHeatAreaPercent=66.5626327424,newProductionSolvedAreaM2=0,physicalSolution=False,
        latestFormalVersion=82,version83IsolatedCompleteEngineeringDelivery=False,goalStatus='active',wholeGoalComplete=False),
    delivery=dict(gitHEAD='42601c9',gitCommitCreated=False,githubPushed=False,smallVersionLocalOnly=True),
    nextAction='Use actual Q24 and its same-problem 3-accepted/13-main-NS history. Do not rerun the Q20 coordinate switch or the failed old-seed 256-step warm solve. Separate current-state construction from direction solution: replay NS and rebuild the Q24 operator/source binding, use the latest direction only as a metric reference, and start the bounded fixed-objective solve directly from zero. Q22 Jacobians/probes cannot qualify the changed Q24. Inspect the new 2000-K first-exit at CV38959/layer2/sourceFacet19207 and the still-dominant layer4/5 nonlinear remainder, without relaxing any original state/temperature/RMS/maximum/space/steady-boundary gates. Keep regional grid/near-wall/exterior/boundary/physics/public validation and scientific UI as highest priorities; no qualified area was added.',
    recorderSHA256=sha(__file__))
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8');print(json.dumps(dict(status=out['status'],path=str(target),sha256=sha(target),selectedEndpoint=selected)))
