# -*- coding: utf-8 -*-
"""绑定一致Schur前置拒绝、诊断修正及工作台，保持物理解状态与正式交付身份。"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,subprocess
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=ROOT/'docs/validation/v83/progress_consistent_schur.json';assert not target.exists()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for b in iter(lambda:stream.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def identity(path):
 path=Path(path);return dict(path=str(path.relative_to(ROOT)),sha256=sha(path),bytes=path.stat().st_size)
names=['consistentschurunit01','consistentschurunit02','consistentschurunit03','consistentschurreal01','consistentschurreal02','consistentdiagrefine01','consistentschurui01'];processes={}
for name in names:
 file=VAL/(name+'_process.json');p=read(file);assert p['status']=='terminal' and p['exitCode']==(1 if name=='consistentschurreal01' else 0) and all(p['sourceStableAfter'].values())
 processes[name]=dict(exitCode=p['exitCode'],engineProcessIDs=p['engineProcessIDs'],sourceStable=True,receipt=identity(file))
audits=[read(VAL/(name+'_independent_consistent_schur.json')) for name in ('consistentschurunit01','consistentschurunit02','consistentschurunit03')]
assert all(a['status'].startswith('passed-') for a in audits)
real=read(VAL/'consistentdiagrefine01_independent_consistent_schur.json');assert real['status'].startswith('passed-') and real['diagnosticOnlyRefinementKeepsAllOtherPayloadValues']
failure=real['fixedMomentumFailure'];assert failure['independentFirstNonfiniteStageAndRowExactlyMatch'] and failure['firstAlgebraicRow']==167894 and failure['originalCV']==4570
assert real['consistentSchurActionsCompleted']==0 and not real['completeGlobalNumericalRankOrUsableInverseProven']
view=VAL/'workbench21';data=read(view/'wall-velocity-study-data.json');research=read(VAL/'block_linear_research_summary_v5.json');assert data['linearResearch']==research and len(research['studies'])==8
row=research['studies'][-1];assert row['displayState']=='动量回代溢出（已定位）' and not row['candidateUsable']
for key in ('innerRelativeCResidual','basisOrthogonalityError','modelRelationRelativeDifference'):assert row[key] is None
assert research['fixedMomentumInverseDiagnostic']['oldNativeScalarZeroPivotLocationStillUnknown']
node=read(view/'node_wall_velocity_summary.json');assert node['status'].startswith('passed-') and node['CSVrows']==592560 and node['sectionSamples']==6624 and node['sameHostEvents']==15
assert not node['browserVisualOrWebGLRuntimeVerified'];payload=[]
for case in data['cases']:
 rel=Path(case['resultsJSONURL'])
 for path in (rel,rel.parent/'case-data.js'):
  assert sha(view/path)==sha(VAL/'workbench20'/path);payload.append(identity(view/path))
assert len(payload)==6
ui=dict(status='passed-independent-eight-source-bound-research-rows-overflow-location-and-unchanged-facet-payload',processID=os.getpid(),workbench='outputs/v83_validation/workbench21/index.html',
 originalSRejectionLocationNotOverwritten=True,nativeScalarZeroPivotLocationStillUnknown=True,sourceFacetPayload=payload,
 previousFullMATFieldAuditReusedByByteIdentity=True,physicalGeometryDefectClaimed=False,browserVisualOrWebGLRuntimeVerified=False)
(view/'independent_consistent_schur_study_view.json').write_text(json.dumps(ui,indent=2)+'\n',encoding='utf8')
old=read(ROOT/'docs/validation/v83/progress_equilibration.json');production=old['production']
for key in ('latestEndpoint','originalNSSource'):assert sha(ROOT/production[key]['path'])==production[key]['sha256']
assert subprocess.check_output(['git','rev-parse','--short','HEAD'],cwd=ROOT,text=True).strip()=='42601c9'
paths=[ROOT/'internalFlow/cutCell'/name for name in ('applyRocketConsistentSchur.m','observeRocketFieldILU0.m','rocketTwoFieldSingularRatio.m','solveRocketFieldILU0.m')]
paths += [ROOT/'tests'/name for name in ('testRocketConsistentSchur.m','studyRocketConsistentSchurActions.m','refineRocketSchurDiagonalDiagnostics.m','auditRocketConsistentSchurUnit.py','auditRocketConsistentSchurReal.py','exportRocketConsistentSchurDiagnostics.py')]
out=dict(status='tested-consistent-Schur-prerequisite-rejected-independent-first-momentum-overflow-and-stable-local-diagonal-audit',version=83,
 recordedUTC=datetime.now(timezone.utc).isoformat(),processID=os.getpid(),processes=processes,smallDomainAudits=audits,realAudit=real,
 repairHistory=dict(firstUnobservedRealFailureRetained=True,firstRealMissingIntermediateEvidenceNotClaimedComplete=True,
 diagonalFormulaFailure=read(VAL/'consistent_schur_diagonal_formula_failure.json'),stableHypotRepairPassedOriginal1eMinus14DiagnosticCheck=True,
 oldAndRefinedMATAllOtherUserPayloadExactlyEqual=True,MATContainerCreationHeaderExcludedFromPayloadComparison=True,
 originalModelGateChanged=False,originalMomentumSolverSourceChanged=False),
 workAccounting=dict(observedRealRunFixedMomentumSolveAttempts=2,secondSolveStoppedAtFirstNonfinite=True,consistentSchurActionsCompleted=0,
 diagnosticRefinementNewMomentumSolves=0,newSelectedJApplications=0,newActualNSEvaluations=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0,
 unobservedFailedRealAttemptInternalCallCountFullySaved=False),
 userInterface=ui,nodeDataVerification=node,codeSourceIdentity=[identity(p) for p in paths],production=production,delivery=old['delivery'],
 nextAction='Do not repeat the unusable global three-momentum ILU, restart the stopped scalar ILUTP, rebuild S, or scan fill/drop thresholds. First inspect a bounded stable inverse/action for the raw momentum A and its actual residual, separating ILU propagation growth from A weak modes. Compare justified local/normal blocks before deciding whether to retain this field Schur split. Preserve all CVs, original C/NS, direction/temperature/descent/convergence gates and regional-validation/scientific-UI priorities.',
 goalStatus='active',wholeGoalComplete=False,recorderSHA256=sha(__file__))
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8');print(json.dumps(dict(status=out['status'],manifest=identity(target),researchRows=8)))
