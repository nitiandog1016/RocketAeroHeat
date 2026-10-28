# -*- coding: utf-8 -*-
"""登记完整父C的可变内层试算与真实全外层方向核对，不晋升为原NS或区域热解。"""
from pathlib import Path
from datetime import datetime,timezone
import json,os,subprocess
from screenRocketFullFLeastSquares import ROOT,VAL,sha


def read(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def receipt(path): return dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),bytes=path.stat().st_size)


target=ROOT/'docs/validation/v83/progress_flexible_reduction.json';assert not target.exists()
prior_file=ROOT/'docs/validation/v83/progress_full_F_image_scale.json'
assert sha(prior_file)=='341bdac948947d00d3365ae75f3bb79114486cdde13235d961538d7013f9d2ea'
prior=read(prior_file);processes={}
for name,code in [('flexiblereductionunit01',1),('flexiblereductionunit02',0),('flexiblereductionreal01',0),('flexiblereductionui01',0)]:
    file=VAL/(name+'_process.json');p=read(file)
    assert p['status']=='terminal' and p['exitCode']==code and all(p['sourceStableAfter'].values())
    processes[name]=dict(exitCode=code,sourceStable=True,engineProcessIDs=p['engineProcessIDs'],receipt=receipt(file))
checks=[]
for name in ('flexiblereductionunit02_independent_flexible_reduction.json','flexiblereductionreal01_independent_flexible_reduction.json',
    'workbench30/independent_flexible_reduction_study_view.json','workbench30/node_wall_velocity_summary.json'):
    file=VAL/name;a=read(file);assert a['status'].startswith('passed-')
    for path,digest in a.get('inputSHA256',{}).items():
        path=Path(path);assert sha(path if path.is_absolute() else ROOT/path)==digest
    checks.append(dict(report=a,receipt=receipt(file)))
failures=[]
for name in ('flexiblereductionunit01_reference_failure.json','flexiblereductionunit02_audit_key_failure.json','parentCimagescalescreen01_reference_failure.json'):
    file=VAL/name;d=read(file);assert d['driverExitCode']==1 and sha(ROOT/d['sourceFile'])==d['sourceSHA256']
    failures.append(dict(report=d,receipt=receipt(file)))
real_file=VAL/'flexiblereductionreal01/flexible_reduction_summary.json';r=read(real_file)
assert r['innerRunComplete'] and not r['candidateRejected'] and not r['physicalSolution']
assert r['newActualNSEvaluations']==r['newFlowAcceptances']==r['newProductionSolvedAreaM2']==0
for key,digest in r['sourceIdentity'].items(): assert sha(ROOT/'internalFlow/cutCell'/(key+'.m'))==digest
diagnostic_file=VAL/'flexiblereductionreal01/flexible_schur_defect_diagnostic.json';diagnostic=read(diagnostic_file)
assert diagnostic['status'].startswith('passed-') and diagnostic['newActualNSEvaluations']==0
for path,digest in diagnostic['inputSHA256'].items(): assert sha(ROOT/path)==digest
screen_file=VAL/'parentCimagescalescreen02/parent_image_scale_screen_summary.json';screen=read(screen_file)
assert screen['status'].startswith('tested-') and screen['allOriginalPhysicalCVs']==217272
for path,digest in screen['inputSHA256'].items(): assert sha(ROOT/path)==digest
assert sha(ROOT/'tests/screenRocketParentImageScale.py')==screen['scriptSHA256']
for name,digest in screen['helperSourceSHA256'].items(): assert sha(ROOT/'tests'/name)==digest
production=prior['production']
assert sha(ROOT/production['latestEndpoint']['path'])==production['latestEndpoint']['sha256']
assert sha(ROOT/production['originalNSSource']['path'])==production['originalNSSource']['sha256']
head=subprocess.check_output(['git','rev-parse','--short','HEAD'],cwd=ROOT,text=True).strip();assert head==prior['delivery']['gitHEAD']=='42601c9'
sources=['internalFlow/cutCell/'+name+'.m' for name in ('applyRocketFlexibleReduction','solveRocketFlexibleReductionKrylov','compactRocketFlexibleReductionReport')]
sources+=['tests/'+name for name in ('testRocketFlexibleReduction.m','studyRocketFlexibleReduction.m','auditRocketFlexibleReduction.py',
    'diagnoseRocketFlexibleSchurDefect.py','screenRocketParentImageScale.py','buildRocketFlexibleReductionSummary.py','auditRocketFlexibleReductionWorkbench.py','recordRocketFlexibleReductionProgress.py','testRocketWallVelocityWorkbench.js')]
sources+=['web/rocket-wall-velocity.js']
next_action=('The complete original C internal target is met for this frozen RHS, but not the original selected J or NS direction. Preserve this explicit qualification and the stored full outer Z evidence, then validate composition into the complete selected-J solve and original-NS total direction under unchanged gates.' if r['internalCTargetMet'] else
    'Do not repeat the failed eight-direction/two-32-step-F configuration or just improve the old coarse inverse: for all eight saved vectors the changed-F-map contribution dominates the old-S local-equation defect. The independent direct full-parent image-column-normalized LSQR screen now improves the same saved original-C endpoint from 0.983921 to 0.354065 in 32 steps versus unscaled 0.916554. It retains all original CVs/five fields and original norm without old S or fixed C representatives, but is not yet a MATLAB delivered solver or a qualified result. Implement and validate that direct whole-parent route in MATLAB, preserving native/default compatibility, all-column scales/adjoint identities and actual full-parent images before selected-J/NS total-direction checks. Original C 1e-3, Arnoldi 1e-8 where applicable, original NS/physics/public validation and scientific UI goals remain unchanged.')
out=dict(version=83,status='tested-complete-parent-C-flexible-variable-inner-preconditioning-and-full-outer-column-audit',
    recordedUTC=datetime.now(timezone.utc).isoformat(),processID=os.getpid(),priorProgress=receipt(prior_file),
    previousGoalTurnClassification='progress-draft-flexible-block-map-and-terminal-dense-absolute-reference-failure',
    currentGoalTurnClassification='progress-propagated-roundoff-reference-fix-dense-flexible-proof-and-bounded-real-parent-C-test',
    processes=processes,independentChecks=checks,retainedFailures=failures,currentRealSummary=r,currentRealSummaryReceipt=receipt(real_file),
    savedDirectionSchurDefectDiagnostic=diagnostic,savedDirectionSchurDefectDiagnosticReceipt=receipt(diagnostic_file),
    directParentImageScalePrototype=screen,directParentImageScalePrototypeReceipt=receipt(screen_file),
    codeSourceIdentity=[receipt(ROOT/path) for path in sources],methodDocumentation=receipt(ROOT/'docs/FLEXIBLE_REDUCTION_V83.md'),
    goalAndReadinessSnapshots=[receipt(ROOT/path) for path in ('docs/PROJECT_GOALS.md','docs/GITHUB_STAGE_READINESS.md')],
    userInterface=dict(workbench='outputs/v83_validation/workbench30/index.html',researchRows=20,oldNineteenRowsAndSixFacetFilesExactlyPreserved=True,
        originalCInternalQualificationNotPromotedToNSOrRegion=True,browserVisualOrWebGLRuntimeVerified=False),
    production=production,delivery=prior['delivery'],newSelectedJApplications=0,newActualNSEvaluations=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0,
    nextAction=next_action,goalStatus='active',wholeGoalComplete=False,recorderSHA256=sha(__file__))
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status=out['status'],manifest=str(target),sha256=sha(target),internalCTargetMet=r['internalCTargetMet'])))
