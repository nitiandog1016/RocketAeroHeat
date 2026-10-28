# -*- coding: utf-8 -*-
"""登记完整像列尺度构造、独立参考、真实有界F解及工作台，保持物理资格分离。"""
from pathlib import Path
from datetime import datetime,timezone
import json
import os
import subprocess
from screenRocketFullFLeastSquares import ROOT,VAL,sha


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def receipt(path):
    return dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),bytes=path.stat().st_size)


target=ROOT/'docs/validation/v83/progress_full_F_image_scale.json'
assert not target.exists()
prior_file=ROOT/'docs/validation/v83/progress_full_F_least_squares.json'
assert sha(prior_file)=='841b4f41a9644d4da9b02ef7b7311c15b7ca35a7211c0a461a45ed3d6a0e49c8'
prior=read(prior_file)
processes={}
for name in ('fullFimagescaleunit01','fullFimagescalereal01','fullFimagescaleui01'):
    file=VAL/(name+'_process.json');p=read(file)
    assert p['status']=='terminal' and p['exitCode']==0 and all(p['sourceStableAfter'].values())
    processes[name]=dict(exitCode=0,sourceStable=True,engineProcessIDs=p['engineProcessIDs'],receipt=receipt(file))
checks=[]
for name in ('fullFimagescaleunit01_independent_full_F_image_scale.json','fullFimagescalereal01_independent_full_F_image_scale.json',
    'workbench29/independent_full_F_image_scale_study_view.json','workbench29/node_wall_velocity_summary.json'):
    file=VAL/name;a=read(file)
    assert a['status'].startswith('passed-')
    for path,digest in a.get('inputSHA256',{}).items():
        path=Path(path);assert sha(path if path.is_absolute() else ROOT/path)==digest
    checks.append(dict(report=a,receipt=receipt(file)))
screen_file=VAL/'fullFimagescalescreen01/full_F_image_scaled_screen_summary.json';screen=read(screen_file)
assert screen['status'].startswith('tested-') and screen['newActualNSEvaluations']==screen['newProductionSolvedAreaM2']==0
for path,digest in screen['inputSHA256'].items(): assert sha(ROOT/path)==digest
assert sha(ROOT/'tests/screenRocketImageScaledFullF.py')==screen['scriptSHA256']
assert sha(ROOT/'tests/screenRocketFullFLeastSquares.py')==screen['helperSourceSHA256']
file=VAL/'fullFimagescalereal01/full_F_image_scaled_summary.json';r=read(file)
assert not r['physicalSolution'] and r['newActualNSEvaluations']==r['newProductionSolvedAreaM2']==r['newFlowAcceptances']==0
for key,digest in r['sourceIdentity'].items(): assert sha(ROOT/'internalFlow/cutCell'/(key+'.m'))==digest
assert sha(VAL/'source_epochs/pre_image_scale_solveRocketFullFLeastSquares.m')=='b90a446df82b795adde1ec4f8dc58a6567228f6b1136f37235c168b52a841d8d'
production=prior['production']
assert sha(ROOT/production['latestEndpoint']['path'])==production['latestEndpoint']['sha256']
assert sha(ROOT/production['originalNSSource']['path'])==production['originalNSSource']['sha256']
head=subprocess.check_output(['git','rev-parse','--short','HEAD'],cwd=ROOT,text=True).strip()
assert head==prior['delivery']['gitHEAD']=='42601c9'
sources=['internalFlow/cutCell/'+name+'.m' for name in ('buildRocketFullFImageScale','solveRocketFullFLeastSquares')]
sources+=['tests/'+name for name in ('screenRocketImageScaledFullF.py','testRocketFullFImageScale.m','studyRocketFullFImageScale.m',
    'auditRocketFullFImageScale.py','buildRocketFullFImageScaleSummary.py','auditRocketFullFImageScaleWorkbench.py','recordRocketFullFImageScaleProgress.py',
    'testRocketWallVelocityWorkbench.js')]
sources+=['web/rocket-wall-velocity.js']
if r['fullFReport']['fullFInternalTargetMet']:
    next_action='The frozen complete F test now meets its original 1e-3 target; this is not original C or NS convergence. Move the main numerical work back to the complete original C using an explicitly flexible right preconditioner with the new controlled F solve. First prove small dense block/restriction/back-substitution identities, variable-inner-work accounting and flexible Arnoldi Z/image combinations. The old S may only be an explicitly approximate preconditioning object; do not claim the old S identity or use a right-hand-side-dependent F inner solve as a fixed Schur operator. Keep the complete original C action fixed, bound inner work/storage, then do a bounded real trial. Original selected-J/NS total-direction and regional physical/public validation gates remain required.'
else:
    next_action='The original F target remains unmet; keep that identity explicit. For the saved new vector the fixed-C parent remaining is now dominated by C rows, so do not merely add F-only LSQR steps. Reintroduce all original C unknowns through a complete-parent fixed-operator flexible right-preconditioning study. The bounded scaled F inner solves remain explicitly approximate and may be nonlinear in the RHS; prove small dense block/back-substitution identities, variable-inner accounting and actual flexible Arnoldi Z/image combinations before a bounded real run. Old S can only be an explicitly approximate preconditioning object, not an exact identity or an operator formed with a RHS-dependent F inverse. Preserve the full-parent 1e-3 and original selected-J/NS total-direction/physics gates. Scientific UI and actual regional heat solutions remain highest priorities.'
out=dict(version=83,status='tested-complete-F-image-column-normalization-bounded-real-LSQR-and-workbench',
    recordedUTC=datetime.now(timezone.utc).isoformat(),processID=os.getpid(),priorProgress=receipt(prior_file),
    previousGoalTurnClassification='progress-complete-forward-adjoint-LSQR-and-native-stop-repair',
    currentGoalTurnClassification='progress-complete-image-column-conditioning-independent-scaled-solver-and-workbench',
    processes=processes,independentChecks=checks,independentPrototype=screen,independentPrototypeReceipt=receipt(screen_file),
    currentRealSummary=r,currentRealSummaryReceipt=receipt(file),codeSourceIdentity=[receipt(ROOT/path) for path in sources],
    savedParentCRowsRemainingSquaredNormFraction=(r['parentDiagnostic']['CRowsRelativeToOriginalRHS']/r['parentDiagnostic']['relativeParentOriginalCResidual'])**2,
    preservedOldLeastSquaresSource=receipt(VAL/'source_epochs/pre_image_scale_solveRocketFullFLeastSquares.m'),
    methodDocumentation=receipt(ROOT/'docs/FULL_F_IMAGE_SCALE_V83.md'),
    goalAndReadinessSnapshots=[receipt(ROOT/path) for path in ('docs/PROJECT_GOALS.md','docs/GITHUB_STAGE_READINESS.md')],
    userInterface=dict(workbench='outputs/v83_validation/workbench29/index.html',researchRows=19,
        oldEighteenRowsAndSixFacetFilesExactlyPreserved=True,fullFParentAndNSQualificationKeptSeparate=True,browserVisualOrWebGLRuntimeVerified=False),
    production=production,delivery=prior['delivery'],newSelectedJApplications=0,newActualNSEvaluations=0,
    newFlowAcceptances=0,newProductionSolvedAreaM2=0,nextAction=next_action,goalStatus='active',wholeGoalComplete=False,recorderSHA256=sha(__file__))
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status=out['status'],manifest=str(target),sha256=sha(target),fullFInternalTargetMet=r['fullFReport']['fullFInternalTargetMet'])))
