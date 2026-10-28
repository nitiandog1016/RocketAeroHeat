"""记录本地有限步诊断、实际终态与交付身份；未完成项保持明确。"""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,subprocess
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v84_validation'
SOURCES=['internalFlow/cutCell/assertRocketNSReplay.m','internalFlow/cutCell/prepareRocketCurrentPeakSystem.m',
 'internalFlow/cutCell/rocketPeakSystemSourceIdentity.m','internalFlow/cutCell/rocketFiniteStepAgreement.m',
 'internalFlow/cutCell/solveRocketMinimumStateCoefficients.m','internalFlow/cutCell/invertRocketPrimitiveUpdateMetric.m',
 'internalFlow/cutCell/selectRocketMinimumStateDirection.m','tests/continueRocketWallVelocityKrylov.m',
 'tests/diagnoseRocketFiniteStep.m','tests/rocketFiniteStepPatchResidual.m','tests/studyRocketFiniteStepPatch.m',
 'tests/diagnoseRocketPatchConditioning.m','tests/studyRocketMinimumStateFull.m','tests/testRocketMinimumStateDirection.m','tests/testRocketMinimumStateReturnedModel.m',
 'tests/auditRocketFiniteStep.py','tests/auditRocketMinimumStateFull.py','tests/diagnoseRocketFiniteStepMechanism.py',
 'tests/plotRocketFiniteStep.py','tests/buildRocketFiniteStepWorkbench.py','tests/testRocketFiniteStepPackage.m',
 'tests/testRocketWallVelocityWorkbench.js','tests/prepareRocketFiniteStepViewSource.py',
 'tests/recordRocketFiniteStepProgress.py','tests/packageRocketFiniteStep.py','tests/verifyRocketFiniteStepIsolation.py','runRocketFiniteStepResearch.m',
 'runRocketFullMinimumStateComparison.m','web/rocket-wall-velocity.js','docs/FINITE_STEP_DIAGNOSIS_V84.md','docs/PROJECT_GOALS.md']
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()
def read(name):
 p=OUT/name
 return json.loads(p.read_text('utf-8-sig')) if p.exists() else None
def main():
 full=read('full02/full_comparison.json');audit=read('full02/independent_audit.json')
 p=read('q24system/peak_system.json')
 report={'version':84,'status':'tested-local-reproduction-full-comparison-'+('terminal' if full else 'incomplete'),
 'recordedUTC':datetime.now(timezone.utc).isoformat(),'sourceSHA256':{n:sha(ROOT/n) for n in SOURCES if (ROOT/n).exists()},
 'inputState':p['sourceBinding'],'inputGeometry':p['geometryBinding'],'mechanism':read('mechanism.json'),
 'originalPhysicsSourcesUnchanged':read('unchanged_physics_sources.json'),
 'diagnosis':read('diagnosis01/diagnosis.json'),'unit':read('unit01/unit.json'),
 'returnedModelRegression':read('returned_model01/returned_model.json'),
 'failedAgreementGuard':read('patch01/patch_comparison.json'),'localSixStepComparison':read('patch02/patch_comparison.json'),
 'independentLocalAudit':read('independent_audit.json'),'failed256StepFullComparison':read('full01/full_comparison.json'),
 'matchedConstraintLocalComparison':read('patch03/patch_comparison.json'),'independentMatchedAudit':read('independent_audit_patch03.json'),
 'fullComparison':full,'independentFullAudit':audit,'initialCurrentSystemConstruction':p,
 'sharedLinearSegments':[read(f'q24direction{k}/peak_zero_summary.json') for k in [1,2,3]],
 'interruptedRun':{'oldSession':1639,'oldEngine':53160,'exitCodeUnknown':True,'completeCheckpointsPreserved':True,
    'unsavedInterruptedWorkUnknown':True,'resumedFrom':'baseline/system2 + direction2_2; first accepted state preserved'},
 'userInterface':read('workbench01/finite_step_view_binding.json'),
 'workbenchNodeAudit':read('workbench01/node_wall_velocity_summary.json'),
 'workbenchIndependentAudit':read('workbench01/independent_wall_velocity_view.json'),
 'scientificVisualReview':read('visual_review.json'),
 'isolatedPackageRun':read('delivery/isolated/RocketAeroHeat/outputs/isolated_validation/isolated_package.json'),
 'isolatedIndependentComparison':read('delivery/isolated_comparison.json'),
 'delivery':read('delivery/package_verification.json'),
 'production':{'physicalSolution':False,'newProductionSolvedAreaM2':0,'formalHeatAreaPercent':66.5626327424,
    'steadyOrExperimentalValidationClaimed':False,'wholeRocketBottleneckResolved':False,'defaultSolverReplaced':False},
 'git':{'HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'committed':False,'pushed':False},
 'nextAction':'Use actual full02 terminal records; do not repeat the failed 0.5-vector-remainder guard or 15-column infeasible direction. Investigate unresolved state-space modes before claiming whole-rocket acceleration; preserve all original gates.'}
 target=ROOT/'docs/validation/v84/progress_finite_step.json';target.parent.mkdir(parents=True,exist_ok=True)
 target.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
 print(json.dumps({'status':report['status'],'path':str(target),'fullAuditPresent':audit is not None}))
if __name__=='__main__':main()
