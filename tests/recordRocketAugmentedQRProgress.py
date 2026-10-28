"""绑定v85真实终态、独立参考、原面UI和研究交付，不提升区域面积资格。"""
from pathlib import Path
from datetime import datetime,timezone
import json,subprocess
from packageRocketAugmentedQR import ROOT,OUT,CHANGED
from screenRocketFullFLeastSquares import sha


def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def main():
    previous=ROOT/'docs/validation/v84/progress_finite_step.json'
    p=read(ROOT/'outputs/v84_validation/full02/baseline/system2/peak_system.json')
    processes={}
    names=['qrunit85','qrreal85','qrview85']
    if (ROOT/'outputs/v83_validation/qrisolated85_process.json').exists():names.append('qrisolated85')
    for name in names:
        file=ROOT/'outputs/v83_validation'/(name+'_process.json');r=read(file)
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['sourceStableAfter'].values())
        processes[name]=dict(file=str(file.relative_to(ROOT)),sha256=sha(file),status=r['status'],exitCode=r['exitCode'],allLockedSourcesStable=True,engineProcessIDs=r['engineProcessIDs'])
    original={}
    for name,digest in p['sourceIdentity'].items():
        path=ROOT/('internalFlow/' if name=='resolveNSCharacteristicBoundary3D' else 'internalFlow/cutCell/')/(name+'.m')
        assert sha(path)==digest;original[name]=digest
    files={
        'independentSameObjectiveScreen':OUT/'qr_screen02/qr_screen.json',
        'unitIndependentAudit':OUT/'qr_unit01/independent_qr_audit.json',
        'realStudy':OUT/'qr_real01/augmented_qr_study.json',
        'realIndependentQRAudit':OUT/'qr_real01/independent_qr_audit.json',
        'wallAndUpdateAudit':ROOT/'outputs/v83_validation/krylov26_independent_audit.json',
        'freshNSDirectionAudit':ROOT/'outputs/v83_validation/krylov26_independent_direct_ns_probe.json',
        'publicOpenFaceAndHistoryAudit':ROOT/'outputs/v83_validation/krylov26_independent_riemann_reservoir.json',
        'rejectedInterface':OUT/'rejected01/rejected_interface_summary.json',
        'independentTemperatureLimit':OUT/'qr_temperature_limit.json',
        'displaySourceBinding':OUT/'qr_view_source_binding.json',
        'workbenchBinding':OUT/'workbench01/qr_view_binding.json',
        'workbenchNodeAudit':OUT/'workbench01/node_wall_velocity_summary.json',
        'workbenchIndependentAudit':OUT/'workbench01/independent_wall_velocity_view.json',
        'sourceRejectionAudit':OUT/'qr_view_binding_rejections.json',
        'scientificVisualReview':OUT/'visual_review.json'}
    evidence={key:dict(file=str(path.relative_to(ROOT)),sha256=sha(path),result=read(path)) for key,path in files.items()}
    status='tested-current-QR-real-NS-step-and-original-facet-workbench'
    delivery=OUT/'delivery/package_verification.json';isolated=OUT/'delivery/isolated_comparison.json'
    r=dict(version=85,status=status,recordedUTC=datetime.now(timezone.utc).isoformat(),
           previousProgress=dict(file=str(previous.relative_to(ROOT)),sha256=sha(previous)),
           previousGoalTurnClassification='no-progress-in-physics-or-UI-version-estimate-and-document-correction-only',
           currentGoalTurnClassification='progress-implemented-independent-validated-new-linear-map-actual-NS-step-and-bound-UI',
           inputState=p['sourceBinding'],inputPreparation=dict(file='outputs/v84_validation/full02/baseline/system2/peak_system.json',
                                                            sha256=sha(ROOT/'outputs/v84_validation/full02/baseline/system2/peak_system.json')),
           originalPhysicsSourcesUnchanged=original,actualTerminalProcesses=processes,evidence=evidence,
           sourceSHA256={name:sha(ROOT/name) for name in CHANGED if (ROOT/name).exists() and name!='docs/validation/v85/progress_augmented_qr.json'},
           isolatedUnitReproduction=read(isolated) if isolated.exists() else None,
           delivery=read(delivery) if delivery.exists() else None,
           production=dict(formalHeatAreaPercent=66.5626327424,newProductionSolvedAreaM2=0,physicalSolution=False,
                           fullAreaMilestoneCompleted=False,wholeRocketNonlinearBottleneckResolved=False,defaultSolverReplaced=False),
           git=dict(HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),committed=False,pushed=False),
           nextAction='Use actual Q26 and first/second-layer face19207 temperature margin. Inspect the interface-temperature chain derivative and constrained-direction/energy-transport feasibility at the current state, preserving original NS/temperature/steady gates. Do not repeat the old 15-column selection or assume a better linear residual resolves nonlinear convergence. Continue regional/near-wall/grid/outer-domain/property/turbulence/public physical validation before adding qualified area.')
    target=ROOT/'docs/validation/v85/progress_augmented_qr.json';target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(status=status,path=str(target),isolatedUnitVerified=isolated.exists(),deliveryPresent=delivery.exists())))


if __name__=='__main__':main()
