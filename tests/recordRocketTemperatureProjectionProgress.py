"""记录v86实际当前态、温度约束/同窗UI与交付，保留未收敛和面积身份。"""
from pathlib import Path
from datetime import datetime,timezone
import json,subprocess
from packageRocketTemperatureProjection import ROOT,OUT,CHANGED
from screenRocketFullFLeastSquares import sha


def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def main():
    prior=ROOT/'docs/validation/v85/progress_augmented_qr.json';p=read(OUT/'current01/peak_system.json');processes={}
    for name in ('temperatureunit86','temperatureunit86b','temperaturefull86','temperatureview86','temperatureisolated86'):
        path=ROOT/'outputs/v83_validation'/(name+'_process.json')
        if not path.exists():continue
        r=read(path);assert r['status']=='terminal' and all(r['sourceStableAfter'].values())
        assert r['exitCode']==(1 if name=='temperatureunit86' else 0)
        processes[name]=dict(file=str(path.relative_to(ROOT)),sha256=sha(path),exitCode=r['exitCode'],sourceStable=True,engineProcessIDs=r['engineProcessIDs'])
    files={'unitReference':OUT/'unit02/independent_temperature_projection.json',
        'realComparison':OUT/'comparison01/temperature_projection_comparison.json',
        'independentDirectionAudit':OUT/'comparison01/independent_temperature_projection.json',
        'scalarKernelSource':OUT/'scalar_guard_source.json','rejectedInterface':OUT/'rejected01/rejected_interface_summary.json',
        'acceptedCurrentTemperatureMargins':OUT/'current_temperature_margins.json','sourceBinding':OUT/'view_source_binding.json',
        'viewBinding':OUT/'workbench01/temperature_view_binding.json','sourceRejections':OUT/'view_binding_rejections.json',
        'nodeViewAudit':OUT/'workbench01/node_wall_velocity_summary.json','independentViewAudit':OUT/'workbench01/independent_wall_velocity_view.json',
        'scientificPlot':OUT/'scientific_plot.json','visualReview':OUT/'visual_review.json'}
    for case in (27,28):
        for suffix in ('independent_audit','independent_direct_ns_probe','independent_riemann_reservoir'):
            files[f'Q{case}_{suffix}']=ROOT/f'outputs/v83_validation/krylov{case}_{suffix}.json'
    if (OUT/'preview_server.json').exists():files['localStaticPreviewResources']=OUT/'preview_server.json'
    evidence={key:dict(file=str(path.relative_to(ROOT)),sha256=sha(path),result=read(path)) for key,path in files.items()}
    old=read(prior)['originalPhysicsSourcesUnchanged'];changed=[]
    for name,h in old.items():
        path=ROOT/('internalFlow/' if name=='resolveNSCharacteristicBoundary3D' else 'internalFlow/cutCell/')/(name+'.m')
        if sha(path)!=h:changed.append(name)
    assert changed==['cutCellExactEulerBranchJacobian']
    delivery=OUT/'delivery/package_verification.json';isolated=OUT/'delivery/isolated_comparison.json'
    out=dict(version=86,status='tested-actual-temperature-projected-QR-step-and-same-host-workbench',recordedUTC=datetime.now(timezone.utc).isoformat(),
        previousProgress=dict(file=str(prior.relative_to(ROOT)),sha256=sha(prior)),previousGoalTurnClassification='progress',currentGoalTurnClassification='progress',
        actualTerminalProcesses=processes,sourceState=p['sourceBinding'],currentPreparation=p,
        originalNSResidualAndBoundarySourcesUnchanged=True,changedSelectedDerivativeSource=changed,
        derivativeChangeScope='scalar-empty-right-selection-guard; old-full-multiface-outputs-bitwise-replay-and-single-vs-duplicate-batch-tested',
        evidence=evidence,sourceSHA256={n:sha(ROOT/n) for n in CHANGED if (ROOT/n).exists() and n!='docs/validation/v86/progress_temperature_projection.json'},
        isolatedUnitReproduction=read(isolated) if isolated.exists() else None,delivery=read(delivery) if delivery.exists() else None,
        production=dict(formalHeatAreaPercent=66.5626327424,newProductionSolvedAreaM2=0,physicalSolution=False,
            fullAreaMilestoneCompleted=False,wholeRocketSteadyConverged=False,defaultSolverReplaced=False,parallelHistoriesMerged=False),
        git=dict(HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),committed=False,pushed=False),
        nextAction='Use actual selected Q28 (6 acceptances/28 main NS) and current 849336-interface margin audit. Rebuild current derivatives before multi-active-temperature-row feasibility/rank/adjoint checks; distinguish heating/cooling and nonlinear curvature. Face19146 is the first rejected face of the Q26 projected 1/64 trial, not the unique or hottest Q28 constraint. Preserve original NS/temperature/steady gates; no area qualification without regional/mesh/near-wall/outer-domain/property/turbulence/public physical validation.')
    target=ROOT/'docs/validation/v86/progress_temperature_projection.json';target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(status=out['status'],file=str(target),isolated=isolated.exists(),delivery=delivery.exists())))


if __name__=='__main__':main()
