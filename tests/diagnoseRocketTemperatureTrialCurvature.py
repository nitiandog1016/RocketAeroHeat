"""用已保存方向检查温度限制与候选路径曲率；不计算新NS或推定候选可接受。"""
from pathlib import Path
import os, sys, json
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,sha

folder=Path(sys.argv[1]).resolve(); target=folder/'temperature_trial_curvature.json'; assert not target.exists()
progress=folder/'progress.json'; digest=sha(progress); steps=json.loads(progress.read_text('utf8'))
steps=[steps] if isinstance(steps,dict) else steps; assert sha(progress)==digest
checks=[]; bindings={}
for step in steps:
    if not step['main'] or not step['main']['acceptedNonlinearIterations']: continue
    file=Path(step['resultFolder'])/'wall_velocity_krylov_evidence.mat'; bindings[str(file.relative_to(ROOT))]=sha(file)
    saved=loadmat(file,simplify_cells=True,variable_names=['result','gas']); result=saved['result']; gas=saved['gas']
    trace=result['trace']; trace=trace if isinstance(trace,dict) else list(trace)[0]
    Q=trace['initialQ']; dQ=trace['linearUpdateU']*result['residualScale']; cv=gas['R']/(gas['gamma']-1)
    rho=Q[:,0]; velocity=Q[:,1:4]/rho[:,None]
    T=(Q[:,4]-.5*rho*np.sum(velocity**2,axis=1))/(rho*cv)
    drho=dQ[:,0]; du=(dQ[:,1:4]-velocity*drho[:,None])/rho[:,None]
    dT=(dQ[:,4]-np.sum(velocity*dQ[:,1:4],axis=1)+(.5*np.sum(velocity**2,axis=1)-cv*T)*drho)/(rho*cv)
    def primitive_path(alpha):
        r=rho+alpha*drho; u=velocity+alpha*du; temperature=T+alpha*dT
        return np.column_stack([r,r[:,None]*u,r*(cv*temperature+.5*np.sum(u*u,axis=1))])
    tiny=1e-6; tangent=(primitive_path(tiny)-primitive_path(-tiny))/(2*tiny)
    tangent_error=float(np.linalg.norm(tangent-dQ)/max(1.,np.linalg.norm(dQ))); assert tangent_error<1e-7
    accepted=Q+float(trace['alpha'])*dQ
    np.testing.assert_array_equal(accepted,result['Q'])
    trials=result['lineLog']; trials=[trials] if isinstance(trials,dict) else list(trials); comparisons=[]
    for trial in trials:
        alpha=float(trial['alpha']); versions=[]
        for label,candidate in [('conservative',Q+alpha*dQ),('primitive-tangent',primitive_path(alpha))]:
            r=candidate[:,0]; u=candidate[:,1:4]/r[:,None]
            internal=candidate[:,4]-.5*r*np.sum(u*u,axis=1); temperature=internal/(r*cv)
            positive_density=bool(np.all(r>0)); positive_internal=bool(np.all(internal>0))
            legal_temperature=bool(np.all((temperature>=gas['temperatureRange'][0])&(temperature<=gas['temperatureRange'][1])))
            if label=='conservative':
                error=trial['errorIdentifier']; error=error if isinstance(error,str) else ''
                if error=='CutNS:State':assert not positive_density
                elif error=='CutNS:InternalEnergy':assert positive_density and not positive_internal
                elif error=='CutNS:Temperature':assert positive_density and positive_internal and not legal_temperature
                else:assert positive_density and positive_internal and legal_temperature
            versions.append(dict(coordinates=label,positiveDensity=positive_density,positiveInternalEnergy=positive_internal,
                allCellTemperaturesInOriginalRange=legal_temperature,
                algebraicTemperatureRangeK=[float(temperature.min()),float(temperature.max())],
                originalCellGuardPassed=positive_density and positive_internal and legal_temperature))
        comparisons.append(dict(alpha=alpha,originalAccepted=bool(trial['accepted']),paths=versions))
    checks.append(dict(step=step['step'],stateFile=str(file),physicalTangentRelativeDifference=tangent_error,
                       comparisons=comparisons,interfacesOrFullNSRecomputed=False,primitiveCandidateAcceptedClaimed=False))
assert all(sha(ROOT/path)==h for path,h in bindings.items())
out=dict(status='passed-saved-conservative-guards-and-same-tangent-primitive-path-cell-only-diagnostic',
         frozenProgressSHA256=digest,completedStepCount=len(checks),steps=checks,inputSHA256=bindings,scriptSHA256=sha(__file__),
         newActualNSEvaluations=0,newFlowAcceptances=0,physicalSolution=False,newProductionSolvedAreaM2=0,
         nextActionConditional='Only if cell feasibility improves, compare the same direction through original complete NS; do not infer interface legality or acceptance.')
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status=out['status'],completedSteps=len(checks),
    cellLegalPaths=[dict(step=c['step'],alpha=t['alpha'],conservative=t['paths'][0]['originalCellGuardPassed'],
        primitive=t['paths'][1]['originalCellGuardPassed']) for c in checks for t in c['comparisons']])))
