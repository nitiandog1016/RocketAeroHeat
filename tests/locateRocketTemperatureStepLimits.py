"""定位已记录的单元温域拒绝，区分当前H支撑和新越界单元，不计算新方向/NS。"""
from pathlib import Path
import os, sys, json
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,sha

snapshot=Path(sys.argv[1]).resolve(); target=Path(sys.argv[2]).resolve(); assert not target.exists()
steps=json.loads(snapshot.read_text('utf8')); steps=steps.get('steps',steps) if isinstance(steps,dict) else steps
steps=[steps] if isinstance(steps,dict) else steps
geometry=ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat'
bindings={str(snapshot.relative_to(ROOT)):sha(snapshot),str(geometry.relative_to(ROOT)):sha(geometry)}
d=loadmat(geometry,simplify_cells=True,variable_names=['domain'])['domain']; nf=int(d['quality']['sourceFaces']); checks=[]
for step in steps:
    file=Path(step['resultFolder'])/'wall_velocity_krylov_evidence.mat'; direction=Path(step['directionFile'])
    for f in (file,direction):bindings[str(f.relative_to(ROOT))]=sha(f)
    s=loadmat(file,simplify_cells=True,variable_names=['result','gas']); r=s['result']; gas=s['gas']; cv=gas['R']/(gas['gamma']-1)
    rows=loadmat(direction,simplify_cells=True,variable_names=['rows'])['rows']; n=rows['H'].shape[1]//5
    support=np.unique(rows['H'].tocsr().indices % n)
    expected=np.unique(np.r_[np.atleast_1d(rows['ownerCV'])[np.any(np.atleast_2d(rows['leftRowsK'])!=0,axis=1)]-1,
                              np.atleast_1d(rows['neighbourCV'])[np.any(np.atleast_2d(rows['rightRowsK'])!=0,axis=1)]-1]).astype(int)
    np.testing.assert_array_equal(support,expected)
    trace=r['trace']; trace=trace if isinstance(trace,dict) else list(trace)[0]
    Q=trace['initialQ']; rho=Q[:,0]; u=Q[:,1:4]/rho[:,None]
    T=(Q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv); delta=trace['linearUpdateU']*r['residualScale']
    dT=(delta[:,4]-np.sum(u*delta[:,1:4],axis=1)+(.5*np.sum(u*u,axis=1)-cv*T)*delta[:,0])/(rho*cv)
    trials=r['lineLog']; trials=[trials] if isinstance(trials,dict) else list(trials)
    rejected=[t for t in trials if isinstance(t['errorIdentifier'],str) and t['errorIdentifier']=='CutNS:Temperature']; assert rejected
    last=rejected[-1]; alpha=float(last['alpha']); q=Q+alpha*trace['linearUpdateU']*r['residualScale']
    rr=q[:,0]; uu=q[:,1:4]/rr[:,None]; temperature=(q[:,4]-.5*rr*np.sum(uu*uu,axis=1))/(rr*cv)
    assert np.all(rr>0) and np.all(temperature>0)
    bad=np.flatnonzero((temperature<gas['temperatureRange'][0])|(temperature>gas['temperatureRange'][1])); assert len(bad)
    sample=[]
    for i in bad:
        sample.append(dict(cellID=int(i+1),sourceFacet=int(i%nf+1),layer=int(i//nf+1),
            centroidM=d['cells']['centroid'][i].tolist(),baseTemperatureK=float(T[i]),
            tangentTemperatureChangeK=float(dT[i]),linearPredictedTemperatureK=float(T[i]+alpha*dT[i]),
            rejectedConservativeTemperatureK=float(temperature[i]),inCurrentHSupport=bool(i in support),
            baseCellWithinOriginal10KBand=bool(min(T[i]-gas['temperatureRange'][0],gas['temperatureRange'][1]-T[i])<=10)))
    checks.append(dict(step=step['step'],rejectedAlpha=alpha,rejectedCellCount=len(bad),
        rejectedCellsInCurrentHSupport=int(np.isin(bad,support).sum()),
        allRejectedCellsOutsideCurrentHSupport=bool(not np.any(np.isin(bad,support))),cells=sample,
        cellOnlyDiagnostic=True,allInterfaceDerivativesRecomputed=False,alternateDirectionOrNSAcceptanceClaimed=False))
assert all(sha(ROOT/path)==value for path,value in bindings.items())
out=dict(status='passed-saved-last-cell-temperature-rejection-locations-and-current-H-support-comparison',
    steps=checks,inputSHA256=bindings,scriptSHA256=sha(__file__),newActualNSEvaluations=0,newFlowAcceptances=0,
    physicalSolution=False,newProductionSolvedAreaM2=0)
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status=out['status'],steps=[{k:v for k,v in c.items() if k!='cells'} for c in checks])))
