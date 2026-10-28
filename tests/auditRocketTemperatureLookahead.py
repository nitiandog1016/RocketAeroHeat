"""独立核对当前候选温度、全部带内面、原H保留及公开Riemann温度导数。"""
from pathlib import Path
import os,sys,json
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,sha
from auditRocketTemperatureProjection import row_reference


def check(Q,info,cut,gas,scale,value,rows,report):
    alpha=report['targetAlpha']; q=Q+alpha*value*scale; rho=q[:,0]; u=q[:,1:4]/rho[:,None]
    T=(q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*gas['R']/(gas['gamma']-1))
    assert np.all(rho>0) and np.all(T>0)
    low,high=gas['temperatureRange'];bad=np.flatnonzero((T<low)|(T>high))
    np.testing.assert_array_equal(bad+1,np.atleast_1d(report['predictedViolatingCellIDs']))
    audit=info['exactEulerFluxAudit'];ids=np.atleast_1d(audit['faceID']).astype(int)
    sample=np.asarray(audit['samplePrimitive']);margin=np.minimum(sample[:,4]-low,high-sample[:,4])
    required=np.atleast_1d(report['requiredOriginalStateFaceIDs']).astype(int)
    slots=np.searchsorted(ids,required);np.testing.assert_array_equal(ids[slots],required)
    region=np.atleast_1d(audit['sampleRegion'])[slots]; assert set(region.tolist())<={1,4}
    minus=np.asarray(cut['openFaces']['cellMinus']).astype(int);plus=np.asarray(cut['openFaces']['cellPlus']).astype(int)
    cv=np.where(region==1,minus[required-1],plus[required-1])-1
    np.testing.assert_array_equal(cv,bad)
    K=max(10.,float(np.max(margin[slots]))); assert K==report['selectedTemperatureBandK']
    np.testing.assert_array_equal(np.atleast_1d(rows['faceIDs']),ids[margin<=K])
    np.testing.assert_array_equal(np.atleast_1d(report['originalBandFaceIDs']).astype(int),ids[margin<=10])
    assert report['originalRowsExactlyRetained'] and not report['physicalBoundaryConditionAdded']
    assert not report['nonlinearTemperatureBoundGuaranteed'] and report['newActualNSEvaluations']==0
    # 两个细化步长都须通过原1e-7门槛；不把首轮粗步长失败改写成通过。
    public=row_reference(rows,3e-7);refined=row_reference(rows,1e-7)
    return dict(violatingCellIDs=(bad+1).tolist(),originalBandFaces=int(np.sum(margin<=10)),
                selectedBandFaces=int(np.sum(margin<=K)),selectedBandK=K,
                allRequiredFacesAreCurrentOriginalStates=True,publicTemperatureDerivative=public,
                finerPublicTemperatureDerivative=refined,finiteDifferenceStepFactors=[3e-7,1e-7],
                originalDerivativeThresholdUnchanged=1e-7)


folder=Path(sys.argv[1]).resolve();file=folder/'temperature_lookahead_unit.mat';before=sha(file)
s=loadmat(file,simplify_cells=True);cases=s['cases']; checks=[]
for c in cases:
    checks.append(check(c['Q'],c['info'],c['cut'],c['gas'],c['scale'],c['value'],c['rows'],c['detail']))
    base=c['base'];rows=c['rows'];old_ids=np.atleast_1d(base['faceIDs']).astype(int);new_ids=np.atleast_1d(rows['faceIDs']).astype(int)
    keep=np.searchsorted(new_ids,old_ids); assert (base['H']-rows['H'][keep]).nnz==0
    np.testing.assert_array_equal(np.atleast_1d(base['temperatureBoundSigns']),np.atleast_1d(rows['temperatureBoundSigns'])[keep])
    assert c['negativeDensityPreviewRetainsOriginalRows'] and c['zeroDirectionRetainsOriginalRows']
source=ROOT/'outputs/v83_validation/krylov35/wall_velocity_krylov_evidence.mat'
next_file=ROOT/'outputs/v83_validation/krylov36/wall_velocity_krylov_evidence.mat'
geometry=ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat'
bindings={str(p.relative_to(ROOT)):sha(p) for p in (source,next_file,geometry)}
a=loadmat(source,simplify_cells=True,variable_names=['result','gas']);b=loadmat(next_file,simplify_cells=True,variable_names=['result'])
d=loadmat(geometry,simplify_cells=True,variable_names=['domain'])['domain'];trace=b['result']['trace']
trace=trace if isinstance(trace,dict) else list(trace)[0]
np.testing.assert_array_equal(a['result']['Q'],trace['initialQ'])
real=check(a['result']['Q'],a['result']['info'],d,a['gas'],b['result']['residualScale'],trace['linearUpdateU'],s['realRows'],s['realReport'])
assert sha(file)==before and all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-temperature-preview-all-rows-pure-face-identity-and-public-derivatives',
    manufacturedCases=checks,realSavedQ35Fixture=real,inputSHA256={**bindings,str(file.relative_to(ROOT)):before},
    scriptSHA256=sha(__file__),newActualNSEvaluations=0,newFlowAcceptances=0,physicalSolution=False,newProductionSolvedAreaM2=0)
out['derivativeReferenceSHA256']=sha(ROOT/'tests/auditRocketTemperatureProjection.py')
out['initialFailedAuditRetained']=str(folder/'independent_audit_initial_failure.json')
out['finiteDifferenceStepDiagnosisSHA256']=sha(folder/'FD_diagnosis.json')
target=folder/'independent_temperature_lookahead.json';assert not target.exists()
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(status=out['status'],manufacturedCases=len(checks),realFixture=real)))
