# -*- coding: utf-8 -*-
"""独立核对父C候选的增量/总方向分母、模型算术及条件原NS差分证据。"""
from pathlib import Path
import os,json,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,VAL,sha

name=sys.argv[1];assert name.replace('_','').isalnum();folder=VAL/name
target=VAL/(name+'_independent_parent_direction_probe.json');assert not target.exists();bindings={}


def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)


def rel(a,b): return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(b)))


s=read(folder/'parent_direction_probe_evidence.mat')
r=s['report'];parent=read(s['parentFile'],['value','rhs','detail'])
source=read(VAL/'innergraphreal01/inner_graph_real_evidence.mat',['Q','scale','rhs','direction'])
assert sha(r['sourceBinding']['sourceFile'])==r['sourceBinding']['sourceFileSHA256']
assert sha(s['parentFile'])==r['parentFileSHA256'] and sha(s['coarseFile'])==r['coarseFileSHA256']
np.testing.assert_array_equal(s['Q'],source['Q']);np.testing.assert_array_equal(s['scale'],source['scale'])
np.testing.assert_array_equal(s['rhs'],source['rhs']);np.testing.assert_array_equal(s['direction0'],source['direction'])
np.testing.assert_array_equal(s['delta'],parent['value']);np.testing.assert_array_equal(s['parentRHS'],parent['rhs'])
assert rel(s['rhs'],s['baseRate']/s['scale'])<1e-15
remaining=s['rhs']+s['initialModel'];assert rel(remaining,s['remaining'])<1e-15
assert rel(remaining,parent['rhs'])<1e-8
np.testing.assert_array_equal(s['rawTotal'],s['direction0']+s['delta'])
assert rel(s['rawModel'],s['initialModel']+s['deltaModel'])<1e-8
assert rel(s['totalDirection'],s['direction0']+r['chosenScalar']*s['delta'])<1e-15
assert rel(s['totalModel'],s['initialModel']+r['chosenScalar']*s['deltaModel'])<1e-8
raw=float(np.linalg.norm(s['rhs']+s['rawModel'])/np.linalg.norm(s['rhs']))
chosen=float(np.linalg.norm(s['rhs']+s['totalModel'])/np.linalg.norm(s['rhs']))
alpha=-float(np.sum(remaining*s['deltaModel'])/np.sum(s['deltaModel']**2))
assert abs(raw-r['rawTotalCurrentModelDefect'])<1e-12 and abs(chosen-r['chosenTotalCurrentModelDefect'])<1e-12
assert abs(alpha-r['optimalScalar'])<1e-12
ns_checked=bool(r['originalNSDirectionChecked']);details={}
if ns_checked:
    e=s['nsEvidence'];h=float(e['actualScalarStep']);scale=s['scale']
    np.testing.assert_array_equal(e['directionU'],s['totalDirection'])
    jv=(e['ratePlus']-e['rateMinusOrBase'])/scale/h
    if e['scheme']=='central': jv/=2
    assert rel(jv,e['value'])<1e-12
    actual=float(np.linalg.norm(s['rhs']+jv)/np.linalg.norm(s['rhs']))
    slope=float(np.sum(s['rhs']*jv)/np.sum(s['rhs']**2))
    assert abs(actual-r['originalNSTotalDirectionDefect'])<1e-12 and abs(slope-r['originalNSNormalizedMeritSlope'])<1e-12
    state=s['Q']+h*s['totalDirection']*scale;rho=state[:,0];velocity=state[:,1:4]/rho[:,None]
    energy=state[:,4]-.5*rho*np.sum(velocity**2,axis=1);temperature=energy/(rho*s['gas']['R']/(s['gas']['gamma']-1))
    assert np.all(rho>0) and np.all(energy>0) and np.all(temperature>=80) and np.all(temperature<=2000)
    assert bool(r['originalNSDirectionTargetMet'])==(actual<=.1 and slope<0)
    assert r['newActualNSEvaluations']==2+(e['scheme']=='central')
    details=dict(actualOriginalNSTotalDirectionDefect=actual,normalizedMeritSlope=slope,
        originalNSDirectionTargetMet=bool(r['originalNSDirectionTargetMet']),plusStateTemperatureRangeK=[float(temperature.min()),float(temperature.max())],
        completeOriginalNSFluxAuditClaimed=False)
else:
    assert not r['originalNSDirectionTargetMet']
    if chosen>.1: assert r['newActualNSEvaluations']==0
assert not r['physicalStateChanged'] and r['newFlowAcceptances']==r['newProductionSolvedAreaM2']==0
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-parent-increment-total-direction-denominator-and-conditional-NS-arithmetic',processID=os.getpid(),
    rawCurrentModelDefect=raw,chosenCurrentModelDefect=chosen,optimalScalar=alpha,
    sameQAndFrozenBestDirectionExactlyVerified=True,parentCStillNotPromotedToNSQualification=True,
    originalNSDirectionChecked=ns_checked,conditionalOriginalNSCheck=details,
    newActualNSEvaluations=int(r['newActualNSEvaluations']),newFlowAcceptances=0,newProductionSolvedAreaM2=0,
    inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
