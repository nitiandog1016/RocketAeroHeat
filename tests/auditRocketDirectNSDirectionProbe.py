# -*- coding: utf-8 -*-
"""独立核对原NS直接初猜的基态/完整几何/来源及两次新鲜方向检查。"""
from pathlib import Path
import json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from riemannReservoirReference import ROOT,VAL,sha

case,kind=sys.argv[1:3];assert case.replace('_','').isalnum() and kind in ('unit','result-unit','real')
target=VAL/(case+'_independent_direct_ns_probe.json');assert not target.exists();bindings={}

def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)

def equal(a,b):
    if isinstance(a,dict) or isinstance(b,dict):return isinstance(a,dict) and isinstance(b,dict) and a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    return np.array_equal(np.asarray(a),np.asarray(b))

file=VAL/case/('direct_ns_probe_unit_evidence.mat' if kind=='unit' else 'result_source_probe_unit_evidence.mat' if kind=='result-unit' else 'wall_velocity_krylov_evidence.mat')
s=read(file);r=s['result'];config=r['options']['steadyInitialGuess'];assert config['variant']=='ns-probe'
p=read(config['file']);source=p['sourceBinding'];geometry=p['geometryBinding']
assert sha(source['file'])==source['sha256'] and sha(geometry['file'])==geometry['sha256']
location=source.get('stateLocation','root');assert location in ('root','result')
if location=='root':b=read(source['file'],['Q','rate','gas','boundary'])
else:
    container=read(source['file'],['result','gas','boundary'])
    b=dict(Q=container['result']['Q'],rate=container['result']['rate'],gas=container['gas'],boundary=container['boundary'])
d=read(geometry['file'],[geometry['variable']])
assert equal(b['Q'],p['Q']) and equal(b['rate'],p['baseRate']) and equal(b['gas'],p['gas']) and equal(b['boundary'],p['boundary'])
assert equal(p['Q'],r['initialQ']) and equal(p['baseRate'],r['initialRate']) and equal(p['scale'],r['residualScale'])
assert equal(p['gas'],s['gas']) and equal(p['boundary'],s['boundary'])
if kind in ('unit','result-unit'):actual_geometry=s['cut']
else:actual_geometry=read(ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat',['domain'])['domain']
assert equal(d[geometry['variable']],actual_geometry)
for name,digest in p['sourceIdentity'].items():
    paths=list((ROOT/'internalFlow').rglob(name+'.m'));assert len(paths)==1 and sha(paths[0])==digest
    bindings[paths[0].relative_to(ROOT).as_posix()]=digest
e=p['nsEvidence'];np.testing.assert_array_equal(p['directionU'],e['directionU']);rhs=p['baseRate']/p['scale']
factor=2 if e['scheme']=='central' else 1
jv=(e['ratePlus']-e['rateMinusOrBase'])/p['scale']/(factor*e['actualScalarStep'])
np.testing.assert_array_equal(jv,e['value'])
probe_defect=float(np.linalg.norm(rhs+jv)/np.linalg.norm(rhs));assert probe_defect<=.1
assert abs(probe_defect-p['report']['originalNSDirectionDefect'])<1e-12
assert not p['report']['physicalStateChanged'] and p['report']['newFlowAcceptances']==0
rows=lambda x:[x] if isinstance(x,dict) else list(x)
a=rows(r['linearLog'])[0];guess=a['initialGuessEvidence'];fresh=a['steadyDirectionEvidence']
np.testing.assert_array_equal(guess['directionU'],p['directionU']*config['scale'])
np.testing.assert_array_equal(fresh['directionU'],guess['directionU']+guess['correctionDirectionU'])
assert a['gmresSkipped'] and np.all(guess['correctionDirectionU']==0)
np.testing.assert_array_equal(guess['currentModelValue'],fresh['value'])
source_report=guess['sourceReport'];assert source_report['sourceStateAndGeometryVerified'] and not source_report['oldDerivativeReused']
assert not source_report['coarseImageRequired'] and not source_report['originalNSAcceptanceClaimed']
assert sha(config['file'])==source_report['sourceFileSHA256']
fresh_factor=2 if fresh['scheme']=='central' else 1
fresh_jv=(fresh['ratePlus']-fresh['rateMinusOrBase'])/r['residualScale']/(fresh_factor*fresh['actualScalarStep'])
np.testing.assert_array_equal(fresh_jv,fresh['value'])
fresh_defect=float(np.linalg.norm(rhs+fresh_jv)/np.linalg.norm(rhs));assert fresh_defect<=a['target']==.1
assert abs(fresh_defect-a['fullSystemRelativeDefect'])<1e-12
actions=[j for j in rows(r['jvpLog']) if j.get('initialGuessModelAction',False)]
assert len(actions)==1 and actions[0]['actualOriginalResidualEvaluations']==1 and actions[0]['scheme']=='forward'
bad_checks=[]
if kind=='unit':
    assert s['report']['oldCompleteDefaultResultReplayed'] and s['report']['sameAcceptedStateAsDefault'] and np.all(s['report']['rejected'])
    assert r['acceptedNonlinearIterations']==1 and r['residualEvaluations']==4
    for i in range(1,9):
        bad=read(VAL/case/('bad_probe_%02d.mat'%i))
        violations=[not equal(bad['Q'],p['Q']),not equal(bad['scale'],p['scale']),not equal(bad['boundary'],p['boundary']),
                    bad['geometryBinding']['sha256']!=p['geometryBinding']['sha256'],not equal(bad['directionU'],p['directionU']),
                    not equal(bad['nsEvidence']['value'],p['nsEvidence']['value']),not equal(bad['sourceIdentity'],p['sourceIdentity']),
                    bad['sourceBinding']['sha256']!=p['sourceBinding']['sha256']]
        assert violations[i-1];bad_checks.append(i)
elif kind=='result-unit':
    assert location=='result' and s['report']['rootAndResultStatesActuallyDiffer'] and s['report']['oldRootProbeExactlyLoaded']
    assert np.all(s['report']['locationRejections']) and s['report']['sameAcceptedStateAsFreshDefault']
    np.testing.assert_array_equal(r['Q'],s['reference']['Q']);np.testing.assert_array_equal(r['rate'],s['reference']['rate'])
    assert r['acceptedNonlinearIterations']==1
    for i in range(1,4):
        bad=read(VAL/case/('bad_location_%d.mat'%i));binding=bad['sourceBinding']
        assert (binding['stateLocation']=='root' if i==1 else binding['stateLocation']=='automatic' if i==2 else binding['file']!=source['file'])
        bad_checks.append(i)
else:
    assert not s['report']['sourceHistoryAppended']
    if not s['report']['sameEquationHistoryContinued']:assert s['report']['totalNewEquationResidualEvaluations']==r['residualEvaluations']
    else:assert s['report']['totalNewEquationResidualEvaluations']>=r['residualEvaluations'] and s['report']['priorAcceptedNonlinearIterations']>=1
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='passed-independent-direct-NS-probe-current-base-geometry-sources-and-two-fresh-checks',case=case,kind=kind,
    processID=os.getpid(),sourceProbeDefect=probe_defect,freshMainNSDefect=fresh_defect,sourceGuardEvidenceChecked=bad_checks,
    coarseImageRequired=False,sourceStateLocation=location,acceptedFlowUpdates=int(r['acceptedNonlinearIterations']),newProductionSolvedAreaM2=0,
    inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
