# -*- coding: utf-8 -*-
"""在同一保存总方向上比较旧图/五带图与原NS差分作用，不执行新模型或NS。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case,parent=sys.argv[1:3]
assert all(x.replace('_','').isalnum() for x in (case,parent)) and case!=parent
target=VAL/case/'wls_normal_action_diagnostic.json';assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
files=[VAL/name/'wall_velocity_krylov_evidence.mat' for name in (case,parent)];bindings={str(p.relative_to(ROOT)):sha(p) for p in files}
r=loadmat(files[0],simplify_cells=True,variable_names=['result'])['result'];old=loadmat(files[1],simplify_cells=True,variable_names=['result'])['result']
np.testing.assert_array_equal(r['initialQ'],old['Q']);assert r['acceptedNonlinearIterations']==0
np.testing.assert_array_equal(r['Q'],old['Q']);np.testing.assert_array_equal(r['rate'],old['rate'])
logs=r['linearLog'];a=logs if isinstance(logs,dict) else list(logs)[-1];e=a['steadyDirectionEvidence'];direction=e['directionU'];full=e['value']
rhs=r['initialRate']/r['residualScale'];remaining=rhs+full;nf,nl=int(r['firstOperator']['sourceFaces']),int(r['firstOperator']['layers'])
def action(op,second):
    out=np.einsum('abn,nb->na',op['diagonal'],direction)
    for offset,key in [(-1,'lower'),(1,'upper')]+([(-2,'secondLower'),(2,'secondUpper')] if second else []):
        if offset<0:
            start=-offset*nf;out[start:]+=np.einsum('abn,nb->na',op[key][:,:,start:],direction[:-start])
        else:
            end=offset*nf;out[:-end]+=np.einsum('abn,nb->na',op[key][:,:,:-end],direction[end:])
    m,p=op['sideMinus'].astype(int)-1,op['sidePlus'].astype(int)-1
    for first in range(0,len(m),8192):
        ids=slice(first,first+8192)
        np.add.at(out,m[ids],np.einsum('abn,nb->na',op['minusCoupling'][:,:,ids],direction[p[ids]]))
        np.add.at(out,p[ids],np.einsum('abn,nb->na',op['plusCoupling'][:,:,ids],direction[m[ids]]))
    return direction/op['pseudoTimeStepSeconds']-out
before,after=action(old['firstOperator'],False),action(r['firstOperator'],True);details=[]
for name,value in [('originalGraph',before),('WlsNormalFiveBand',after)]:
    difference=value-full;den=float(np.sum(difference*difference))
    details.append(dict(name=name,relativeDifferenceToSavedOriginalNSDirection=float(np.linalg.norm(difference)/np.linalg.norm(full)),
        equationDifferenceRelativeToSavedOriginalNSDirection=(np.linalg.norm(difference,axis=0)/np.linalg.norm(full,axis=0)).tolist(),
        layerSquaredDifferenceFractions=(np.sum(difference.reshape((nl,nf,5))**2,axis=(1,2))/den).tolist(),
        approximateSteadyLinearDefect=float(np.linalg.norm(rhs+value)/np.linalg.norm(rhs))))
mass_full=float(np.linalg.norm(full[:,0]));mass_difference=float(np.linalg.norm(after[:,0]-full[:,0]))
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='source-bound-same-direction-original-graph-five-band-and-original-NS-action-diagnostic',processID=os.getpid(),details=details,
    actualNSDefect=float(a['fullSystemRelativeDefect']),remainingSquaredNormEquationFractions=(np.sum(remaining*remaining,axis=0)/np.sum(remaining*remaining)).tolist(),
    massActionRelativeDifference=mass_difference/mass_full,NSDirectionalImageIsFiniteDifference=True,
    approximateActionDifferenceIsNotFullJacobianConditioningOrUniqueStagnationCause=True,newModelApplications=0,newActualNSEvaluations=0,
    newFlowUpdates=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
