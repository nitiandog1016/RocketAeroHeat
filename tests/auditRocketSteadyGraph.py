# -*- coding: utf-8 -*-
"""独立检查整体去时间对角、稠密Schur/五带求解及真实候选的当前模型账。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case,kind=sys.argv[1:3]
assert case.replace('_','').isalnum() and kind in ('unit','real')
target=VAL/(case+'_independent_steady_graph.json');assert not target.exists();bindings={}
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path;bindings[str(path.relative_to(ROOT))]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)
def rows(value):return [value] if isinstance(value,dict) else list(np.atleast_1d(value).ravel())
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
def scan(C,rhs,groups,n):
    result=np.zeros(5*n);order=list(range(len(groups)))+list(range(len(groups)-2,-1,-1))
    for color in order:
        cells=np.atleast_1d(groups[color]['cellIDs']).astype(int)-1;ids=np.concatenate([cells+j*n for j in range(5)]);other=np.setdiff1d(np.arange(5*n),ids)
        result[ids]=np.linalg.solve(C[np.ix_(ids,ids)],rhs[ids]-C[np.ix_(ids,other)]@result[other])
    return result
if kind=='unit':
    s=read(VAL/case/'steady_graph_unit_evidence.mat');C=s['B']-np.eye(len(s['C']))/s['original']['pseudoTimeStepSeconds'];np.testing.assert_array_equal(C,s['C'])
    source=read(Path(s['sourceFile']),['operator','newMatrix']);np.testing.assert_array_equal(source['newMatrix'],s['B'])
    np.testing.assert_array_equal(source['operator']['diagonal'],s['original']['diagonal'])
    op=s['operator'];original=s['original'];tau=original['pseudoTimeStepSeconds'];n=len(s['Q'])
    expectedD=original['diagonal'].copy()
    for i in range(5):expectedD[i,i]-=1/tau
    np.testing.assert_array_equal(op['diagonal'],expectedD)
    for key in ('lower','upper','secondLower','secondUpper','minusCoupling','plusCoupling'):np.testing.assert_array_equal(op[key],original[key])
    groups=rows(op['colorGroups']);rhs=s['rhs'].ravel(order='F');reference=scan(C,rhs,groups,n)
    assert rel(reference,s['value'].ravel(order='F'))<2e-10
    inner=s['innerReport'];kr=inner['krylov'];assert not inner['outerFixedLinearMapClaimed'] and inner['requiresFlexibleOuterMethod']
    assert kr['basisStorageUpperBoundBytes']<=kr['maximumBasisBytes'] and kr['modelApplications']<=kr['maximumModelApplications']
    checks=[]
    for cycle in rows(kr['cycles']):
        j=int(cycle['arnoldiSteps']);V=np.asarray(cycle['V']).reshape((5*n,j+1),order='F');Z=np.asarray(cycle['Z']).reshape((5*n,j),order='F');H=np.asarray(cycle['smallH']).reshape((j+1,j),order='F')
        relation=rel(C@Z,V@H);assert relation<2e-10
        for col in range(j):assert rel(scan(C,V[:,col],groups,n),Z[:,col])<2e-10
        checks.append(dict(arnoldiSteps=j,independentCArnoldiRelationDifference=relation))
    direct=np.linalg.solve(C,rhs);residual=float(np.linalg.norm(C@s['inner'].ravel(order='F')-rhs)/np.linalg.norm(rhs));assert residual<1.1e-10
    assert rel(direct,s['direct'])<2e-12 and rel(direct,s['inner'].ravel(order='F'))<2e-8
    assert s['report']['doubleUnshiftRejected']
    details=dict(scanReferenceDifference=rel(reference,s['value'].ravel(order='F')),innerCRelativeResidual=residual,cycles=checks,
        originalBMatrixKept=True,allFiveFieldsUnshifted=True,newActualNS=0)
else:
    summary=read(VAL/case/'steady_graph_real_evidence.mat');s=summary;r=read(Path(s['sourceFile']),['result'])['result'];op=r['firstOperator'];nf=int(op['sourceFaces']);n=len(r['initialQ'])
    coarse=read(Path(s['coarseFile']),['coarse'])['coarse'];np.testing.assert_array_equal(s['Q'],r['initialQ']);np.testing.assert_array_equal(r['Q'],r['initialQ'])
    np.testing.assert_array_equal(s['remaining'],s['rhs']+s['initialModel']);assert s['report']['factorReport']['allFiveConservativeTimeDiagonalsRemoved']
    def Caction(value):
        # 对角先去时间项再作用，避免两个大向量相减掩盖小项。
        D=op['diagonal'].copy()
        for a in range(5):D[a,a]-=1/op['pseudoTimeStepSeconds']
        out=np.einsum('abn,nb->na',D,value)
        for offset,key in [(-1,'lower'),(1,'upper'),(-2,'secondLower'),(2,'secondUpper')]:
            if offset<0:
                start=-offset*nf;out[start:]+=np.einsum('abn,nb->na',op[key][:,:,start:],value[:-start])
            else:
                end=offset*nf;out[:-end]+=np.einsum('abn,nb->na',op[key][:,:,:-end],value[end:])
        m,p=op['sideMinus'].astype(int)-1,op['sidePlus'].astype(int)-1
        for first in range(0,len(m),8192):
            ids=slice(first,first+8192);np.add.at(out,m[ids],np.einsum('abn,nb->na',op['minusCoupling'][:,:,ids],value[p[ids]]));np.add.at(out,p[ids],np.einsum('abn,nb->na',op['plusCoupling'][:,:,ids],value[m[ids]]))
        return out
    checks=[];remaining=s['remaining'];rhs=s['rhs']
    for item in rows(s['evidence']):
        Cvalue=Caction(item['value']);assert rel(Cvalue,item['Cvalue'])<2e-10
        Avalue=item['Avalue'];model=item['totalAction'];right=coarse['image'].T@(remaining-Avalue).ravel(order='F')
        weights=np.linalg.solve(coarse['upper'],np.linalg.solve(coarse['upper'].T,right));prediction=Avalue+(coarse['image']@weights).reshape((n,5),order='F')
        assert rel(prediction,model)<2e-9;np.testing.assert_array_equal(item['total'],item['value']+item['correction'])
        alpha=float(np.sum(remaining*model)/np.sum(model**2));best=remaining-alpha*model
        values=dict(relativeCResidual=float(np.linalg.norm(remaining-Cvalue)/np.linalg.norm(remaining)),
            fullADirectionRelativeResidual=float(np.linalg.norm(remaining-Avalue)/np.linalg.norm(remaining)),
            afterCoarseRelativeResidual=float(np.linalg.norm(remaining-model)/np.linalg.norm(remaining)),bestScalarAfterCoarse=alpha,
            predictedTotalDefectAfterBestScalar=float(np.linalg.norm(best)/np.linalg.norm(rhs)),
            remainingSquaredNormRemovedByBestScalar=float(1-np.linalg.norm(best)**2/np.linalg.norm(remaining)**2),
            steadyGraphToFullAActionRelativeDifference=rel(Cvalue,Avalue))
        for key,value in values.items():np.testing.assert_allclose(value,item['summary'][key],rtol=3e-10,atol=2e-12)
        assert not item['summary']['originalNSAcceptanceClaimed'];checks.append(dict(name=item['summary']['name'],**values))
    assert s['report']['newSelectedJApplications']==5 and s['report']['coarseReuseProbeModelApplications']==1
    assert s['report']['newActualNSEvaluations']==s['report']['newFlowAcceptances']==s['report']['newCoarseConstructionApplications']==0
    details=dict(candidates=checks,allFiveFieldsUnshifted=True,physicalEquationChanged=False,newActualNS=0)
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-all-five-unshifted-steady-graph-'+kind+'-current-reference-and-work-audit',case=case,kind=kind,
    processID=os.getpid(),details=details,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out if kind=='unit' else dict(status=out['status'],processID=os.getpid(),details=details)))
