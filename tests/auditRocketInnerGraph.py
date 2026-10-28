# -*- coding: utf-8 -*-
"""独立核对内部B求解、当前模型预测、时间项拆分及存储/作用预算。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case,kind=sys.argv[1:3]
assert case.replace('_','').isalnum() and kind in ('unit','real')
target=VAL/(case+'_independent_inner_graph.json');assert not target.exists();bindings={}
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path
    bindings[str(path.relative_to(ROOT))]=sha(path);return loadmat(path,simplify_cells=True,variable_names=variables)
def convert(value):
    if hasattr(value,'_fieldnames'):return {key:convert(getattr(value,key)) for key in value._fieldnames}
    if isinstance(value,dict):return {key:convert(item) for key,item in value.items()}
    if isinstance(value,np.ndarray) and value.dtype==object:
        out=np.empty(value.shape,dtype=object)
        for index in np.ndindex(value.shape):out[index]=convert(value[index])
        return out
    return value
def rows(value):
    value=convert(value)
    if isinstance(value,dict):return [value]
    return list(np.atleast_1d(value).ravel())
def relative(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
def spectrum(matrix):
    singular=np.linalg.svd(matrix,compute_uv=False);threshold=max(matrix.shape)*np.finfo(float).eps*singular[0]
    return dict(minimumSingularValue=float(singular[-1]),maximumSingularValue=float(singular[0]),
        conditionNumber=float(singular[0]/singular[-1]) if singular[-1]>0 else None,
        numericalRank=int(np.count_nonzero(singular>threshold)),rankThreshold=float(threshold))
def check_work(detail):
    r=detail['krylov'];assert not detail['outerFixedLinearMapClaimed'] and detail['requiresFlexibleOuterMethod']
    assert not detail['steadyNewtonEquationSolved'] and not detail['originalNSAcceptanceClaimed']
    assert r['options']['maximumCycles']==1 and r['options']['retainedDirections']==0
    assert r['basisStorageUpperBoundBytes']<=r['maximumBasisBytes'] and r['modelApplications']<=r['maximumModelApplications']
    assert r['preconditionerApplications']==r['arnoldiSteps']
    assert detail['totalBApplicationsIncludingScanDiagnostics']==r['modelApplications']+r['preconditionerApplications']
    assert detail['fullSelectedJApplications']==detail['newActualNSEvaluations']==0
if kind=='unit':
    s=read(VAL/case/'inner_graph_unit_evidence.mat');B=s['B'];P=s['P'];A=s['steadyA'];checks=[]
    source=read(Path(s['sourceFile']),['newMatrix','unitA']);np.testing.assert_array_equal(B,source['newMatrix']);np.testing.assert_array_equal(A,source['unitA'])
    for item in rows(s['evidence']):
        item=convert(item);b=item['rhs'].ravel(order='F');x=item['value'].ravel(order='F');detail=item['detail'];check_work(detail);r=detail['krylov']
        direct=np.linalg.solve(B,b);assert relative(direct,item['directBReference'])<2e-13
        measured=float(np.linalg.norm(b-B@x)/np.linalg.norm(b));assert measured<1.1e-10
        np.testing.assert_allclose(measured,detail['relativeBResidual'],rtol=1e-3,atol=1e-13)
        maximum_relation=0.;maximum_residual_formation=0.
        for c in rows(r['cycles']):
            j=int(c['arnoldiSteps']);V=np.asarray(c['V']).reshape((len(b),j+1),order='F');Z=np.asarray(c['Z']).reshape((len(b),j),order='F')
            H=np.asarray(c['smallH']).reshape((j+1,j),order='F');weights=np.atleast_1d(c['weights']);maximum_relation=max(maximum_relation,relative(B@Z,V@H))
            assert maximum_relation<2e-11 and relative(P@V[:,:j],Z)<2e-11
            np.testing.assert_allclose(c['candidate'],c['xStart']+Z@weights,rtol=1e-12,atol=1e-12)
            # 两个大项相减形成很小的剩余，独立形成误差应按原右端归一化。
            formation=float(np.linalg.norm(c['actualRemaining']-(b-B@c['candidate']))/max(1,np.linalg.norm(b)))
            maximum_residual_formation=max(maximum_residual_formation,formation);assert formation<2e-12
            independent=np.linalg.lstsq(B@Z,c['rStart'],rcond=1e-13)[0]
            assert np.linalg.norm(B@Z@(weights-independent))/max(1,np.linalg.norm(b))<2e-11
        checks.append(dict(BRelativeResidual=measured,directBReferenceDifference=relative(x,direct),maximumArnoldiRelationDifference=maximum_relation,
            maximumResidualFormationDifferenceNormalizedByRhs=maximum_residual_formation,
            completeSteadyARelativeResidualOfBInverse=float(np.linalg.norm(b-A@x)/np.linalg.norm(b)),
            arnoldiSteps=int(r['arnoldiSteps']),BApplicationsIncludingScanDiagnostics=int(detail['totalBApplicationsIncludingScanDiagnostics'])))
    assert s['report']['budgetRejected'] and np.all(s['zero']==0) and s['zeroReport']['totalBApplicationsIncludingScanDiagnostics']==0
    Binverse=np.linalg.solve(B,np.eye(len(B)))
    details=dict(rightHandSides=checks,spectra=dict(steadyA=spectrum(A),approximateB=spectrum(B),
        BPscan=spectrum(B@P),steadyAExactBInverse=spectrum(A@Binverse)),
        oneScanInverseResidualNorm=float(np.linalg.norm(np.eye(len(B))-B@P)),
        exactBInverseIsNotASteadyInverse=True,newActualNS=0)
else:
    s=read(VAL/case/'inner_graph_real_evidence.mat');r=read(Path(s['sourceFile']),['result'])['result'];op=r['firstOperator'];n=len(r['initialQ']);nf=int(op['sourceFaces'])
    np.testing.assert_array_equal(s['Q'],r['initialQ']);np.testing.assert_array_equal(r['Q'],r['initialQ']);np.testing.assert_array_equal(s['scale'],r['residualScale'])
    log=rows(r['linearLog'])[-1];np.testing.assert_array_equal(s['direction'],log['steadyDirectionEvidence']['directionU'])
    np.testing.assert_array_equal(s['savedNSValue'],log['steadyDirectionEvidence']['value']);np.testing.assert_array_equal(s['remaining'],s['rhs']+s['initialModel'])
    assert relative(s['initialModel'],s['savedNSValue'])<2e-6
    coarse=read(Path(s['coarseFile']),['coarse'])['coarse'];assert coarse['imageBytes']<=256*1024**2
    np.testing.assert_array_equal(coarse['baseQ'],s['Q']);check_work(s['innerReport'])
    def Baction(value):
        out=np.einsum('abn,nb->na',op['diagonal'],value)
        for offset,key in [(-1,'lower'),(1,'upper'),(-2,'secondLower'),(2,'secondUpper')]:
            if offset<0:
                start=-offset*nf;out[start:]+=np.einsum('abn,nb->na',op[key][:,:,start:],value[:-start])
            else:
                end=offset*nf;out[:-end]+=np.einsum('abn,nb->na',op[key][:,:,:-end],value[end:])
        m,p=op['sideMinus'].astype(int)-1,op['sidePlus'].astype(int)-1
        for first in range(0,len(m),8192):
            ids=slice(first,first+8192)
            np.add.at(out,m[ids],np.einsum('abn,nb->na',op['minusCoupling'][:,:,ids],value[p[ids]]))
            np.add.at(out,p[ids],np.einsum('abn,nb->na',op['plusCoupling'][:,:,ids],value[m[ids]]))
        return out
    remaining=s['remaining'];rhs=s['rhs'];items=rows(s['evidence']);checks=[]
    for item in items:
        item=convert(item);value=item['value'];Bx=Baction(value);Ax=item['Avalue'];summary=item['summary'];temporal=value/op['pseudoTimeStepSeconds']
        assert relative(Bx,item['Bvalue'])<2e-12
        coarse_rhs=(remaining-Ax).ravel(order='F');right=coarse['image'].T@coarse_rhs
        weights=np.linalg.solve(coarse['upper'],np.linalg.solve(coarse['upper'].T,right))
        np.testing.assert_allclose(weights,item['coarseWeights'],rtol=2e-7,atol=1e-8)
        predicted=Ax+(coarse['image']@weights).reshape((n,5),order='F');assert relative(predicted,item['totalAction'])<1e-9
        np.testing.assert_array_equal(item['total'],value+item['correction'])
        alpha=float(np.sum(remaining*item['totalAction'])/np.sum(item['totalAction']**2));best=remaining-alpha*item['totalAction']
        values=dict(relativeBResidual=float(np.linalg.norm(remaining-Bx)/np.linalg.norm(remaining)),
            fullADirectionRelativeResidual=float(np.linalg.norm(remaining-Ax)/np.linalg.norm(remaining)),
            afterCoarseRelativeResidual=float(np.linalg.norm(remaining-item['totalAction'])/np.linalg.norm(remaining)),
            bestScalarAfterCoarse=alpha,predictedTotalDefectAfterBestScalar=float(np.linalg.norm(best)/np.linalg.norm(rhs)),
            remainingSquaredNormRemovedByBestScalar=float(1-np.linalg.norm(best)**2/np.linalg.norm(remaining)**2),
            temporalActionNormRelativeToRemaining=float(np.linalg.norm(temporal)/np.linalg.norm(remaining)),
            temporalProjectionOnRemaining=float(np.sum(remaining*temporal)/np.sum(remaining**2)),
            graphAToFullAActionRelativeDifference=relative(Bx-temporal,Ax))
        for key,val in values.items():np.testing.assert_allclose(val,summary[key],rtol=3e-10,atol=2e-12)
        assert not summary['originalNSAcceptanceClaimed'];checks.append(dict(name=item['name'],**values))
    report=s['report'];assert report['newSelectedJApplications']==5 and report['coarseReuseProbeModelApplications']==1
    assert report['newActualNSEvaluations']==report['newFlowAcceptances']==report['newCoarseConstructionApplications']==0
    assert not report['newPhysicalState']
    details=dict(directions=checks,innerBApplications=int(s['innerReport']['totalBApplicationsIncludingScanDiagnostics']),
        newSelectedModelApplications=5,coarseReuseProbeApplications=1,newActualNS=0,currentModelDefect=float(report['initialCurrentModelDefect']))
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-inner-B-'+kind+'-Arnoldi-action-projection-and-work-audit',case=case,kind=kind,processID=os.getpid(),details=details,
    fullIndependentNSSolver=False,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out if kind=='unit' else dict(status=out['status'],processID=os.getpid(),details=details)))
