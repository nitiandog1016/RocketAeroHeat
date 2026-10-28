# -*- coding: utf-8 -*-
"""独立核对同线限制、五带矩阵/着色直接解、原NS方向与主链接入身份。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case,kind=sys.argv[1:3]
assert case.replace('_','').isalnum() and kind in ('kernel','integration','real','blocks')
target=VAL/(case+'_independent_wls_normal.json');assert not target.exists();bindings={}
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def read(path,variables=None):
    path=Path(path);bindings[str(path.relative_to(ROOT))]=sha(path);return loadmat(path,simplify_cells=True,variable_names=variables)
def rows(value):return [value] if isinstance(value,dict) else list(np.atleast_1d(value))
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
def matrix(op,second):
    nf,nl=int(op['sourceFaces']),int(op['layers']);n=nf*nl;M=np.zeros((5*n,5*n))
    for cell in range(n):
        ids=cell+np.arange(5)*n;M[np.ix_(ids,ids)]=op['diagonal'][:,:,cell]
        for offset,key in [(-1,'lower'),(1,'upper')]+([(-2,'secondLower'),(2,'secondUpper')] if second else []):
            other=cell+offset*nf
            if 0<=other<n:M[np.ix_(ids,other+np.arange(5)*n)]=op[key][:,:,cell]
    for j,(m,p) in enumerate(zip(op['sideMinus'].astype(int)-1,op['sidePlus'].astype(int)-1)):
        a,b=m+np.arange(5)*n,p+np.arange(5)*n
        M[np.ix_(a,b)]+=op['minusCoupling'][:,:,j];M[np.ix_(b,a)]+=op['plusCoupling'][:,:,j]
    return M
def action(op,value,second):
    nf=int(op['sourceFaces']);out=np.einsum('abn,nb->na',op['diagonal'],value)
    for offset,key in [(-1,'lower'),(1,'upper')]+([(-2,'secondLower'),(2,'secondUpper')] if second else []):
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
details={}
if kind=='kernel':
    s=read(VAL/case/'wls_normal_unit_evidence.mat');n=len(s['Q']);nf=int(s['graph']['sourceFaces']);nl=int(s['graph']['layers'])
    old=matrix(s['graph'],False);new=matrix(s['operator'],True)
    assert rel(old,s['oldMatrix'])<2e-13 and rel(new,s['newMatrix'])<2e-13
    reference=read(Path(s['sourceFile']),['unitA'])['unitA'];np.testing.assert_array_equal(reference,s['unitA'])
    cells=np.tile(np.arange(n),5);same=(cells[:,None]%nf)==(cells[None,:]%nf)
    difference=s['fullViscous']-s['graphViscous'];restricted=difference*same
    extraction=rel(old-new,restricted);assert extraction<2e-10
    far=abs(cells[:,None]//nf-cells[None,:]//nf)>2;assert np.linalg.norm(restricted*far)<2e-12*max(1,np.linalg.norm(restricted))
    expected=np.zeros(5*n);groups=rows(s['operator']['colorGroups']);order=list(range(len(groups)))+list(range(len(groups)-2,-1,-1));rhs=s['rhs'].ravel(order='F')
    for sweep in range(int(s['operator']['fixedSymmetricSweeps'])):
        for color in order:
            selected=np.atleast_1d(groups[color]['cellIDs']).astype(int)-1
            selected=np.concatenate([selected+j*n for j in range(5)]);other=np.setdiff1d(np.arange(5*n),selected)
            expected[selected]=np.linalg.solve(new[np.ix_(selected,selected)],rhs[selected]-new[np.ix_(selected,other)]@expected[other])
    scan=rel(s['update'].ravel(order='F'),expected);assert scan<2e-11
    assert rel(action(s['operator'],s['update'],True).ravel(order='F'),new@s['update'].ravel(order='F'))<2e-13
    assert rel(s['combined'],.7*s['update']-1.3*s['other'])<2e-11
    errors=[]
    for direction,plus,minus,model in zip(s['directions'],s['fdPlus'],s['fdMinus'],s['fdModels']):
        measured=(plus-minus)/(2*s['h'])/s['scale'];errors.append(rel(measured,model));assert errors[-1]<2e-6
        assert rel(model.ravel(order='F'),-reference@direction.ravel(order='F'))<2e-13
    difference2=float(np.sum(difference*difference));retained_fraction=float(np.sum(restricted*restricted)/difference2)
    details=dict(lineExtractionRelativeDifference=extraction,denseColoredScanRelativeDifference=scan,originalNSDirectionFDErrors=errors,
        sameLineFractionOfSmallDomainViscousDifferenceSquaredNorm=retained_fraction,
        droppedCrossLineDifferenceNorm=float(np.linalg.norm(difference*(~same))),completeReferenceApplications=160,actualUnitNS=5)
    assert s['report']['budgetRejected']
    if 'batchSizes' in s['report']:
        np.testing.assert_array_equal(s['report']['batchSizes'],[1,3,17])
        assert np.max(s['report']['batchAssemblyRelativeDifferences'])<2e-12
        details['batchSizes']=[1,3,17];details['maximumBatchAssemblyDifference']=float(np.max(s['report']['batchAssemblyRelativeDifferences']))
elif kind in ('integration','real'):
    file=VAL/case/('steady_newton_unit_evidence.mat' if kind=='integration' else 'wall_velocity_krylov_evidence.mat')
    s=read(file);r=s['result'];op=r['firstOperator'];correction=r['wlsNormalCorrection']
    assert r['options']['wlsNormalCorrection']['maximumAdditionalBytes']==correction['maximumAdditionalBytes']
    assert correction['additionalPersistentBlockUpperBoundBytes']<=correction['maximumAdditionalBytes']
    assert correction['stressHeatAndBothWorkTermsIncluded'] and not correction['crossOriginalFacetWlsBlocksIncluded'] and not correction['physicalEquationChanged']
    assert r['originalNSCombinedDirectionCheckPreserved'] and not r['originalSpatialEquationReplaced']
    assert r['coarseOperatorReuses']==1 and r['coarseOperatorBuilds']==r['coarseConstructionModelApplications']==0
    assert r['coarseImageBytes']<=r['options']['coarseMaximumImageBytes']
    assert r['options']['linearRelativeTolerance']==.1 and r['options']['steadyResidualTolerance']==2e-6
    assert r['options']['steadyFluxTolerance']==5e-10 and r['options']['spatialConservationTolerance']==1e-9
    if kind=='integration':
        assert s['report']['oldSingleGuessCompleteOutputReplayed'] and s['report']['additionalBlockBudgetRejected']
        assert r['acceptedNonlinearIterations']==1 and s['limited']['acceptedNonlinearIterations']==s['limited']['candidateTrials']==0
        assert r['residualEvaluations']==3 and s['limited']['residualEvaluations']==2
    else:
        source=read(Path(r['steadyInitialGuessReport']['sourceFile']),['result'])['result']
        np.testing.assert_array_equal(source['Q'],r['initialQ']);np.testing.assert_array_equal(source['rate'],r['initialRate'])
        np.testing.assert_array_equal(source['info']['wallPressure'],r['initialInfo']['wallPressure'])
        np.testing.assert_array_equal(s['gas']['temperatureRange'],[80,2000]);assert not r['steadyConverged']
    details=dict(acceptedIterations=int(r['acceptedNonlinearIterations']),actualNS=int(r['residualEvaluations']),
        actualDefects=[float(a['fullSystemRelativeDefect']) for a in rows(r['linearLog'])],
        sameLineOnly=True,additionalPersistentBlockUpperBoundBytes=int(correction['additionalPersistentBlockUpperBoundBytes']),
        originalGatesKept=True,coarseImageUnchanged=True)
else:
    s=read(VAL/case/'wls_normal_real_blocks.mat');new=read(Path(s['sourceFile']),['result'])['result'];old=read(Path(s['parentFile']),['result'])['result']
    np.testing.assert_array_equal(new['initialQ'],old['Q']);op=new['firstOperator'];graph=old['firstOperator'];nf=int(op['sourceFaces']);n=len(new['initialQ']);errors=[]
    for e in np.asarray(s['evidence'],dtype=object).ravel():
        # scipy的二维MATLAB cell在simplify_cells下仍可能保留mat_struct成员。
        if not isinstance(e,dict):e={key:getattr(e,key) for key in e._fieldnames}
        direction=np.zeros((n,5));cell,eq=int(e['cellID'])-1,int(e['equationID'])-1;direction[cell,eq]=1
        measured=action(graph,direction,False)-action(op,direction,True);expected=e['fullViscousAction']-e['approxViscousAction'];expected[np.arange(n)%nf!=cell%nf]=0
        errors.append(rel(measured,expected));assert errors[-1]<2e-9 and rel(e['actual'],measured)<2e-12
    assert len(errors)==10 and s['report']['newActualNSEvaluations']==0
    details=dict(directions=10,maximumCurrentReferenceDifference=max(errors),newActualNS=0)
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-WLS-normal-'+kind+'-matrix-reference-and-original-gates-audit',case=case,kind=kind,processID=os.getpid(),details=details,
    fullIndependentNSSolver=False,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(dict(status=out['status'],case=case,processID=os.getpid(),details=details)))
