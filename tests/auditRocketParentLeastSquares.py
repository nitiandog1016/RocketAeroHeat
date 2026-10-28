# -*- coding: utf-8 -*-
"""独立检查完整父C全部小域系数/伴随及真实全部像列尺度、坐标提升和最终像。"""
from pathlib import Path
import os,json,sys,gc
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,VAL,sha,local
from screenRocketParentImageScale import setup

mode,name=sys.argv[1:3];assert mode in ('unit','real') and name.replace('_','').isalnum()
folder=VAL/name;target=VAL/(name+'_independent_parent_least_squares.json');assert not target.exists();bindings={}


def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)


def rel(a,b): return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(b)))


if mode=='unit':
    s=read(folder/'parent_least_squares_unit_evidence.mat');results=[]
    for case in s['cases']:
        A,plan,order=setup(case['packed'],case['factor']);n=len(order);N=5*n;eye=np.eye(N)
        Malg=np.column_stack([local(plan,eye[:,j]) for j in range(N)])
        ids=(order[:,None]+np.arange(5)[None,:]*n).ravel();C=np.zeros((N,N));M=C.copy()
        C[np.ix_(ids,ids)]=A.toarray();M[np.ix_(ids,ids)]=Malg
        errors=dict(C=rel(case['C'],C),CT=rel(case['CT'],C.T),M=rel(case['M'],M),MT=rel(case['MT'],M.T))
        assert max(errors.values())<1e-12
        expected_alg=np.linalg.norm(A.toarray()@Malg,axis=0).reshape(n,5)
        expected=np.empty_like(expected_alg);expected[order]=expected_alg
        norm_error=float(np.max(abs(case['columnNorm']-expected)/expected));assert norm_error<1e-12
        assert np.max(abs(case['columnScale']*case['columnNorm']-1))<1e-12
        residual=float(np.linalg.norm(case['rhs'].ravel(order='F')-C@case['value'].ravel(order='F'))/np.linalg.norm(case['rhs']))
        direct=np.linalg.solve(C,case['rhs'].ravel(order='F'))
        assert residual<1e-10 and rel(case['value'].ravel(order='F'),direct)<1e-8
        correction=local(plan,(case['coordinates'][order]*case['columnScale'][order]).ravel()).reshape(n,5)
        expected_value=case['initial'].copy();expected_value[order]+=correction
        assert rel(case['value'],expected_value)<1e-12
        assert case['detail']['adjointIdentityRelativeDifference']<1e-12
        assert case['earlyDetail']['nativeIterations']==case['zeroDetail']['nativeIterations']==0
        results.append(dict(allDenseMatrixRelativeDifferences=errors,allColumnNormRelativeDifference=norm_error,
            completeParentResidual=residual,coordinateLiftDifference=rel(case['value'],expected_value),
            nativeIterations=int(case['detail']['nativeIterations']),referenceBackslashRoundoffDifference=float(case['referenceRoundoffDifference'])))
    out=dict(status='passed-independent-two-orderings-all-parent-matrices-adjoints-image-columns-and-dense-LSQR',processID=os.getpid(),cases=results)
else:
    s=read(folder/'parent_least_squares_real_evidence.mat')
    summary=json.loads((folder/'parent_least_squares_summary.json').read_text(encoding='utf8'))
    assert not summary['physicalSolution'] and summary['newActualNSEvaluations']==summary['newFlowAcceptances']==summary['newProductionSolvedAreaM2']==0
    for key,digest in s['sourceIdentity'].items(): assert sha(ROOT/'internalFlow/cutCell'/(key+'.m'))==digest
    assert sha(s['sourceBinding']['sourceFile'])==s['sourceBinding']['sourceFileSHA256']
    assert sha(s['factorFile'])=='09b04be0c77f76f8c04bd717745a0dd2fbcef7cdb6809d53bc403a07ed6349e0'
    assert sha(s['packedFile'])=='beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'
    packed=read(s['packedFile'],['packed'])['packed'];factor=read(s['factorFile'],['factor'])['factor']
    A,plan,order=setup(packed,factor);n=len(order);del factor;gc.collect()
    data=read(s['scaleFile'],['columnScale','columnNorm','construction'])
    reference=read(VAL/'parentCimagescalescreen02/parent_image_scale_evidence.mat',['columnNormAlgebraic','originalOrder','packedEvidenceSHA256','coupledFactorSHA256'])
    np.testing.assert_array_equal(order+1,reference['originalOrder'])
    norms=np.empty_like(data['columnNorm']);norms[order]=reference['columnNormAlgebraic']
    norm_error=float(np.max(abs(norms-data['columnNorm'])/norms));assert norm_error<1e-12
    assert np.max(abs(data['columnNorm']*data['columnScale']-1))<1e-12
    base=read(s['baseFile'],['rhs','value']);np.testing.assert_array_equal(s['rhs'],base['rhs']);np.testing.assert_array_equal(s['initial'],base['value'])
    correction=local(plan,(s['coordinates'][order]*data['columnScale'][order]).ravel()).reshape(n,5)
    expected=s['initial'].copy();expected[order]+=correction;lift_error=rel(s['value'],expected);assert lift_error<1e-8
    image_alg=(A@s['value'][order].ravel()).reshape(n,5);image=np.empty_like(image_alg);image[order]=image_alg
    image_error=rel(image,s['image']);assert image_error<1e-12
    residual=float(np.linalg.norm(s['rhs']-image)/np.linalg.norm(s['rhs']))
    assert abs(residual-s['detail']['relativeCResidual'])<1e-12
    assert bool(s['detail']['internalCTargetMet'])==(residual<=1e-3)
    assert not s['detail']['oldSchurUsed'] and not s['detail']['fixedCRepresentativesUsed']
    out=dict(status='passed-independent-all-real-parent-image-columns-coordinate-lift-and-complete-C-final-image',processID=os.getpid(),
        allOriginalCVs=n,allPreconditionedColumns=5*n,allColumnNormMaximumRelativeDifference=norm_error,
        completeCoordinateLiftDifference=lift_error,completeParentImageDifference=image_error,finalRelativeOriginalCResidual=residual,
        originalCTargetMet=residual<=1e-3,actualNativeIterations=int(s['detail']['nativeIterations']),
        originalParentNormAndAllFiveFieldsPreserved=True,oldSchurUsed=False,fixedCRepresentativesUsed=False,
        independentAllRealBidiagonalizationColumnsClaimed=False)
out.update(newActualNSEvaluations=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0,physicalSolution=False,inputSHA256=bindings,scriptSHA256=sha(__file__),
    helperSourceSHA256={name:sha(ROOT/'tests'/name) for name in ('screenRocketFullFLeastSquares.py','screenRocketParentImageScale.py')})
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
