# -*- coding: utf-8 -*-
"""独立核对完整像列尺度、小域缩放解和真实全部尺度/完整F及父C最终像。"""
from pathlib import Path
import json
import os
import sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import bsr_matrix
from screenRocketFullFLeastSquares import ROOT, VAL, sha, arrays

mode, name = sys.argv[1:3]
assert mode in ('unit','real') and name.replace('_','').isalnum()
folder=VAL/name
target=VAL/(name+'_independent_full_F_image_scale.json')
assert not target.exists()
bindings={}


def read(path, variables=None):
    path=Path(path)
    if not path.is_absolute(): path=ROOT/path
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)


def relative(a,b):
    return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(b)))


def action(g,value):
    n=int(g['cellCount'])
    ptr=np.atleast_1d(g['rowPointer']).astype(np.int64)-1
    col=np.atleast_1d(g['column']).astype(np.int64)-1
    order=np.atleast_1d(g['order']).astype(np.int64)-1
    data=np.asarray(g['values']).reshape(5,5,len(col),order='F').transpose(2,0,1)
    A=bsr_matrix((data,col,ptr),shape=(5*n,5*n))
    result=(A@value[order].ravel()).reshape(n,5)
    output=np.empty_like(result)
    output[order]=result
    return output


def dense(g):
    n=int(g['cellCount'])
    ptr=np.atleast_1d(g['rowPointer']).astype(np.int64)-1
    col=np.atleast_1d(g['column']).astype(np.int64)-1
    order=np.atleast_1d(g['order']).astype(np.int64)-1
    blocks=np.asarray(g['values']).reshape(5,5,len(col),order='F')
    A=np.zeros((5*n,5*n))
    for row in range(n):
        for edge in range(ptr[row],ptr[row+1]):
            rr=5*order[row]+np.arange(5);cc=5*order[col[edge]]+np.arange(5)
            A[np.ix_(rr,cc)]=blocks[:,:,edge]
    return A


if mode=='unit':
    s=read(folder/'full_F_image_scale_unit_evidence.mat')
    old=read(s['sourceFile'],['cases'])
    results=[]
    for case,prior in zip(s['cases'],old['cases']):
        n=int(case['graph']['cellCount']);A=dense(case['graph']);M=np.zeros_like(A)
        for members,inverses in arrays(case['factor']):
            for ids,inverse in zip(members,inverses):
                dofs=(5*ids[:,None]+np.arange(5)[None,:]).ravel()
                M[np.ix_(dofs,dofs)]=inverse
        norms=np.linalg.norm(A@M,axis=0).reshape(n,5)
        permuted=np.linalg.norm(dense(case['permutedGraph'])@M,axis=0).reshape(n,5)
        error=float(np.max(abs(norms-case['columnNorm'])/norms))
        permuted_error=float(np.max(abs(permuted-case['permutedNorms'])/permuted))
        assert error<1e-12 and permuted_error<1e-12
        assert np.max(abs(case['columnScale']*norms-1))<1e-12
        residual=float(np.linalg.norm(case['rhs'].ravel()-A@case['value'].ravel())/np.linalg.norm(case['rhs']))
        direct=np.linalg.solve(A,case['rhs'].ravel())
        assert residual<1e-10 and relative(case['value'].ravel(),direct)<1e-8
        np.testing.assert_array_equal(case['plainValue'],prior['value'])
        assert case['detail']['adjointIdentityRelativeDifference']<1e-12
        results.append(dict(allColumnNormRelativeDifference=error,permutedColumnNormRelativeDifference=permuted_error,
            fullFRelativeResidual=residual,oldDefaultValueExactlyReproduced=True))
    out=dict(status='passed-independent-all-small-complete-image-columns-permutations-scaled-LSQR-and-old-default',
        processID=os.getpid(),cases=results,newActualNSEvaluations=0,newProductionSolvedAreaM2=0)
else:
    s=read(folder/'full_F_image_scaled_real_evidence.mat')
    report=json.loads((folder/'full_F_image_scaled_summary.json').read_text(encoding='utf8'))
    for source,digest in s['sourceIdentity'].items():
        assert sha(ROOT/'internalFlow/cutCell'/(source+'.m'))==digest
    assert sha(s['sourceBinding']['sourceFile'])==s['sourceBinding']['sourceFileSHA256']
    scale=read(s['scaleFile'],['columnNorm','columnScale','construction'])
    expected=read(VAL/'fullFimagescalescreen01/full_F_image_scale_evidence.mat',['columnNorm','columnScale','graphEvidenceSHA256','originalFactorAssemblySHA256'])
    graph_path=Path(s['graphFile']);graph_path=graph_path if graph_path.is_absolute() else ROOT/graph_path
    assembly_path=Path(s['assemblyFile']);assembly_path=assembly_path if assembly_path.is_absolute() else ROOT/assembly_path
    assert sha(graph_path)==expected['graphEvidenceSHA256'] and sha(assembly_path)==expected['originalFactorAssemblySHA256']
    scale_error=float(np.max(abs(scale['columnNorm']-expected['columnNorm'])/expected['columnNorm']))
    assert scale_error<1e-12 and np.max(abs(scale['columnScale']*scale['columnNorm']-1))<1e-12
    graph=read(graph_path,['graph'])['graph']
    base=read(s['baseFile'],['rhsF','valueF','value'])
    np.testing.assert_array_equal(base['rhsF'],s['rhsF']);np.testing.assert_array_equal(base['valueF'],s['initial'])
    image=action(graph,s['valueF']);image_error=relative(image,s['Fimage'])
    residual=float(np.linalg.norm(s['rhsF']-image)/np.linalg.norm(s['rhsF']))
    assert image_error<1e-12 and abs(residual-s['detail']['relativeFullFResidual'])<1e-12
    f=np.atleast_1d(graph['originalPhysicalCVs']).astype(np.int64)-1
    value=base['value'].copy();value[f]=s['valueF'];np.testing.assert_array_equal(value,s['value'])
    packed=read(s['packedFile'],['packed'])['packed'];parent_image=action(packed,s['value'])
    parent_error=relative(parent_image,s['image']);remaining=s['originalRHS']-parent_image
    parent_residual=float(np.linalg.norm(remaining)/np.linalg.norm(s['originalRHS']))
    consistency=float(np.linalg.norm(remaining[f]-(s['rhsF']-image))/max(1.,np.linalg.norm(s['originalRHS'])))
    assert parent_error<1e-12 and consistency<1e-8
    assert abs(parent_residual-s['diagnostic']['relativeParentOriginalCResidual'])<1e-10
    assert bool(s['detail']['fullFInternalTargetMet'])==(residual<=1e-3)
    assert not report['physicalSolution'] and report['newActualNSEvaluations']==report['newProductionSolvedAreaM2']==0
    out=dict(status='passed-independent-all-real-complete-image-column-norms-scaled-F-and-fixed-C-parent-diagnostic',
        processID=os.getpid(),allPreconditionedColumns=int(scale['columnNorm'].size),maximumColumnNormRelativeDifference=scale_error,
        scaledFImageDifference=image_error,parentOriginalCImageDifference=parent_error,finalFullFResidual=residual,
        fullFInternalTargetMet=residual<=1e-3,fixedCParentOriginalCResidual=parent_residual,FRowsSubmatrixConsistency=consistency,
        oldApproximateSchurIdentityClaimed=False,independentAllRealBidiagonalizationColumnsClaimed=False,
        originalFResidualNormUnchanged=True,newActualNSEvaluations=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0,physicalSolution=False)
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out.update(inputSHA256=bindings,scriptSHA256=sha(__file__),helperSourceSHA256=sha(ROOT/'tests/screenRocketFullFLeastSquares.py'))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
