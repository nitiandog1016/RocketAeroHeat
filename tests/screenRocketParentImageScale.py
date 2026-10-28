# -*- coding: utf-8 -*-
"""直接完整父C的五量局部映射/像列归一对照，不使用旧Schur或固定C代表值。"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,gc
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from scipy.io import loadmat,savemat
from scipy.sparse import bsr_matrix
from scipy.sparse.linalg import LinearOperator,lsqr
from screenRocketFullFLeastSquares import ROOT,VAL,sha,local
from screenRocketImageScaledFullF import column_norms


def setup(packed,factor):
    n=int(packed['cellCount']);order=np.atleast_1d(packed['order']).astype(np.int64)-1
    inverse_order=np.empty(n,np.int64);inverse_order[order]=np.arange(n)
    ptr=np.atleast_1d(packed['rowPointer']).astype(np.int64)-1;col=np.atleast_1d(packed['column']).astype(np.int64)-1
    blocks=np.asarray(packed['values']).reshape(5,5,len(col),order='F').transpose(2,0,1)
    A=bsr_matrix((blocks,col,ptr),shape=(5*n,5*n));plan=[]
    buckets=[factor['buckets']] if isinstance(factor['buckets'],dict) else list(np.atleast_1d(factor['buckets']).ravel())
    for k,b in enumerate(buckets,1):
        count=np.size(b['patchIDs'])
        if not count: continue
        physical=np.asarray(b['members']).reshape(k,count,order='F').T.astype(np.int64)-1
        members=inverse_order[physical]
        inv=np.asarray(b['inverse']).reshape(5*k,5*k,count,order='F').transpose(2,0,1)
        plan.append((members,inv))
    return A,plan,order


def make_operator(A,plan,scale):
    d=np.ones(A.shape[1]) if scale is None else scale.ravel()
    count=dict(C=0,CT=0,M=0,MT=0)
    def forward(v):
        count['C']+=1;count['M']+=1
        return A@local(plan,d*v)
    def transpose(v):
        count['CT']+=1;count['MT']+=1
        return d*local(plan,A.T@v,transpose=True)
    return LinearOperator(A.shape,matvec=forward,rmatvec=transpose,dtype=np.float64),count


def adjoint(B):
    x=np.sin(np.arange(B.shape[0])*.17);y=np.cos(np.arange(B.shape[0])*.23)
    bx=B@x;bty=B.rmatvec(y)
    return abs(float(np.dot(bx,y)-np.dot(x,bty)))/max(1.,np.linalg.norm(bx)*np.linalg.norm(y),np.linalg.norm(bty)*np.linalg.norm(x))


def main():
    folder=VAL/'parentCimagescalescreen02';assert not folder.exists();folder.mkdir();bindings={}
    def read(path,variables=None):
        path=Path(path);path=path if path.is_absolute() else ROOT/path
        bindings[path.relative_to(ROOT).as_posix()]=sha(path)
        return loadmat(path,simplify_cells=True,variable_names=variables)
    s=read(VAL/'coupledpatchunit03/coupled_patches_unit_evidence.mat',['cases']);small=[]
    for case in s['cases']:
        A,plan,order=setup(case['packed'],case['factor']);norms,construction=column_norms(A,plan);scale=1/norms
        B,counts=make_operator(A,plan,scale)
        eye=np.eye(A.shape[0]);dense=np.column_stack([A@local(plan,eye[:,j]) for j in range(A.shape[0])])
        expected=np.linalg.norm(dense,axis=0).reshape(norms.shape)
        error=float(np.max(abs(expected-norms)/expected));assert error<1e-12 and adjoint(B)<1e-12
        b=case['rhs'][order].ravel();result=lsqr(B,b,atol=0.,btol=1e-11,conlim=1e12,iter_lim=600)
        value=local(plan,scale.ravel()*result[0]);residual=float(np.linalg.norm(b-A@value)/np.linalg.norm(b));assert residual<1e-10
        physical=np.empty_like(case['rhs']);physical[order]=value.reshape(len(order),5)
        direct=np.asarray(case['direct']).reshape(physical.shape,order='F')
        assert np.linalg.norm(physical-direct)/max(1.,np.linalg.norm(direct))<1e-8
        small.append(dict(allColumnNormMaximumRelativeDifference=error,completeParentCResidual=residual,adjointDifference=adjoint(B)))
    del s,case,A,plan,B,eye,dense,result
    gc.collect();print('SMALL_ORIGINAL_C_ORDERING_IMAGE_NORMS_AND_DENSE_SOLUTION_PASSED',flush=True)
    factor_file=VAL/'coupledpatchreal01/coupled_patch_factor_evidence.mat';assert sha(factor_file)=='09b04be0c77f76f8c04bd717745a0dd2fbcef7cdb6809d53bc403a07ed6349e0'
    raw=read(factor_file,['factor','packedFile','sourceBinding']);factor=raw['factor'];packed_file=Path(raw['packedFile']);binding=raw['sourceBinding'];del raw
    assert sha(packed_file)=='beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'
    packed=read(packed_file,['packed'])['packed'];A,plan,order=setup(packed,factor);del factor;gc.collect()
    norms,construction=column_norms(A,plan);scale=1/norms
    assert construction['groups']==29523 and construction['originalBlocksUsedExactlyOnce']==1619664
    baseline=read(VAL/'spatialreductionreal01/spatial_reduction_real_evidence.mat',['rhs','value','detail'])
    b=baseline['rhs'][order].ravel();initial=baseline['value'][order].ravel();image0=A@initial;r0=b-image0
    initial_residual=float(np.linalg.norm(r0)/np.linalg.norm(b));assert abs(initial_residual-baseline['detail']['relativeCResidual'])<1e-12
    comparisons=[]
    for label,d in [('unscaled-original-all-five-field-local-map',None),('complete-parent-image-column-normalization',scale)]:
        B,counts=make_operator(A,plan,d);adj=adjoint(B);assert adj<1e-12;before=counts.copy()
        result=lsqr(B,r0,atol=0.,btol=1e-3*np.linalg.norm(b)/np.linalg.norm(r0),conlim=1e12,iter_lim=32)
        actual_counts={key:counts[key]-before[key] for key in counts};coordinates=result[0] if d is None else d.ravel()*result[0]
        value=initial+local(plan,coordinates);image=A@value;residual=float(np.linalg.norm(b-image)/np.linalg.norm(b))
        physical=np.empty_like(baseline['rhs']);physical[order]=value.reshape(len(order),5)
        physical_image=np.empty_like(physical);physical_image[order]=image.reshape(len(order),5)
        out=dict(coordinates=label,initialOriginalCResidual=initial_residual,finalOriginalCResidual=residual,actualIterations=int(result[2]),nativeStopCode=int(result[1]),
            nativeEstimatedCondition=float(result[6]),actualSolverOperatorApplications=actual_counts,adjointIdentityDifference=adj,originalCTargetMet=residual<=1e-3)
        comparisons.append(out);print(json.dumps(out),flush=True)
        savemat(folder/('unscaled_evidence.mat' if d is None else 'image_scaled_evidence.mat'),dict(value=physical,image=physical_image,rhs=baseline['rhs'],
            sourceBinding=binding,originalOrder=order+1,packedEvidenceSHA256=sha(packed_file),coupledFactorSHA256=sha(factor_file)),do_compression=True)
    savemat(folder/'parent_image_scale_evidence.mat',dict(columnNormAlgebraic=norms,columnScaleAlgebraic=scale,originalOrder=order+1,
        sourceBinding=binding,packedEvidenceSHA256=sha(packed_file),coupledFactorSHA256=sha(factor_file)),do_compression=True)
    assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
    out=dict(status='tested-direct-complete-parent-C-forward-adjoint-image-column-normalization-and-same-budget-LSQR-screen',
        processID=os.getpid(),recordedUTC=datetime.now(timezone.utc).isoformat(),smallReferences=small,construction=construction,
        allOriginalPhysicalCVs=len(order),allFiveFieldsAndOriginalNormPreserved=True,columnNormRange=[float(norms.min()),float(norms.max())],
        columnNormSpan=float(norms.max()/norms.min()),comparisons=comparisons,oldSchurUsed=False,fixedCRepresentativesUsed=False,
        matlabProductionSolverModified=False,physicalStateChanged=False,newActualNSEvaluations=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0,
        processPeakMemoryClaimed=False,inputSHA256=bindings,scriptSHA256=sha(__file__),
        helperSourceSHA256={name:sha(ROOT/'tests'/name) for name in ('screenRocketFullFLeastSquares.py','screenRocketImageScaledFullF.py')})
    (folder/'parent_image_scale_screen_summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256','smallReferences','comparisons')}))


if __name__=='__main__': main()
