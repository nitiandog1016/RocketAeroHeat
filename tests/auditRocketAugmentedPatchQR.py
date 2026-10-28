# -*- coding: utf-8 -*-
"""独立核对MATLAB的完整QR右映射、方向提升、固定目标及原始状态幅度。"""
from pathlib import Path
from datetime import datetime, timezone
import os, sys, json, gc
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT, sha, local
from screenRocketParentImageScale import setup
from screenRocketAugmentedPatchQR import metric_blocks, block_weight


def relative(a, b):
    return float(np.linalg.norm(a-b)/max(1., np.linalg.norm(b)))


def unit(folder, read):
    s = read(folder/'augmented_qr_unit_evidence.mat', ['cases'])
    checks = []
    for c in s['cases']:
        n = int(c['packed']['cellCount']);N = n*5;b = c['source'];penalty = c['factor']['penaltyInOriginalResidualUnits']
        A = np.row_stack([b['base']['C'], penalty*b['W']])
        rhs = np.r_[c['rhs'].ravel(order='F'), np.zeros(N)]
        dense = np.linalg.lstsq(A, rhs, rcond=None)[0]
        errs = [relative(c[k].ravel(order='F'), dense) for k in ('value', 'warm')]
        assert max(errs) < 1e-8
        assert relative(c['P']@c['coordinates'].ravel(order='F'), c['value'].ravel(order='F')) < 1e-12
        gram_max = 0.
        buckets = c['factor']['buckets'];buckets = [buckets] if isinstance(buckets, dict) else np.atleast_1d(buckets)
        for k, bucket in enumerate(buckets, 1):
            count = np.size(bucket['patchIDs'])
            if not count:continue
            members = np.asarray(bucket['members']).reshape(k,count,order='F').T.astype(np.int64)-1
            maps = np.asarray(bucket['map']).reshape(5*k,5*k,count,order='F').transpose(2,0,1)
            for group, mapping in zip(members,maps):
                ids = (group[:,None]+n*np.arange(5)[None,:]).ravel()
                image = A[:,ids]@mapping
                error = float(np.linalg.norm(image.T@image-np.eye(5*k),ord=np.inf))
                assert error < 1e-8;gram_max=max(gram_max,error)
        checks.append(dict(zeroAndWarmDenseDifferences=errs,maximumFullAugmentedLocalGramDefect=gram_max))
    return checks


def real(folder, read):
    s=read(folder/'augmented_qr_direction.mat', None)
    pfile=Path(s['preparationFile']);p=json.loads(pfile.read_text(encoding='utf-8-sig'))
    assert sha(pfile)==s['preparationSHA256']
    for a in p['artifacts'].values():assert sha(a['file'])==a['sha256']
    assert sha(p['sourceBinding']['sourceFile'])==s['sourceBinding']['sourceFileSHA256']
    packed=read(p['artifacts']['packed']['file'], ['packed'])['packed']
    factor=read(folder/'augmented_qr_factor.mat',['factor'])['factor']
    buckets=[factor['buckets']] if isinstance(factor['buckets'],dict) else list(np.atleast_1d(factor['buckets']))
    for b in buckets:b['inverse']=b['map']
    factor['buckets']=buckets
    A,plan,order=setup(packed,factor);del packed,factor,buckets;gc.collect()
    value=s['value'][order].ravel();lift=local(plan,s['coordinates'][order].ravel())
    lift_error=relative(lift,value);assert lift_error<1e-12
    image=A@value;rhs=s['rhs'][order].ravel();relative_c=float(np.linalg.norm(rhs-image)/np.linalg.norm(rhs))
    np.testing.assert_allclose(relative_c,s['detail']['relativeCResidual'],rtol=1e-11)
    c=read(p['artifacts']['cache']['file'],['Q','gas','scale']);q=c['Q'];gas=c['gas'];rho=q[:,0];u=q[:,1:4]/rho[:,None]
    cv=gas['R']/(gas['gamma']-1);T=(q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv)
    dq=s['value']*c['scale'];dr=dq[:,0]
    du=(dq[:,1:4]-u*dr[:,None])/rho[:,None]
    dt=(dq[:,4]-np.sum(u*dq[:,1:4],axis=1)+(.5*np.sum(u*u,axis=1)-cv*T)*dr)/(rho*cv)
    raw=np.column_stack([dr/rho,du/np.sqrt(gas['gamma']*gas['R']*T)[:,None],dt/T])
    raw_error=relative(raw,s['raw']);assert raw_error<1e-12
    W=metric_blocks(s['metric'],order);weighted=block_weight(W,value)
    np.testing.assert_allclose(weighted.reshape(-1,5), (raw*s['metric']['peakRowWeights'])[order], rtol=1e-10,atol=1e-10)
    penalty=s['detail']['penaltyInOriginalResidualUnits']
    np.testing.assert_allclose(penalty,.1*np.linalg.norm(rhs)/p['metricNormalization'],rtol=1e-12)
    objective=relative_c**2+(penalty*np.linalg.norm(weighted)/np.linalg.norm(rhs))**2
    np.testing.assert_allclose(objective,s['detail']['normalizedAugmentedObjective'],rtol=1e-11)
    # 核对每个局部映射对完整原C像和W像的Gram；不能只抽查对角块。
    counts=np.bincount(A.indices,minlength=len(order));bycol=np.argsort(A.indices,kind='stable')
    ptr=np.r_[0,np.cumsum(counts)];owners=np.repeat(np.arange(len(order)),np.diff(A.indptr))
    maximum=0.;used=0;groups=0
    for members,maps in plan:
        for group,P in zip(members,maps):
            k=len(group);ids=np.concatenate([bycol[ptr[j]:ptr[j+1]] for j in group])
            cols=np.concatenate([np.full(counts[j],a) for a,j in enumerate(group)])
            rows,rr=np.unique(owners[ids],return_inverse=True)
            panel=np.zeros((len(rows),5,k,5));panel[rr,:,cols,:]=A.data[ids]
            bottom=np.zeros((k,5,k,5))
            for j,cell in enumerate(group):bottom[j,:,j,:]=penalty*W[cell]
            aug=np.row_stack([panel.reshape(-1,5*k),bottom.reshape(5*k,5*k)])
            im=aug@P;error=float(np.linalg.norm(im.T@im-np.eye(5*k),ord=np.inf))
            assert error<1e-8;maximum=max(maximum,error);used+=len(ids);groups+=1
    assert used==len(A.indices)==s['construction']['originalBlocksUsedExactlyOnce']
    model=read(folder/'augmented_qr_model.mat',['value','fullModel','report'])
    np.testing.assert_array_equal(model['value'],s['value'])
    full=float(np.linalg.norm(s['rhs']+model['fullModel'])/np.linalg.norm(rhs))
    np.testing.assert_allclose(full,model['report']['completeSelectedModelRelativeDefect'],rtol=1e-11)
    return [dict(allOriginalCVs=len(order),allOriginalBlocks=used,allLocalMaps=groups,
                 maximumFullAugmentedLocalGramDefect=maximum,coordinateLiftRelativeDifference=lift_error,
                 rawStateMetricRelativeDifference=raw_error,originalCRelativeResidual=relative_c,
                 normalizedAugmentedObjective=float(objective),completeSelectedModelRelativeDefect=full,
                 originalCTargetMet=bool(s['detail']['internalCTargetMet']),
                 allNativeBidiagonalizationVectorsIndependentlyVerified=False)]


def main(mode,folder):
    assert mode in ('unit','real');folder=Path(folder).resolve();bindings={}
    def read(path,variables):
        path=Path(path);bindings[str(path)]=sha(path)
        return loadmat(path,simplify_cells=True,variable_names=variables)
    checks=unit(folder,read) if mode=='unit' else real(folder,read)
    assert all(sha(path)==h for path,h in bindings.items())
    out=dict(status='passed-independent-complete-augmented-QR-map-objective-and-state-metric',mode=mode,checks=checks,
             processID=os.getpid(),recordedUTC=datetime.now(timezone.utc).isoformat(),inputSHA256=bindings,
             scriptSHA256=sha(__file__),physicalSolution=False,newProductionSolvedAreaM2=0)
    target=folder/'independent_qr_audit.json';assert not target.exists()
    target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))


if __name__=='__main__':main(*sys.argv[1:3])
