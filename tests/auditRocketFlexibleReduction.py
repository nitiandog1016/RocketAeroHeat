# -*- coding: utf-8 -*-
"""独立核对全部可变块回代中间量和真实存储Z的完整父C Arnoldi关系。"""
from pathlib import Path
import os,json,sys,gc
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,VAL,sha

mode,name=sys.argv[1:3]
assert mode in ('unit','real') and name.replace('_','').isalnum()
folder=VAL/name;target=VAL/(name+'_independent_flexible_reduction.json');assert not target.exists()
bindings={}


def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)


def rows(value): return [value] if isinstance(value,dict) else list(np.atleast_1d(value).ravel())
def rel(a,b): return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(b)))


class Packed:
    def __init__(self,g):
        self.n=int(g['cellCount']);self.order=np.atleast_1d(g['order']).astype(np.int64)-1
        ptr=np.atleast_1d(g['rowPointer']).astype(np.int64)-1
        self.col=np.atleast_1d(g['column']).astype(np.int64)-1
        self.owner=np.repeat(np.arange(self.n),np.diff(ptr))
        self.blocks=np.asarray(g['values']).reshape(5,5,len(self.col),order='F').transpose(2,0,1)

    def __call__(self,x):
        result=np.zeros((self.n,5));ordered=x[self.order]
        for first in range(0,len(self.col),4096):
            ii=slice(first,first+4096)
            np.add.at(result,self.owner[ii],np.einsum('kab,kb->ka',self.blocks[ii],ordered[self.col[ii]]))
        value=np.empty_like(result);value[self.order]=result;return value


def rect(g,x):
    owner=np.atleast_1d(g['owner']).astype(np.int64)-1;col=np.atleast_1d(g['column']).astype(np.int64)-1
    blocks=np.asarray(g['values']).reshape(5,5,len(col),order='F').transpose(2,0,1)
    result=np.zeros((int(g['rowCount']),5))
    for first in range(0,len(col),4096):
        ii=slice(first,first+4096)
        np.add.at(result,owner[ii],np.einsum('kab,kb->ka',blocks[ii],x[col[ii]]))
    return result


def coarse(factor,rhs):
    result=np.zeros_like(rhs)
    for k,b in enumerate(rows(factor['buckets']),1):
        count=np.size(b['patchIDs'])
        if not count: continue
        members=np.asarray(b['members']).reshape(k,count,order='F').T.astype(np.int64)-1
        inverse=np.asarray(b['inverse']).reshape(5*k,5*k,count,order='F').transpose(2,0,1)
        result[members]=np.einsum('pij,pj->pi',inverse,rhs[members].reshape(count,5*k)).reshape(count,k,5)
    return result


def check(case):
    C=Packed(case['packed']);F=Packed(case['Fgraph']);r=case['reduction'];split=r['split']
    f=np.atleast_1d(split['fineOriginalCVs']).astype(np.int64)-1;c=np.atleast_1d(split['coarseOriginalCVs']).astype(np.int64)-1
    detail=case['detail'];k=detail['krylov'];cycle=rows(k['cycles'])[0]
    V=np.asarray(cycle['V']);Z=np.asarray(cycle['Z']);H=np.asarray(cycle['smallH'])
    used=int(cycle['arnoldiSteps']);assert Z.shape==(5*C.n,used)
    reports=rows(detail['preconditionerReports']);assert len(reports)==used
    scales=np.atleast_1d(cycle['imageNormalizationScales']);sequence=np.atleast_1d(detail['options']['innerIterationSequence']).astype(int)
    total=dict(F=0,FT=0,M=0,MT=0);maximum_block_error=maximum_inner_residual_error=maximum_z_error=0.;unmet=0
    for j,p in enumerate(reports):
        d=read(p['options']['evidenceFile']);b=np.asarray(d['rhs'])
        np.testing.assert_array_equal(b.ravel(order='F'),V[:,j])
        expected_coarse=b[c]-rect(r['E'],d['firstF']);xc=coarse(case['coarseFactor'],expected_coarse)
        expected_second=b[f]-rect(r['B'],xc)
        block_error=max(rel(d['coarseRHS'],expected_coarse),rel(d['xC'],xc),rel(d['secondRHS'],expected_second))
        assert block_error<1e-8
        raw=np.zeros_like(b);raw[f]=d['secondF'];raw[c]=d['xC']
        z_error=rel(Z[:,j],raw.ravel(order='F')/scales[j]);assert z_error<1e-8
        for label,rhs,x in [('firstReport',b[f],d['firstF']),('secondReport',d['secondRHS'],d['secondF'])]:
            q=d[label];actual=float(np.linalg.norm(rhs-F(x))/max(np.finfo(float).tiny,np.linalg.norm(rhs)))
            discrepancy=abs(actual-q['relativeFullFResidual']);assert discrepancy<1e-8
            assert bool(q['fullFInternalTargetMet'])==(actual<=q['options']['relativeTolerance'])
            assert q['nativeIterations']<=sequence[j%len(sequence)]
            maximum_inner_residual_error=max(maximum_inner_residual_error,discrepancy)
            unmet+=int(not q['fullFInternalTargetMet'])
            for key in total: total[key]+=int(q['operatorApplications'][key])
        assert not p['fixedLinearMapClaimed'] and not p['oldSchurIdentityClaimedForNewMap']
        maximum_block_error=max(maximum_block_error,block_error);maximum_z_error=max(maximum_z_error,z_error)
        print('INDEPENDENT_FLEXIBLE_APPLICATION='+str(j+1),flush=True)
    for key in total: assert total[key]==detail['innerOperatorApplications'][key]
    assert unmet==detail['FInnerTargetsNotMet']
    assert total['F']+total['FT']==detail['actualInnerFAndFTApplications']<=detail['options']['maximumInnerOperatorApplications']
    images=np.column_stack([C(Z[:,j].reshape(C.n,5,order='F')).ravel(order='F') for j in range(used)])
    relation=rel(images,V@H);assert relation<1e-8
    active=used+1
    if np.all(H[-1]==0): active=used
    orth=float(np.linalg.norm(V[:,:active].T@V[:,:active]-np.eye(active)));assert orth<1e-8
    candidate=np.asarray(cycle['xStart']).ravel()+Z@np.atleast_1d(cycle['weights'])
    assert rel(candidate,np.asarray(cycle['candidate']).ravel())<1e-8
    image=C(case['value']);residual=float(np.linalg.norm(case['rhs']-image)/np.linalg.norm(case['rhs']))
    assert abs(residual-detail['relativeCResidual'])<1e-8
    return dict(allActualPreconditionerApplications=used,maximumBlockBackSubstitutionDifference=maximum_block_error,
        maximumStoredZFromActualVariableMapDifference=maximum_z_error,maximumInnerFResidualDifference=maximum_inner_residual_error,
        allActualOuterArnoldiColumnsRelativeDifference=relation,outerOrthogonalityError=orth,
        finalOriginalCRelativeResidual=residual,innerOperatorApplications=total,FInnerTargetsNotMet=unmet,
        approximateOldSNotClaimedExact=True,variableInnerMapNotClaimedLinear=True,
        independentAllRealInnerLSQRColumnsClaimed=False)


if mode=='unit':
    s=read(folder/'flexible_reduction_unit_evidence.mat');results=[]
    for case in rows(s['cases']):
        result=check(case)
        expected=np.linalg.solve(case['C'],case['rhs'].ravel(order='F'))
        assert rel(case['value'].ravel(order='F'),expected)<1e-8 and result['finalOriginalCRelativeResidual']<1e-10
        assert case['nonlinearityDifference']>1e-8
        result['measuredInnerMapNonlinearity']=float(case['nonlinearityDifference']);results.append(result)
    out=dict(status='passed-independent-two-small-cases-variable-block-maps-all-flexible-Arnoldi-columns-and-dense-C-solution',processID=os.getpid(),cases=results)
else:
    s=read(folder/'flexible_reduction_real_evidence.mat');summary=json.loads((folder/'flexible_reduction_summary.json').read_text())
    assert summary['innerRunComplete'] and not summary['candidateRejected']
    for name,digest in s['sourceIdentity'].items(): assert sha(ROOT/'internalFlow/cutCell'/(name+'.m'))==digest
    assert sha(s['sourceBinding']['sourceFile'])==s['sourceBinding']['sourceFileSHA256']
    file=Path(s['assemblyFile']);assert sha(file)=='b04c5c8d604c52c9289ce5d08127270a41c31322eef5ac44cd0fda4ebc0c76ca'
    r=read(file,['reduction'])['reduction'];cf=read(s['coarseFactorFile'],['coarseFactor'])['coarseFactor']
    assert sha(s['coarseFactorFile'])=='fc7fcdf76174f97829eab02cc7b0405f42cdef327abeac2cb9ea691be076efaf'
    assert sha(s['packedFile'])=='beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'
    assert sha(s['graphFile'])=='b8ac74781499c0054ea49af77f3bad41f6ae6600dea9045264bba0230c5e4a04'
    packed=read(s['packedFile'],['packed'])['packed'];graph=read(s['graphFile'],['graph'])['graph']
    c=dict(packed=packed,Fgraph=graph,reduction=r,coarseFactor=cf,rhs=s['rhs'],value=s['value'],detail=s['detail'])
    result=check(c);image=Packed(packed)(s['value']);assert rel(image,s['image'])<1e-8
    out=dict(status='passed-independent-real-variable-block-maps-all-stored-outer-Arnoldi-columns-and-original-C-residual',
        processID=os.getpid(),result=result,internalCTargetMet=result['finalOriginalCRelativeResidual']<=1e-3,
        allOriginalCVs=int(packed['cellCount']),allFiveConservativeFieldsRetained=True,
        independentAllRealOuterArnoldiColumnsVerified=True,independentAllRealInnerLSQRColumnsClaimed=False)
out.update(newActualNSEvaluations=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0,physicalSolution=False,inputSHA256=bindings,scriptSHA256=sha(__file__))
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
