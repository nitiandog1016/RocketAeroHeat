"""模块职责：只读定位差分缺口，检验它能否由实际面通量传播解释。"""
from pathlib import Path
import json,time
import numpy as np
from readEilmerNative import read_grid,sha
from auditFrozenMappedGeometryV213 import FACES
from auditFrozenOperatorV213 import scaled
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v213_validation'

def main():
    start=time.monotonic();case=BASE/'control03';grid=read_grid(case/'lmrsim/grid/grid-0000')
    g=np.genfromtxt(str(case/'geometry.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
    c=g[g['kind']=='cell'];f=g[g['kind']=='face'];vol=c['volumeOrArea']
    a=np.genfromtxt(str(case/'samples.csv'),delimiter=',',names=True,dtype=None,encoding='utf8')
    keys={tuple(sorted(row[row>=0])):i for i,row in enumerate(grid['faces'])}
    cell_faces=np.array([[keys[tuple(sorted(row[list(local)]))] for local in FACES] for row in grid['cells'][:,:6]])
    xyz=lambda t:np.stack([t[n] for n in ['x','y','z']],axis=1)
    n=np.stack([f[k] for k in ['nx','ny','nz']],axis=1)
    signs=np.sign(np.einsum('ijk,ijk->ij',xyz(f)[cell_faces]-xyz(c)[:,None,:],n[cell_faces])).astype(int)
    assert np.all(np.abs(signs)==1)
    tags=np.full(1197,'internal',dtype='U16')
    for tag,b in grid['boundaries'].items():
        tags[b['ids']]=tag
        for fid,s,owner in zip(b['ids'],b['signs'],b['owners']):assert signs[owner][cell_faces[owner]==fid][0]==s
    records=[]
    for regime in range(3):
        r=a[a['regime']==regime];rr=r[r['kind']=='residual'];ff=r[r['kind']=='integratedFlux']
        residual=rr['frozen'].reshape(416,6);flux=ff['frozen'].reshape(1197,6)
        direct=-np.sum(flux[cell_faces]*signs[:,:,None],axis=1)/vol[:,None]
        # 仅前五式无体源；第六SA式含生成/耗散，不能套无源账。
        balance=np.abs(direct-residual)[:,:5]/np.maximum(1,np.sum(np.abs(flux[cell_faces]),axis=1)[:,:5]/vol[:,None])
        assert balance.max()<1e-11
        fd=(ff['plus']-ff['minus'])/2e-4;diff=fd-ff['frozenDerivative'];err=scaled(fd,ff['frozenDerivative'])
        worst=int(np.argmax(err));fid=int(ff[worst]['id']);comp=int(ff[worst]['component'])
        base=float(ff[worst]['frozen']);der=float(ff[worst]['frozenDerivative'])
        failures={t:dict(failedComponents=int(((tags[ff['id']]==t)&(err>1e-5)).sum()),
            maximum=float(err[tags[ff['id']]==t].max())) for t in ['wall','outer','exchange','internal']}
        records.append(dict(regime=regime,cellBalanceScaled=float(balance.max()),
            faceFailureByBoundary=failures,worstFaceID=fid,worstFaceBoundary=str(tags[fid]),component=comp,
            complexDerivative=der,forwardDifference=(float(ff[worst]['plus'])-base)/1e-4,
            backwardDifference=(base-float(ff[worst]['minus']))/1e-4,
            cancellationFreeUniqueCauseProven=False))
    result=dict(passed=True,records=records,sourceSHA256=sha(case/'samples.csv'),
        note='Base and complex response match resident donors; real plus/minus resident samples were not collected. Gradient finite differences pass. Root cause is unresolved.',
        actualNativeNKMatVecInvoked=False,nativeEntryDifference='probe directly calls evalResidualWorker with primitive perturbations; native NK calls setJacobianEvalSettings then perturbs U[1] and decode_conserved.',
        seconds=time.monotonic()-start)
    out=BASE/'derivative_mismatch_corrected.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
