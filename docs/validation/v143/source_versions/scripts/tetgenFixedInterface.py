"""TetGen固定表面接口：准确连接导入、显式点标记恢复与完整域前置检查。"""
from pathlib import Path
import itertools,json
import numpy as np
from prepareTetGenRuntime import ROOT,DEST,sha
from buildCGALRemeshingRuntime import run
from cgalTetrahedralContract import faces_and_neighbors
from auditSU2HybridTopology import paired,polyvolume

def sorted_cells(T):
    s=np.sort(T,axis=1);return s[np.lexsort(s.T[::-1])]

def read_native(path):
    with path.open() as f:
        nv,nt,nf=map(int,f.readline().split())
        assert 4<=nv<=410000 and 1<=nt<=1000000 and 4<=nf<=900000
        P=np.loadtxt(itertools.islice(f,nv),ndmin=2)
        T=np.loadtxt(itertools.islice(f,nt),dtype=np.int64,ndmin=2)
        F=np.loadtxt(itertools.islice(f,nf),dtype=np.int64,ndmin=2)
        assert P.shape==(nv,4) and T.shape==(nt,4) and F.shape==(nf,4) and not f.read().strip()
    assert np.isfinite(P).all() and np.array_equal(P[:,3],P[:,3].astype(np.int64))
    assert np.min(T)>=0 and np.max(T)<nv and np.min(F[:,:3])>=0 and np.max(F[:,:3])<nv
    return P[:,:3],T,F[:,:3],P[:,3].astype(np.int64),F[:,3]

def tet_audit(P,T):
    assert len(np.unique(np.sort(T,axis=1),axis=0))==len(T),'duplicate tetrahedra'
    q=P[T];vol=np.einsum('ij,ij->i',q[:,1]-q[:,0],np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0]))/6
    assert np.all(vol>0),'native tetrahedra are not directly positive'
    F,_=faces_and_neighbors(T)
    pair=paired(T[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3),F)
    total=float(np.sum(vol));boundary=polyvolume(P,F);err=abs(total-boundary)/total
    assert err<=2e-10
    return F,dict(pair,volume=total,boundaryVolume=boundary,volumeClosureRelative=err,minimumVolume=float(min(vol)))

def execute(folder,P,T,timeout,sourcePaths):
    folder.mkdir(parents=True,exist_ok=False);F,before=tet_audit(P,T);fixed=np.unique(F)
    refs=np.zeros(len(P),np.int64);refs[fixed]=fixed+1
    inp=folder/'input.dat'
    with inp.open('w') as f:
        f.write('%d %d %d\n'%(len(P),len(T),len(F)))
        np.savetxt(f,np.c_[P,refs],fmt=['%.17g']*3+['%d'])
        np.savetxt(f,T,fmt='%d');np.savetxt(f,np.c_[F,np.ones(len(F),np.int64)],fmt='%d')
    runtime=json.loads((DEST/'runtime.json').read_text());exe=Path(runtime['executable'])
    assert sha(exe)==runtime['executableSHA256']
    paths=sourcePaths+[inp,exe,DEST/'runtime.json',ROOT/'scripts/tetgenFixedInterface.py',
        ROOT/'scripts/buildCGALRemeshingRuntime.py',ROOT/'scripts/cgalTetrahedralContract.py',
        ROOT/'scripts/auditSU2HybridTopology.py',ROOT/'third_party/tetgen/fixed_tetra_repair.cpp']
    locked={str(p.resolve()):sha(p) for p in paths}
    receipt=run([str(exe),str(inp),str(folder/'native')],folder/'process',timeout,locked)
    assert receipt['exitCode']==0 and all(receipt['sourceStableAfter'].values()),'native TetGen failed'
    ip,it,iff,ir,im=read_native(folder/'native_import.dat')
    assert np.array_equal(ip,P) and np.array_equal(it,T) and np.array_equal(iff,F) and np.array_equal(ir,refs) and np.all(im==1)
    Q,U,G,ref,mark=read_native(folder/'native_result.dat')
    actualF,after=tet_audit(Q,U)
    assert np.array_equal(sorted_cells(actualF),sorted_cells(G)) and np.all(mark==1)
    assert len(Q)<=len(P)+10000
    protected=ref>0;assert np.all(ref>=0) and np.all(ref<=len(P))
    assert np.array_equal(np.sort(ref[protected]-1),fixed),'fixed point identities missing or duplicated'
    assert np.all(protected[np.unique(G)]),'new or unmarked boundary point'
    delta=float(np.max(np.abs(Q[protected]-P[ref[protected]-1])))
    tol=16*np.finfo(float).eps*max(1.,float(np.ptp(P,axis=0).max()))
    assert delta<=tol,'fixed point moved beyond roundtrip tolerance'
    mapping=np.full(len(Q),-1,np.int64);mapping[protected]=ref[protected]-1
    free=np.flatnonzero(~protected);mapping[free]=np.arange(len(P),len(P)+len(free),dtype=np.int64)
    R=np.r_[P,Q[free]];newT=mapping[U];newF=mapping[G]
    assert np.array_equal(sorted_cells(newF),sorted_cells(F)),'surface triangles changed'
    assert abs(after['volume']/before['volume']-1)<=2e-10
    report=dict(exactImport=True,boundaryExactlyRestored=True,fixedNodeCount=len(fixed),
        maximumNativeFixedNodeDifferenceM=delta,roundtripToleranceM=tol,nativeNodes=len(Q),nativeTets=len(U),
        netAdditionalNativeNodes=len(Q)-len(P),freeOutputNodes=len(free),nativeSeconds=receipt['seconds'],
        nativeExitCode=receipt['exitCode'],nativeBefore=before,nativeAfter=after,
        switches='rY1q2/5O2/7S10000zVC',physicsQualified=False,newQualifiedArea=0)
    (folder/'interface.json').write_text(json.dumps(report,indent=2)+'\n')
    return R,newT,report
