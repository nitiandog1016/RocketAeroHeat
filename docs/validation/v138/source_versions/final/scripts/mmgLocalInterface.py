"""MMG原始连接输入/输出与固定区域恢复；仅按显式reference匹配，不做最近邻补值。"""
import json
import numpy as np
from cgalTetrahedralContract import faces_and_neighbors
from prepareCGALRemeshingRuntime import ROOT,sha
from prepareMMGLocalRuntime import DEST
from buildCGALRemeshingRuntime import run

def write_input(path,P,T,selected):
    F,_=faces_and_neighbors(T);fixed=np.unique(np.r_[F.ravel(),T[~selected].ravel()])
    refs=np.zeros(len(P),int);refs[fixed]=fixed+1;required=np.isin(np.arange(len(P)),fixed)
    with path.open('w') as out:
        out.write('%d %d %d\n'%(len(P),len(T),len(F)))
        np.savetxt(out,np.c_[P,refs,required],fmt=['%.17g','%.17g','%.17g','%d','%d'])
        np.savetxt(out,np.c_[T+1,~selected],fmt='%d');np.savetxt(out,F+1,fmt='%d')
    return F,fixed,refs

def read_output(path):
    with path.open() as f:
        nv,nt,nf=map(int,f.readline().split())
        P=np.array([list(map(float,f.readline().split())) for _ in range(nv)])
        T=np.array([list(map(int,f.readline().split())) for _ in range(nt)])
        F=np.array([list(map(int,f.readline().split())) for _ in range(nf)])
        assert P.shape==(nv,5) and T.shape==(nt,5) and F.shape==(nf,5) and not f.read().strip()
    return P[:,:3],T[:,:4]-1,F[:,:3]-1,P[:,3].astype(int),T[:,4].astype(bool)

def sorted_cells(T):
    s=np.sort(T,axis=1);return s[np.lexsort(s.T[::-1])]

def restore(P,T,F,selected,fixed,raw):
    Q,U,G,ref,req=raw;protected=ref>0
    assert np.array_equal(np.sort(ref[protected]-1),fixed),'fixed node references changed or duplicated'
    assert np.all(ref>=0) and np.all(ref<=len(P))
    delta=np.max(np.abs(Q[protected]-P[ref[protected]-1]));tol=16*np.finfo(float).eps*max(1.,float(np.ptp(P,axis=0).max()))
    assert delta<=tol,'required node moved beyond floating point roundtrip'
    mapping=np.full(len(Q),-1,int);mapping[protected]=ref[protected]-1
    # 先保留原数组作为身份基准；旧自由点未必还存在，最后统一压缩未引用节点。
    free=np.flatnonzero(~protected);mapping[free]=np.arange(len(P),len(P)+len(free))
    reconstructed=np.r_[P,Q[free]];newT=mapping[U];newF=mapping[G]
    assert np.array_equal(sorted_cells(F),sorted_cells(newF)),'boundary connectivity changed'
    assert np.array_equal(sorted_cells(T[~selected]),sorted_cells(newT[req])),'required tetrahedra changed'
    nativeCoordinates=Q.copy();nativeCoordinates[protected]=P[ref[protected]-1]
    return reconstructed,newT,dict(fixedNodeCount=len(fixed),maximumNativeFixedNodeDifference=float(delta),
        roundtripTolerance=float(tol),boundaryExactlyRestored=True,requiredTetsExactlyPreserved=True,
        originalRequiredTets=int(sum(~selected)),nativeRequiredTets=int(sum(req)),
        nativeNodes=len(Q),nativeTets=len(U),freeOutputNodes=len(free),physicsQualified=False,newQualifiedArea=0)

def execute(folder,P,T,selected,timeout,sourcePaths):
    folder.mkdir(parents=True,exist_ok=False);inp=folder/'input.dat';F,fixed,refs=write_input(inp,P,T,selected)
    runtime=json.loads((DEST/'runtime.json').read_text());exe=DEST/'local_tetra_remesh.exe'
    assert sha(exe)==runtime['executableSHA256']
    paths=sourcePaths+[inp,exe,DEST/'runtime.json',ROOT/'scripts/mmgLocalInterface.py',ROOT/'scripts/buildCGALRemeshingRuntime.py',
        ROOT/'scripts/cgalTetrahedralContract.py',ROOT/'third_party/mmg/local_tetra_remesh.cpp']
    locked={p.resolve().as_posix():sha(p) for p in paths}
    receipt=run([str(exe),str(inp),str(folder/'native')],folder/'process',timeout,locked)
    assert receipt['exitCode']==0 and all(receipt['sourceStableAfter'].values()),'native MMG failed'
    ip,it,iff,ir,iq=read_output(folder/'native_import.dat')
    assert np.array_equal(ip,P) and np.array_equal(it,T) and np.array_equal(iff,F) and np.array_equal(ir,refs) and np.array_equal(iq,~selected)
    raw=read_output(folder/'native_result.dat');PP,TT,report=restore(P,T,F,selected,fixed,raw)
    report.update(exactImport=True,nativeSeconds=receipt['seconds'],nativeExitCode=receipt['exitCode'],selectedInputTets=int(sum(selected)))
    (folder/'interface.json').write_text(json.dumps(report,indent=2)+'\n')
    return PP,TT,report
