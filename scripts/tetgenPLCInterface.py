"""封闭双表面PLC、显式孔种子和原边界身份合同；不使用旧体内点或最近邻。"""
from pathlib import Path
import itertools,json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from prepareTetGenRuntime import ROOT,sha
from buildTetGenPLCRuntime import PLC
from buildCGALRemeshingRuntime import run
from cgalTetrahedralContract import faces_and_neighbors
from tetgenFixedInterface import sorted_cells,read_native,tet_audit

def winding(P,F,point):
    a,b,c=(P[F]-point).transpose(1,0,2);la=np.linalg.norm(a,axis=1);lb=np.linalg.norm(b,axis=1);lc=np.linalg.norm(c,axis=1)
    assert min(np.r_[la,lb,lc])>1e-14
    numerator=np.einsum('ij,ij->i',a,np.cross(b,c))
    denominator=la*lb*lc+np.einsum('ij,ij->i',a,b)*lc+np.einsum('ij,ij->i',b,c)*la+np.einsum('ij,ij->i',c,a)*lb
    return float(np.sum(2*np.arctan2(numerator,denominator))/(4*np.pi))

def plc_contract(P,F,hole,fluidPoint):
    assert np.array_equal(np.unique(F),np.arange(len(P))) and len(np.unique(np.sort(F,axis=1),axis=0))==len(F)
    edges=np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]);key,inv,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    assert np.all(count==2) and np.all(np.bincount(inv,weights=np.sign(edges[:,1]-edges[:,0]))==0),'surface edge pairing failure'
    graph=coo_matrix((np.ones(2*len(edges)),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(P),len(P))).tocsr()
    nc,label=connected_components(graph,directed=False);assert nc==2
    components=[]
    for i in range(nc):
        faces=F[label[F[:,0]]==i];assert np.all(label[faces]==i)
        signed=float(np.sum(np.einsum('ij,ij->i',P[faces[:,0]],np.cross(P[faces[:,1]],P[faces[:,2]])))/6)
        wh=winding(P,faces,hole);wf=winding(P,faces,fluidPoint)
        assert abs(abs(wh)-1)<1e-8,'hole seed is not strictly inside both shells'
        assert abs(wf-(1 if signed>0 else 0))<1e-8,'fluid witness is not in intended domain'
        components.append(dict(nodes=int(sum(label==i)),faces=len(faces),signedVolume=signed,holeWinding=wh,fluidWinding=wf))
    assert sum(c['signedVolume']>0 for c in components)==1 and abs(sum(c['holeWinding'] for c in components))<1e-8
    volume=sum(c['signedVolume'] for c in components);assert volume>0
    return dict(components=components,expectedFluidVolume=volume,holeSeed=hole.tolist(),fluidWitness=fluidPoint.tolist(),
        closedOrientedTwoShells=True,sourceBoundaryNodes=len(P),sourceBoundaryFaces=len(F))

def read_plc(path):
    with path.open() as f:
        nv,nf,nh=map(int,f.readline().split())
        P=np.loadtxt(itertools.islice(f,nv),ndmin=2);F=np.loadtxt(itertools.islice(f,nf),dtype=np.int64,ndmin=2)
        H=np.loadtxt(itertools.islice(f,nh),ndmin=2);assert not f.read().strip()
    assert P.shape==(nv,4) and F.shape==(nf,4) and H.shape==(1,3)
    return P,F,H

def execute(folder,P,F,hole,fluidPoint,timeout,sourcePaths):
    folder.mkdir(parents=True,exist_ok=False);contract=plc_contract(P,F,hole,fluidPoint)
    (folder/'preflight.json').write_text(json.dumps(contract,indent=2)+'\n')
    inp=folder/'input.dat';refs=np.arange(1,len(P)+1,dtype=np.int64)
    with inp.open('w') as f:
        f.write('%d %d 1\n'%(len(P),len(F)))
        np.savetxt(f,np.c_[P,refs],fmt=['%.17g']*3+['%d'])
        np.savetxt(f,np.c_[F,np.ones(len(F),np.int64)],fmt='%d');np.savetxt(f,hole[None],fmt='%.17g')
    runtime=json.loads((PLC/'runtime.json').read_text());exe=Path(runtime['executable']);assert sha(exe)==runtime['executableSHA256']
    paths=sourcePaths+[inp,PLC/'runtime.json',exe,Path(__file__).resolve(),ROOT/'third_party/tetgen/fixed_surface_plc.cpp',
        ROOT/'scripts/buildCGALRemeshingRuntime.py',ROOT/'scripts/tetgenFixedInterface.py',ROOT/'scripts/cgalTetrahedralContract.py',
        ROOT/'scripts/auditSU2HybridTopology.py',ROOT/'scripts/buildTetGenPLCRuntime.py']
    locked={str(p.resolve()):sha(p) for p in paths}
    receipt=run([str(exe),str(inp),str(folder/'native')],folder/'process',timeout,locked)
    assert receipt['exitCode']==0 and all(receipt['sourceStableAfter'].values()),'native PLC failed'
    ip,iff,ih=read_plc(folder/'native_import.dat')
    assert np.array_equal(ip[:,:3],P) and np.array_equal(ip[:,3],refs) and np.array_equal(iff[:,:3],F) and np.all(iff[:,3]==1) and np.array_equal(ih,hole[None])
    Q,T,G,ref,mark=read_native(folder/'native_result.dat');actualF,ga=tet_audit(Q,T)
    assert np.array_equal(sorted_cells(actualF),sorted_cells(G)) and np.all(mark==1)
    fixed=ref>0;assert np.all(ref>=0) and np.array_equal(np.sort(ref[fixed]),refs)
    assert np.all(fixed[np.unique(G)]) and sum(~fixed)<=30000
    delta=float(np.max(np.abs(Q[fixed]-P[ref[fixed]-1])));tol=16*np.finfo(float).eps*max(1.,float(np.ptp(P,axis=0).max()))
    assert delta<=tol
    mapping=np.full(len(Q),-1,np.int64);mapping[fixed]=ref[fixed]-1;mapping[~fixed]=np.arange(len(P),len(P)+sum(~fixed),dtype=np.int64)
    R=np.r_[P,Q[~fixed]];U=mapping[T];newF=mapping[G]
    assert np.array_equal(sorted_cells(F),sorted_cells(newF)),'fixed surface changed'
    _,graph=faces_and_neighbors(T);nc,_=connected_components(graph,directed=False);assert nc==1
    J=np.stack([Q[T[:,1]]-Q[T[:,0]],Q[T[:,2]]-Q[T[:,0]],Q[T[:,3]]-Q[T[:,0]]],axis=2)
    bary=np.linalg.solve(J,hole-Q[T[:,0]])
    contains=np.all(bary>=-1e-12,axis=1)&(np.sum(bary,axis=1)<=1+1e-12);assert not np.any(contains),'hole seed filled'
    err=abs(ga['volume']/contract['expectedFluidVolume']-1);assert err<=2e-10
    report=dict(contract,exactPLCImport=True,fixedBoundaryExactlyRestored=True,maximumNativeFixedNodeDifferenceM=delta,
        roundtripToleranceM=tol,nativeNodes=len(Q),nativeTets=len(T),newInteriorNodes=int(sum(~fixed)),
        nativeSeconds=receipt['seconds'],nativeExitCode=receipt['exitCode'],geometry=ga,
        connectedFluidComponents=int(nc),holeSeedInFluid=False,volumeRelativeDifference=err,
        switches='pY1q2/5O2/7S30000zVC',physicsQualified=False,newQualifiedArea=0)
    (folder/'interface.json').write_text(json.dumps(report,indent=2)+'\n');return R,U,report
