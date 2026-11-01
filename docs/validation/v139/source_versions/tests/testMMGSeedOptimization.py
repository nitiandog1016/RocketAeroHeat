"""固定外边界的偏心双锥：验证完整边星细分和保留新点的MMG优化。"""
from pathlib import Path
import json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from splitClosedTetEdges import split,volumes
from mmgLocalInterface import execute,sorted_cells,read_output
from cgalTetrahedralContract import faces_and_neighbors
from reviewSU2ExplicitQuality import quality
from auditSU2HybridTopology import paired
from prepareCGALRemeshingRuntime import ROOT,sha

def positive(P,T):
    T=T.copy();neg=volumes(P,T)<0;T[neg]=T[neg][:,[1,0,2,3]]
    assert np.all(volumes(P,T)>0)
    return T

def fixture():
    # 双锥体积8/3；偏心内边靠近右壁，固定表面不含退化三角形。
    P=np.array([[.999,0,-1],[.999,0,1],[1,1,0],[-1,1,0],[-1,-1,0],[1,-1,0]])
    T=positive(P,np.array([[0,1,2+i,2+(i+1)%4] for i in range(4)],dtype=np.int64))
    F,_=faces_and_neighbors(T);ids=np.flatnonzero(np.all(np.sort(F,axis=1)==[0,3,4],axis=1));assert len(ids)==1
    tri=P[F[ids[0]]];n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);area=np.linalg.norm(n)/2
    # 由扇区中心辨别外侧，再接一个冻结四面体；高0.3，解析体积area*0.3/3。
    n/=np.linalg.norm(n);center=tri.mean(axis=0)
    if np.dot(n,P.mean(axis=0)-center)>0:n=-n
    P=np.r_[P,(center+.3*n)[None]];T=np.r_[T,positive(P,np.array([[0,3,4,6]]))]
    return P,T,np.array([True]*4+[False]),8/3+area*.3/3

def topology(P,T):
    F,_=faces_and_neighbors(T)
    return paired(T[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3),F)

def main():
    folder=ROOT/'docs/validation/v139/seed_test';folder.mkdir(parents=True,exist_ok=False)
    P,T,selected,known=fixture();before=topology(P,T)
    assert abs(sum(volumes(P,T))-known)<2e-12
    S,U,mask,parent,sr=split(P,T,selected,[(0,1)])
    topology(S,U);assert len(S)==len(P)+1 and len(U)==9
    assert np.array_equal(np.bincount(parent,weights=volumes(S,U)).shape,volumes(P,T).shape)
    assert np.allclose(np.bincount(parent,weights=volumes(S,U)),volumes(P,T),rtol=2e-12,atol=0)
    np.savez_compressed(str(folder/'fixture.npz'),vertices=P,tetrahedra=T,selected=selected,
        seededVertices=S,seededTetrahedra=U,seededSelected=mask,parentTetRows=parent)
    (folder/'split_report.json').write_text(json.dumps(sr,indent=2)+'\n')
    paths=[Path(__file__).resolve()]+[ROOT/'scripts'/name for name in ['splitClosedTetEdges.py','reviewSU2ExplicitQuality.py','auditSU2HybridTopology.py']]
    Q,V,r=execute(folder/'native',S,U,mask,60,paths,ROOT/'tools/runtime/mmg/5.8.0-seeds/runtime.json')
    topo=topology(Q,V);vp=volumes(Q,V);q0=float(min(quality(P[T])));qs=float(min(quality(S[U])));q1=float(min(quality(Q[V])))
    raw=read_output(folder/'native/native_result.dat');free=raw[0][raw[3]==0]
    movement=float(np.linalg.norm(free[0]-S[-1])) if len(free)==1 else None
    passed=len(free)==1 and movement>1e-6 and q1>q0 and q1>=.1 and min(vp)>0 and abs(sum(vp)-known)<2e-12
    r.update(passed=bool(passed),knownVolume=float(known),actualVolume=float(sum(vp)),minimumSICNBefore=q0,
        minimumSICNSeeded=qs,minimumSICNAfter=q1,seedBefore=S[-1].tolist(),freePointsAfter=free.tolist(),
        seedMovement=movement,inputFreeNodes=1,seedCountPreserved=len(free)==1,topology=topo,
        positiveVolumes=bool(min(vp)>0),lineageConservation=True,physicsQualified=False,newQualifiedArea=0,
        filesSHA256={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file()})
    (folder/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
    assert passed,'fixed-boundary seed test did not meet predeclared gates'

if __name__=='__main__':main()
