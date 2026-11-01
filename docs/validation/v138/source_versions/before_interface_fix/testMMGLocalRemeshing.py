"""带孔域上核对MMG固定边界/冻结邻域与真实局部质量改善。"""
from pathlib import Path
import json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from testCGALLocalRemeshing import hollow,orient
from mmgLocalInterface import execute
from cgalTetrahedralContract import faces_and_neighbors
from reviewSU2ExplicitQuality import quality
from auditSU2HybridTopology import paired
from prepareCGALRemeshingRuntime import ROOT,sha

def main():
    P,T=hollow();seed=13;tet=T[seed];point=P[tet].T@np.array([.333,.333,.333,.001]);nid=len(P);P=np.r_[P,point[None]]
    split=np.tile(tet,(4,1));split[np.arange(4),np.arange(4)]=nid
    T=np.r_[np.delete(T,seed,axis=0),orient(P,split)];selected=np.zeros(len(T),bool);selected[-4:]=True
    _,graph=faces_and_neighbors(T);selected|=np.asarray(graph@selected.astype(int)>0).ravel()
    folder=ROOT/'docs/validation/v138/hollow_test'
    paths=[Path(__file__).resolve(),ROOT/'tests/testCGALLocalRemeshing.py',ROOT/'scripts/reviewSU2ExplicitQuality.py',ROOT/'scripts/auditSU2HybridTopology.py']
    Q,U,r=execute(folder,P,T,selected,60,paths)
    boundary,_=faces_and_neighbors(U);allF=U[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3);paired(allF,boundary)
    q=Q[U];vol=np.linalg.det(np.stack([q[:,1]-q[:,0],q[:,2]-q[:,0],q[:,3]-q[:,0]],axis=2))/6
    assert min(vol)>0 and abs(sum(vol)-26)<2e-12
    q0=float(min(quality(P[T])));q1=float(min(quality(Q[U])));assert q1>q0 and q1>=.1
    r.update(passed=True,originalNodes=len(P),originalTets=len(T),selectedInputTets=int(sum(selected)),
        minimumSICNBefore=q0,minimumSICNAfter=q1,knownVolume=26,actualVolume=float(sum(vol)),knownHoleVolume=1,
        positiveVolumes=True,oppositeSharedFaces=True,filesSHA256={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file()})
    (folder/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
