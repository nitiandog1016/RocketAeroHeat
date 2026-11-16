"""独立SVD复核显式过渡的低质单元及边界约束；单面上界不冒充全域可网格性。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
from prepareSU2CheckpointContinuation import ROOT,sha,read

REFERENCE=np.array([[1,.5,.5],[0,np.sqrt(3)/2,np.sqrt(3)/6],[0,0,np.sqrt(2/3)]])
def quality(q):
    J=np.stack([q[:,1]-q[:,0],q[:,2]-q[:,0],q[:,3]-q[:,0]],axis=2)@np.linalg.inv(REFERENCE)
    s=np.linalg.svd(J,compute_uv=False)
    return 3*np.sign(np.linalg.det(J))/np.sqrt(np.sum(s*s,axis=1)*np.sum(s**-2,axis=1))
def bound(tri):
    a=tri[:,1]-tri[:,0];b=(tri[:,2]-tri[:,0]-a/2)*2/np.sqrt(3)
    s=np.linalg.svd(np.stack([a,b],axis=2),compute_uv=False);ratio=s[:,0]/s[:,1]
    return 3/(1+ratio+1/ratio),a,b,s
def main():
    start=time.monotonic();folder=ROOT/'outputs/v136_validation/explicit';d=loadmat(str(folder/'hybrid_volume.mat'));prepared=loadmat(str(folder/'prepared_transition.mat'))
    # 两个固定三角面的最优法向顶点达到解析上界；随机顶点仅用于反例筛查。
    rng=np.random.RandomState(136);tests=[]
    for tri in [np.array([[0.,0,0],[1,0,0],[.5,np.sqrt(3)/2,0]]),np.array([[0.,0,0],[1,0,0],[.2,.03,0]])]:
        upper,a,b,s=bound(tri[None]);normal=np.cross(a[0],b[0]);normal/=np.linalg.norm(normal)
        height=np.sqrt(np.prod(s[0]));apex=tri[0]+.5*a[0]+np.sqrt(3)/6*b[0]+np.sqrt(2/3)*height*normal
        optimum=quality(np.r_[tri,apex[None]][None])[0];assert abs(optimum-upper[0])<1e-12
        arbitrary=tri[0]+rng.normal(size=(100,3));q=quality(np.concatenate([np.tile(tri,(100,1,1)),arbitrary[:,None]],axis=1))
        assert max(abs(q))<=upper[0]+1e-12;tests.append(dict(bound=float(upper[0]),attained=float(optimum),randomApexCount=100))
    P=d['vertices'];T=d['tetrahedra'].astype(int)-1;Q=d['pyramids'].astype(int)-1;W=d['wedges'].astype(int)-1
    values=quality(P[T]);bad=np.flatnonzero(values<.1);fixed=np.unique(np.r_[W.ravel(),Q.ravel(),d['farfield'].ravel()-1])
    count=np.isin(T[bad],fixed).sum(axis=1);shellP=prepared['vertices'];F=prepared['exteriorInnerFaces'].astype(int)-1;upper,_,_,_=bound(shellP[F])
    r=read(folder/'candidate_report.json');assert len(bad)==r['tetsBelowSICN01'] and abs(min(values)-r['minimumTetSICN'])<1e-12
    report=dict(minimumTetSICN=float(min(values)),badTets=len(bad),badTetRowsOneBased=(bad+1).tolist(),badTetSICN=values[bad].tolist(),
        badTetCentersM=P[T[bad]].mean(axis=1).tolist(),fixedVertexCountHistogram=np.bincount(count,minlength=5).tolist(),
        shellFaces=len(F),minimumSingleFaceOptimalTetSICN=float(min(upper)),shellFacesWithBoundBelow01=int(sum(upper<.1)),
        boundFormula='q_max=3/(1+kappa+1/kappa), kappa=equilateral-reference face Jacobian singular-value ratio',
        boundScope='one isolated fixed triangular face only; not a global mesh existence or quality proof',analyticFixtures=tests,
        seconds=time.monotonic()-start,physicsQualified=False,newQualifiedArea=0,
        sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in [folder/'hybrid_volume.mat',folder/'prepared_transition.mat',Path(__file__).resolve()]})
    target=folder/'independent_quality.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['badTetRowsOneBased','badTetSICN','badTetCentersM','sourceSHA256']},indent=2))
if __name__=='__main__':main()
