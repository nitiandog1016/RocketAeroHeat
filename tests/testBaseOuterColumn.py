"""固定近壁及外环的底面内部列重分布：解析体积与全参考正性。"""
from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from recenterBaseOuterColumn import recenter
from meshSU2HybridVolume import wedge_jacobian_min
from auditSU2HybridTopology import paired
from prepareTetGenRuntime import ROOT,sha

def main():
    start=time.monotonic();folder=ROOT/'docs/validation/v144/fixture';folder.mkdir(parents=True,exist_ok=False)
    # 有向截面法向-X，沿-X外推；内部点偏心但真实截面始终同一正方形。
    V=np.array([[0,-1.,-1],[0,-1,1],[0,1,1],[0,1,-1],[0,.75,.2]])
    F=np.array([[4,0,1],[4,1,2],[4,2,3],[4,3,0]],np.int64);h=np.array([0,.1,.25,.45,.75]);n=len(V)
    P=np.concatenate([V+[-z,0,0] for z in h]);W=np.concatenate([np.c_[F+j*n,F+(j+1)*n] for j in range(len(h)-1)])
    E=np.array([[0,1],[1,2],[2,3],[3,0]])+(len(h)-1)*n;ids=np.arange(len(h))*n+4
    Q,offset,r=recenter(P,ids,E,h,2);jac=float(min(wedge_jacobian_min(Q[W])));assert jac>0
    assert np.array_equal(Q[:3*n],P[:3*n]) and np.array_equal(Q[E],P[E])
    assert np.max(abs(Q[ids[-1],1:]))<1e-14 and r['maximumDisplacementM']>0
    tri=np.r_[W[:,[0,2,1]],W[:,[3,4,5]]];quad=np.r_[W[:,[0,1,4,3]],W[:,[1,2,5,4]],W[:,[2,0,3,5]]]
    bquad=np.concatenate([np.c_[E-4*n+j*n,E[:,::-1]-4*n+(j+1)*n] for j in range(4)])
    tr=paired(tri,np.r_[F[:,[0,2,1]],F+4*n]);qr=paired(quad,bquad)
    # 三角形一阶+轴向两点Gauss对线性棱柱detJ足够精确，独立于最小Jacobian函数。
    volume=0.;low=Q[W[:,:3]];d=Q[W[:,3:]]-low
    for z in [.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)]:
        a=low[:,1]-low[:,0]+z*(d[:,1]-d[:,0]);b=low[:,2]-low[:,0]+z*(d[:,2]-d[:,0]);c=d.mean(axis=1)
        volume+=float(np.sum(np.linalg.det(np.stack([a,b,c],axis=2))))*.25
    assert abs(volume-3)<1e-12
    np.savez_compressed(str(folder/'fixture.npz'),originalVertices=P,vertices=Q,wedges=W,offsets=offset,columnIDs=ids,rimEdgesTop=E,heights=h)
    paths=[Path(__file__).resolve(),ROOT/'scripts/recenterBaseOuterColumn.py',ROOT/'scripts/meshSU2HybridVolume.py',ROOT/'scripts/auditSU2HybridTopology.py',ROOT/'docs/validation/v144/plan.md']
    r.update(passed=True,knownVolume=3,actualVolume=volume,minimumPrismJacobian=jac,triangles=tr,quadrilaterals=qr,
        seconds=time.monotonic()-start,sourceSHA256={str(p):sha(p) for p in paths},fixtureSHA256=sha(folder/'fixture.npz'))
    (folder/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
