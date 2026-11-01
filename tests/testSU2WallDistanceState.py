"""一般壁距初场的解析/穷举/置换检查；没有SU2或几何重生成。"""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2WallDistanceState import project_triangle,surface_distance,initial_values
from su2SSTInterface import triangle_distances
from prepareSU2CheckpointContinuation import sha

def main():
    t=np.array([[0.,0,0],[2,0,0],[0,1,0]])
    p=np.array([[.5,.25,3],[1,-2,0],[-1,-1,0],[2,1,0]])
    d,w=project_triangle(p,t);assert np.allclose(d,[9,4,2,.8],rtol=0,atol=1e-14)
    assert np.all(w>=0) and np.allclose(w.sum(axis=1),1)
    rng=np.random.RandomState(145);v=rng.normal(size=(153,3));F=np.arange(153).reshape(-1,3);q=rng.normal(size=(181,3))
    distance,face,weights,pairs=surface_distance(q,v,F,1.7)
    brute=np.array([min(triangle_distances(x,v[F[:,0]],v[F[:,1]],v[F[:,2]])) for x in q])
    error=float(max(abs(distance-np.minimum(brute,1.7))));assert error<2e-14
    # 单个斜三角壁面，面内投影点离全部顶点很远；保证不是最近顶点或层号距离。
    V=np.array([[0.,0,0],[2,0,.4],[0,1,.3]]);normal=np.cross(V[1]-V[0],V[2]-V[0]);normal/=np.linalg.norm(normal)
    P=np.r_[V,[(V.mean(axis=0)+z*normal) for z in [.01,.13,.9]]]
    wall=np.arange(3);far=np.array([5]);lengths=np.array([.001,.002,.003]);U=np.array([-1000.,20,400])
    args=(.2,U,1500.,230.,300.,287.,1.4,25.,180.,.0008)
    field,info=initial_values(P,P,np.array([[0,1,2]]),wall,lengths,far,*args)
    assert np.allclose(info['distance'][3:],[.01,.13,.2],rtol=0,atol=1e-15)
    rho=field[:,0];vel=field[:,1:4]/rho[:,None];temp=(field[:,4]/rho-.5*np.sum(vel**2,axis=1)-field[:,5])*.4/287.
    assert np.allclose(rho*287*temp,1500,rtol=1e-14) and np.allclose(temp[wall],300,rtol=1e-14)
    assert np.all(vel[wall]==0) and np.allclose(vel[far],U,rtol=1e-14) and np.allclose(temp[far],230,rtol=1e-14)
    assert np.allclose(field[wall,6],60*.0008/(.075*lengths**2),rtol=1e-14)
    assert np.all(field[3:,5]>0) and np.all(field[:,6]>0)
    perm=np.array([4,0,5,2,1,3]);inv=np.argsort(perm)
    permuted,_=initial_values(P[perm],P[perm],inv[np.array([[0,1,2]])],inv[wall],lengths,inv[far],*args)
    assert np.allclose(permuted[inv],field,rtol=1e-14,atol=1e-14)
    rejected=False
    try:initial_values(P,P,np.array([[0,1,2]]),wall,lengths,np.array([3]),*args)
    except AssertionError:rejected=True
    assert rejected
    out=ROOT/'docs/validation/v145/unit';out.mkdir(parents=True,exist_ok=False)
    report=dict(passed=True,analyticFaceEdgeVertex=True,bruteDistanceMaxErrorM=error,permutationInvariant=True,
        exactFarfield=True,invalidFarfieldRejected=True,positiveThermalEnergyIncludingK=True,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__),ROOT/'scripts/prepareSU2WallDistanceState.py',ROOT/'scripts/su2SSTInterface.py']})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
