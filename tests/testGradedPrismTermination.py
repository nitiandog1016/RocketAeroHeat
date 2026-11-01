"""已知体积平面渐变补片：独立边界、两种终止模板与近壁原列保持。"""
from pathlib import Path
import json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from terminatePrismLayers import terminate,mixed_faces,graph_layers
from meshSU2HybridVolume import wedge_jacobian_min
from su2PyramidGeometry import pyramid_geometry
from auditSU2HybridTopology import paired,polyvolume
from reviewSU2ExplicitQuality import quality
from prepareCGALRemeshingRuntime import ROOT,sha

def main():
    out=ROOT/'docs/validation/v140/fixture';out.mkdir(parents=True,exist_ok=False)
    V=np.array([[float(x),float(y),0.] for x in range(4) for y in range(2)])
    F=np.array([[2*x,2*x+2,2*x+3] for x in range(3)]+[[2*x,2*x+3,2*x+1] for x in range(3)],np.int64)
    h=np.arange(4);P=np.concatenate([V+[0,0,float(z)] for z in h]);nv=len(V)
    W=np.concatenate([np.c_[F+j*nv,F+(j+1)*nv] for j in range(3)])
    K,d=graph_layers(F,nv,V[:,0]<=1,1,3);assert np.array_equal(K,np.repeat([1,1,2,3],2))
    r=terminate(F,W,K,3);w,q,t=[r[n] for n in ['wedges','pyramids','tetrahedra']]
    assert (len(w),len(q),len(t))==(8,2,2)
    np.savez_compressed(str(out/'fixture.npz'),vertices=P,sourceVertices=V,sourceFaces=F,originalWedges=W,**r)
    # 独立地按补片边及原高度构造封闭外边界，不从体单元计数推断预期答案。
    E=F[:,[[0,1],[1,2],[2,0]]].reshape(-1,2);_,i,c=np.unique(np.sort(E,axis=1),axis=0,return_index=True,return_counts=True)
    tris=list(F[:,[0,2,1]])+list(F+K[F]*nv);quads=[]
    for a,b in E[i[c==1]]:
        for j in range(max(K[a],K[b])):
            cycle=[a+min(j,K[a])*nv,b+min(j,K[b])*nv,b+min(j+1,K[b])*nv,a+min(j+1,K[a])*nv]
            unique=list(dict.fromkeys(cycle))
            if len(unique)==3:tris.append(unique)
            elif len(unique)==4:quads.append(unique)
    tri,quad=mixed_faces(w,q,t);tr=paired(tri,np.array(tris));qr=paired(quad,np.array(quads))
    jac=float(min(wedge_jacobian_min(P[w])));pv,pg=pyramid_geometry(P,q);tv=np.linalg.det(np.stack([P[t[:,j]]-P[t[:,0]] for j in [1,2,3]],axis=2))/6
    wp=P[w];wv=np.cross(wp[:,1]-wp[:,0],wp[:,2]-wp[:,0])[:,2]*(wp[:,3,2]-wp[:,0,2])/2
    total=float(sum(wv)+sum(pv)+sum(tv));assert abs(total-5)<2e-12 and jac>0 and min(tv)>0
    sicn=float(min(quality(P[t])));assert sicn>=.1
    assert np.array_equal(w[:len(F)],W[:len(F)])
    report=dict(passed=True,knownVolume=5,computedVolume=total,prisms=len(w),pyramids=len(q),tetrahedra=len(t),
        minimumPrismJacobian=jac,pyramidGeometry=pg,minimumTetSICN=sicn,triangular=tr,quadrilateral=qr,
        originalFirstLayerExactlyPreserved=True,bothTerminationTemplatesTested=True,newVertices=0,
        physicsQualified=False,newQualifiedArea=0,fixtureSHA256=sha(out/'fixture.npz'))
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
