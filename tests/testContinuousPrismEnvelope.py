"""用非对称正八面体包络的解析体积，核对连续外层变形及完整棱柱正性。"""
from pathlib import Path
import json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from scipy.spatial import ConvexHull
from redistributePrismEnvelope import redistribute
from meshSU2HybridVolume import prepare_layers,wedge_jacobian_min
from auditSU2HybridTopology import audit
from prepareCGALRemeshingRuntime import ROOT,sha

def main():
    folder=ROOT/'docs/validation/v141/fixture';folder.mkdir(parents=True,exist_ok=False)
    V=np.array([[1.,0,0],[-1.,0,0],[0,1.,0],[0,-1.,0],[0,0,1.],[0,0,-1.]])
    F=ConvexHull(V).simplices.copy();p=V[F];flip=np.einsum('ij,ij->i',np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),p.mean(axis=1))<0
    F[flip]=F[flip][:,[0,2,1]];h=np.array([0,.1,.25,.45,.75]);P,W,_,_,_,_=prepare_layers(V,F,dict(layerHeights=h))
    Q,H,r=redistribute(P,V,h,2,-.5,.5);top=F+4*len(V)
    known=(1.55+1.75)*(1.65+1.65)**2/6-4/3
    assert np.allclose(Q[-6:],V*np.array([1.55,1.75,1.65,1.65,1.65,1.65])[:,None],rtol=0,atol=1e-14)
    data=dict(vertices=Q,wedges=W+1,tetrahedra=np.empty((0,4),int),farfield=top+1,sourceVertices=V,sourceFaces=F+1)
    geom=audit(data);jac=float(min(wedge_jacobian_min(Q[W])));assert jac>0 and abs(geom['totalFluidVolume']-known)<2e-12
    np.savez_compressed(str(folder/'fixture.npz'),originalVertices=P,newHeights=H,**data)
    r.update(passed=True,knownVolume=known,actualVolume=geom['totalFluidVolume'],geometry=geom,minimumPrismJacobian=jac,
        nonzeroEnvelopeChange=True,fixtureSHA256=sha(folder/'fixture.npz'))
    (folder/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
