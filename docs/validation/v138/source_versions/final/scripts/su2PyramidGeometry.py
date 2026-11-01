"""五节点金字塔的全参考Jacobian与中位面体积，支持非平面四边底。"""
import numpy as np

FACES=[[0,3,2,1],[0,1,4],[1,2,4],[2,3,4],[3,0,4]]

def pyramid_geometry(P,Q):
    """去除顶点自然退化因子(1-t)^2后，Jacobian双线性，极值在四角。"""
    p=P[Q];values=[]
    a=p[:,1]-p[:,0];b=p[:,3]-p[:,0];c=p[:,2]-p[:,1]-p[:,3]+p[:,0]
    for r,s in [(0.,0.),(1.,0.),(1.,1.),(0.,1.)]:
        point=p[:,0]+r*a+s*b+r*s*c
        values.append(np.einsum('ij,ij->i',np.cross(a+s*c,b+r*c),p[:,4]-point))
    values=np.array(values).T;minimum=np.min(values,axis=1);volume=np.sum(values,axis=1)/12
    assert np.all(minimum>0),'金字塔全参考归一化Jacobian非正'
    center=p.mean(axis=1);integral=np.zeros(len(Q));closure=np.zeros((len(Q),3));area=np.zeros(len(Q));submin=np.inf
    for face in FACES:
        fc=p[:,face].mean(axis=1)-center
        for j in range(len(face)):
            x=p[:,face[j]]-center;y=p[:,face[(j+1)%len(face)]]-center
            av=np.cross(y-x,fc-x)/2;sub=np.einsum('ij,ij->i',x,av)/3
            assert np.all(sub>0),'中位金字塔子体积非正'
            if len(sub):submin=min(submin,float(min(sub)))
            integral+=sub;closure+=av;area+=np.linalg.norm(av,axis=1)
    assert np.all(abs(integral-volume)<=2e-12*volume)
    assert np.all(np.linalg.norm(closure,axis=1)<=2e-12*area)
    return volume,dict(minimumScaledJacobian=float(min(minimum)) if len(Q) else None,
        minimumMedianSubVolume=submin if len(Q) else None,
        maximumVolumeFormulaRelativeError=float(max(abs(integral-volume)/volume)) if len(Q) else 0.)
