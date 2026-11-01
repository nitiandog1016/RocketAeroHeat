"""保留原棱柱层连接和受保护层，只沿原法向列构造连续外包络。"""
import numpy as np

def height_table(V,h,protected,join,cut):
    h=np.asarray(h,float);assert 1<=protected<len(h)-1 and join<cut and np.all(np.diff(h)>0) and h[0]==0
    u=np.clip((V[:,0]-join)/(cut-join),0,1);weight=u**3*(10+u*(-15+6*u))
    cap=h[protected]+np.arange(1,len(h)-protected)*(h[protected]-h[protected-1])
    target=np.tile(h[:,None],(1,len(V)))
    target[protected+1:]=(1-weight[None])*h[protected+1:,None]+weight[None]*cap[:,None]
    assert np.all(np.diff(target,axis=0)>0) and np.all(cap<=h[protected+1:]+1e-14)
    return target,weight,cap

def redistribute(P,V,h,protected,join,cut):
    """输入须具有层序编号；只修改非零权重列，原尾缘扇及外域节点均保持原字节。"""
    target,weight,cap=height_table(V,h,protected,join,cut);h=np.asarray(h,float);nv=len(V);n=len(h)-1
    assert len(P)>=(n+1)*nv and np.array_equal(P[:nv],V)
    selected=np.flatnonzero(weight>0);N=(P[n*nv+selected]-V[selected])/h[-1]
    assert np.all(abs(np.linalg.norm(N,axis=1)-1)<1e-10),'modified columns must be original unit-normal columns'
    Q=P.copy()
    for j in range(protected+1,n+1):
        ids=j*nv+selected
        assert np.max(abs(P[ids]-(V[selected]+h[j]*N)))<2e-13
        Q[ids]=V[selected]+target[j,selected,None]*N
    assert np.array_equal(Q[:(protected+1)*nv],P[:(protected+1)*nv])
    moved=np.flatnonzero(np.any(Q!=P,axis=1));delta=np.linalg.norm(Q[moved]-P[moved],axis=1)
    bound=float(max(h[protected+1:]-cap));assert max(delta)<=bound+2e-13
    return Q,target,dict(protectedLayers=protected,normalLayers=n,joinSourceXM=float(join),cutSourceXM=float(cut),
        modifiedSourceColumns=len(selected),movedOriginalNodes=len(moved),maximumDisplacementM=float(max(delta)),
        analyticalDisplacementBoundM=bound,capOuterHeightM=float(cap[-1]),capSpacingM=float(h[protected]-h[protected-1]),
        minimumNewSpacingM=float(np.diff(target,axis=0).min()),originalProtectedPointsExactlyPreserved=True,
        newNodes=0,weightRule='quintic smoothstep of original wall X',physicsQualified=False,newQualifiedArea=0)
