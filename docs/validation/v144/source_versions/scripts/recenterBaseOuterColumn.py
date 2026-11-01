"""仅重分布底面内部列的受保护层外切向位置；原法向参考坐标和连接保留。"""
import numpy as np

def recenter(P,columnIDs,rimEdgesTop,h,protected):
    h=np.asarray(h,float);columnIDs=np.asarray(columnIDs,dtype=np.int64);E=np.asarray(rimEdgesTop,dtype=np.int64)
    assert len(columnIDs)==len(h) and 0<protected<len(h)-1 and h[0]==0 and np.all(np.diff(h)>0)
    assert not np.any(np.isin(columnIDs,E)),'interior column must not share rim points'
    a=P[E[:,0],1:];b=P[E[:,1],1:];cross=a[:,0]*b[:,1]-a[:,1]*b[:,0]
    area2=float(sum(cross));assert abs(area2)>1e-14
    centroid=np.sum((a+b)*cross[:,None],axis=0)/(3*area2)
    N=(P[columnIDs[-1]]-P[columnIDs[0]])/h[-1];assert abs(np.linalg.norm(N)-1)<1e-10
    assert np.max(abs(P[columnIDs]-(P[columnIDs[0]]+h[:,None]*N)))<2e-13
    target=np.r_[0.,centroid-P[columnIDs[-1],1:]];delta=target-N*np.dot(target,N)
    u=np.clip((h-h[protected])/(h[-1]-h[protected]),0,1);w=u**3*(10+u*(-15+6*u))
    offsets=w[:,None]*delta;Q=P.copy();active=np.flatnonzero(w>0);Q[columnIDs[active]]+=offsets[active]
    bound=float(np.linalg.norm(delta));actual=float(np.max(np.linalg.norm(Q-P,axis=1)))
    assert actual<=bound+2e-13 and np.array_equal(Q[columnIDs[:protected+1]],P[columnIDs[:protected+1]])
    ratios=np.linalg.norm(np.diff(offsets,axis=0),axis=1)/np.diff(h)
    return Q,offsets,dict(projectedPolygonCentroidYZM=centroid.tolist(),originalColumnNormal=N.tolist(),
        targetTangentialDisplacementM=delta.tolist(),maximumDisplacementM=actual,analyticalDisplacementBoundM=bound,
        maximumTangentialStepOverReferenceNormalStep=float(max(ratios)),
        analyticalSlopeBound=1.875*bound/(h[-1]-h[protected]),
        protectedLayers=protected,movedPoints=len(active),columnCoordinateMode='reference-normal-height-plus-tangential-offset',
        normalReferenceProjectionDifferenceM=float(max(abs(offsets@N))),physicsQualified=False,newQualifiedArea=0)
