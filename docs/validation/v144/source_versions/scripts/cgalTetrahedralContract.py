"""检查CGAL直接导入的必要凸包条件；不把通过此检查当成全几何或物理资格。"""
import numpy as np
from scipy.spatial import ConvexHull
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

def faces_and_neighbors(T):
    faces=T[:,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]]].reshape(-1,3)
    key,inv,count=np.unique(np.sort(faces,axis=1),axis=0,return_inverse=True,return_counts=True)
    assert np.all(count<=2),'nonmanifold face'
    order=np.argsort(inv,kind='stable');start=np.r_[0,np.cumsum(count)[:-1]];pairs=count==2
    a=order[start[pairs]]//4;b=order[start[pairs]+1]//4
    graph=coo_matrix((np.ones(2*len(a)),(np.r_[a,b],np.r_[b,a])),shape=(len(T),len(T))).tocsr()
    return faces[count[inv]==1],graph

def condition(P,T):
    q=P[T];v=np.einsum('ij,ij->i',q[:,1]-q[:,0],np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0]))/6
    assert np.all(v>0) and len(np.unique(np.sort(T,axis=1),axis=0))==len(T)
    boundary,graph=faces_and_neighbors(T);nc,labels=connected_components(graph,directed=False)
    total=float(sum(v));hull=ConvexHull(P[np.unique(T)]);hvol=float(hull.volume)
    gap=abs(hvol-total)/hvol
    return dict(nodes=len(np.unique(T)),tets=len(T),boundaryTriangles=len(boundary),components=int(nc),
        tetVolume=total,convexHullVolume=hvol,relativeConvexHullGap=gap,
        directImportNecessaryConditions=bool(nc==1 and gap<=2e-10),
        scope='necessary only; not a sufficient intersection/topology/quality/physics test'),graph,labels
