"""完整内部边邻域中点细分：保存父单元身份，保持边界与冻结单元，无悬挂面。"""
import numpy as np
from cgalTetrahedralContract import faces_and_neighbors
from mmgLocalInterface import sorted_cells

def volumes(P,T):
    q=P[T];return np.einsum('ij,ij->i',q[:,1]-q[:,0],np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0]))/6

def link_cycle(cells,a,b):
    adj={}
    for cell in cells:
        opposite=[int(v) for v in cell if v not in (a,b)];assert len(opposite)==2
        c,d=opposite;adj.setdefault(c,set()).add(d);adj.setdefault(d,set()).add(c)
    assert all(len(v)==2 for v in adj.values()),'internal edge link must have degree two'
    seen=set();queue=[next(iter(adj))]
    while queue:
        v=queue.pop()
        if v not in seen:seen.add(v);queue.extend(adj[v]-seen)
    assert seen==set(adj),'edge link must be one connected cycle'

def split(P,T,selected,edges):
    P=np.asarray(P,float).copy();T=np.asarray(T,np.int64).copy();selected=selected.copy();parent=np.arange(len(T),dtype=np.int64)
    originalP=P.copy();originalT=T.copy();frozen=T[~selected].copy();before,_=faces_and_neighbors(T)
    boundaryEdges=set(tuple(sorted(map(int,e))) for e in before[:,[[0,1],[1,2],[2,0]]].reshape(-1,2))
    order=sorted(set(tuple(sorted(map(int,e))) for e in edges));assert len(order)==len(edges)
    events=[]
    for a,b in order:
        assert a!=b and 0<=a<b<len(originalP) and (a,b) not in boundaryEdges
        rows=np.flatnonzero(np.any(T==a,axis=1)&np.any(T==b,axis=1));assert len(rows)>=3 and np.all(selected[rows])
        old=T[rows].copy();link_cycle(old,a,b);mid=(P[a]+P[b])/2;node=len(P);P=np.r_[P,mid[None]]
        left=old.copy();right=old.copy();left[left==b]=node;right[right==a]=node;children=np.r_[left,right]
        vo=volumes(P,old);vc=volumes(P,children);assert min(vo)>0 and min(vc)>0
        err=float(max(abs(vc[:len(vo)]+vc[len(vo):]-vo)/vo));assert err<2e-10
        keep=np.ones(len(T),bool);keep[rows]=False;orig=parent[rows]
        T=np.r_[T[keep],children];selected=np.r_[selected[keep],np.ones(len(children),bool)];parent=np.r_[parent[keep],orig,orig]
        events.append(dict(edgeIDsOneBased=[a+1,b+1],newVertexIDOneBased=node+1,splitCells=len(rows),localVolumeRelative=err))
    after,_=faces_and_neighbors(T)
    assert np.array_equal(sorted_cells(before),sorted_cells(after))
    assert np.array_equal(sorted_cells(frozen),sorted_cells(T[~selected])) and np.array_equal(P[:len(originalP)],originalP)
    return P,T,selected,parent,dict(newVertices=len(order),addedTets=len(T)-len(originalT),
        originalTets=len(originalT),finalTets=len(T),fixedTets=len(frozen),boundaryExactlyPreserved=True,
        originalCoordinatesExactlyPreserved=True,allSplitLinksConnectedCycles=True,events=events)
