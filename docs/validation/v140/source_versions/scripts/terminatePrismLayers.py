"""按顶点层数连续终止原棱柱列：输出真实棱柱/金字塔/四面体及原面层身份。"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import breadth_first_order
from collections import deque

def graph_layers(F,nv,seeds,minimum,maximum):
    """多源无权表面边距离，每跨一条原边最多增加一层，不平移切面搜索。"""
    E=np.unique(np.sort(F[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1),axis=0)
    G=coo_matrix((np.ones(2*len(E)),(np.r_[E[:,0],E[:,1]],np.r_[E[:,1],E[:,0]])),shape=(nv,nv)).tocsr()
    distance=np.full(nv,-1,np.int64);distance[seeds]=0;todo=deque(map(int,np.flatnonzero(seeds)))
    assert todo and 1<=minimum<=maximum
    while todo:
        v=todo.popleft()
        for w in G.indices[G.indptr[v]:G.indptr[v+1]]:
            if distance[w]<0:distance[w]=distance[v]+1;todo.append(int(w))
    assert np.all(distance>=0)
    K=minimum+np.minimum(distance,maximum-minimum)
    assert np.all(abs(K[E[:,1]]-K[E[:,0]])<=1)
    return K,distance

def terminate(F,W,K,layers):
    """保留原列节点；终止模板的活动边数决定单元拓扑，不生成悬挂面或新坐标。"""
    nf=len(F);W=np.asarray(W,np.int64);regular=W[:nf*layers].reshape(layers,nf,6)
    columns=np.concatenate([regular[:,:,:3],regular[-1:,:,3:]],axis=0)
    assert np.array_equal(regular[:-1,:,3:],regular[1:,:,:3])
    k=K[F];assert np.all(k>=1) and np.all(k<=layers) and np.all(np.ptp(k,axis=1)<=1)
    level=np.minimum(np.arange(layers+1)[:,None,None],k[None])
    ids=columns[level,np.arange(nf)[None,:,None],np.arange(3)[None,None,:]]
    wedges=[];pyramids=[];tets=[];wi=[];qi=[];ti=[]
    for j in range(layers):
        lo=ids[j];hi=ids[j+1];active=hi!=lo;count=active.sum(axis=1);source=np.arange(nf)
        take=count==3;wedges.append(np.c_[lo[take],hi[take]]);wi.append(np.c_[source[take],np.full(sum(take),j)])
        take=count==2;rows=np.flatnonzero(take);stop=np.argmin(active[take],axis=1);a=(stop+1)%3;b=(stop+2)%3
        pyramids.append(np.c_[lo[rows,a],hi[rows,a],hi[rows,b],lo[rows,b],lo[rows,stop]])
        qi.append(np.c_[rows,np.full(len(rows),j)])
        take=count==1;rows=np.flatnonzero(take);moving=np.argmax(active[take],axis=1)
        tets.append(np.c_[lo[take],hi[rows,moving]]);ti.append(np.c_[rows,np.full(len(rows),j)])
    # 尾缘扇不属于逐原面列，全部按原连接保留，并用-1标明不伪造原面谱系。
    tail=W[nf*layers:];wedges.append(tail);wi.append(np.full((len(tail),2),-1,np.int64))
    result=dict(wedges=np.concatenate(wedges),pyramids=np.concatenate(pyramids),tetrahedra=np.concatenate(tets),
        outerFaces=ids[-1],wedgeSourceFaceLayer=np.concatenate(wi),pyramidSourceFaceLayer=np.concatenate(qi),
        tetSourceFaceLayer=np.concatenate(ti),columnLayers=K)
    assert np.array_equal(result['wedges'][:nf*min(K)],W[:nf*min(K)])
    return result

def mixed_faces(W,Q,T):
    """按标准单元方向列出真实三角/四边面，供独立配对审核。"""
    tri=np.concatenate([W[:,[0,2,1]],W[:,[3,4,5]],T[:,[0,2,1]],T[:,[0,1,3]],T[:,[1,2,3]],T[:,[2,0,3]],
        Q[:,[0,1,4]],Q[:,[1,2,4]],Q[:,[2,3,4]],Q[:,[3,0,4]]])
    quad=np.concatenate([W[:,[0,1,4,3]],W[:,[1,2,5,4]],W[:,[2,0,3,5]],Q[:,[0,3,2,1]]])
    return tri,quad
