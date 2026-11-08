"""模块职责：补齐一次已定位的母单元图连接；保留原失败设计，不改网格/原壁责任。"""
import json,time
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components,breadth_first_order
from eilmerNKSAResourceBudget import ROOT,sha,save,lower_bound,EXPECTED_SIZES,VOLUME_SHA
BASE=ROOT/'outputs/v216_validation'
def main():
    start=time.monotonic();out=BASE/'connected01';out.mkdir(exist_ok=False)
    old=json.loads((BASE/'design03/report.json').read_text());source=BASE/'design03/partition.npz'
    assert sha(source)==old['partitionSHA256'] and not old['allRegionsConnected']
    z=np.load(str(source));a={k:z[k] for k in z.files};n=old['cells'];owner=a['faceOwnerCell'];other=a['faceNeighbourCell'];internal=other>=0
    graph=coo_matrix((np.ones(2*np.sum(internal),np.int8),(np.r_[owner[internal],other[internal]],np.r_[other[internal],owner[internal]])),shape=(n,n)).tocsr()
    i=2;name='afterbody';selected=a['membership'][i].copy();ids=np.flatnonzero(selected)
    number,labels=connected_components(graph[ids][:,ids],directed=False);assert number==2
    count=np.bincount(labels);large=int(np.argmax(count));small=ids[labels!=large];mainmask=np.zeros(n,bool);mainmask[ids[labels==large]]=True
    # 唯一虚拟节点同时连接所有小分量节点，得到真正到主分量的最短无权路径。
    edges=coo_matrix(graph);row=np.r_[edges.row,np.full(len(small),n),small];col=np.r_[edges.col,small,np.full(len(small),n)]
    extended=coo_matrix((np.ones(len(row),np.int8),(row,col)),shape=(n+1,n+1)).tocsr();extended.sort_indices()
    order,pred=breadth_first_order(extended,n,directed=False,return_predecessors=True)
    target=next(int(c) for c in order if c<n and mainmask[c]);path=[target]
    while pred[path[-1]]!=n:path.append(int(pred[path[-1]]))
    path=np.array(path[::-1],np.int32);added=path[~selected[path]];assert len(added)>0
    oldwall=selected[owner[a['originalWallMotherFaceIDs']]];selected[added]=True
    assert np.array_equal(oldwall,selected[owner[a['originalWallMotherFaceIDs']]])
    ids=np.flatnonzero(selected);assert connected_components(graph[ids][:,ids],directed=False,return_labels=False)==1
    # 完整计数重新按原整数面图计算，只更新这个区域；其余区域结果直接复用。
    left=selected[owner];right=np.zeros(len(owner),bool);right[internal]=selected[other[internal]]
    used=left|right;inner=left&right;cut=left^right;interface=cut&internal
    receiver=np.where(left[interface],owner[interface],other[interface]);donor=np.where(left[interface],other[interface],owner[interface])
    mat=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat';assert sha(mat)==VOLUME_SHA
    d=loadmat(str(mat),variable_names=['vertices','wedges','tetrahedra']);w=d['wedges'].astype(np.int64)-1;t=d['tetrahedra'].astype(np.int64)-1
    nodes=np.unique(np.r_[w[ids[ids<len(w)]].reshape(-1),t[ids[ids>=len(w)]-len(w)].reshape(-1)])
    sizes=np.sum(a['motherFaces']>=0,axis=1)
    counts=dict(cells=len(ids),ghosts=int(cut.sum()),faces=int(used.sum()),vertices=len(nodes),internalFaces=int(inner.sum()),
        cellVertexLinks=int(6*np.sum(ids<len(w))+4*np.sum(ids>=len(w))),cellFaceLinks=int(used.sum()+inner.sum()),faceVertexLinks=int(sizes[used].sum()))
    lower=lower_bound(counts,EXPECTED_SIZES);peers,pc=np.unique(a['cellOwnerRegion'][donor],return_counts=True);names=[r['name'] for r in old['regions']]
    r=old['regions'][i];r.update(counts=counts,cellComponents=1,interfaceFaces=int(interface.sum()),
        payloadLowerBoundBytes=lower['payloadLowerBoundBytes'],payloadLowerBoundGiB=lower['payloadLowerBoundGiB'],
        interfaceHexPayloadLowerBoundBytes=int(interface.sum()*140*34),oneV2FileImpossible=int(interface.sum()*140*34)>16*1024**2,
        neighborPatches=[dict(donorRegion=names[j],faces=int(k),hexPayloadLowerBoundBytes=int(k*140*34),compactFileBoundUnder64KiBHeader=int(k*(140*35+256)+65536)) for j,k in zip(peers,pc)])
    a['membership'][i]=selected;a[name+'_cells']=ids.astype(np.int32);a[name+'_interfaceMotherFaces']=np.flatnonzero(interface).astype(np.int32)
    a[name+'_interfaceReceiverCells']=receiver.astype(np.int32);a[name+'_interfaceDonorCells']=donor.astype(np.int32)
    a[name+'_physicalMotherFaces']=np.flatnonzero(cut&~internal).astype(np.int32)
    np.savez_compressed(out/'partition.npz',**a)
    old.update(allRegionsConnected=True,partitionSHA256=sha(out/'partition.npz'),originalDisconnectedDesignAccepted=False,
        connectivityCorrection=dict(region=name,originalComponentSizes=count.tolist(),smallComponentMotherCellIDs=small.tolist(),
            shortestPathMotherCellIDs=path.tolist(),addedMotherCellIDs=added.tolist(),addedCellTypes=a['originalCellTypes'][added].tolist(),
            sourcePartitionSHA256=sha(source),originalWallCoverageUnchanged=True,ownershipUnchanged=True,
            addedCellCentersM=d['vertices'][t[added-len(w)]].mean(axis=1).tolist() if np.all(added>=len(w)) else None),
        allNaturalPatchesWithin16MiBUnderHeaderContract=all(p['compactFileBoundUnder64KiBHeader']<=16*1024**2 for r in old['regions'] for p in r['neighborPatches']),
        designStage='exact-mother-cell-indices-with-one-graph-closure; no physical qualification',seconds=time.monotonic()-start)
    save(out/'report.json',old);print(json.dumps({k:old[k] for k in ['passed','connectivityCorrection','allRegionsConnected','seconds']}))
if __name__=='__main__':main()
