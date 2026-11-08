"""模块职责：完整母域的索引分区、物理边界和准确资源核算；不生成/改写网格。"""
from pathlib import Path
import json,time
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from eilmerNKSAResourceBudget import ROOT,VOLUME_SHA,sha,save,lower_bound,EXPECTED_SIZES
BASE=ROOT/'outputs/v216_validation'

def main():
    start=time.monotonic();out=BASE/'design02';out.mkdir(exist_ok=False)
    source=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat';assert sha(source)==VOLUME_SHA
    d=loadmat(str(source));p=d['vertices'];sv=d['sourceVertices'];f=d['sourceFaces'].astype(np.int64)-1
    w=d['wedges'].astype(np.int64)-1;t=d['tetrahedra'].astype(np.int64)-1;far=d['farfield'].astype(np.int64)-1
    assert len(f)==19752 and len(w)+len(t)==705490 and 'pyramids' not in d
    nw,nt,nv=len(w),len(t),len(sv);n=nw+nt
    cells=np.full((n,6),-1,np.int32);cells[:nw]=w;cells[nw:,:4]=t
    center=np.concatenate([p[w].mean(axis=1),p[t].mean(axis=1)])
    # 模块一：从原单元面建立全域整数图。此处只做组合映射，不重新认证网格质量。
    face_rows=[];owner_rows=[]
    for cell,offset,local in [(w,0,[[0,2,1],[3,4,5],[0,1,4,3],[1,2,5,4],[2,0,3,5]]),
                              (t,nw,[[0,2,1],[0,1,3],[1,2,3],[2,0,3]])]:
        for ids in local:
            face=np.full((len(cell),4),-1,np.int32);face[:,:len(ids)]=cell[:,ids];face_rows.append(face)
            owner_rows.append(np.arange(offset,offset+len(cell),dtype=np.int32))
    half=np.concatenate(face_rows);half_owner=np.concatenate(owner_rows);del face_rows,owner_rows
    keys=np.sort(half,axis=1);unique,first,inv,count=np.unique(keys,axis=0,return_index=True,return_inverse=True,return_counts=True)
    assert np.all((count==1)|(count==2));faces=half[first];owner=half_owner[first]
    second=np.full(len(unique),-1,np.int32);dup=np.arange(len(inv))!=first[inv];second[inv[dup]]=half_owner[dup]
    del half,half_owner,keys,first,inv,dup
    internal=second>=0;face_size=np.sum(faces>=0,axis=1).astype(np.int8)
    # 壁面和远场标签通过精确规范键复用，不凭位置猜测。
    boundary=np.concatenate([f,far]);bk=np.sort(np.column_stack((np.full(len(boundary),-1),boundary)),axis=1).astype(np.int32)
    dtype=np.dtype([('n'+str(i),np.int32) for i in range(4)])
    lookup=np.searchsorted(np.ascontiguousarray(unique).view(dtype).reshape(-1),np.ascontiguousarray(bk).view(dtype).reshape(-1))
    assert np.array_equal(unique[lookup],bk) and np.all(~internal[lookup]) and np.sum(~internal)==len(boundary)
    wall_faces=lookup[:len(f)];outer_faces=lookup[len(f):];assert len(np.unique(lookup))==len(lookup)
    graph=coo_matrix((np.ones(2*np.sum(internal),np.int8),
        (np.r_[owner[internal],second[internal]],np.r_[second[internal],owner[internal]])),shape=(n,n)).tocsr()
    assert connected_components(graph,directed=False,return_labels=False)==1
    # 原规整壁柱按母面编号绑定；完整尾缘扇保留为非规整原单元，不强行解释成32层柱。
    regular=(w<33*nv).all(axis=1)&np.all(np.sort(w[:,:3]%nv,axis=1)==np.sort(w[:,3:]%nv,axis=1),axis=1)
    wallkey=lambda x:np.sort(x,axis=1)@np.array([1,nv,nv*nv],np.int64)
    fk=wallkey(f);order=np.argsort(fk);qkey=wallkey(w[regular,:3]%nv)
    wall_for_wedge=np.full(nw,-1,np.int32);idx=order[np.searchsorted(fk[order],qkey)]
    assert np.array_equal(fk[idx],qkey);wall_for_wedge[regular]=idx
    assert np.all(np.bincount(idx,minlength=len(f))==32)
    old=json.loads((ROOT/'outputs/v212_validation/design01/regions.json').read_text())
    xmin=sv[f,0].min(axis=1);xmax=sv[f,0].max(axis=1);rings=np.unique(sv[:,0]);ring=lambda x:float(rings[np.argmin(abs(rings-x))])
    area=np.linalg.norm(np.cross(sv[f[:,1]]-sv[f[:,0]],sv[f[:,2]]-sv[f[:,0]]),axis=1)/2
    cuts=old['ownershipCutsXM']+[old['extensionEndXM']];end=old['extensionEndXM'];tail=float(sv[:,0].min())
    intervals=int(np.ceil((end-tail)/.5));downcuts=[ring(end+(tail-end)*j/intervals) for j in range(1,intervals)]
    cuts+=downcuts;names=['nose','contraction','afterbody']+['aft_'+str(i+1) for i in range(intervals)]
    owned=np.sum(center[:,0,None]<np.array(cuts),axis=1).astype(np.int8)
    wall_owned=np.sum(sv[f].mean(axis=1)[:,0,None]<np.array(cuts),axis=1).astype(np.int8)
    for i,name in enumerate(names[:3]):wall_owned[np.array(old['regions'][name]['ownedFaceIDs'])-1]=i
    owned[np.flatnonzero(regular)]=wall_owned[wall_for_wedge[regular]]
    # 同一壁面的直接owner也必须归该壁面的责任域（底面/扇采用实际拓扑）。
    owned[owner[wall_faces]]=wall_owned
    membership=np.zeros((len(names),n),bool);records=[];arrays=dict(cellOwnerRegion=owned,wallOwnerRegion=wall_owned,
        motherFaces=faces,faceOwnerCell=owner,faceNeighbourCell=second,originalWallMotherFaceIDs=wall_faces,
        originalFarfieldMotherFaceIDs=outer_faces,originalCellTypes=np.r_[np.full(nw,13,np.int8),np.full(nt,10,np.int8)])
    for i,name in enumerate(names):
        if i<3:
            oldr=old['regions'][name];wallmask=np.zeros(len(f),bool);wallmask[np.array(oldr['originalFaceIDs'])-1]=True
            low=oldr['boundaryRingXM'][0];high=np.inf if i==0 else oldr['boundaryRingXM'][-1]
        else:
            high=ring(cuts[i-1]+.09);low=-np.inf if i==len(names)-1 else ring(cuts[i]-.09)
            wallmask=(xmin>=low-1e-10)&(xmax<=high+1e-10)
        selected=(center[:,0]>=low)&(center[:,0]<=high)
        selected[np.flatnonzero(regular)]=wallmask[wall_for_wedge[regular]]
        selected[owner[wall_faces[wallmask]]]=True;selected|=owned==i
        # 原责任是必须保留量，不能让几何外扩舍入删除其原单元。
        membership[i]=selected;ids=np.flatnonzero(selected)
        adjacent=selected[owner];other=np.zeros(len(faces),bool);other[internal]=selected[second[internal]]
        used=adjacent|other;inner=adjacent&other;cut=adjacent^other
        interface=cut&internal;physical=cut&~internal
        own=np.where(adjacent[interface],owner[interface],second[interface]);donor=np.where(adjacent[interface],second[interface],owner[interface])
        nodes=np.unique(cells[ids]);nodes=nodes[nodes>=0]
        counts=dict(cells=len(ids),ghosts=int(cut.sum()),faces=int(used.sum()),vertices=len(nodes),
            internalFaces=int(inner.sum()),cellVertexLinks=int(np.sum(cells[ids]>=0)),
            cellFaceLinks=int(np.sum(used)+np.sum(inner)),faceVertexLinks=int(face_size[used].sum()))
        bound=lower_bound(counts,EXPECTED_SIZES)
        independent=(13217*(counts['cells']+counts['ghosts'])+3003*counts['faces']+745*counts['vertices']+3170*counts['cells']+
            32*counts['cellVertexLinks']+124*counts['cellFaceLinks']+24*counts['faceVertexLinks']+
            584*(counts['cells']+2*counts['internalFaces'])+344464)
        assert independent==bound['payloadLowerBoundBytes']
        connected=connected_components(graph[ids][:,ids],directed=False,return_labels=False)
        fids=np.flatnonzero(interface);file_min=int(len(fids)*140*34) # 每32位字符串至少2引号，不含逗号/元数据
        expected_wall=np.flatnonzero(selected[owner[wall_faces]])
        assert np.all(np.isin(np.flatnonzero(wall_owned==i),expected_wall))
        for key,value in dict(cells=ids,ownedCells=np.flatnonzero(owned==i),wallFaces=expected_wall+1,
            ownedWallFaces=np.flatnonzero(wall_owned==i)+1,interfaceMotherFaces=fids,
            interfaceReceiverCells=own,interfaceDonorCells=donor,physicalMotherFaces=np.flatnonzero(physical)).items():
            arrays[name+'_'+key]=value.astype(np.int32)
        peers,peer_counts=np.unique(owned[donor],return_counts=True)
        patches=[dict(donorRegion=names[j],faces=int(k),hexPayloadLowerBoundBytes=int(k*140*34),
            compactFileBoundUnder64KiBHeader=int(k*(140*35+256)+65536)) for j,k in zip(peers,peer_counts)]
        assert i not in peers
        record=dict(name=name,axialBoundsM=[None if not np.isfinite(low) else low,None if not np.isfinite(high) else high],
            counts=counts,cellComponents=int(connected),wallFaceCount=len(expected_wall),ownedWallFaceCount=int(np.sum(wall_owned==i)),
            ownedWallAreaM2=float(area[wall_owned==i].sum()),interfaceFaces=len(fids),physicalOuterFaces=int(np.sum(selected[owner[outer_faces]])),
            payloadLowerBoundBytes=independent,payloadLowerBoundGiB=independent/1024**3,
            payloadExcludesExterior=False,exclusions=[s for s in bound['exclusions'] if not s.startswith('all exterior')],
            interfaceHexPayloadLowerBoundBytes=file_min,oneV2FileImpossible=file_min>16*1024**2,neighborPatches=patches,
            fitProven=False,physicalQualification=False,runAuthorized=False)
        records.append(record)
    assert np.all(membership.any(axis=0)) and np.sum(np.bincount(wall_owned))==19752
    arrays['membership']=membership
    assert all(np.array_equal(arrays[n+'_ownedWallFaces'],np.array(old['regions'][n]['ownedFaceIDs'])) for n in names[:3])
    np.savez(out/'partition.npz',**arrays)
    # 模块四：全原外边界的均匀初态传播标签，不能代替发展后的适用性。
    normal=np.cross(p[far[:,1]]-p[far[:,0]],p[far[:,2]]-p[far[:,0]])
    normal/=np.linalg.norm(normal,axis=1)[:,None]
    mach=5*normal@np.array([-np.cos(np.deg2rad(20)),0,np.sin(np.deg2rad(20))])
    outer=dict(count=len(far),normalMachMin=float(mach.min()),normalMachMax=float(mach.max()),
        supersonicInflow=int(np.sum(mach<-1)),supersonicOutflow=int(np.sum(mach>1)),mixed=int(np.sum(abs(mach)<=1)),
        initialOnly=True,physicalBoundaryQualified=False)
    report=dict(passed=True,schema='RAH.CompleteRegionDesign.v1',motherSHA256=VOLUME_SHA,
        originalFaces=19752,originalAreaM2=float(area.sum()),cells=n,faces=len(faces),vertices=len(p),
        regularLayerCells=int(regular.sum()),preservedNonRegularPrisms=int((~regular).sum()),regions=records,
        ownershipCutsXM=cuts,allMotherCellsCovered=True,allOriginalWallsOwnedExactlyOnce=True,
        oldThreeWallResponsibilitiesUnchanged=True,outerBoundary=outer,
        allDomainPayloadsBelow8GiB=all(r['payloadLowerBoundBytes']<8*1024**3 for r in records),
        allSingleV2FilesFeasible=all(not r['oneV2FileImpossible'] for r in records),
        allRegionsConnected=all(r['cellComponents']==1 for r in records),
        completeCandidateAccepted=False,newQualifiedAreaM2=0,CFDRuns=0,meshGenerated=False,
        partitionSHA256=sha(out/'partition.npz'),seconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(report))
if __name__=='__main__':main()
