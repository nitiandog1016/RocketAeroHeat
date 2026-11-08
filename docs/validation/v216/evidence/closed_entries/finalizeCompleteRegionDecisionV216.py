"""模块职责：独立审核七域整数责任/连接/资源，画真实箭体分区，并形成明确的入口决定。"""
import json,time,csv,hashlib
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components,shortest_path
from eilmerNKSAResourceBudget import ROOT,VOLUME_SHA,sha,save
BASE=ROOT/'outputs/v216_validation'
def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False)
    r=json.loads((BASE/'connected01/report.json').read_text());path=BASE/'connected01/partition.npz'
    assert sha(path)==r['partitionSHA256'];z=np.load(str(path));owner=z['faceOwnerCell'];other=z['faceNeighbourCell'];internal=other>=0
    membership=z['membership'];cell_owner=z['cellOwnerRegion'];wall_owner=z['wallOwnerRegion'];wall=z['originalWallMotherFaceIDs']
    assert membership.shape==(7,705490) and membership.any(axis=0).all() and len(wall)==19752
    assert np.all(membership[cell_owner,np.arange(len(cell_owner))])
    assert np.all(membership[wall_owner,owner[wall]])
    source=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat';assert sha(source)==VOLUME_SHA
    d=loadmat(str(source));v=d['sourceVertices'];f=d['sourceFaces'].astype(np.int64)-1;w=d['wedges'].astype(np.int64)-1;t=d['tetrahedra'].astype(np.int64)-1
    area=np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)/2
    assert len(w)==634916 and len(t)==70574 and abs(float(area.sum())-1.9646292565715553)<1e-14
    graph=coo_matrix((np.ones(2*int(internal.sum()),np.int8),(np.r_[owner[internal],other[internal]],np.r_[other[internal],owner[internal]])),shape=(705490,705490)).tocsr()
    previous=json.loads((ROOT/'outputs/v212_validation/design01/regions.json').read_text());checks=[]
    for i,region in enumerate(r['regions']):
        mask=membership[i];ids=np.flatnonzero(mask);name=region['name'];assert np.array_equal(z[name+'_cells'],ids)
        assert connected_components(graph[ids][:,ids],directed=False,return_labels=False)==1
        a=mask[owner];b=np.zeros(len(owner),bool);b[internal]=mask[other[internal]];present=a|b;both=a&b
        nodes=np.unique(np.r_[w[ids[ids<len(w)]].reshape(-1),t[ids[ids>=len(w)]-len(w)].reshape(-1)])
        count=dict(cells=len(ids),ghosts=int(np.count_nonzero(a!=b)),faces=int(present.sum()),vertices=len(nodes),
            internalFaces=int(both.sum()),cellVertexLinks=int(6*np.sum(ids<len(w))+4*np.sum(ids>=len(w))),
            cellFaceLinks=int(5*np.sum(ids<len(w))+4*np.sum(ids>=len(w))),faceVertexLinks=int(np.sum(z['motherFaces'][present]>=0)))
        assert count==region['counts']
        value=13217*(count['cells']+count['ghosts'])+3003*count['faces']+745*count['vertices']+3170*count['cells']+32*count['cellVertexLinks']+124*count['cellFaceLinks']+24*count['faceVertexLinks']+584*(count['cells']+2*count['internalFaces'])+344464
        assert value==region['payloadLowerBoundBytes']
        owned=np.flatnonzero(wall_owner==i)+1;assert np.array_equal(owned,z[name+'_ownedWallFaces'])
        if i<3:assert np.array_equal(owned,np.asarray(previous['regions'][name]['ownedFaceIDs']))
        assert abs(float(area[owned-1].sum())-region['ownedWallAreaM2'])<1e-14
        ifids=z[name+'_interfaceMotherFaces'];assert np.array_equal(ifids,np.flatnonzero((a^b)&internal))
        donors=np.where(a[ifids],other[ifids],owner[ifids]);assert np.array_equal(donors,z[name+'_interfaceDonorCells'])
        assert np.all(~mask[donors]) and np.all(cell_owner[donors]!=i)
        checks.append(dict(name=name,counts=count,payloadLowerBoundBytes=value))
    corr=r['connectivityCorrection'];chain=np.array(corr['shortestPathMotherCellIDs']);assert len(chain)==3
    assert all(graph[int(a),int(b)]==1 for a,b in zip(chain[:-1],chain[1:]))
    old=np.load(str(BASE/'design03/partition.npz'));changed=np.flatnonzero(membership[2]!=old['membership'][2])
    assert changed.tolist()==corr['addedMotherCellIDs']==[675194]
    assert np.array_equal(cell_owner,old['cellOwnerRegion']) and np.array_equal(wall_owner,old['wallOwnerRegion'])
    assert np.array_equal(membership[[0,1,3,4,5,6]],old['membership'][[0,1,3,4,5,6]])
    assert not old['membership'][2,675194] and not old['membership'][2,635870-0:635870+1].size==0
    distance=shortest_path(graph,directed=False,indices=635870,unweighted=True)
    large=old['membership'][2].copy();large[635870]=False
    assert distance[np.flatnonzero(large)].min()==2
    num=json.loads((BASE/'numerical02/report.json').read_text());assert not num['fullNumericalAcceptance']
    assert all(x['assembly']['firstFiveSavedStateCoveredByGamma7'] for x in num['records'])
    assert all(x['assembly']['firstFiveFailedAfter70DigitAssembly']>0 for x in num['records'])
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text());assert len(protected)==29
    assert all(sha(ROOT/k.replace('\\','/'))==h for k,h in protected.items())
    # 模块二：物理外边界的下一选择。此处仅连续几何设计，没有生成表面/体网格。
    inv=json.loads((BASE/'inventory.json').read_text());lo,hi=np.array(inv['domainBounds']);length=float(v[:,0].max()-v[:,0].min())
    front=float(hi[0]+.1*length);back=float(lo[0]);rf=float(1.05*np.linalg.norm(np.maximum(abs(lo[1:]),abs(hi[1:]))))
    half=np.deg2rad(40);rb=float(rf+(front-back)*np.tan(half))
    cone=dict(status='continuous-geometry-design-only',purpose='enclose existing full cube; keep all original wall/layer/fan/cells',
        frontXM=front,backXM=back,frontRadiusM=rf,backRadiusM=rb,halfAngleDeg=40,alphaDeg=20,Mach=5,
        sideInitialNormalMachBounds=[float(-5*np.sin(np.deg2rad(60))),float(-5*np.sin(np.deg2rad(20)))],
        frontInitialNormalMach=float(-5*np.cos(np.deg2rad(20))),backInitialNormalMach=float(5*np.cos(np.deg2rad(20))),
        meshGenerated=False,developmentBoundaryQualified=False,
        admissionConditions=['actual discrete normals/closed outer shell/preserve all old cells',
            'supersonic upstream and side inflow stays outside disturbed region',
            'downstream backflow/subsonic state explicitly rejects simple supersonic extrapolation',
            'complete extra-cell/IO/allocator peak budget and independent outer-boundary sensitivity'])
    decision=dict(schema='RAH.CompleteRegionDecision.v1',geometryPartitionAccepted=True,resourcePeakAccepted=False,
        filePatchDesignAcceptedUnderHeaderContract=True,originalCubeBoundaryAccepted=False,normalMixedOuterFaces=2604,
        originalFiniteDifferenceAccepted=False,productionSchedulerAccepted=False,CFDAdmitted=False,
        preferredBoundaryDesign=cone,nextNumericalEvidence='capture exact real primitive face/cell stencil inputs and native weights at the already located worst cases; independent high-precision local evaluation; no repeated three-group Jv/FD sweep',
        noNewQualifiedArea=True,goalComplete=False)
    save(out/'decision.json',decision)
    # 模块三：科研图以原19752三角形责任和实际下界为主体，不伪装成压力或热流云图。
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    colors=['#446e91','#d28b49','#689679','#9b7aa0','#647f98','#ae9b60','#a87373']
    labels=['Nose','Contraction','Short aft','Aft 1','Aft 2','Aft 3','Base / wake']
    fig=plt.figure(figsize=(12,5.5),facecolor='white');ax=fig.add_subplot(121,projection='3d')
    collection=Poly3DCollection(v[f],facecolors=np.array(colors)[wall_owner],edgecolors='none',linewidths=0,rasterized=True)
    ax.add_collection3d(collection);ax.set_xlim(v[:,0].min(),v[:,0].max());ax.set_ylim(-.17,.17);ax.set_zlim(-.17,.17)
    ax.set_box_aspect((5.7,1,1));ax.view_init(elev=16,azim=-68);ax.set_xlabel('X (m)');ax.set_ylabel('Y (m)');ax.set_zlabel('Z (m)')
    ax.set_title('Original rocket: unique wall responsibility\n19,752 faces / 1.9646292566 m2',fontsize=11,pad=18)
    bx=fig.add_subplot(122);values=[x['payloadLowerBoundGiB'] for x in r['regions']]
    bars=bx.barh(labels,values,color=colors,height=.65);bx.invert_yaxis();bx.set_xlim(0,8.6)
    bx.axvline(8,color='#555555',linestyle='--',linewidth=1);bx.text(8.08,6.6,'WSL 8 GiB',rotation=90,ha='left',va='bottom',fontsize=9)
    for bar,val in zip(bars,values):bx.text(val+.07,bar.get_y()+bar.get_height()/2,'%.3f'%val,va='center',fontsize=10)
    bx.set_xlabel('Necessary solver payload (GiB)');bx.set_title('Full original domain included\nPeak memory / physical validity not established',fontsize=11)
    bx.spines['top'].set_visible(False);bx.spines['right'].set_visible(False);bx.grid(axis='x',alpha=.15);bx.set_axisbelow(True)
    fig.text(.5,.025,'Design and resource accounting only. No new solved pressure or heat flux.',ha='center',fontsize=10,color='#444444')
    fig.subplots_adjust(left=.01,right=.95,bottom=.15,top=.85,wspace=.35);fig.savefig(out/'regions.png',dpi=150);plt.close(fig)
    report=dict(passed=True,regionChecks=checks,allOriginalCellsAndWallsPreserved=True,protectedFilesUnchanged=29,
        oneShortestPathCellAdded=675194,allRegionsConnected=True,originalDisconnectedDesignRemainsRejected=True,
        numericalWholeOperatorAccepted=False,sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in [path,source,BASE/'numerical02/report.json',Path(__file__).resolve()]},
        decisionSHA256=sha(out/'decision.json'),seconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(dict(passed=True,protected=29,seconds=report['seconds'])))
if __name__=='__main__':main()
