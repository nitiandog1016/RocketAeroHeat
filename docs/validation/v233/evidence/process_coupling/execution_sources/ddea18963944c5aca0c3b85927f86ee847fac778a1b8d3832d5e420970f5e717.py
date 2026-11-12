"""模块职责：监督一次准确源和真实区域边界检查，输出可复用母网格索引及实施依据。"""
from pathlib import Path
import json,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseEulerRegionalVolume import read_tetra_mesh,topology,extract_region,interface_characteristics
BASE=ROOT/'outputs/v233_validation';DOC=ROOT/'docs/validation/v233'
OLD=ROOT/'tools/runtime/su2-thermal-build/source';PREV=ROOT/'outputs/v232_validation'
MESH=PREV/'rocket01/mesh.su2';FRAME=PREV/'rocket01/state_000800.csv'
REGION=ROOT/'outputs/v228_validation/prepare03/region.json'
GEOM=ROOT/'outputs/v230_validation/audit02/geometry.npz'
SPECS=[
 (ROOT/'outputs/v231_validation/recovery_prepare04/CEulerVariable.cpp','bool CEulerVariable::SetPrimVar',72),
 (ROOT/'outputs/v229_validation/source/SU2_CFD/src/solvers/CEulerSolver.cpp','unsigned long CEulerSolver::SetPrimitive_Variables',42),
 (OLD/'SU2_CFD/include/solvers/CFVMFlowSolverBase.hpp','inline void BC_Euler_Wall',10),
 (OLD/'SU2_CFD/include/solvers/CFVMFlowSolverBase.inl','::BC_Sym_Plane(',200),
 (OLD/'SU2_CFD/include/solvers/CFVMFlowSolverBase.inl','::CompleteImplicitIteration(',57),
 (OLD/'SU2_CFD/src/integration/CMultiGridIntegration.cpp','void CMultiGridIntegration::PreSmoothing',125),
 (OLD/'SU2_CFD/src/integration/CIntegration.cpp','void CIntegration::Space_Integration',255)]


def sources():
    """模块一：只读绑定旧对象实际源与新修复源，保存带行号的必要调用段。"""
    out=BASE/'source01';out.mkdir(exist_ok=False)
    known={}
    for p in [ROOT/'outputs/v167_validation/build01/source_tree.json',ROOT/'docs/validation/v231/evidence/actual_dependencies.json',
      ROOT/'outputs/v229_validation/process_build/report.json',ROOT/'outputs/v231_validation/recovery_prepare04/object_dependencies.json']:
        d=read(p);known.update(d.get('inputSHA256',d))
    report=[]
    for k,(p,needle,count) in enumerate(SPECS):
        h=sha(p);assert known[str(p)]==h,('SOURCE_NOT_BOUND',str(p))
        lines=p.read_text(encoding='utf8').splitlines();matches=[i for i,s in enumerate(lines) if needle in s]
        assert len(matches)==1,(needle,matches)
        start=max(0,matches[0]-2);end=min(len(lines),matches[0]+count)
        excerpt='\n'.join('%6d %s'%(i+1,lines[i]) for i in range(start,end))+'\n'
        target=out/('%02d_%s.txt'%(k,p.name));target.write_text(excerpt,encoding='utf8')
        report.append(dict(path=str(p),sha256=h,firstLine=start+1,lastLine=end,anchor=needle,excerptSHA256=sha(target)))
    save(out/'report.json',dict(passed=True,sources=report,uniqueSources=len(set(str(s[0]) for s in SPECS)),
      flowUpdates=0,nativeQueries=0,builds=0,observationImplemented=False))
    print(json.dumps(dict(passed=True,sourceExcerpts=len(report))))


def regions():
    """模块二：原切口的前/后互补控制体，精确原面归属及双向特征需求。"""
    out=BASE/'region01';out.mkdir(exist_ok=False)
    previous=read(PREV/'final01/report.json');assert previous['passed'] and not previous['scientificAcceptance']
    assert sha(MESH)==previous['assets']['outputs/v194_validation/rocket01/mesh.su2']
    frame_info=next(r for r in read(PREV/'audit01/report.json')['savedFrameChecks'] if r['update']==800)
    assert sha(FRAME)==frame_info['sourceSHA256']
    xyz,cells,markers=read_tetra_mesh(MESH);assert cells.shape==(64982,4) and xyz.shape==(18027,3)
    topo=topology(xyz,cells,markers);region=read(REGION);g=np.load(str(GEOM))
    original_nodes=g['volumeNodeIDs'][g['faces']];area=g['area'];assert len(area)==19752
    mapping={tuple(sorted(row)):int(i) for row,i in zip(markers['wall'],topo['markers']['wall'])}
    wall=np.array([mapping[tuple(sorted(f))] for f in original_nodes]);assert len(np.unique(wall))==19752
    fore_ids=np.array(region['originalFaceIDs'],int)-1;nose_ids=np.array(region['retainedNoseCandidateIDs'],int)-1
    front_wall=np.zeros(19752,bool);front_wall[np.r_[fore_ids,nose_ids]]=True
    assert len(fore_ids)==3128 and len(nose_ids)==1306 and sum(front_wall)==4434
    cut=float(region['originalOwnershipCutXM']);selected=topo['center'][:,0]>=cut
    # 壁面owner必须完整保留原责任；若同一个母单元跨冲突责任，则拒绝，不能偷偷删除原面。
    wall_owner=topo['owner'][wall]
    if np.intersect1d(wall_owner[front_wall],wall_owner[~front_wall]).size:raise ValueError('WALL_OWNER_RESPONSIBILITY_CONFLICT')
    initial=selected.copy();selected[wall_owner[front_wall]]=True;selected[wall_owner[~front_wall]]=False
    front=extract_region(topo,selected);aft=extract_region(topo,~selected)
    assert np.array_equal(np.sort(np.r_[front['cells'],aft['cells']]),np.arange(len(cells)))
    assert np.array_equal(front['interface'],aft['interface'])
    assert np.array_equal(selected[wall_owner],front_wall)
    fids=front['interface'];sign=np.where(selected[topo['owner'][fids]],1.,-1.);normal=topo['normal'][fids]*sign[:,None]
    with FRAME.open(encoding='utf8') as f:columns=f.readline().strip().split(',')
    assert columns==['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy','Volume']+['RHS%d'%i for i in range(5)]
    raw=np.loadtxt(str(FRAME),delimiter=',',skiprows=1);raw=raw[np.argsort(raw[:,0])]
    assert np.array_equal(raw[:,0],np.arange(len(xyz))) and np.array_equal(raw[:,1:4],xyz)
    result=interface_characteristics(raw[:,4:9],topo['faces'][fids],normal)
    reverse=interface_characteristics(raw[:,4:9],topo['faces'][fids],-normal)
    assert not np.any(result['sonic']) and np.array_equal(result['incoming']+reverse['incoming'],np.full(len(fids),5))
    assert np.array_equal(result['p'],reverse['p']) and np.array_equal(result['T'],reverse['T'])
    receivers=np.where(sign>0,topo['owner'][fids],topo['other'][fids]);donors=np.where(sign>0,topo['other'][fids],topo['owner'][fids])
    columns=[fids,receivers,donors]+[topo['faces'][fids,j] for j in range(3)]+[normal[:,j] for j in range(3)]
    columns += [result[k] for k in ('p','T','H','S','normalMach','incoming')]+[reverse['incoming']]
    np.savetxt(str(out/'interface.csv'),np.column_stack(columns),delimiter=',',fmt='%.17g',comments='',
      header='motherFaceID,frontCellID,aftCellID,node0,node1,node2,outwardAx,outwardAy,outwardAz,pPa,TK,HJkg,SJkgK,normalMach,incomingFront,incomingAft')
    np.savetxt(str(out/'wall_responsibility.csv'),np.column_stack([np.arange(1,19753),area,front_wall.astype(int),wall_owner]),
      delimiter=',',fmt='%.17g',comments='',header='originalFaceID,areaM2,frontRegion,motherWallOwnerCellID')
    np.savez_compressed(str(out/'regional_volume_index.npz'),frontCells=front['cells'],aftCells=aft['cells'],
      interfaceFaces=fids,interfaceNodes=topo['faces'][fids],interfaceNormalOutOfFront=normal,frontCell=receivers,aftCell=donors,
      originalWallMotherFaces=wall,originalWallOwnerCells=wall_owner,frontOriginalFaces=np.flatnonzero(front_wall)+1,
      aftOriginalFaces=np.flatnonzero(~front_wall)+1,frontBoundary=front['boundary'],aftBoundary=aft['boundary'])
    h0=read(PREV/'regional01/report.json')['referenceFreestreamTotalEnthalpyJkg']
    report=dict(passed=True,method='original tet control-volume split / fixed v228 ownership cut',cutXM=cut,
      originalFaces=19752,originalAreaM2=float(sum(area)),wallForcedCellChanges=int(sum(initial!=selected)),
      front=front['summary'],aft=aft['summary'],frontWallFaces=int(sum(front_wall)),aftWallFaces=int(sum(~front_wall)),
      frontWallAreaM2=float(sum(area[front_wall])),aftWallAreaM2=float(sum(area[~front_wall])),
      targetEngineeringFaces=len(fore_ids),targetEngineeringAreaM2=float(sum(area[fore_ids])),
      exactPartition=True,originalWallOwnershipExact=True,frontInterface=result['summary'],aftInterface=reverse['summary'],
      complementCharacteristicCountsExact=True,oneWaySupersonicFrontOutletPossible=bool(np.all(result['incoming']==0)),
      interfaceTotalHRelativeToFreestreamRange=[float(min((result['H']-h0)/h0)),float(max((result['H']-h0)/h0))],
      nativeSameStageMomentumObserved=False,actualRejectedIDsKnown=False,inletProfilesAvailable=False,
      spatialFarfieldQualified=False,nativeRegionalFluxBalance=False,physicalQualification=False,newQualifiedAreaM2=0,
      sourceSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in (MESH,FRAME,REGION,GEOM)},
      sourceStateUpdate=800,physicalTime=None,flowUpdates=0,newMeshGenerated=False,
      fullVolumeIndexSHA256=sha(out/'regional_volume_index.npz'),interfaceSHA256=sha(out/'interface.csv'))
    save(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('sourceSHA256',)}))


def main(stage):
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original_size=guard.size_bytes;guard.size_bytes=lambda _:original_size(BASE)
    inputs=[Path(__file__).resolve(),ROOT/'scripts/phaseEulerRegionalVolume.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
    if stage=='sources':
        inputs += [p for p,_,_ in SPECS]+[ROOT/'outputs/v167_validation/build01/source_tree.json',
          ROOT/'docs/validation/v231/evidence/actual_dependencies.json',ROOT/'outputs/v229_validation/process_build/report.json',
          ROOT/'outputs/v231_validation/recovery_prepare04/object_dependencies.json']
    elif stage=='regions':
        inputs += [MESH,FRAME,REGION,GEOM,PREV/'final01/report.json',PREV/'audit01/report.json',PREV/'regional01/report.json',
          ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',ROOT/'data/thermo/phase_aware_dry_air_v1.json']
    elif stage=='coupling':
        inputs += [ROOT/'outputs/v167_validation/build01/source_tree.json',ROOT/'outputs/v216_validation/connected01/partition.npz',
          ROOT/'outputs/v216_validation/final01/report.json',GEOM,REGION]
        inputs += [OLD/p for p in ('SU2_CFD/src/drivers/CDriver.cpp','SU2_CFD/src/drivers/CMultizoneDriver.cpp',
          'SU2_CFD/src/interfaces/cfd/CSlidingInterface.cpp','SU2_CFD/include/solvers/CFVMFlowSolverBase.inl')]
    else:raise ValueError('STAGE_NOT_ADMITTED')
    return guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_'+stage],120,128,list(dict.fromkeys(inputs)),BASE)


def coupling():
    """模块三：完整黏性母域双区设计与原生传递来源；不把软件接口当守恒或收敛证明。"""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    out=BASE/'coupling01';out.mkdir(exist_ok=False);known=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    chunks=[('SU2_CFD/src/drivers/CDriver.cpp','else if (fluid_donor && fluid_target)',53),
      ('SU2_CFD/src/drivers/CMultizoneDriver.cpp','void CMultizoneDriver::RunGaussSeidel',70),
      ('SU2_CFD/src/drivers/CMultizoneDriver.cpp','bool CMultizoneDriver::TransferData',52),
      ('SU2_CFD/src/interfaces/cfd/CSlidingInterface.cpp','void CSlidingInterface::GetDonor_Variable',61),
      ('SU2_CFD/include/solvers/CFVMFlowSolverBase.inl','::BC_Fluid_Interface(',180)]
    sources=[]
    for k,(rel,needle,count) in enumerate(chunks):
        p=OLD/rel;h=sha(p);assert known[str(p)]==h
        lines=p.read_text(encoding='utf8').splitlines();idx=[i for i,s in enumerate(lines) if needle in s];assert len(idx)==1
        start=idx[0];end=min(len(lines),start+count);target=out/('%02d_%s.txt'%(k,p.name))
        target.write_text('\n'.join('%6d %s'%(i+1,lines[i]) for i in range(start,end))+'\n',encoding='utf8')
        sources.append(dict(path=str(p),sha256=h,firstLine=start+1,lastLine=end,excerptSHA256=sha(target)))
    p=ROOT/'outputs/v216_validation/connected01/partition.npz';previous=read(ROOT/'outputs/v216_validation/final01/report.json')
    assert sha(p)==previous['sourceSHA256']['outputs\\v216_validation\\connected01\\partition.npz']
    z=np.load(str(p));owned=z['cellOwnerRegion'];wall_owned=z['wallOwnerRegion'];front=owned==0
    owner=z['faceOwnerCell'];other=z['faceNeighbourCell'];internal=other>=0;n=len(owned);assert n==705490
    graph=coo_matrix((np.ones(2*sum(internal)),(np.r_[owner[internal],other[internal]],np.r_[other[internal],owner[internal]])),shape=(n,n)).tocsr()
    split=np.zeros(len(owner),bool);split[internal]=front[owner[internal]]!=front[other[internal]];interface=np.flatnonzero(split)
    wall=z['originalWallMotherFaceIDs'];assert np.array_equal(front[owner[wall]],wall_owned==0)
    g=np.load(str(GEOM));area=g['area'];old_region=read(REGION)
    expected=np.sort(np.array(old_region['originalFaceIDs']+old_region['retainedNoseCandidateIDs']))
    assert np.array_equal(np.flatnonzero(wall_owned==0)+1,expected)
    info=[]
    for selected,label in ((front,'upstream-viscous-provider'),(~front,'aft-viscous-feedback')):
        ids=np.flatnonzero(selected);mask=(wall_owned==0) if label.startswith('upstream') else (wall_owned!=0)
        info.append(dict(name=label,cells=len(ids),components=int(connected_components(graph[ids][:,ids],directed=False,return_labels=False)),
          wallFaces=int(sum(mask)),wallAreaM2=float(sum(area[mask])),
          prismCells=int(sum(z['originalCellTypes'][ids]==13)),tetCells=int(sum(z['originalCellTypes'][ids]==10))))
    receiver=np.where(front[owner[interface]],owner[interface],other[interface]);donor=np.where(front[owner[interface]],other[interface],owner[interface])
    nodes=z['motherFaces'][interface];assert nodes.shape[1]==4
    np.savez_compressed(str(out/'viscous_two_zone_index.npz'),upstreamCells=np.flatnonzero(front),aftCells=np.flatnonzero(~front),
      interfaceMotherFaces=interface,interfaceMotherNodes=nodes,upstreamInterfaceCells=receiver,aftInterfaceCells=donor,
      upstreamOriginalFaces=np.flatnonzero(wall_owned==0)+1,aftOriginalFaces=np.flatnonzero(wall_owned!=0)+1)
    save(out/'report.json',dict(passed=True,sourceSHA256=sources,parentPartitionSHA256=sha(p),regions=info,interfaceFaces=len(interface),
      originalCells=n,originalFaces=19752,originalAreaM2=float(sum(area)),allOriginalCellsAndWallsPreserved=True,
      exactOriginalResponsibility=True,allRegionsConnected=all(r['components']==1 for r in info),
      interfaceNodeIDsExact=True,old32LayerAssetsUnchanged=True,nativeEulerToRANSTurbulenceDonorAbsent=True,
      directEulerRANSReady=False,chosenCandidate='same-phase-Air RANS/SA-neg two zones / mature Gauss-Seidel + sliding fluid interface',
      artificialInternalFreestreamBoundary=False,requiresInventedInletProfile=False,nativeMultiZoneTested=False,
      interfaceViscousGradientUsesReceiverBothSides=True,pairFluxConservationProven=False,wallDistanceScopeVerified=False,
      fullPeakBytes=None,oldSingleDomainMeasuredRSSBytes=2153095168,
      newMeshGenerated=False,nativeQueries=0,builds=0,flowUpdates=0,newQualifiedAreaM2=0,physicalQualification=False))
    print(json.dumps(dict(passed=True,regions=info,interfaceFaces=len(interface),directEulerRANSReady=False,nativeMultiZoneTested=False)))


if __name__=='__main__':
    if sys.argv[1]=='task_sources':sources()
    elif sys.argv[1]=='task_regions':regions()
    elif sys.argv[1]=='task_coupling':coupling()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
