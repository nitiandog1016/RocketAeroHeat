"""模块职责：导出真实前体黏性区域及双侧接口证据，拒绝缺物理供体的启动。"""
from pathlib import Path
import sys
import numpy as np
from scipy.io import loadmat
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from singleDomainFrame import read_frame,COLUMNS
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy

BASE=ROOT/'outputs/v244_validation';DOC=ROOT/'docs/validation/v244';CLOSED=False
VOLUME=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat'
PART=ROOT/'outputs/v216_validation/connected01/partition.npz'
INDEX=ROOT/'outputs/v233_validation/coupling01/viscous_two_zone_index.npz'
FRAME=ROOT/'outputs/v236_validation/rocket01/state_000300.bin'
SOURCES=ROOT/'outputs/v235_validation/source01'


def wall_contact_stats(interface,wall,table):
    """壁面接触按拓扑识别；浮点壁距的严格零值另列，不能代替接触数。"""
    nodes=np.unique(wall[wall>=0]);contact=np.isin(interface,nodes).any(axis=1)
    return dict(interfaceFacesTouchingWall=int(sum(contact)),
        interfaceFacesWithExactlyZeroStoredWallDistance=int(sum(table[:,17]==0)))


def sources(out):
    """模块一：绑定真实构建来源；源码存在不等于当前配置允许该边界。"""
    specs=[('SU2_CFD/src/solvers/CEulerSolver.cpp',[
        ('case STATIC_PRESSURE:',23,5290),('/*--- Compute the characteristic jumps ---*/',24,5290),
        ('/*--- Viscous contribution ---*/',79,5630)]),
        ('SU2_CFD/src/solvers/CTurbSASolver.cpp',[('void CTurbSASolver::BC_Outlet',27,0)]),
        ('SU2_CFD/src/integration/CIntegration.cpp',[('case RIEMANN_BOUNDARY:',13,0)])]
    bound={r['relative']:r for r in read(SOURCES/'bindings.json')};rows=[]
    for rel,anchors in specs:
        p=SOURCES/rel;assert sha(p)==bound[rel]['sha256'];lines=p.read_text(encoding='utf8').splitlines()
        for k,(anchor,count,offset) in enumerate(anchors):
            first=next(i for i,s in enumerate(lines) if i>=offset and anchor in s);last=min(len(lines),first+count)
            dest=out/(p.stem+'_%d.txt'%k)
            dest.write_text('\n'.join('%6d %s'%(i+1,lines[i]) for i in range(first,last))+'\n',encoding='utf8')
            rows.append(dict(path=str(p),sha256=sha(p),firstLine=first+1,lastLine=last,excerpt=dest.name))
    p=ROOT/'outputs/v220_validation/source/Common/src/CConfig.cpp'
    known=read(ROOT/'outputs/v220_validation/build02/dependencies.json');assert sha(p)==known[str(p)]
    lines=p.read_text(encoding='utf8').splitlines();first=next(i for i,s in enumerate(lines) if 'if (Kind_FluidModel == RAH_PHASE_AIR) {' in s)
    (out/'phase_configuration.txt').write_text('\n'.join('%6d %s'%(i+1,lines[i]) for i in range(first,first+13))+'\n',encoding='utf8')
    rows.append(dict(path=str(p),sha256=sha(p),firstLine=first+1,lastLine=first+13,excerpt='phase_configuration.txt'))
    return rows


def prepare():
    """模块二：从原母单元精确导出区域，不重新网格化、不改变近壁结构。"""
    out=BASE/'region02';out.mkdir(exist_ok=False)
    known=read(ROOT/'docs/validation/v243/evidence/final_report.json')['assets']
    assert sha(VOLUME)==known['outputs/v144_validation/actual/hybrid_volume.mat']
    assert sha(PART)==read(ROOT/'outputs/v233_validation/coupling01/report.json')['parentPartitionSHA256']
    assert sha(INDEX)==read(ROOT/'docs/validation/v233/evidence/local_asset_index.json')['files']['coupling01/viscous_two_zone_index.npz']['sha256']
    receipt=read(ROOT/'outputs/v236_validation/process_audit3/report.json')
    assert receipt['passed'] and sha(FRAME)==receipt['inputSHA256'][str(FRAME)]
    d=loadmat(str(VOLUME),variable_names=['vertices','wedges','tetrahedra','sourceFaces'])
    xyz=d['vertices'];w=d['wedges'].astype(np.int64)-1;t=d['tetrahedra'].astype(np.int64)-1
    cells=np.full((len(w)+len(t),6),-1,dtype=np.int64);cells[:len(w)]=w;cells[len(w):,:4]=t
    z=np.load(str(PART));idx=np.load(str(INDEX));selected=idx['upstreamCells'].astype(np.int64)
    assert len(selected)==160304 and np.all(selected[:-1]<selected[1:])
    in_region=np.zeros(len(cells),bool);in_region[selected]=True
    regioncells=cells[selected];nodeids=np.unique(regioncells[regioncells>=0])
    old_to_local=np.full(len(xyz),-1,dtype=np.int64);old_to_local[nodeids]=np.arange(len(nodeids))
    owner=z['faceOwnerCell'];other=z['faceNeighbourCell'];f=z['motherFaces'].astype(np.int64)
    internal=other>=0;left=in_region[owner];right=np.zeros(len(f),bool);right[internal]=in_region[other[internal]]
    boundary=np.flatnonzero(left^right)
    wall_all=z['originalWallMotherFaceIDs'];far_all=z['originalFarfieldMotherFaceIDs']
    wall_ids=idx['upstreamOriginalFaces'].astype(np.int64)-1
    markers=dict(wall=wall_all[wall_ids],farfield=far_all[in_region[owner[far_all]]],
        unclosed_internal_interface=idx['interfaceMotherFaces'])
    assert np.array_equal(np.sort(np.concatenate(list(markers.values()))),boundary)
    assert len(wall_ids)==4434 and len(markers['unclosed_internal_interface'])==1984
    # 单独保存原号，负填充不会经由负下标偷偷变成末节点。
    localcells=regioncells.copy();mask=localcells>=0;localcells[mask]=old_to_local[localcells[mask]]
    localfaces={}
    for name,ids in markers.items():
        q=f[ids].copy();valid=q>=0;q[valid]=old_to_local[q[valid]]
        assert np.all(q[valid]>=0);localfaces[name]=q
    path=out/'region_geometry.su2'
    with path.open('w',encoding='utf8',newline='\n') as stream:
        stream.write('NDIME= 3\nNELEM= %d\n'%len(selected))
        for i,(cid,row) in enumerate(zip(selected,localcells)):
            stream.write('%d '%(13 if cid<len(w) else 10)+' '.join(map(str,row[row>=0]))+' %d\n'%i)
        stream.write('NPOIN= %d\n'%len(nodeids))
        for i,point in enumerate(xyz[nodeids]):stream.write('%.17g %.17g %.17g %d\n'%(*point,i))
        stream.write('NMARK= 3\n')
        for name,rows in localfaces.items():
            stream.write('MARKER_TAG= %s\nMARKER_ELEMS= %d\n'%(name,len(rows)))
            for row in rows:
                row=row[row>=0];stream.write(('%d '%(5 if len(row)==3 else 9))+' '.join(map(str,row))+'\n')
    state=read_frame(FRAME,335726)['values'];state=state[np.argsort(state[:,0])]
    assert np.array_equal(state[:,1:4],xyz)
    # 原存场直接逐位提取；物理能量换算另列，禁止将其当收敛边界条件。
    np.savez_compressed(str(out/'region.npz'),globalNodeIDs=nodeids,motherCellIDs=selected,
        originalFaceIDs=wall_ids+1,localCells=localcells,vertices=xyz[nodeids],
        regionState=state[nodeids],interfaceMotherFaces=markers['unclosed_internal_interface'],
        **{k+'LocalFaces':v for k,v in localfaces.items()})
    """模块三：实际接口状态与法向特征；黏性扩散不等同五条对流特征。"""
    faces=f[markers['unclosed_internal_interface']];air=PhaseAwareViscousAir();records=[]
    for fid,row,up,down in zip(markers['unclosed_internal_interface'],faces,
            idx['upstreamInterfaceCells'],idx['aftInterfaceCells']):
        nodes=row[row>=0];points=xyz[nodes];center=points.mean(axis=0)
        normal=.5*np.cross(points-center,np.roll(points,-1,axis=0)-center).sum(axis=0)
        cc=cells[up];upcenter=xyz[cc[cc>=0]].mean(axis=0)
        if np.dot(normal,center-upcenter)<0:normal=-normal
        area=np.linalg.norm(normal);assert area>0;normal/=area
        q=state[nodes,5:10].mean(axis=0);rho=q[0];v=q[1:4]/rho
        e=q[4]/rho-.5*np.dot(v,v)-300000
        temperature=float(inverse_energy(air,np.array([e]))[0]);pressure=rho*air.R*temperature
        thermo=air.evaluate(pressure,temperature);sound=float(thermo['soundSpeedMps'])
        un=float(np.dot(v,normal));lambdas=np.r_[np.full(3,un),un+sound,un-sound]
        incoming=int(sum(lambdas<0));reverse=int(sum(-lambdas<0));sonic=int(sum(abs(lambdas)<1e-10))
        nodeMn=(state[nodes,12:15]@normal)/state[nodes,18]
        records.append([fid,up,down,len(nodes),*center,*normal,area,pressure,temperature,un/sound,
            incoming,reverse,sonic,state[nodes,24].min(),state[nodes,24].max(),
            state[nodes,10].min(),state[nodes,10].max(),nodeMn.min(),nodeMn.max()])
    table=np.asarray(records)
    header=('motherFaceID,upstreamCellID,downstreamCellID,vertices,x,y,z,nx,ny,nz,areaM2,'
        'pPa,TK,normalMach,incomingUpstream,incomingDownstream,nearZeroEigenvalues,'
        'minimumFullDomainWallDistanceM,maximumFullDomainWallDistanceM,minimumNuTilde,maximumNuTilde,'
        'minimumVertexNormalMach,maximumVertexNormalMach')
    np.savetxt(str(out/'interface.csv'),table,delimiter=',',header=header,comments='',fmt='%.17g')
    source_rows=sources(out)
    classification=dict(allOutgoing=int(sum(table[:,14]==0)),oneIncoming=int(sum(table[:,14]==1)),
        fourIncoming=int(sum(table[:,14]==4)),fiveIncoming=int(sum(table[:,14]==5)),
        other=int(sum(~np.isin(table[:,14],[0,1,4,5]))),
        normalMachRange=[float(table[:,13].min()),float(table[:,13].max())])
    classification.update(wall_contact_stats(localfaces['unclosed_internal_interface'],localfaces['wall'],table))
    requires=int(sum(table[:,14]>0));assert requires>0
    proposal=dict(schema='RAH.RegionalBoundaryPlan.v1',regionGeometry=str(path),
        geometryReady=True,sourceState='v236 update 300 / research only',
        productionAdmitted=False,restartAdmitted=False,canBeUsedAsAcceptedEdge=False,
        missing=['current downstream acoustic/thermodynamic donor',
            'viscous and SA interface consistency / feedback',
            'native regional wall-distance scope and original dual-volume equivalence'],
        rejected=['all-outgoing Euler boundary on viscous cut','ambient p infinity at internal cut',
            'old unqualified state frozen as new physical boundary','MZ17 failed sliding interface'],
        adoptedNextCandidate='mature single-zone internal domain coupling with downstream retained',
        candidateMeaning='internal partition is not a physical boundary; no new CFD admission',
        stateArchiveIsNotBoundarySupply=True,sourceBinding=source_rows)
    save(out/'boundary_plan.json',proposal)
    report=dict(passed=True,geometryExported=True,meshSHA256=sha(path),regionalNodes=len(nodeids),
        regionalCells=len(selected),prisms=int(sum(selected<len(w))),tetrahedra=int(sum(selected>=len(w))),
        wallParents=4434,farfieldFaces=len(markers['farfield']),interfaceFaces=1984,
        originalFaces=19752,originalAreaM2=1.9646292565715553,fullOriginalStateCoordinatesBitEqual=True,
        sourceFrameSHA256=sha(FRAME),sourceFrameUpdate=300,sourceScientificQualification=False,
        sourceFullDomainWallDistancePreserved=True,regionalWallDistanceRecomputed=False,
        primitiveFaceCharacteristics=classification,facesRequiringDownstreamInformation=requires,
        actualNativeDualBoundaryCharacteristics=False,viscousDiffusionAccountedByEigenvalues=False,
        fixedCutProductionRejected=True,sourceRows=source_rows,
        flowUpdates=0,nativeQueries=0,builds=0,newQualifiedAreaM2=0,physicalQualification=False)
    save(out/'report.json',report)
    print({k:report[k] for k in ('regionalNodes','regionalCells','farfieldFaces',
        'primitiveFaceCharacteristics','facesRequiringDownstreamInformation')})


def main():
    if CLOSED:raise SystemExit('v244 preparation closed; use NEXT_DECISION')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;size=guard.size_bytes
    guard.size_bytes=lambda p:size(BASE) if Path(p)==BASE or Path(p).name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
        VOLUME,PART,INDEX,FRAME,ROOT/'scripts/singleDomainFrame.py',ROOT/'scripts/phaseAwareViscousAir.py',
        ROOT/'scripts/su2PhaseAirFields.py',ROOT/'scripts/phaseAwareDryAir.py',
        ROOT/'data/thermo/phase_aware_dry_air_v1.json',
        ROOT/'docs/validation/v243/evidence/final_report.json',
        ROOT/'outputs/v233_validation/coupling01/report.json',
        ROOT/'docs/validation/v233/evidence/local_asset_index.json',
        ROOT/'outputs/v236_validation/process_audit3/report.json',SOURCES/'bindings.json',
        ROOT/'outputs/v220_validation/build02/dependencies.json',
        ROOT/'outputs/v220_validation/source/Common/src/CConfig.cpp']
    locks += [SOURCES/rel for rel in ('SU2_CFD/src/solvers/CEulerSolver.cpp',
        'SU2_CFD/src/solvers/CTurbSASolver.cpp','SU2_CFD/src/integration/CIntegration.cpp')]
    return guard.supervise('prepare2',[sys.executable,str(Path(__file__).resolve()),'task'],120,256,locks,BASE,'offline')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v244 preparation closed; use NEXT_DECISION')
    if sys.argv[1]=='task':prepare()
    else:raise SystemExit(0 if main() else 1)
