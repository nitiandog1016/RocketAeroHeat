"""模块职责：监督原前体的平面裁切、原面映射、同气体出口检查和实际求解配置导出。"""
from pathlib import Path
import json, sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT, read, save, sha
from clipEulerForebody import clip_tetrahedra, interpolate_saved_field, write_su2
from phaseEulerRegionalVolume import read_tetra_mesh, topology, extract_region, interface_characteristics

BASE=ROOT/'outputs/v240_validation'; DOC=ROOT/'docs/validation/v240'
MESH=ROOT/'outputs/v232_validation/rocket01/mesh.su2'
STATE=ROOT/'outputs/v232_validation/rocket01/state_000800.csv'
REGION=ROOT/'outputs/v228_validation/prepare03/region.json'
GEOM=ROOT/'outputs/v230_validation/audit02/geometry.npz'
SOURCE=ROOT/'outputs/v229_validation/source/SU2_CFD/src/solvers/CEulerSolver.cpp'
CLOSED=True


def sources():
    """模块一：绑定真实已编译的一般EOS边界源，保留采用/排除理由。"""
    out=BASE/'source01';out.mkdir(exist_ok=False)
    known=read(ROOT/'outputs/v229_validation/process_build/report.json')['inputSHA256']
    assert sha(SOURCE)==known[str(SOURCE)]
    rows=SOURCE.read_text(encoding='utf8').splitlines()
    excerpts=[]
    for name, needle, count in [('riemann','void CEulerSolver::BC_Riemann(',620),
                                ('supersonic_outlet','void CEulerSolver::BC_Supersonic_Outlet(',123)]:
        indices=[i for i,s in enumerate(rows) if needle in s];assert len(indices)==1
        start=indices[0];end=min(len(rows),start+count)
        path=out/(name+'.txt');path.write_text('\n'.join('%d %s'%(i+1,rows[i]) for i in range(start,end))+'\n',encoding='utf8')
        excerpts.append(dict(source=str(SOURCE),sha256=sha(SOURCE),firstLine=start+1,lastLine=end,file=path.name))
    save(out/'report.json',dict(passed=True,excerpts=excerpts,adopted='mature single-region phase-Air Euler / characteristic outlet conditional on all outgoing',
      rejected='direct SUPERSONIC_OUTLET lacks SetSecondary for general EOS in exact bound source',
      pnsAdopted=False,flowUpdates=0,nativeQueries=0,builds=0,physicalQualification=False))


def produce():
    """模块二：从原面边界环生成实际区域网格，而非质心截取或替代几何。"""
    out=BASE/'region01';out.mkdir(exist_ok=False)
    previous=read(ROOT/'outputs/v233_validation/region01/report.json')
    for p in (MESH,STATE,REGION,GEOM):assert sha(p)==previous['sourceSHA256'][p.relative_to(ROOT).as_posix()]
    xyz,cells,markers=read_tetra_mesh(MESH);g=np.load(str(GEOM));region=read(REGION)
    original= g['volumeNodeIDs'][g['faces']]
    responsibility=np.r_[region['retainedNoseCandidateIDs'],region['originalFaceIDs']].astype(int)-1
    assert len(responsibility)==4434 and len(np.unique(responsibility))==4434
    selected=original[responsibility]
    edge=np.sort(np.concatenate([selected[:,[0,1]],selected[:,[1,2]],selected[:,[2,0]]]),axis=1)
    keys,count=np.unique(edge,axis=0,return_counts=True);assert np.all(count<=2)
    ring_edges=keys[count==1];ring_nodes=np.unique(ring_edges);ring_x=xyz[ring_nodes,0]
    assert len(ring_nodes)>3 and np.all(np.unique(ring_edges,return_counts=True)[1]==2)
    cut=float(np.median(ring_x));spread=float(np.ptp(ring_x))
    assert spread<=64*np.finfo(float).eps*max(1.,abs(cut)),('ORIGINAL_RING_NOT_PLANAR',spread)
    result=clip_tetrahedra(xyz,cells,markers,cut)
    topo=topology(result['xyz'],result['cells'],result['markers'])
    closure=extract_region(topo,np.ones(len(result['cells']),bool));assert closure['summary']['components']==1
    wall=result['markers']['wall'];oldwall=result['originalNode'][wall];assert np.all(oldwall>=0)
    lookup={tuple(sorted(row)):int(i) for i,row in enumerate(original)}
    ids=np.array([lookup[tuple(sorted(row))] for row in oldwall]);assert np.array_equal(np.sort(ids),np.sort(responsibility))
    assert np.array_equal(result['xyz'][wall],xyz[oldwall])
    normals=topo['normal'][topo['markers']['outlet']]
    unit=normals/np.linalg.norm(normals,axis=1)[:,None]
    assert np.max(abs(unit-[-1.,0.,0.]))<1e-11
    write_su2(out/'mesh.su2',result)
    xx,cc,mm=read_tetra_mesh(out/'mesh.su2')
    assert np.array_equal(xx,result['xyz']) and np.array_equal(cc,result['cells'])
    assert all(np.array_equal(mm[k],v) for k,v in result['markers'].items())
    np.savez_compressed(str(out/'provenance.npz'),**{k:v for k,v in result.items() if k not in ('markers',)})
    np.savetxt(str(out/'original_faces.csv'),np.column_stack([ids+1,wall,g['area'][ids]]),delimiter=',',fmt='%.17g',comments='',
      header='originalFaceID,node0,node1,node2,originalAreaM2')
    # 只使用已有800场判断候选出口；没有旧场重启或新的流动推进。
    raw=np.loadtxt(str(STATE),delimiter=',',skiprows=1);raw=raw[np.argsort(raw[:,0])]
    assert np.array_equal(raw[:,0],np.arange(len(xyz))) and np.array_equal(raw[:,1:4],xyz)
    u=interpolate_saved_field(result,raw[:,4:9])
    faces=result['markers']['outlet'];nodes=np.unique(faces)
    tri=interface_characteristics(u,faces,normals)
    point=interface_characteristics(u,np.repeat(nodes[:,None],3,axis=1),np.tile([-1.,0.,0.],(len(nodes),1)))
    np.savetxt(str(out/'outlet_nodes.csv'),np.column_stack([nodes,result['xyz'][nodes],u[nodes],point['p'],point['T'],point['H'],point['S'],point['normalMach'],point['incoming']]),
      delimiter=',',fmt='%.17g',comments='',header='nodeID,x,y,z,rho,rhou,rhov,rhow,rhoEsolver,pPa,TK,HJkg,SJkgK,normalMach,incoming')
    np.savetxt(str(out/'outlet_faces.csv'),np.column_stack([faces,normals,tri['normalMach'],tri['incoming']]),delimiter=',',fmt='%.17g',comments='',
      header='node0,node1,node2,Ax,Ay,Az,normalMach,incoming')
    potential=bool(np.all(point['normalMach']>1) and np.all(tri['normalMach']>1))
    config=(ROOT/'outputs/v232_validation/rocket01/case.cfg').read_text(encoding='utf8')
    old=[s for s in config.splitlines() if s.startswith('MARKER_RIEMANN=')];assert len(old)==1
    reservoir=old[0].split('(farfield,')[1].strip().rstrip(')')
    config=config.replace(old[0],old[0][:-1]+', outlet, '+reservoir+')')
    # ITER只是未来候选上限；实际执行仍必须有另预算准入和动态特征守卫。
    config=config.replace('ITER= 800','ITER= 400')
    (out/'case.cfg').write_text(config,encoding='ascii')
    volume=topo['volume'];parent=result['motherCell'];cut_mask=np.isin(parent,result['partialMotherCells'])
    local_volume=volume[cut_mask]
    report=dict(passed=True,method='original tetrahedra clipped at exact original wall-ring plane',
      originalFaces=19752,originalAreaM2=float(np.sum(g['area'])),regionWallFaces=len(ids),regionWallAreaM2=float(np.sum(g['area'][ids])),
      engineeringTargetFaces=3128,engineeringTargetAreaM2=region['areaM2'],cutXM=cut,ringNodes=len(ring_nodes),ringXSpreadM=spread,
      originalWallFacesUnsplit=True,originalWallCoordinatesBitEqual=True,su2RoundTripBitEqual=True,
      cells=len(result['cells']),nodes=len(result['xyz']),partialMotherCells=len(result['partialMotherCells']),
      retainedMotherCells=len(np.unique(parent)),newCellsFromCutMothers=int(sum(cut_mask)),newNodes=int(sum(result['originalNode']<0)),
      closure=closure['summary'],minimumCutCellVolumeM3=float(min(local_volume)),maximumCutCellVolumeM3=float(max(local_volume)),
      outletFaces=tri['summary'],outletNodes=dict(nodes=len(nodes),normalMachRange=point['summary']['normalMachRange'],
        incomingNodeCounts=point['summary']['characteristicFaceCounts'],method='same phase Air at individual outlet vertices; no nodal area assigned'),candidateOneWayOutletPassed=potential,
      sourceStateUpdate=800,sourceStateQualified=False,interpolatedStateUsedAsRestart=False,
      originalViscousMeshChanged=False,wholeAreaResponsibilityPreserved=True,
      configurationExported=True,actualNativeBoundaryTested=False,CFDAdmitted=False,
      physicalQualification=False,newQualifiedAreaM2=0,flowUpdates=0,
      sourceSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in (MESH,STATE,REGION,GEOM)},
      meshSHA256=sha(out/'mesh.su2'),configSHA256=sha(out/'case.cfg'))
    save(out/'report.json',report);print(json.dumps(report))


def main(stage):
    if CLOSED:raise SystemExit('v240 preparation closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original_size=guard.size_bytes;guard.size_bytes=lambda _:original_size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/clipEulerForebody.py',ROOT/'scripts/phaseEulerRegionalVolume.py',
      ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',
      ROOT/'data/thermo/phase_aware_dry_air_v1.json',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
    if stage=='sources':locks += [SOURCE,ROOT/'outputs/v229_validation/process_build/report.json']
    elif stage=='produce':locks += [MESH,STATE,REGION,GEOM,ROOT/'outputs/v233_validation/region01/report.json',ROOT/'outputs/v232_validation/rocket01/case.cfg']
    elif stage=='test':locks += [ROOT/'tests/testForebodyClipV240.py']
    else:raise ValueError('STAGE_NOT_ADMITTED')
    argv=[sys.executable,str(Path(__file__).resolve()),'task_'+stage]
    if stage=='test':argv=[sys.executable,str(ROOT/'tests/testForebodyClipV240.py'),str(BASE/'test01')]
    return guard.supervise(stage,argv,120,256,locks,BASE)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v240 preparation closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_sources':sources()
    elif sys.argv[1]=='task_produce':produce()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
