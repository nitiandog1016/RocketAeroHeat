"""模块职责：原完整壁面相态Euler的离线准备/独立映射审查；不运行任何引擎。"""
from pathlib import Path
import sys,json,re,copy
import numpy as np
from scipy.spatial import cKDTree
from phaseAwareViscousAir import PhaseAwareViscousAir
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v228_validation';DOC=ROOT/'docs/validation/v228'
MESH=ROOT/'outputs/v194_validation/rocket01/mesh.su2'
ENGINE=ROOT/'outputs/v220_validation/build02/SU2_CFD.exe'
IDENTITY='rah-su2-phase-air-euler-eplus300k-v1'
SOURCE_PATHS=['Common/src/CConfig.cpp','Common/include/option_structure.hpp',
 'SU2_CFD/src/solvers/CEulerSolver.cpp','SU2_CFD/src/variables/CEulerVariable.cpp',
 'SU2_CFD/include/solvers/CFVMFlowSolverBase.hpp','SU2_CFD/include/solvers/CFVMFlowSolverBase.inl',
 'SU2_CFD/src/solvers/CSolverFactory.cpp','SU2_CFD/src/integration/CIntegration.cpp',
 'SU2_CFD/src/output/CFlowCompOutput.cpp','SU2_CFD/src/numerics/flow/convection/roe.cpp']

def read_mesh():
    # 只解码现成SU2连接及壁面ID；不生成网格或重复旧质量检查。
    lines=MESH.read_text(encoding='utf8').splitlines();i=1
    assert lines[0].strip()=='NDIME= 3';nc=int(lines[i].split('=')[1]);i+=1
    cells=np.array([[int(z) for z in l.split()[1:5]] for l in lines[i:i+nc]],int);i+=nc
    npnt=int(lines[i].split('=')[1].split()[0]);i+=1
    V=np.array([[float(z) for z in l.split()[:3]] for l in lines[i:i+npnt]]);i+=npnt
    nm=int(lines[i].split('=')[1]);i+=1;markers={}
    for _ in range(nm):
        name=lines[i].split('=')[1].strip();nt=int(lines[i+1].split('=')[1]);i+=2
        markers[name]=np.array([[int(z) for z in l.split()[1:4]] for l in lines[i:i+nt]],int);i+=nt
    assert (nc,npnt)==(64982,18027) and set(markers)=={'wall','farfield'}
    return V,cells,markers

def bound_sources():
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    old=read(ROOT/'outputs/v167_validation/build01/source_tree.json');result={};text={}
    for rel in SOURCE_PATHS:
        p=ROOT/'outputs/v220_validation/source'/rel;baseline='v220-build02-dependency'
        if p.is_file():expected=deps[str(p)]
        else:
            p=ROOT/'tools/runtime/su2-thermal-build/source'/rel;expected=old[str(p)];baseline='v167-object-source-reused-by-v220'
        assert sha(p)==expected,rel
        text[rel]=p.read_text(encoding='utf8');result[rel]=dict(path=str(p),sha256=expected,binding=baseline)
    return result,text

def check_design(c):
    """本地设计合同守卫，测试通过不表示原生新身份已经注册。"""
    assert c['identity']==IDENTITY and c['nativeFluidModel']=='RAH_PHASE_AIR_EULER'
    assert c['solver']=='EULER' and c['energyShiftJkg']==300000 and c['physicalEnthalpyReference']=='hs(298.15K)=0'
    assert c['meshSHA256']==sha(MESH) and c['originalFaces']==19752
    assert c['wallBoundary']=='MARKER_EULER' and c['outerBoundary']=='general-EOS-ambient-PT-characteristics'
    assert c['turbulenceModel']=='NONE' and c['heatFluxAvailable'] is False
    assert c['engineRegistered'] is False and c['runAuthorized'] is False
    assert c['fullSurfaceResponsibility'] is True and c['engineeringHeatReady'] is False

def prepare():
    out=BASE/'prepare01';out.mkdir(exist_ok=False)
    origin=read(ROOT/'outputs/v194_validation/rocket01/manifest.json');assert sha(MESH)==origin['filesSHA256']['mesh.su2']
    assert read(ROOT/'outputs/v194_validation/native_audit01/report.json')['passed']
    V,cells,markers=read_mesh();surface=read(ROOT/'outputs/v194_validation/rocket01/surface.json')
    W=np.array(surface['vertices']);F=np.array(surface['facesOneBased'],int)-1
    distance,point_map=cKDTree(V).query(W);assert np.max(distance)==0 and len(set(point_map))==len(W)
    original={tuple(sorted(row)):i+1 for i,row in enumerate(point_map[F])}
    face_map=np.array([original[tuple(sorted(row))] for row in markers['wall']],int)
    assert np.array_equal(np.sort(face_map),np.arange(1,19752+1))
    area=np.linalg.norm(np.cross(W[F[:,1]]-W[F[:,0]],W[F[:,2]]-W[F[:,0]]),axis=1)/2
    assert abs(area.sum()-1.9646292565715553)<2e-13
    design=read(ROOT/'docs/validation/v212/evidence/design01/regions.json')
    owned=set(design['regions']['nose']['ownedFaceIDs'])
    oldfaces=np.loadtxt(str(ROOT/'outputs/v226_validation/matlab02/rocket_faces.csv'),delimiter=',',skiprows=1)
    existing=set(oldfaces[oldfaces[:,13]==1,0].astype(int));assert len(existing)==1306 and existing<=owned
    target=sorted(owned-existing);remainder=sorted(set(range(1,19753))-owned)
    centers=W[F].mean(axis=1)
    region=dict(name='forebody-after-fitted-cap-before-contraction',status='selected-responsibility-not-qualified',
      originalFaceIDs=target,faces=len(target),areaM2=float(area[np.array(target)-1].sum()),
      minimumCenterXM=float(centers[np.array(target)-1,0].min()),maximumCenterXM=float(centers[np.array(target)-1,0].max()),
      originalOwnershipCutXM=design['ownershipCutsXM'][0],retainedNoseCandidateIDs=sorted(existing),
      remainderFaceIDs=remainder,remainderAreaM2=float(area[np.array(remainder)-1].sum()),
      totalOriginalAreaM2=float(area.sum()),physicalApplicabilityAssumed=False)
    sources,text=bound_sources();cfg=text['Common/src/CConfig.cpp'];enum=text['Common/include/option_structure.hpp']
    assert '(Kind_Solver != MAIN_SOLVER::NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::RANS)' in cfg
    assert 'RAH_PHASE_AIR_EULER' not in enum
    assert 'BC_Sym_Plane(geometry, solver_container' in text['SU2_CFD/include/solvers/CFVMFlowSolverBase.hpp']
    assert 'FluidModel->SetTDState_rhoe(density, staticEnergy);' in text['SU2_CFD/src/variables/CEulerVariable.cpp']
    # 对准确源只提取新Euler分派差异；不重新查询整套物性、热壁、Roe控制。
    needles={
      'Common/src/CConfig.cpp':['if (Kind_FluidModel == RAH_PHASE_AIR)', 'requires direct SI steady single-zone NS/SA'],
      'Common/include/option_structure.hpp':['RAH_PHASE_AIR =','MakePair("RAH_PHASE_AIR"'],
      'SU2_CFD/src/variables/CEulerVariable.cpp':['bool CEulerVariable::SetPrimVar'],
      'SU2_CFD/include/solvers/CFVMFlowSolverBase.hpp':['inline void BC_Euler_Wall'],
      'SU2_CFD/include/solvers/CFVMFlowSolverBase.inl':['void CFVMFlowSolverBase<V, FlowRegime>::BC_Sym_Plane'],
      'SU2_CFD/src/solvers/CSolverFactory.cpp':['MAIN_SOLVER::EULER'],
      'SU2_CFD/src/integration/CIntegration.cpp':['BC_Euler_Wall'],
      'SU2_CFD/src/solvers/CEulerSolver.cpp':['case RAH_PHASE_AIR:', 'void CEulerSolver::BC_Riemann'],
      'SU2_CFD/src/numerics/flow/convection/roe.cpp':['RAH_PHASE_ROE_AVERAGE_REJECTED'],
      'SU2_CFD/src/output/CFlowCompOutput.cpp':['LoadVolumeData']}
    excerpts={}
    for rel,terms in needles.items():
        lines=text[rel].splitlines();selected=set()
        for term in terms:
            hits=[i for i,l in enumerate(lines) if term in l];assert hits,(rel,term)
            for i in hits[:2]:selected.update(range(max(0,i-3),min(len(lines),i+24)))
        excerpts[rel]=[dict(line=i+1,text=lines[i]) for i in sorted(selected)]
    condition=read(ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json')['condition'];air=PhaseAwareViscousAir()
    state=air.evaluate(condition['pressurePa'],condition['temperatureK']);U=condition['speedMps'];rho=float(state['densityKgm3'])
    velocity=np.array(condition['flowDirection'])*U
    shifted=float(state['enthalpyJkg'])-air.R*condition['temperatureK']+U*U/2+300000
    initial=dict(densityKgm3=rho,velocityMps=velocity.tolist(),physicalTotalEnthalpyJkg=float(state['enthalpyJkg'])+U*U/2,
      shiftedTotalEnergyJkg=shifted,conservedDensityMomentumEnergy=[rho]+(rho*velocity).tolist()+[rho*shifted],
      flowTimeSeconds=0,solutionStatus='uncomputed-compatible-initial-state')
    # 仅核算五变量稀疏矩阵/一个状态的明确子项，不把它当进程峰值。
    ed=set()
    for a,b in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)):
        ed.update(map(tuple,np.sort(cells[:,[a,b]],axis=1)))
    blocks=len(V)+2*len(ed);resources=dict(vertices=len(V),tetrahedra=len(cells),uniquePrimalEdges=len(ed),variables=5,
      singleStateBytes=len(V)*5*8,blockMatrixScalarPayloadBytes=blocks*25*8,blockCount=blocks,
      matrixScope='edge block graph payload only; solver actual allocation/halos/index/gradient/preconditioner not measured',
      peakBytes=None,fitProven=False,proposedObservedGroupLimitBytes=4*1024**3,proposedCpuLimit=1,
      previousRANSObservedPeakBytes=read(ROOT/'outputs/v222_validation/process_engine/report.json')['peakObservedGroupRSSBytes'])
    contract=dict(schema='RAH.PhaseEulerDesign.v1',identity=IDENTITY,nativeFluidModel='RAH_PHASE_AIR_EULER',solver='EULER',
      energyShiftJkg=300000,physicalEnthalpyReference='hs(298.15K)=0',meshPath=str(MESH),meshSHA256=sha(MESH),
      originalFaces=19752,originalAreaM2=float(area.sum()),wallBoundary='MARKER_EULER',
      outerBoundary='general-EOS-ambient-PT-characteristics',turbulenceModel='NONE',condition=condition,
      heatFluxAvailable=False,wallTemperatureSolved=False,engineRegistered=False,runAuthorized=False,
      fullSurfaceResponsibility=True,engineeringHeatReady=False,physicalQualification=False,
      currentEnginePath=str(ENGINE),currentEngineSHA256=sha(ENGINE),currentEngineRejectedByExistingNSGuard=True,
      pending=['new isolated registration and all gas dispatch branches','native Euler slip/primitive/characteristic/save-restart compatibility',
               'accurate actual process resource and full development budget','steady residual/stability/balance and spatial/farfield error',
               'wall-adjacent vs BL-edge sampling and original-face qualification',
               'labelled three-dimensional paths/entropy/total-enthalpy/metric and nonsimilar inlet',
               'separation/base/strong crossflow mature viscous feedback'])
    check_design(contract)
    # 安全准备为ITER=0的提案。尚未登记实际发展；不提供可自动启动CFD的入口。
    native={}
    for line in (ROOT/'outputs/v222_validation/rocket01/case.cfg').read_text(encoding='utf8').splitlines():
        if '=' in line:k,v=line.split('=',1);native[k.strip()]=v.strip()
    remove=['KIND_TURB_MODEL','VISCOSITY_MODEL','CONDUCTIVITY_MODEL','MARKER_ISOTHERMAL','TURBULENT_CONDUCTIVITY_MODEL',
      'PRANDTL_TURB','CONV_NUM_METHOD_TURB','MUSCL_TURB','TIME_DISCRE_TURB','CFL_REDUCTION_TURB','SA_OPTIONS',
      'FREESTREAM_NU_FACTOR','KIND_TRANS_MODEL','USE_ACCURATE_TURB_JACOBIANS']
    for key in remove:native.pop(key,None)
    native.update(SOLVER='EULER',FLUID_MODEL='RAH_PHASE_AIR_EULER',MARKER_EULER='(wall)',ITER='0',
      MESH_FILENAME=MESH.as_posix(),HISTORY_OUTPUT='(ITER, RMS_RES, AERO_COEFF, LINSOL)',
      CONV_FIELD='(RMS_DENSITY, RMS_MOMENTUM-X, RMS_MOMENTUM-Y, RMS_MOMENTUM-Z, RMS_ENERGY)')
    (out/'case.proposed.cfg').write_text('\n'.join(k+'= '+v for k,v in native.items())+'\n',encoding='utf8')
    np.savetxt(str(out/'original_face_map.csv'),np.c_[np.arange(len(face_map)),face_map],fmt='%d',delimiter=',',header='wallMarkerRowZeroBased,originalFaceID',comments='')
    for name,value in [('contract.json',contract),('region.json',region),('sources.json',sources),('source_excerpts.json',excerpts),
      ('resources.json',resources),('initial_state.json',initial)]:save(out/name,value)
    save(out/'report.json',dict(passed=True,designOnly=True,originalFaceMappingExact=True,newMeshGenerated=False,
      targetFaces=region['faces'],targetAreaM2=region['areaM2'],totalAreaM2=float(area.sum()),
      currentEngineCanRunCandidate=False,sourceFiles=len(sources),CFDRuns=0,newHeatAreaM2=0,newQualifiedAreaM2=0))
    print(json.dumps(read(out/'report.json')))

def audit():
    out=BASE/'audit01';out.mkdir(exist_ok=False);p=BASE/'prepare01';c=read(p/'contract.json');check_design(c)
    # 独立从坐标三角多重集核对映射；不复用prepare的KD树映射。
    V,_,markers=read_mesh();s=read(ROOT/'outputs/v194_validation/rocket01/surface.json')
    W=np.array(s['vertices']);F=np.array(s['facesOneBased'],int)-1
    key=lambda tri:tuple(sorted(tuple(row) for row in tri))
    lookup={key(W[f]):i+1 for i,f in enumerate(F)}
    expected=np.array([lookup[key(V[f])] for f in markers['wall']],int)
    actual=np.loadtxt(str(p/'original_face_map.csv'),delimiter=',',skiprows=1,dtype=int)
    assert np.array_equal(expected,actual[:,1]) and len(set(expected))==19752
    r=read(p/'region.json');groups=[set(r[k]) for k in ('originalFaceIDs','retainedNoseCandidateIDs','remainderFaceIDs')]
    assert not(groups[0]&groups[1] or groups[0]&groups[2] or groups[1]&groups[2])
    assert set.union(*groups)==set(range(1,19753))
    lengths=np.stack([np.linalg.norm(W[F[:,j]]-W[F[:,(j+1)%3]],axis=1) for j in range(3)])
    half=lengths.sum(axis=0)/2;area=np.sqrt(half*(half-lengths[0])*(half-lengths[1])*(half-lengths[2]))
    assert abs(area[np.array(r['originalFaceIDs'])-1].sum()-r['areaM2'])<1e-14
    rejected=0
    for key,value in [('identity','old'),('nativeFluidModel','RAH_PHASE_AIR'),('solver','RANS'),('energyShiftJkg',0),
      ('wallBoundary','MARKER_ISOTHERMAL'),('outerBoundary','fixed-pinf-internal'),('heatFluxAvailable',True),
      ('fullSurfaceResponsibility',False),('runAuthorized',True),('engineeringHeatReady',True)]:
        bad=copy.deepcopy(c);bad[key]=value
        try:check_design(bad)
        except AssertionError:rejected+=1
    assert rejected==10
    for rel,v in read(p/'sources.json').items():assert sha(v['path'])==v['sha256'],rel
    report=dict(passed=True,designOnly=True,independentCoordinateTriangleMap=True,originalFaces=19752,
      responsibilityPartitionComplete=True,targetFaces=r['faces'],targetAreaM2=r['areaM2'],rejections=rejected,
      currentNativeEngineAdmission=False,newQualifiedAreaM2=0,numericalRuns=0)
    save(out/'report.json',report);print(json.dumps(report))
if __name__=='__main__':
    if sys.argv[1]=='prepare01':prepare()
    elif sys.argv[1]=='audit01':audit()
    else:raise SystemExit('No native or CFD execution entry')
