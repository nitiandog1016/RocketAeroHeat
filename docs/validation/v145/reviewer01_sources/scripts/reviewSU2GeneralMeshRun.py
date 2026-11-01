"""只读核对v145初场、原面导出、旧下限位置及偏移列离散；不重跑流场。"""
from pathlib import Path
import csv,json,time
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from su2SSTInterface import state,triangle_distances
from auditSU2SSTState import vtk_tokens
from prepareSU2SSTPocketStudy import floor_nodes
from reviewSU2SSTLocalDiffusion import local_geometry

def main():
    start=time.monotonic();base=ROOT/'outputs/v145_validation';target=base/'independent_review.json';assert not target.exists()
    case=Path(read(base/'main01/engine_folder.json')['folder']);r=read(case/'engine_process.json');assert r['exitCode']==0 and all(r['inputsStable'].values())
    assert read(base/'main01_process.json')['status']=='terminal' and all(read(base/'main01_process.json')['sourceStableAfter'].values())
    for n in ['restart_flow.csv','surface_flow.vtk','history.csv']:assert sha(case/n)==r['outputSHA256'][n]
    mesh=loadmat(str(base/'mesh/hybrid_volume.mat'));P=mesh['vertices'];V=mesh['sourceVertices'];F=mesh['sourceFaces'].astype(np.int64)-1
    W=mesh['wedges'].astype(np.int64)-1;Tet=mesh['tetrahedra'].astype(np.int64)-1
    a,keys,_,rho,vel,energy=state(case/'restart_flow.csv',True,True);assert np.array_equal(a[:,1:4],P)
    initial,ki,_,ri,ui,ei=state(case/'initial_flow.csv',True);assert np.array_equal(initial[:,1:4],P)
    geo=np.load(str(case/'wall_distance_initialization.npz'));b=geo['distance']/.001;b=b*b*(3-2*b)
    cfg=config(case/'case.cfg');R=float(cfg['GAS_CONSTANT']);gamma=float(cfg['GAMMA_VALUE']);cp=gamma*R/(gamma-1)
    Ti=ei*(gamma-1)/R;U=np.asarray(read(case/'case_manifest.json')['expectedVelocity']);pinf=float(cfg['FREESTREAM_PRESSURE']);Tinf=float(cfg['FREESTREAM_TEMPERATURE'])
    hInitial=cp*Ti+.5*np.sum(ui**2,axis=1);hExpected=cp*300*(1-b)+(cp*Tinf+.5*(U@U))*b
    initError=float(max(abs(hInitial/hExpected-1)));assert initError<1e-12 and max(abs(ri*R*Ti/pinf-1))<1e-12
    # 不同算法穷举全部原三角形，对均匀抽样和全部13个偏移点核对实际壁距。
    moved=mesh['baseColumnNodeIDsOneBased'].ravel().astype(np.int64)[20:]-1
    sample=np.unique(np.r_[np.linspace(0,len(P)-1,181).astype(np.int64),moved,geo['wallIDs'][::211]])
    exact=np.array([min(triangle_distances(p,V[F[:,0]],V[F[:,1]],V[F[:,2]])) for p in P[sample]])
    distanceError=float(max(abs(geo['distance'][sample]-np.minimum(exact,.001))));assert distanceError<2e-13
    # 原7点属于尾缘；用明确列/层及双法向扇映射到两侧，不做最近邻旧场填值。
    old=read(ROOT/'docs/validation/v132/evidence/gradient_comparison.json')['variants']['gg']['snapshot']['sourceFloorNodesRemaining']
    fan=read(base/'mesh/mesh_report.json');rim=list(fan['rimNodeIDs']);mapped=[]
    oldP=loadmat(str(ROOT/'outputs/v120_validation/mesh32/hybrid_volume.mat'),variable_names=['vertices'])['vertices']
    omega=a[:,keys.index('Omega')];k=a[:,keys.index('Turb_Kin_Energy')];oldCandidates=[]
    for node in old:
        layer,col=divmod(node,len(V));pos=rim.index(col);sides={}
        for label,key in [('side','sideNodeIDs'),('base','baseNodeIDs')]:
            new=int(fan[key][layer][pos]);oldCandidates.append(new)
            sides[label]=dict(nodeIDZeroBased=new,xyzM=P[new].tolist(),oldCoordinateDisplacementM=float(np.linalg.norm(P[new]-oldP[node])),
                omega=float(omega[new]),k=float(k[new]),omegaAtFloor=bool(new in set(floor_nodes(omega))))
        mapped.append(dict(oldNodeIDZeroBased=node,sourceColumnOneBased=col+1,layer=layer,branches=sides))
    # 只读局部几何算子：是仿射再现/条件数诊断，不冒充原生RHS或网格敏感性证明。
    local=local_geometry(P,[W,Tet],set(map(int,moved)))
    distortion=[dict(nodeIDZeroBased=int(i),xyzM=P[i].tolist(),weightedLeastSquaresCondition=g['condition'],
        weightedLeastSquaresAffineError=g['linearError'],greenGaussAffineError=g['ggLinearError'],
        volumeClosure=g['closure'],temperatureK=float(a[i,keys.index('Temperature')])) for i,g in local.items()]
    with (case/'history.csv').open() as f:rows=list(csv.reader(f,skipinitialspace=True))
    h=dict(zip([n.strip() for n in rows[0]],np.asarray(rows[1:],float).T));assert np.array_equal(h['Inner_Iter'],np.arange(150))
    residual={name:dict(final=float(values[-1]),minimum=float(min(values)),last25Mean=float(np.mean(values[-25:]))) for name,values in h.items() if name.startswith('rms[')}
    assert len(residual)==7
    xyz,fields=vtk_tokens(case/'surface_flow.vtk');d,ids=cKDTree(V).query(xyz);assert max(d)<1e-6 and len(set(ids))==len(V)
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    with (base/'main01/fields.csv').open() as f:export=list(csv.DictReader(f))
    col=lambda n:np.array([float(row[n]) for row in export]);assert len(export)==len(F)
    assert np.array_equal(col('originalFaceID'),np.arange(1,len(F)+1)) and max(abs(col('area_m2')/area-1))<1e-12
    face={};errors={}
    for name,key in [('pressureRaw_Pa','Pressure'),('heatFluxRaw_W_m2','Heat_Flux'),('fluidNodeTemperatureMean_K','Temperature')]:
        nodes=np.empty(len(V));nodes[ids]=fields[key][:,0];face[key]=nodes[F].mean(axis=1)
        errors[name]=float(max(abs(col(name)-face[key])/np.maximum(abs(face[key]),1)));assert errors[name]<1e-12
    assert np.all(col('givenWallTemperature_K')==300) and np.all(np.isnan(col('solvedWallTemperature_K')))
    assert not np.any(col('productionQualified')) and not np.any(col('heatFluxIsBoundaryInput'))
    midxyz,mid=vtk_tokens(case/'surface_flow_000075.vtk');assert np.array_equal(midxyz,xyz)
    nodes=np.empty(len(V));nodes[ids]=mid['Heat_Flux'][:,0];midface=nodes[F].mean(axis=1)
    delta=float(area@abs(face['Heat_Flux']-midface)/(area@abs(midface)))
    with (base/'main01/section.csv').open(encoding='utf-8-sig') as f:section=list(csv.reader(f))
    assert len(section)==93;endpointError=0.;qerr=0.
    for row in section[1:]:
        x=float(row[0]);i=int(row[2])-1;assert x==-1.5;tri=V[F[i]]
        for point in [np.array([x,float(row[5]),float(row[6])]),np.array([x,float(row[7]),float(row[8])])]:
            uv=np.linalg.lstsq((tri[1:]-tri[0]).T,point-tri[0],rcond=None)[0]
            err=float(np.linalg.norm((tri[1:]-tri[0]).T@uv-(point-tri[0])));endpointError=max(endpointError,err)
            assert err<1e-10 and min(uv)>-1e-10 and sum(uv)<1+1e-10
        qerr=max(qerr,abs(float(row[10])-face['Heat_Flux'][i])/max(abs(float(row[10])),1));assert 'SST-2003m' in row[12]
    assert qerr<1e-12
    # 同阶段流动primitive总焓不含更新后的k；只作为诊断，不能混阶段解释因果。
    T=a[:,keys.index('Temperature')];up=np.column_stack([a[:,keys.index('Velocity_'+q)] for q in 'xyz']);ht=cp*T+.5*np.sum(up*up,axis=1)
    extrem=int(np.argmax(ht));hInf=cp*Tinf+.5*(U@U);sa=read(case/'sst_state_audit.json')
    result=dict(passed=True,seconds=time.monotonic()-start,initialTotalEnthalpyRelativeError=initError,exactDistanceSampleCount=len(sample),
        bruteDistanceErrorM=distanceError,oldSevenNodeBranches=mapped,finalOmegaFloorCount=len(floor_nodes(omega)),
        finalKFloorCount=int(sum(k<=1e-10*(1+8*np.finfo(float).eps))),sevenResiduals=residual,
        outerTangentialColumnDiagnostics=distortion,maximumPrimitiveTotalEnthalpyToFreestream=float(ht[extrem]/hInf),
        maximumPrimitiveTotalEnthalpyNodeZeroBased=extrem,maximumPrimitiveTotalEnthalpyXYZM=P[extrem].tolist(),
        originalFaces=len(F),originalAreaM2=float(sum(area)),fieldRelativeErrors=errors,sectionHeatRelativeError=qerr,
        maximumSectionEndpointDistanceM=endpointError,midpoint75ToFinal149HeatRelativeL1Change=delta,
        originalFaceHeatPowerW=float(area@face['Heat_Flux']),
        screeningGates=dict(originalResiduals=sa['originalResidualThresholdPassed'],heatL1Below10Percent=delta<.1,
            last25HeatPowerBelow2Percent=abs(sa['last25HeatPowerRelativeChange'])<.02,yPlusP99Below1=sa['yPlusPercentiles']['p99']<=1),
        physicsQualified=False,newQualifiedArea=0,causalityOfGridVersusInitializationNotIsolated=True,
        localGeometryIsNotNativeResidualReplay=True,solverRuns=0,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'scripts/reviewSU2SSTLocalDiffusion.py',case/'engine_process.json',
            case/'sst_state_audit.json',case/'restart_flow.csv',case/'initial_flow.csv',case/'surface_flow_000075.vtk',base/'mesh/hybrid_volume.mat',base/'main01/fields.csv',base/'main01/section.csv']})
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['passed','seconds','finalOmegaFloorCount','midpoint75ToFinal149HeatRelativeL1Change','screeningGates','maximumPrimitiveTotalEnthalpyToFreestream']}))
if __name__=='__main__':main()
