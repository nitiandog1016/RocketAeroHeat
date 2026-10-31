"""独立从原生场审核SST能量、热流、y+和具名迭代；不将启动场判作物理解。"""
from pathlib import Path
import argparse,csv,json,re,sys
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2CheckpointContinuation import sha,read,config
from su2SSTInterface import wall_geometry,state,PHASE_OUTPUT
def vtk_tokens(path):
    """只读本项目要求的原生ASCII节点场，不按文件名猜物理含义。"""
    with path.open() as f:
        assert f.readline().startswith('# vtk');f.readline();assert f.readline().strip()=='ASCII'
        assert f.readline().strip()=='DATASET UNSTRUCTURED_GRID';words=iter(f.read().split())
    assert next(words)=='POINTS';n=int(next(words));next(words)
    xyz=np.fromiter((float(next(words)) for _ in range(3*n)),float).reshape(n,3)
    assert next(words)=='CELLS';nc=int(next(words));size=int(next(words))
    for _ in range(size):next(words)
    assert next(words)=='CELL_TYPES' and int(next(words))==nc
    for _ in range(nc):next(words)
    assert next(words)=='POINT_DATA' and int(next(words))==n;fields={}
    for tag in words:
        name=next(words);next(words)
        if tag=='SCALARS':
            width=int(next(words));assert next(words)=='LOOKUP_TABLE' and next(words)=='default'
        else:assert tag=='VECTORS';width=3
        fields[name]=np.fromiter((float(next(words)) for _ in range(n*width)),float).reshape(n,width)
    return xyz,fields


def audit(folder,actual_surface=False):
    folder=Path(folder).resolve();r=read(folder/'engine_process.json');c=config(folder/'case.cfg')
    assert r['status']=='terminal' and r['exitCode']==0 and not r['timedOut'] and not r['cancelled']
    assert all(r['inputsStable'].values()) and r['physicsModel']['model']=='SST-2003m'
    for name in ('case.cfg','mesh.su2','sst_contract.json','case_manifest.json'):
        assert sha(folder/name)==r['inputSHA256'][str(folder/name)]
    for name in ('restart_flow.csv','surface_flow.vtk','history.csv'):
        assert sha(folder/name)==r['outputSHA256'][name]
    extended=r['physicsModel'].get('outputStateContract')==PHASE_OUTPUT
    a,cols,nd,_,_,_=state(folder/'restart_flow.csv',True,extended)
    rho=a[:,cols.index('Density')];E=a[:,cols.index('Energy')]/rho
    uv=np.column_stack([a[:,cols.index('Momentum_'+v)]/rho for v in ['x','y','z'][:nd]])
    k=a[:,cols.index('Turb_Kin_Energy')];w=a[:,cols.index('Omega')]
    thermal=E-np.sum(uv**2,axis=1)/2-k;assert min(rho)>0 and min(thermal)>0 and min(k)>=0 and min(w)>0
    gamma=float(c['GAMMA_VALUE']);R=float(c['GAS_CONSTANT'])
    phase='saved-current-k';stageThermal=thermal
    initialization=read(folder/'sst_initialization.json') if (folder/'sst_initialization.json').exists() else {}
    if initialization.get('kind')=='same-mesh-SST-single-step-phase-diagnostic':
        # 单步的输入k是已知的上一个耦合阶段。普通稀疏保存的多步结果不能套用此分支。
        from prepareSU2SSTPhaseProbe import validate_source
        validate_source(initialization,folder,c)
        for path,h in initialization['sourceSHA256'].items():assert sha(Path(path))==h
        assert sha(folder/'initial_flow.csv')==r['inputSHA256'][str(folder/'initial_flow.csv')]
        initial,oldcols,oldnd,_,_,_=state(folder/'initial_flow.csv',True)
        assert oldcols==cols and oldnd==nd and np.array_equal(initial[:,:nd+1],a[:,:nd+1])
        stageThermal=thermal+k-initial[:,cols.index('Turb_Kin_Energy')]
        assert min(stageThermal)>0
        phase='updated-flow-with-input-k; single diagnostic update only'
    T=stageThermal*(gamma-1)/R;p=rho*R*T
    primitiveEOS=None
    if extended:
        assert c['WRT_RESTART_COMPACT']=='NO' and c['VOLUME_OUTPUT']=='(COORDINATES, SOLUTION, PRIMITIVE)'
        T=a[:,cols.index('Temperature')];p=a[:,cols.index('Pressure')]
        assert min(T)>0 and min(p)>0
        primitiveEOS=float(max(abs(p/(rho*R*T)-1)));assert primitiveEOS<2e-12
        phase=PHASE_OUTPUT
    g=wall_geometry(folder/'mesh.su2');xyz,f=vtk_tokens(folder/'surface_flow.vtk')
    points=np.column_stack([a[:,1:nd+1],np.zeros((len(a),3-nd))]);d,ids=cKDTree(points).query(xyz)
    assert max(d)<1e-6 and set(ids)==set(g['wallIDs']);order={v:i for i,v in enumerate(g['wallIDs'])};ix=np.array([order[j] for j in ids])
    eos=max(float(max(abs(p[ids]/f['Pressure'][:,0]-1))),float(max(abs(T[ids]/f['Temperature'][:,0]-1))))
    # CTurbSSTSolver.cpp:139 / 原生隐式更新的固定下限为1e-10；边界赋零后保存值可为该下限。
    # 仅接受原生零或准确下限，不把未知大k裁剪成零；EOS始终扣实际保存的k。
    assert eos<2e-6 and np.all((k[ids]==0)|(k[ids]==1e-10))
    def mu(t):return float(c['MU_REF'])*(t/float(c['MU_T_REF']))**1.5*(float(c['MU_T_REF'])+float(c['SUTHERLAND_CONSTANT']))/(t+float(c['SUTHERLAND_CONSTANT']))
    viscosity=mu(T[ids]);assert max(abs(viscosity/f['Laminar_Viscosity'][:,0]-1))<2e-6
    assert 'Eddy_Viscosity' in f
    eddyFraction=float(max(f['Eddy_Viscosity'][:,0]/viscosity));assert eddyFraction<1e-12
    Tw=float(re.fullmatch(r'\(\s*wall\s*,\s*([-+\deE.]+)\s*\)',c['MARKER_ISOTHERMAL'])[1])
    cp=gamma*R/(gamma-1)
    q=viscosity*cp/float(c['PRANDTL_LAM'])*(T[g['normalNeighbor'][ix]]-Tw)/g['thermalDistance'][ix]
    heatError=float(max(abs(q-f['Heat_Flux'][:,0])/np.maximum(abs(q),1)));assert heatError<2e-6
    highPrecisionHeatError=None
    if extended:
        nativeq=a[ids,cols.index('Heat_Flux')]
        highPrecisionHeatError=float(max(abs(q-nativeq)/np.maximum(abs(q),1)));assert highPrecisionHeatError<2e-6
    assert max(abs(q))>1,'必须有非零壁热'
    rhoi=float(c['FREESTREAM_PRESSURE'])/(R*float(c['FREESTREAM_TEMPERATURE']))
    speed=float(c['MACH_NUMBER'])*np.sqrt(gamma*R*float(c['FREESTREAM_TEMPERATURE']))
    shear=np.linalg.norm(f['Skin_Friction_Coefficient'],axis=1)*.5*rhoi*speed**2
    yp=g['nearestWallDistance'][ix]*np.sqrt(rho[ids]*shear)/viscosity
    ypError=float(max(abs(yp-f['Y_Plus'][:,0])/np.maximum(abs(yp),1e-20)));assert ypError<2e-6 and max(yp)>0
    # 壁ω由该迭代更新前流动状态施加，终态ρ/μ可有差别，记录差值而不裁剪。
    expectedOmega=60*viscosity/(rho[ids]*.075*g['nearestWallDistance'][ix]**2)
    omegaLag=float(max(abs(w[ids]/expectedOmega-1)))
    with (folder/'history.csv').open() as h:
        rows=list(csv.reader(h,skipinitialspace=True));keys=[v.strip() for v in rows[0]];hist=np.array(rows[1:],float)
    iterations=hist[:,keys.index('Inner_Iter')];assert np.array_equal(iterations,np.arange(len(hist)))
    step=int(iterations[-1]);suffix='_%06d'%step
    for stem,ext in [('restart_flow','csv'),('surface_flow','vtk')]:
        name=stem+suffix+'.'+ext;assert sha(folder/name)==r['outputSHA256'][name]==sha(folder/(stem+'.'+ext))
    residual=[float(hist[-1,keys.index(name)]) for name in ('rms[Rho]','rms[RhoE]','rms[k]','rms[w]')]
    hf=hist[:,keys.index('HF')];native=float(g['nodeArea'][ix]@f['Heat_Flux'][:,0])
    assert abs(native-hf[-1])/max(1,abs(native))<2e-6
    report=dict(interfacePassed=True,nodes=len(a),wallNodes=len(ids),dimension=nd,model=r['physicsModel'],
        thermalStatePhase=phase,diagnosticOnly=bool(initialization.get('diagnosticOnly',False)),
        fullPrimitiveEOSRelativeError=primitiveEOS,nativePrecisionHeatGradientError=highPrecisionHeatError,
        savedCurrentKTemperatureRangeK=[float(min(thermal*(gamma-1)/R)),float(max(thermal*(gamma-1)/R))],
        elapsedSeconds=r['elapsedSeconds'],checkpointLocalIteration=step,finalResiduals=residual,
        originalResidualThresholdPassed=all(v<-8 for v in residual) and step>=10,
        minimumDensity=float(min(rho)),minimumThermalInternalEnergy=float(min(thermal)),kRange=[float(min(k)),float(max(k))],
        omegaRange=[float(min(w)),float(max(w))],temperatureRangeK=[float(min(T)),float(max(T))],
        EOSRelativeError=eos,heatGradientRelativeError=heatError,yPlusRelativeError=ypError,
        yPlusPercentiles=dict(zip(['p50','p95','p99','max'],np.percentile(yp,[50,95,99,100]).tolist())),
        omegaBoundaryLagRelativeDifference=omegaLag,wallKRange=[float(min(k[ids])),float(max(k[ids]))],
        nativeWallKFloor=1e-10,wallEddyToLaminarMaximum=eddyFraction,
        last25HeatPowerRelativeChange=(float((hf[-1]-hf[-26])/max(abs(hf[-26]),1)) if len(hf)>=26 else None),
        nativeHeatPowerW=native,newQualifiedArea=0,physicsQualified=False,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__),ROOT/'scripts/su2SSTInterface.py',folder/'engine_process.json']})
    if actual_surface:
        assert nd==3
        from scipy.io import loadmat
        s=loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'),struct_as_record=False,squeeze_me=True)['s']
        V=s.V;F=s.F.astype(int)-1;area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
        assert len(F)==19752 and abs(sum(area)-1.9646292565715553)<1e-12
        _,vi=cKDTree(V).query(xyz);assert len(set(vi))==len(V)
        nodeq=np.zeros(len(V));nodeq[vi]=f['Heat_Flux'][:,0];face=nodeq[F].mean(axis=1)
        y=np.zeros(len(V));y[vi]=yp;above=float(sum(area[y[F].mean(axis=1)>1])/sum(area))
        checkpoints=sorted(p for p in folder.glob('surface_flow_*.vtk') if int(p.stem.rsplit('_',1)[1])<step)
        assert checkpoints or initialization.get('diagnosticOnly',False),'没有较早保存场供本段热流变化核对'
        midpath=checkpoints[-1] if checkpoints else None
        if initialization.get('diagnosticOnly',False):
            midpath=Path(initialization['sourceDirectory'])/('surface_flow_%06d.vtk'%initialization['sourceLocalIteration'])
            upstream=read(Path(initialization['sourceDirectory'])/'engine_process.json')
            assert sha(midpath)==upstream['outputSHA256'][midpath.name]
        else:assert sha(midpath)==r['outputSHA256'][midpath.name]
        midxyz,mid=vtk_tokens(midpath);assert np.array_equal(midxyz,xyz)
        nq=np.zeros(len(V));nq[vi]=mid['Heat_Flux'][:,0];mf=nq[F].mean(axis=1)
        l1=float(area@abs(face-mf)/(area@abs(mf)))
        report.update(originalFaces=len(F),originalAreaM2=float(sum(area)),originalFaceHeatPowerW=float(area@face),
                      heatComparisonSource=str(midpath),faceMeanYPlusAbove1AreaFraction=above)
        if initialization.get('diagnosticOnly',False):report['diagnosticSourceToOutputHeatRelativeL1Change']=l1
        else:report['withinRunHeatRelativeL1Change']=l1
        report['screeningGates']=dict(heatL1Below10Percent=l1<.1,last25HeatPowerBelow2Percent=report['last25HeatPowerRelativeChange'] is not None and abs(report['last25HeatPowerRelativeChange'])<.02,
                                     yPlusP99Below1=report['yPlusPercentiles']['p99']<=1)
        if initialization.get('diagnosticOnly',False):
            # 一个更新没有25步时间窗，不能用0变化自动通过稳定性筛选。
            report['last25HeatPowerRelativeChange']=None
            report['screeningGates']['last25HeatPowerBelow2Percent']=False
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',required=True);p.add_argument('--output',required=True);p.add_argument('--actual-surface',action='store_true');a=p.parse_args()
    target=Path(a.output);assert not target.exists();r=audit(a.folder,a.actual_surface)
    target.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(r,ensure_ascii=False,indent=2))
