"""模块职责：独立反算真实前体场、物理能量残差和停止规则，映射回完整原面责任。"""
from pathlib import Path
import csv,json,math,sys
import numpy as np
import runForebodyEulerV240 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseEulerRegionalVolume import read_tetra_mesh
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy
BASE=run.BASE;CASE=run.CASE;DOC=run.DOC


def field(path,xyz):
    a=np.loadtxt(str(path),delimiter=',',skiprows=1);a=a[np.argsort(a[:,0])]
    assert a.shape==(len(xyz),15) and np.array_equal(a[:,0],np.arange(len(xyz)))
    assert np.array_equal(a[:,1:4],xyz) and np.isfinite(a[:,:10]).all()
    rho=a[:,4];assert np.all(rho>0) and np.all(a[:,9]>0)
    v=a[:,5:8]/rho[:,None];energy=a[:,8]/rho-.5*np.sum(v*v,axis=1)-300000
    air=PhaseAwareViscousAir();T=inverse_energy(air,energy);p=rho*air.R*T;s=air.evaluate(p,T)
    return dict(raw=a,p=p,T=T,H=np.asarray(s['enthalpyJkg'])+.5*np.sum(v*v,axis=1),S=np.asarray(s['entropyJkgK']),v=v,
      M=np.linalg.norm(v,axis=1)/np.asarray(s['soundSpeedMps']))


def audit():
    out=BASE/'audit01';out.mkdir(exist_ok=False);receipt=read(BASE/'process_run/report.json');end=read(CASE/'native_completion.json')
    assert receipt['passed'] and receipt['inputsStable'] and not receipt['liveGroupAfterCleanup']
    assert end['normalReturn'] and end['forcedFinalNativeOutput'];n=int(end['updates'])
    xyz,_,markers=read_tetra_mesh(CASE/'mesh.su2');mapping=np.loadtxt(str(CASE/'original_faces.csv'),delimiter=',',skiprows=1)
    ids=mapping[:,0].astype(int)-1;tri=mapping[:,1:4].astype(int);area=mapping[:,4]
    assert len(ids)==4434 and len(np.unique(ids))==4434
    h=np.loadtxt(str(CASE/'physical_history.csv'),delimiter=',',skiprows=1)
    oh=np.loadtxt(str(CASE/'outlet_history.csv'),delimiter=',',skiprows=1)
    assert h.shape==(n,22) and oh.shape==(n,6) and np.array_equal(h[:,0],np.arange(1,n+1))
    assert np.all(h[:,1]==6111) and np.all(h[:,2:4]==0)
    assert np.array_equal(oh[:,0],h[:,0]) and np.all(oh[:,1]==807) and np.all(oh[:,4]==4035) and np.all(oh[:,2]>1)
    # 原生末预处理拒绝数只代表该阶段，不推断不存在任何隐藏子阶段事件。
    initial=field(CASE/'state_000000.csv',xyz)
    baseline=np.loadtxt(str(ROOT/'outputs/v232_validation/rocket01/state_000000.csv'),delimiter=',',skiprows=1)
    assert np.all(initial['raw'][:,4:9]==baseline[0,4:9])
    assert np.isnan(initial['raw'][:,10:]).all()
    air=PhaseAwareViscousAir();freestream=air.evaluate(1197.0279007383137,226.50908361133003)
    U=initial['v'][0];rho=float(freestream['densityKgm3']);H0=float(freestream['enthalpyJkg'])+.5*np.sum(U*U)
    outer=np.concatenate([markers['farfield'],markers['outlet']]);points=xyz[outer]
    normals=np.cross(points[:,1]-points[:,0],points[:,2]-points[:,0])/2
    mass_scale=float(np.sum(abs(normals@U))*rho);energy_scale=mass_scale*abs(H0)
    scales=np.array([mass_scale,energy_scale]);assert np.all(scales>0)
    assembly=np.abs(h[:,[9,13]]-h[:,[14,18]])/scales
    assert float(assembly.max())<1e-8
    steps=sorted(int(p.stem.split('_')[1]) for p in CASE.glob('state_*.csv'))
    assert steps==sorted(set([0,1]+list(range(50,n+1,50))+[n]))
    stats=[];last=None;rmsmax=0.;prev_p=None;pressure_changes={}
    for step in steps:
        value=initial if step==0 else field(CASE/('state_%06d.csv'%step),xyz);a=value['raw']
        if step:
            assert np.isfinite(a).all()
            residual=-a[:,10:15].copy();residual[:,4]-=300000*residual[:,0]
            rms=np.sqrt(np.mean(residual*residual,axis=0));error=float(np.max(abs(rms-h[step-1,4:9])/np.maximum(1,abs(h[step-1,4:9]))))
            rmsmax=max(rmsmax,error);assert error<2e-9
            # 接近抵消的总量用明确累加舍入界核对，不能用净值作分母。
            sums=np.array([math.fsum(residual[:,j]) for j in range(5)])
            mag=np.sum(abs(residual),axis=0);mag[4]=np.sum(abs(a[:,14])+300000*abs(a[:,10]))
            u=np.finfo(float).eps/2;k=4*len(a)+8;bound=k*u/(1-k*u)*mag
            assert np.all(abs(sums-h[step-1,9:14])<=bound)
        data=np.column_stack([ids+1,area]+[value[k][tri].mean(axis=1) for k in ('p','T','H','S')]+
          [value['v'][tri,j].mean(axis=1) for j in range(3)]+[value['M'][tri].mean(axis=1)])
        np.savetxt(str(out/('original_faces_%06d.csv'%step)),data,delimiter=',',fmt='%.17g',comments='',
          header='originalFaceID,areaM2,pPa,fluidTK,totalHJkg,entropyJkgK,Ux,Uy,Uz,Mach')
        if step and step%50==0:
            if prev_p is not None:pressure_changes[step]=float(np.sum(area*abs(data[:,2]-prev_p))/np.sum(area*abs(prev_p)))
            prev_p=data[:,2]
        stats.append(dict(update=step,sourceSHA256=sha(CASE/('state_%06d.csv'%step)),jointDomainAllNodes=True,
          range={k:[float(value[k].min()),float(value[k].max())] for k in ('p','T','H','S','M')}))
        last=value
    trend=np.loadtxt(str(CASE/'pressure_trend.csv'),delimiter=',',skiprows=1)
    for step,value in pressure_changes.items():
        native=float(trend[trend[:,0]==step,1][0]);assert abs(value-native)<2e-12,(step,value,native)
    decisions=np.atleast_2d(np.loadtxt(str(CASE/'benefit.csv'),delimiter=',',skiprows=1));stop=False
    for row in decisions:
        step=int(row[0]);ratio=h[step-50:step,4:9].mean(axis=0)/h[step-100:step-50,4:9].mean(axis=0)
        gm=float(np.exp(np.log(ratio).mean()));stop=bool(np.any(ratio>1.25) or (gm>=.95 and pressure_changes[step]>.01))
        assert np.max(abs(row[1:6]-ratio))<2e-12 and abs(row[6]-gm)<2e-12 and bool(row[8])==stop
    assert end['benefitStop']==stop
    # 正常强制写出的原生末CSV必须与完整末帧坐标/五守恒量逐位一致。
    matches=[]
    for p in sorted(CASE.glob('restart_flow*.csv')):
        data=np.loadtxt(str(p),delimiter=',',skiprows=1);data=data[np.argsort(data[:,0])]
        if np.array_equal(data[:,:9],last['raw'][:,:9]):matches.append(p.name)
    assert matches,'NO_NATIVE_FINAL_CSV_MATCH'
    p=last['p'][tri].mean(axis=1);hh=last['H'][tri].mean(axis=1)
    report=dict(passed=True,normalReturn=True,completedUpdates=n,benefitStop=stop,independentStopPredicatePassed=True,
      savedUpdates=steps,uniformInitialFiveBitEqualToV232=True,nativeFinalCSVBitEqual=matches,
      actualBoundaryCalls=int(oh[:,1].sum()),actualBoundaryFiveComponentsBitEqual=int(oh[:,4].sum()),
      minimumActualOutletNormalMach=float(oh[:,2].min()),maximumActualOutletNormalMach=float(oh[:,3].max()),
      lastPreprocessingRejectionMax=int(oh[:,5].max()),savedFrameChecks=stats,rmsIndependentMaxScaledDifference=rmsmax,
      physicalFinalLogRMS=np.log10(h[-1,4:9]).tolist(),physicalMinus8Passed=bool(np.all(np.log10(h[-1,4:9])<=-8)),
      massEnergyBoundaryAssemblyMaxScaled=float(assembly.max()),massEnergyReferenceScales=scales.tolist(),
      finalNetMassEnergyScaled=(abs(h[-1,[14,18]])/scales).tolist(),finalPressureAreaChange=pressure_changes[n],
      wallMomentumHistorySameStage=False,fullDomainMomentumBalanceVerified=False,
      originalFaces=19752,originalAreaM2=1.9646292565715553,regionFaces=len(ids),regionAreaM2=float(area.sum()),
      regionalPressurePaRange=[float(p.min()),float(p.max())],freestreamHJkg=H0,
      regionalHAbsoluteRelativeAreaMean=float(np.sum(area*abs(hh/H0-1))/sum(area)),regionalHAbsoluteRelativeMax=float(max(abs(hh/H0-1))),
      pressureStabilityPassed=pressure_changes[n]<=.01,spatialFarfieldQualified=False,edgeStreamtubeQualified=False,
      physicalQualification=False,newQualifiedAreaM2=0,wallHeatFluxComputed=False,sourceRunSHA256=sha(BASE/'process_run/report.json'))
    save(out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k!='savedFrameChecks'}))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':audit()
    else:
        run.guard.BASE=BASE;run.guard.CASE=BASE;run.guard.DOC=DOC
        size=run.guard.size_bytes;run.guard.size_bytes=lambda _:size(BASE)
        inputs=[Path(__file__).resolve(),ROOT/'scripts/runForebodyEulerV240.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',DOC/'REGIONAL_PLAN.md',
          ROOT/'scripts/phaseEulerRegionalVolume.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',
          ROOT/'data/thermo/phase_aware_dry_air_v1.json',BASE/'process_run/report.json',ROOT/'outputs/v232_validation/rocket01/state_000000.csv']+[p for p in CASE.iterdir() if p.is_file()]
        raise SystemExit(0 if run.guard.supervise('audit',[sys.executable,str(Path(__file__).resolve()),'task'],120,256,inputs,BASE,'audit') else 1)
