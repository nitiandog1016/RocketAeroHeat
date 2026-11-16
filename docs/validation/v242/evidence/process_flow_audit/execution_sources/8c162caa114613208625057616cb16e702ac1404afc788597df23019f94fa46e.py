"""模块职责：独立核验新网格真实发展、物理残差/守恒账及原面映射，不放宽门槛。"""
from pathlib import Path
import math,sys
import numpy as np
import runForebodyEulerV242 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseEulerRegionalVolume import read_tetra_mesh
from phaseAwareViscousAir import PhaseAwareViscousAir
from auditForebodyEulerV240 import field
BASE=run.BASE;CASE=run.CASE;DOC=run.DOC
CLOSED=False


def audit():
    """模块一：原始完整场独立反算，再核对原生历史和停止规则。"""
    out=BASE/'flow_audit01';out.mkdir(exist_ok=False)
    receipt=read(BASE/'process_flow_run/report.json');end=read(CASE/'native_completion.json');mesh=read(BASE/'mesh_audit01/report.json')
    assert receipt['passed'] and receipt['inputsStable'] and not receipt['liveGroupAfterCleanup']
    assert end['normalReturn'] and end['forcedFinalNativeOutput'];n=int(end['updates'])
    xyz,_,markers=read_tetra_mesh(CASE/'mesh.su2');mapping=np.loadtxt(str(CASE/'original_faces.csv'),delimiter=',',skiprows=1)
    ids=mapping[:,0].astype(int)-1;tri=mapping[:,1:4].astype(int);area=mapping[:,4]
    assert len(ids)==4434 and len(np.unique(ids))==4434 and np.array_equal(tri,markers['wall'])
    h=np.loadtxt(str(CASE/'physical_history.csv'),delimiter=',',skiprows=1)
    oh=np.loadtxt(str(CASE/'outlet_history.csv'),delimiter=',',skiprows=1)
    assert h.shape==(n,22) and oh.shape==(n,6) and np.array_equal(h[:,0],np.arange(1,n+1))
    assert np.all(h[:,1]==mesh['nodes']) and np.all(h[:,2:4]==0)
    assert np.array_equal(oh[:,0],h[:,0]) and np.all(oh[:,1]==mesh['outletNodes']) and np.all(oh[:,4]==5*mesh['outletNodes']) and np.all(oh[:,2]>1)
    initial=field(CASE/'state_000000.csv',xyz)
    baseline=np.loadtxt(str(ROOT/'outputs/v240_validation/run01/state_000000.csv'),delimiter=',',skiprows=1)
    assert np.all(initial['raw'][:,4:9]==baseline[0,4:9]) and np.isnan(initial['raw'][:,10:]).all()
    air=PhaseAwareViscousAir();free=air.evaluate(1197.0279007383137,226.50908361133003)
    velocity=initial['v'][0];rho=float(free['densityKgm3']);H0=float(free['enthalpyJkg'])+.5*np.sum(velocity**2)
    points=xyz[np.concatenate([markers['farfield'],markers['outlet']])]
    normals=np.cross(points[:,1]-points[:,0],points[:,2]-points[:,0])/2
    mass_scale=float(np.sum(abs(normals@velocity))*rho);scales=np.array([mass_scale,mass_scale*abs(H0)])
    assembly=abs(h[:,[9,13]]-h[:,[14,18]])/scales;assert assembly.max()<1e-8
    steps=sorted(int(p.stem.split('_')[1]) for p in CASE.glob('state_*.csv'))
    assert steps==sorted(set([0,1]+list(range(50,n+1,50))+[n]))
    stats=[];rmsmax=0.;previous=None;changes={};last=None
    for step in steps:
        value=initial if step==0 else field(CASE/('state_%06d.csv'%step),xyz);a=value['raw']
        if step:
            assert np.isfinite(a).all()
            residual=-a[:,10:15].copy();residual[:,4]-=300000*residual[:,0]
            rms=np.sqrt(np.mean(residual**2,axis=0));err=float(np.max(abs(rms-h[step-1,4:9])/np.maximum(1,abs(h[step-1,4:9]))))
            rmsmax=max(rmsmax,err);assert err<2e-9
            sums=np.array([math.fsum(residual[:,j]) for j in range(5)])
            magnitude=np.sum(abs(residual),axis=0);magnitude[4]=np.sum(abs(a[:,14])+300000*abs(a[:,10]))
            eps=np.finfo(float).eps/2;k=4*len(a)+8;bound=k*eps/(1-k*eps)*magnitude
            assert np.all(abs(sums-h[step-1,9:14])<=bound)
        data=np.column_stack([ids+1,area]+[value[k][tri].mean(axis=1) for k in ('p','T','H','S')]+
          [value['v'][tri,j].mean(axis=1) for j in range(3)]+[value['M'][tri].mean(axis=1)])
        np.savetxt(str(out/('original_faces_%06d.csv'%step)),data,delimiter=',',fmt='%.17g',comments='',
          header='originalFaceID,areaM2,pPa,fluidTK,totalHJkg,entropyJkgK,Ux,Uy,Uz,Mach')
        if step and step%50==0:
            if previous is not None:changes[step]=float(np.sum(area*abs(data[:,2]-previous))/np.sum(area*abs(previous)))
            previous=data[:,2]
        stats.append(dict(update=step,sourceSHA256=sha(CASE/('state_%06d.csv'%step)),jointDomainAllNodes=True,
          range={k:[float(value[k].min()),float(value[k].max())] for k in ('p','T','H','S','M')}));last=value
    trend=np.loadtxt(str(CASE/'pressure_trend.csv'),delimiter=',',skiprows=1);trenderr=0.
    for step,value in changes.items():
        native=float(trend[trend[:,0]==step,1][0]);err=abs(native-value);trenderr=max(trenderr,err);assert err<2e-12
    decisions=np.atleast_2d(np.loadtxt(str(CASE/'benefit.csv'),delimiter=',',skiprows=1));stop=False;decision_rows=[]
    for row in decisions:
        step=int(row[0]);ratio=h[step-50:step,4:9].mean(axis=0)/h[step-100:step-50,4:9].mean(axis=0)
        gm=float(np.exp(np.log(ratio).mean()));stop=bool(np.any(ratio>1.25) or (gm>=.95 and changes[step]>.01))
        assert np.max(abs(row[1:6]-ratio))<2e-12 and abs(row[6]-gm)<2e-12 and bool(row[8])==stop
        decision_rows.append(dict(update=step,ratio=ratio.tolist(),geometricMean=gm,pressureChange=changes[step],stop=stop))
    assert end['benefitStop']==stop
    matches=[]
    for path in sorted(CASE.glob('restart_flow*.csv')):
        data=np.loadtxt(str(path),delimiter=',',skiprows=1);data=data[np.argsort(data[:,0])]
        if np.array_equal(data[:,:9],last['raw'][:,:9]):matches.append(path.name)
    assert matches,'NO_NATIVE_FINAL_CSV_MATCH'
    p=last['p'][tri].mean(axis=1);hh=last['H'][tri].mean(axis=1);logrms=np.log10(h[-1,4:9]);net=abs(h[-1,[14,18]])/scales
    old=read(ROOT/'outputs/v240_validation/audit02/report.json')
    report=dict(passed=True,normalReturn=True,completedUpdates=n,benefitStop=stop,independentStopPredicatePassed=True,
      stopReason='benefit-stop' if stop else 'registered-400-update-budget',benefitChecks=decision_rows,
      savedUpdates=steps,uniformInitialFiveBitEqualToV240=True,nativeFinalCSVBitEqual=matches,
      actualBoundaryCalls=int(oh[:,1].sum()),actualBoundaryFiveComponentsBitEqual=int(oh[:,4].sum()),
      minimumActualOutletNormalMach=float(oh[:,2].min()),maximumActualOutletNormalMach=float(oh[:,3].max()),
      lastPreprocessingRejectionMax=int(oh[:,5].max()),lastPreprocessingRejectionFinal=int(oh[-1,5]),
      savedFrameChecks=stats,rmsIndependentMaxScaledDifference=rmsmax,pressureTrendMaxAbsoluteDifference=trenderr,
      physicalFinalLogRMS=logrms.tolist(),physicalMinus8ByEquation=(logrms<=-8).tolist(),physicalMinus8Passed=bool(np.all(logrms<=-8)),
      massEnergyBoundaryAssemblyMaxScaled=float(assembly.max()),massEnergyReferenceScales=scales.tolist(),
      finalNetMassEnergyScaled=net.tolist(),netMassEnergyPassed=bool(np.all(net<=1e-8)),finalPressureAreaChange=changes[n],
      originalFaces=19752,originalAreaM2=1.9646292565715553,regionFaces=len(ids),regionAreaM2=float(area.sum()),
      regionalPressurePaRange=[float(p.min()),float(p.max())],freestreamHJkg=H0,
      regionalHAbsoluteRelativeAreaMean=float(np.sum(area*abs(hh/H0-1))/sum(area)),regionalHAbsoluteRelativeMax=float(max(abs(hh/H0-1))),
      pressureStabilityPassed=changes[n]<=.01,spatialFarfieldQualified=False,edgeStreamtubeQualified=False,
      wallMomentumHistorySameStage=False,fullDomainMomentumBalanceVerified=False,
      physicalQualification=False,newQualifiedAreaM2=0,wallHeatFluxComputed=False,
      v240Comparison=dict(finalLogRMS=old['physicalFinalLogRMS'],netMassEnergy=old['finalNetMassEnergyScaled'],
        pressureChange=old['finalPressureAreaChange'],HMean=old['regionalHAbsoluteRelativeAreaMean'],HMax=old['regionalHAbsoluteRelativeMax']),
      sourceRunSHA256=sha(BASE/'process_flow_run/report.json'))
    save(out/'report.json',report)
    print({k:report[k] for k in ('completedUpdates','physicalFinalLogRMS','finalNetMassEnergyScaled','finalPressureAreaChange','regionalHAbsoluteRelativeAreaMean','regionalHAbsoluteRelativeMax','lastPreprocessingRejectionMax')})


def main():
    if CLOSED:raise SystemExit('v242 flow audit closed; use NEXT_DECISION')
    run.guard.BASE=BASE;run.guard.CASE=CASE;run.guard.DOC=DOC;size=run.guard.size_bytes
    run.guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),ROOT/'scripts/runForebodyEulerV242.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',DOC/'CFD_PLAN.md',
      ROOT/'scripts/auditForebodyEulerV240.py',ROOT/'scripts/runForebodyEulerV240.py',ROOT/'scripts/prepareSU2PhaseAirV220.py',
      ROOT/'scripts/phaseEulerRegionalVolume.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',
      ROOT/'data/thermo/phase_aware_dry_air_v1.json',BASE/'process_flow_run/report.json',BASE/'mesh_audit01/report.json',
      ROOT/'outputs/v240_validation/run01/state_000000.csv',ROOT/'outputs/v240_validation/audit02/report.json']+[p for p in CASE.iterdir() if p.is_file()]
    return run.guard.supervise('flow_audit',[sys.executable,str(Path(__file__).resolve()),'task'],120,256,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v242 flow audit closed; use NEXT_DECISION')
    if len(sys.argv)>1 and sys.argv[1]=='task':audit()
    else:raise SystemExit(0 if main() else 1)
