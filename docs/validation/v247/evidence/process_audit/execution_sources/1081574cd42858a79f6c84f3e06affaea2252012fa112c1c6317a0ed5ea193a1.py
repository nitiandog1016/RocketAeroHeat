"""模块职责：审核真实发展、原面结果和停止规则，独立核对完整终态空间残差的物理能量。"""
from pathlib import Path
import csv,subprocess,sys
import numpy as np
from scipy.optimize import brentq
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy
from singleDomainFrame import read_frame
from fasDevelopmentFieldsV247 import state,residual
import runFASDevelopmentV247 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;CASE=run.CASE;OUT=BASE/'audit01';CLOSED=False
GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz'


def rows(path):
    with Path(path).open(encoding='utf8') as f:return list(csv.DictReader(f))


def audit():
    OUT.mkdir(exist_ok=False);process=read(BASE/'process_rocket/report.json');assert process['inputsStable'] and not process['liveGroupAfterCleanup']
    success=process['passed'] and (CASE/'native_completion.json').exists();steps=rows(CASE/'steps.csv');completed=len(steps)
    assert [int(r['update']) for r in steps]==list(range(1,completed+1))
    air=PhaseAwareViscousAir();errors={};records=[]
    def check(name,a,b,tol=2e-9):
        a=np.asarray(a);b=np.asarray(b);assert np.isfinite(a).all() and np.isfinite(b).all(),name
        e=float(np.max(abs(a-b)/np.maximum(1,abs(b))));errors[name]=max(errors.get(name,0),e);assert e<tol,(name,e)
    g={k:v for k,v in np.load(str(GEOMETRY)).items()};tri=g['triangles'];wallids=g['wallIDs'];area=g['faceArea'];N=335726
    assert sha(GEOMETRY)==read(ROOT/'outputs/v222_validation/geometry01/report.json')['cacheSHA256']
    saved=sorted(int(p.stem.split('_')[1]) for p in CASE.glob('state_*.bin'))
    stage=rows(CASE/'stages.csv');complete_stage=[r for r in stage if int(r['update'])<=completed]
    assert len(complete_stage)==16*completed
    template=[(r['system'],r['stage']) for r in rows(run.PREV/'control_observe01/stages.csv')[:16]]
    for step in range(1,completed+1):assert [(r['system'],r['stage']) for r in complete_stage if int(r['update'])==step]==template
    after={(int(r['update']),int(r['system'])):r for r in complete_stage if r['stage']=='after_time'}
    check('initialElevenBitHeader',read(CASE/'initial_check.json')['initialElevenColumnBitMismatches'],0)
    states_checked=0;brent_samples=0
    for step in saved:
        raw=state(CASE/('state_%06d.bin'%step),initial=step==0);a=raw[np.argsort(raw[:,0])];states_checked+=N
        assert np.array_equal(a[:,1:4],g['points'])
        if step in (0,1):
            ref=read_frame(run.PREV/('control_observe01/state_%06d.bin'%step),N,step==0)['values'];ref=ref[np.argsort(ref[:,0])]
            assert np.array_equal(a[:,:11],ref[:,:11]),('initialTrajectoryBitMismatch',step)
        vel=a[:,6:9]/a[:,5,None];ep=a[:,9]/a[:,5]-.5*np.sum(vel*vel,axis=1)-300000;t=inverse_energy(air,ep)
        phys=air.evaluate(a[:,5]*air.R*t,t);check('currentT',a[:,11],t);check('currentP',a[:,15],phys['pressurePa']);check('currentU',a[:,12:15],vel)
        check('physicalH',a[:,17]-300000,np.asarray(phys['enthalpyJkg'])+.5*np.sum(vel*vel,axis=1))
        transport=air.interpolated_transport(a[:,15],a[:,11]);mu=np.asarray(transport['dynamicViscosityPas']);k=np.asarray(transport['conductivityWmK'])
        check('mu',a[:,19],mu);check('k',a[:,21],k);chi=a[:,10]/(mu/a[:,5]);mut=np.where(a[:,10]<0,0,a[:,5]*a[:,10]*chi**3/(chi**3+7.1**3));check('latestSAmut',a[:,23],mut)
        selected=np.unique(np.r_[np.linspace(0,N-1,12).astype(int),np.argmin(ep),np.argmax(ep)])
        tb=[brentq(lambda x:float(air.thermal(x)[1]-air.R*x)-ep[i],60,2000,xtol=1e-10) for i in selected];check('BrentCurrentT',t[selected],tb);brent_samples+=len(selected)
        q=np.full(N,np.nan);q[wallids]=k[wallids]*(a[g['normalNeighbor'],11]-300)/g['thermalDistance']
        f=np.column_stack((np.arange(1,19753),area,a[tri,15].mean(axis=1),q[tri].mean(axis=1),a[tri,11].mean(axis=1),np.full(19752,300.)))
        np.savetxt(str(OUT/('original_faces_%06d.csv'%step)),f,delimiter=',',header='faceID,areaM2,pressurePa,heatFluxIntoWallWm2,wallFluidTemperatureK,givenWallTemperatureK',comments='',fmt='%.17g')
        records.append(dict(update=step,pRangePa=[float(min(f[:,2])),float(max(f[:,2]))],qRangeWm2=[float(min(f[:,3])),float(max(f[:,3]))],wallFluidRangeK=[float(min(f[:,4])),float(max(f[:,4]))],powerIntoWallW=float(area@f[:,3])))
        if success and step==completed:
            native=np.loadtxt(str(CASE/'restart_flow.csv'),delimiter=',',skiprows=1,usecols=range(10));native=native[np.argsort(native[:,0])]
            assert np.array_equal(native[:,:4],a[:,:4]) and np.array_equal(native[:,4:10],a[:,5:11])
    # 50步原壁小表独立映射与面积变化：不能用壁节点数替代原面分母。
    face_history={};previous=None;changes={}
    for path in sorted(CASE.glob('wall_*.csv')):
        step=int(path.stem.split('_')[1]);w=np.loadtxt(str(path),delimiter=',',skiprows=1);assert len(w)==9878
        w=w[np.argsort(w[:,0])];assert np.array_equal(w[:,0],np.sort(wallids));assert np.isfinite(w[:,:4]).all()
        check('wallGradientFormula',w[:,3],w[:,9]*(w[:,8]-300)/w[:,7]);check('wallThermalDistance',w[:,7],g['thermalDistance'][np.argsort(wallids)])
        values=np.full((N,3),np.nan);values[w[:,0].astype(int)]=w[:,[1,3,2]];f=np.column_stack((np.arange(1,19753),area,values[tri].mean(axis=1),np.full(19752,300.)))
        if step in saved:check('savedAndWallFaces',f,np.loadtxt(str(OUT/('original_faces_%06d.csv'%step)),delimiter=',',skiprows=1))
        else:np.savetxt(str(OUT/('original_faces_%06d.csv'%step)),f,delimiter=',',header='faceID,areaM2,pressurePa,heatFluxIntoWallWm2,wallFluidTemperatureK,givenWallTemperatureK',comments='',fmt='%.17g')
        face_history[step]=f
        if step==0 or step%50==0:
            changes[step]=np.zeros(2) if previous is None else np.sum(area[:,None]*abs(f[:,2:4]-previous[:,2:4]),axis=0)/np.maximum(1,np.sum(area[:,None]*abs(previous[:,2:4]),axis=0))
            previous=f
    for r in rows(CASE/'surface_trend.csv'):
        step=int(r['update'])
        if step in changes:check('originalAreaTrend',[float(r['pressureAreaL1Change']),float(r['heatAreaL1Change'])],changes[step])
    decisions=rows(CASE/'benefit_decisions.csv');stop_expected=False
    if completed>=200:
        assert len(decisions)==2;log=np.array([[np.log10(max(1e-300,float(after[(i,0 if j<5 else 2)]['rms%d'%j]))) for j in range(6)] for i in range(1,completed+1)])
        goods=[];worses=[]
        for row,end in zip(decisions,(150,200)):
            ratio=10**(log[end-50:end].mean(axis=0)-log[end-100:end-50].mean(axis=0));gm=float(np.exp(np.mean(np.log(ratio))))
            old,now=changes[end-50],changes[end];notworse=bool(np.all(ratio<=1.25) and np.all(now<=1.25*np.maximum(old,[.01,.05])))
            area_good=bool(np.all((now<=[.01,.05])|(now<=.95*old)));allsmall=bool(np.all(log[end-50:end].mean(axis=0)<=-8))
            good=notworse and area_good and (gm<=.95 or allsmall);goods.append(good);worses.append(not notworse)
            check('benefitRatios',[float(row['ratio%d'%j]) for j in range(6)],ratio);check('benefitGM',float(row['geometricMeanRatio']),gm)
            assert bool(int(row['noSignificantDeterioration']))==notworse and bool(int(row['benefit']))==good
        stop_expected=any(worses) or not any(goods)
        reason='significant-deterioration' if any(worses) else 'two-windows-no-benefit' if not any(goods) else 'maximum-updates'
        if success:assert read(CASE/'flow_completion.json')['stopReason']==reason
    # 实际拒绝独立物性与逐位恢复核对；old并不等于延拓前状态，不能反推出真实粗层修正。
    rejection=[]
    initial=state(CASE/'state_000000.bin',initial=True);initial=initial[np.argsort(initial[:,0])]
    for event in rows(CASE/'rejections.csv'):
        trial=np.array([float(event['trial%d'%j]) for j in range(5)]);old=np.array([float(event['old%d'%j]) for j in range(5)]);returned=np.array([float(event['returned%d'%j]) for j in range(5)])
        assert np.array_equal(old.view(np.uint64),returned.view(np.uint64)) and int(event['recovered'])==1
        def energy(u):return u[4]/u[0]-.5*float(np.sum((u[1:4]/u[0])**2))-300000
        et,eo=energy(trial),energy(old);low=float(air.thermal(60)[1]-air.R*60)
        oldT=brentq(lambda x:float(air.thermal(x)[1]-air.R*x)-eo,60,2000,xtol=1e-10);oldp=old[0]*air.R*oldT;air.evaluate(oldp,oldT)
        assert et<low and eo>=low
        gid=int(event['globalID']);idx=int(event['localID']);assert int(state(CASE/'state_000000.bin',initial=True)[idx,0])==gid
        rejection.append(dict(update=int(event['update']),lastObservedPhase=event['phase'],level=int(event['level']),globalID=gid,localID=idx,
          xyzM=g['points'][gid].tolist(),wallDistanceM=float(initial[gid,24]),nodeVolumeM3=float(initial[gid,4]),
          trialPhysicalEnergyJkg=et,oldPhysicalEnergyJkg=eo,minimumPhysicalEnergyJkg=low,trialBelowMinimumJkg=low-et,
          oldTemperatureK=oldT,oldPressurePa=oldp,returnedFiveStateBitEqualOld=True,nativeRecoverySucceeded=True,
          observerStoppedByContract=True,oldIsNotProlongationInput=True,actualFinePreCorrectionAndCorrectionNotSaved=True))
    save(OUT/'actual_rejection.json',dict(passed=True,events=rejection,uniqueRootCauseProved=False))
    terminal=None
    if success:
        done=read(CASE/'native_completion.json');assert done['updates']==completed and completed in saved
        summary=read(CASE/'terminal_summary.json');assert summary['completedUpdates']==completed and summary['stateAndSAmutBitChanges']==0 and summary['additionalTimeUpdates']==0
        raw,r=residual(CASE/'terminal_spatial.bin');ext=r.astype(np.longdouble);rms=np.sqrt(np.mean(ext*ext,axis=0)).astype(float);summed=np.sum(ext,axis=0).astype(float);ab=np.sum(abs(ext),axis=0).astype(float)
        check('terminalSixRMS',summary['physicalRMS'],rms);check('terminalSixSum',summary['physicalSum'],summed);check('terminalSixSumAbs',summary['physicalSumAbs'],ab)
        terminal=dict(physicalLogRMS=np.log10(np.maximum(rms,1e-300)).tolist(),allSixMinus8=bool(np.all(rms<=1e-8)),
          netResidualToL1Ratio=(abs(summed)/np.maximum(ab,1e-300)).tolist(),ratioIsNotBoundaryFluxBalance=True,nativeSevenFieldBitIdentity=True,noTimeAdvance=True)
        for rstage in rows(CASE/'terminal_stages.csv'):
            if rstage['stage']=='after_strong':
                sysid=int(rstage['system']);cols=list(range(5)) if sysid==0 else [5];check('terminalStageRMS',[float(rstage['rms%d'%j]) for j in cols],rms[cols])
        assert read(CASE/'flow_completion.json')['transfersObserved']==9*completed
    tests=subprocess.run([sys.executable,str(ROOT/'tests/testFASDevelopmentFramesV247.py')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (OUT/'reader_test.txt').write_bytes(tests.stdout);assert tests.returncode==0
    report=dict(passed=True,developmentSucceeded=success,completedUpdates=completed,completeSavedUpdates=saved,originalFaceUpdates=sorted(face_history),
      originalFaces=19752,originalAreaM2=float(sum(area)),independentSavedStateChecks=states_checked,brentCrosschecks=brent_samples,
      initialZeroAndOneConservedBitEqual=True,nativeFinalSixStateBitEqual=success,stageRows=len(stage),completeStageRows=len(complete_stage),
      rejectionRecords=len(rejection),actualRejection=rejection,benefitStopExpected=stop_expected,benefitDecisionReached=completed>=200,terminal=terminal,frames=records,errorMetrics=errors,
      physicalQualification=False,newQualifiedAreaM2=0,wholeSurfaceGoalComplete=False,spaceOuterNearwallExperimentalPassed=False,fullSteadyBalanceVerified=False)
    save(OUT/'report.json',report);print({k:report[k] for k in ('passed','developmentSucceeded','completedUpdates','terminal')})


def main():
    if CLOSED:raise SystemExit('v247 audit closed; follow NEXT_DECISION')
    g=run.setup_guard();locks=[Path(__file__).resolve(),Path(run.__file__),DOC/'PLAN.md',GEOMETRY,ROOT/'scripts/fasDevelopmentFieldsV247.py',
      ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',ROOT/'scripts/singleDomainFrame.py',
      ROOT/'data/thermo/phase_aware_dry_air_v1.json',ROOT/'tests/testFASDevelopmentFramesV247.py',BASE/'process_rocket/report.json']
    locks += [p for p in CASE.iterdir() if p.is_file()]+[run.PREV/'control_observe01/state_000000.bin',run.PREV/'control_observe01/state_000001.bin']
    return g.supervise('audit',[sys.executable,str(Path(__file__).resolve()),'task'],120,1536,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v247 audit closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':audit()
    else:raise SystemExit(0 if main() else 1)
