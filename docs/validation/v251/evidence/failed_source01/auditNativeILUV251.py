"""模块职责：独立核对原生NK阶段、实际恢复方向、保存场及原面研究结果；不授热资格。"""
from pathlib import Path
import csv
import sys
import numpy as np
from scipy.optimize import brentq
import runNativeILUV251 as run
from terminalCapacityV249 import supervise
from singleDomainFrame import read_frame
from fasRecoveryFieldsV248 import binary,residual
from phaseAwareViscousAir import PhaseAwareViscousAir
from phaseAwareDryAir import pressure_ceiling
from su2PhaseAirFields import inverse_energy
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;CASE=run.CASE;CLOSED=False
GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz'


def rows(path):
    with Path(path).open(encoding='utf8') as f:return list(csv.DictReader(f))


def audit():
    out=BASE/'audit01';out.mkdir(exist_ok=False);errors={}
    def check(key,a,b,tol=2e-9):
        a=np.asarray(a);b=np.asarray(b);assert np.isfinite(a).all() and np.isfinite(b).all(),key
        err=float(np.max(abs(a-b)/np.maximum(1,abs(b))));errors[key]=err;assert err<=tol,(key,err)
    process=read(BASE/'process_rocket/report.json');capacity=read(BASE/'process_rocket/terminal_capacity.json')
    assert process['status']=='terminal' and process['inputsStable'] and not process['liveGroupAfterCleanup']
    assert capacity['passed']
    steps=rows(CASE/'steps.csv');completed=len(steps);assert [int(r['update']) for r in steps]==list(range(1,completed+1))
    normal=(CASE/'native_completion.json').exists()
    if normal:
        done=read(CASE/'native_completion.json');assert process['passed'] and done['updates']==completed and done['actualNewtonFactory']
    else:
        done=None;assert not process['passed']
        assert (CASE/'failure.txt').exists() or (CASE/'caught_failure.txt').exists()
    # 模块一：复核可移除观察代码和实际阶段，不能将最后Jv的RHS当基准趋势。
    original=(BASE/'stage01/CNewtonIntegration_original.cpp').read_text(encoding='utf8')
    observed=(BASE/'stage01/CNewtonIntegration_observed.cpp').read_text(encoding='utf8')
    observed=observed.replace((run.GAS/'source/SU2_CFD/include/integration/CNewtonIntegration.hpp').as_posix(),'../../include/integration/CNewtonIntegration.hpp')
    observed=observed.replace((run.GAS/'source/Common/include/linear_algebra/CMatrixVectorProduct.hpp').as_posix(),'../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp')
    stripped='\n'.join(x.strip() for x in observed.splitlines() if x.strip() and 'rah_newton_observer::event' not in x and 'NewtonHooksV250.hpp' not in x)
    assert stripped=='\n'.join(x.strip() for x in original.splitlines() if x.strip())
    phases=rows(CASE/'newton_phases.csv');contexts=rows(CASE/'evaluation_context.csv');stages=rows(CASE/'stages.csv')
    regular=[c for c in contexts if c['role']!='terminal'];assert len(regular)==len(stages)
    trend={};flow_patterns={};sa_patterns={}
    for c,s in zip(regular,stages):
        assert (c['update'],c['system'],c['stage'])==(s['update'],s['system'],s['stage'])
        step=int(c['update']);system=int(c['system']);label=c['stage']
        if system==0:flow_patterns.setdefault(step,[]).append(label)
        else:sa_patterns.setdefault(step,[]).append(label)
        if system==0 and c['role']=='baseline' and c['evaluation']=='DEFAULT' and label=='after_strong':
            assert step not in trend;trend[step]=[float(s['rms%d'%j]) for j in range(5)]+[None]
        if system==2 and label=='after_time':trend[step][5]=float(s['rms5'])
    pattern=['after_sources','after_weak','before_strong','after_strong']
    for step in range(1,completed+1):
        assert len(flow_patterns[step])%4==0 and flow_patterns[step]==pattern*(len(flow_patterns[step])//4)
        assert sa_patterns[step]==pattern+['before_time','after_time'] and all(x is not None for x in trend[step])
        tags=[p['event'] for p in phases if int(p['update'])==step]
        assert tags[0]=='default_begin' and tags[-1]=='flow_complete'
        assert tags.count('default_begin')==tags.count('default_end')==tags.count('update_begin')==tags.count('update_end')==1
        assert tags.count('jv_begin')==tags.count('jv_end')==tags.count('explicit_begin')==tags.count('explicit_end')
        assert len(flow_patterns[step])//4==1+tags.count('jv_begin')
    begins=[p for p in phases if p['event']=='jv_begin'];ends=[p for p in phases if p['event']=='jv_end']
    assert [int(p['productCount']) for p in begins]==list(range(1,len(begins)+1))
    if normal:assert len(begins)==len(ends)==done['matrixFreeProducts']
    # 原生线性历史只证明引擎报告，不替代独立Jv/线性残差检验；达到上限与达到容差分开。
    history=[{k.strip().strip('"'):v.strip() for k,v in r.items()} for r in rows(CASE/'history.csv')]
    assert len(history) in (completed,completed-1)
    if len(history)!=completed:assert not normal and 'SA_NATIVE_LINEAR_TOLERANCE' in (CASE/'failure.txt').read_text(encoding='utf8')
    cap_fail=[int(r['Inner_Iter'])+1 for r in history if int(r['Linear_Solver_Iterations'])==20 and float(r['LinSolRes'])>-5]
    linear=dict(nativeReportedOnly=True,iterationLimit=20,reportedTargetLog10=-5,
      capReachedAboveReportedTolerance=len(cap_fail),firstCapFailureUpdate=cap_fail[0] if cap_fail else None,
      worstReportedLog10Residual=max([float(r['LinSolRes']) for r in history] or [-300]),
      lastReportedLog10Residual=float(history[-1]['LinSolRes']) if history else None,independentLinearAccuracyProved=False)
    np.savetxt(str(out/'baseline_trends.csv'),[[s]+trend[s] for s in range(1,completed+1)],delimiter=',',comments='',
      header='update,Rmass,Rmx,Rmy,Rmz,REphysical,RSAafterTime',fmt='%.17g')
    # 模块二：恢复事件用实际方向/倍率重建；合法旧态和越域试探独立核对。
    events=rows(CASE/'rejections.csv');domains=rows(CASE/'domain_context.csv');assert len(events)==len(domains)
    air=PhaseAwareViscousAir();elow=float(air.thermal(60)[1]-air.R*60);ehigh=float(air.thermal(2000)[1]-air.R*2000)
    rejection_data=[];roles={};trial_low=0;pressure_rejections=0
    for n,(e,d) in enumerate(zip(events,domains)):
        assert (e['update'],e['localID'],e['globalID'])==(d['update'],d['localID'],d['globalID'])
        u={word:np.array([float(e[word+str(j)]) for j in range(5)]) for word in ('trial','old','returned')}
        assert e['recovered']=='1' and np.array_equal(u['old'].view(np.uint64),u['returned'].view(np.uint64))
        assert d['directionAvailable']=='1'
        direction=np.array([float(d['dir%d'%j]) for j in range(5)])
        factor=float(d['magnitude']) if d['role']=='jv' else float(d['localRelaxation'])
        check('actualTrialDirection%d'%n,u['trial'],u['old']+factor*direction,2e-13)
        energies={k:float(v[4]/v[0]-.5*np.sum((v[1:4]/v[0])**2)-300000) for k,v in u.items()}
        assert u['old'][0]>0 and elow<=energies['old']<=ehigh
        t=brentq(lambda x:float(air.thermal(x)[1]-air.R*x)-energies['old'],60,2000,xtol=1e-10)
        air.evaluate(u['old'][0]*air.R*t,t)
        trial_t=None;trial_p=None;trial_ceiling=None
        if elow<=energies['trial']<=ehigh:
            trial_t=brentq(lambda x:float(air.thermal(x)[1]-air.R*x)-energies['trial'],60,2000,xtol=1e-10)
            trial_p=u['trial'][0]*air.R*trial_t;trial_ceiling=float(pressure_ceiling(trial_t))
            assert e['reason']=='RAH_PHASE_PT_DOMAIN' and not (0<trial_p<=trial_ceiling)
            pressure_rejections+=1
        else:assert e['reason']=='RAH_PHASE_ENERGY_ENTROPY_DOMAIN'
        roles[d['role']]=roles.get(d['role'],0)+1;trial_low+=energies['trial']<elow
        rejection_data.append(dict(update=int(e['update']),globalID=int(e['globalID']),role=d['role'],
          oldTemperatureK=t,trialPhysicalEnergyJkg=energies['trial'],physicalLowerEnergyJkg=elow,
          trialTemperatureK=trial_t,trialPressurePa=trial_p,trialPressureCeilingPa=trial_ceiling,
          nativeReason=e['reason'],returnedBitEqualOld=True,directionAvailable=True,magnitude=factor))
    if roles.get('jv',0):assert not normal and roles['jv']==1
    save(out/'recovery_events.json',rejection_data)
    # 模块三：全部完整保存帧、原面量和正常末空间/CSV独立复算。
    g={k:v for k,v in np.load(str(GEOMETRY)).items()};tri=g['triangles'];wall=g['wallIDs'];area=g['faceArea']
    frame_paths=sorted(CASE.glob('state_*.bin'));summaries=[];first_bit=None;first_error=None;states=0;last=None;face_by_step={}
    for path in frame_paths:
        step=int(path.stem.split('_')[1]);a=read_frame(path,335726,step==0)['values'];a=a[np.argsort(a[:,0])]
        assert np.array_equal(a[:,1:4],g['points']);states+=len(a)
        if step in (0,1):
            ref=read_frame(run.PREV/'rocket01'/path.name,335726,step==0)['values'];ref=ref[np.argsort(ref[:,0])]
            if step==0:assert np.array_equal(a[:,:27].copy().view(np.uint64),ref[:,:27].copy().view(np.uint64))
            else:
                first_bit=np.array_equal(a.copy().view(np.uint64),ref.copy().view(np.uint64))
                first_error=float(np.max(abs(a-ref)/np.maximum(1,abs(ref)))) # 新预条件允许真实轨迹不同，不能要求首步相同。
        vel=a[:,6:9]/a[:,5,None];energy=a[:,9]/a[:,5]-.5*np.sum(vel**2,axis=1)-300000;t=inverse_energy(air,energy)
        phys=air.evaluate(a[:,5]*air.R*t,t);transport=air.interpolated_transport(a[:,15],a[:,11])
        check('temperature%d'%step,a[:,11],t);check('pressure%d'%step,a[:,15],phys['pressurePa']);check('velocity%d'%step,a[:,12:15],vel)
        check('conductivity%d'%step,a[:,21],transport['conductivityWmK']);check('H%d'%step,a[:,17]-300000,np.asarray(phys['enthalpyJkg'])+.5*np.sum(vel**2,axis=1))
        chi=a[:,10]/(np.asarray(transport['dynamicViscosityPas'])/a[:,5]);mut=np.where(a[:,10]<0,0,a[:,5]*a[:,10]*chi**3/(chi**3+7.1**3))
        check('SAmut%d'%step,a[:,23],mut)
        ids=np.unique(np.r_[np.linspace(0,len(a)-1,11).astype(int),np.argmin(t),np.argmax(t)])
        independent=[brentq(lambda x:float(air.thermal(x)[1]-air.R*x)-energy[i],60,2000,xtol=1e-10) for i in ids]
        check('Brent%d'%step,t[ids],independent)
        q=np.full(len(a),np.nan);q[wall]=a[wall,21]*(a[g['normalNeighbor'],11]-300)/g['thermalDistance']
        faces=np.column_stack((np.arange(1,19753),area,a[tri,15].mean(axis=1),q[tri].mean(axis=1),a[tri,11].mean(axis=1),np.full(19752,300)))
        np.savetxt(str(out/('original_faces_%06d.csv'%step)),faces,delimiter=',',header='faceID,areaM2,pressurePa,heatFluxIntoWallWm2,wallFluidTemperatureK,givenWallTemperatureK',comments='',fmt='%.17g')
        face_by_step[step]=faces;summaries.append(dict(update=step,Tmin=float(min(t)),Tmax=float(max(t)),researchPowerIntoWallW=float(area@faces[:,3])))
        last=a
    assert summaries[0]['update']==0 and (completed<1 or summaries[1]['update']==1)
    terminal_result=None
    if normal:
        assert summaries[-1]['update']==completed
        native=np.loadtxt(str(CASE/'restart_flow.csv'),delimiter=',',skiprows=1,usecols=range(10));native=native[np.argsort(native[:,0])]
        assert np.array_equal(native[:,:4],last[:,:4]) and np.array_equal(native[:,4:10].copy().view(np.uint64),last[:,5:11].copy().view(np.uint64))
        post=binary(CASE/'terminal_post_state.bin',b'RAH248PS',335726,8);post=post[np.argsort(post[:,0])]
        expected=np.column_stack((last[:,:1],last[:,5:11],last[:,23]))
        assert np.array_equal(post.copy().view(np.uint64),expected.copy().view(np.uint64))
        raw,physical=residual(CASE/'terminal_spatial.bin');rms=np.sqrt(np.mean(physical.astype(np.longdouble)**2,axis=0));summ=read(CASE/'terminal_summary.json')
        check('terminalRMS',rms,np.array(summ['physicalRMS']),2e-12)
        scaled=np.abs(physical.astype(np.longdouble).sum(axis=0))/np.maximum(1,np.abs(physical.astype(np.longdouble)).sum(axis=0))
        terminal_result=dict(stateAndSAmutBitUnchanged=True,nativeFinalSixCoordinateBitEqual=True,physicalRMS=np.asarray(rms,dtype=float).tolist(),
          scaledNetMassEnergy=np.asarray(scaled[[0,4]],dtype=float).tolist(),allSixLogRMSBelowMinus8=bool(np.all(rms<=1e-8)),fullSteadyBalanceVerified=False)
    # 独立核对在线原面积变化，第一步对照不推进50步窗口的参照场。
    previous=face_by_step[0]
    for row in rows(CASE/'surface_trend.csv'):
        step=int(row['update'])
        if step not in face_by_step:continue
        cur=face_by_step[step];value=np.sum(area[:,None]*abs(cur[:,2:4]-previous[:,2:4]),axis=0)/np.maximum(1,np.sum(area[:,None]*abs(previous[:,2:4]),axis=0))
        check('surfaceTrend%d'%step,value,[float(row['pressureAreaL1Change']),float(row['heatAreaL1Change'])],2e-12)
        if step!=1:previous=cur
    # 只有完整两窗才复核收益；未达到检查点明确为未执行。
    benefit=rows(CASE/'benefit_decisions.csv')
    for row in benefit:
        end=int(row['windowEnd']);values=np.log10(np.maximum(1e-300,np.array([trend[s] for s in range(end-99,end+1)])))
        ratio=10**(values[50:].mean(axis=0)-values[:50].mean(axis=0))
        check('benefitRatios%d'%end,ratio,[float(row['ratio%d'%j]) for j in range(6)],2e-12)
        gm=float(10**np.mean(values[50:].mean(axis=0)-values[:50].mean(axis=0)))
        check('benefitGeometricMean%d'%end,gm,float(row['geometricMeanRatio']),2e-12)
        before=np.array([float(row['oldPressureChange']),float(row['oldHeatChange'])])
        after=np.array([float(row['newPressureChange']),float(row['newHeatChange'])]);threshold=np.array([.01,.05])
        no_worse=bool(np.all(ratio<=1.25) and np.all(after<=1.25*np.maximum(before,threshold)))
        good=bool(no_worse and np.all((after<=threshold)|(after<=.95*before)) and (gm<=.95 or np.all(values[50:].mean(axis=0)<=-8)))
        assert int(row['noSignificantDeterioration'])==int(no_worse) and int(row['benefit'])==int(good)
    failure_text=''
    if not normal:
        failure_path=CASE/('failure.txt' if (CASE/'failure.txt').exists() else 'caught_failure.txt')
        failure_text=failure_path.read_text(encoding='utf8')
    # 模块四：精度门逐条独立判定；原生报告/有效容差与独立Jacobian证明分开。
    native=read(CASE/'native_resources.json');prepared=read(BASE/'stage01/report.json')
    assert native['points']==prepared['points'] and native['lowerBlocks']==native['upperBlocks']==prepared['edges']
    assert native['outerScalarBytes']==native['preconditionerScalarBytes']==8
    gate_rows=rows(CASE/'linear_acceptance.csv');flow_gate=[a for a in gate_rows if a['system']=='flow'];sa_gate=[a for a in gate_rows if a['system']=='SA']
    rejected=[]
    for a in gate_rows:
        value=float(a['scaledReport']);factor=float(a['toleranceFactor']);target=1e-5
        accepted=np.isfinite(value) and value>=0 and value<=target
        assert int(a['accepted'])==int(accepted) and 1<=factor<=2
        check('gateRaw'+a['system']+a['update'],float(a['rawEstimatedResidual']),factor*value,2e-15)
        check('gateTolerance'+a['system']+a['update'],float(a['effectiveTolerance']),factor*target,2e-15)
        if not accepted:rejected.append(dict(update=int(a['update']),system=a['system'],report=value,iterations=int(a['iterations'])))
    assert len(sa_gate)==completed
    assert all(a['accepted']=='1' for a in flow_gate[:completed])
    first_native=-20.
    for h,a in zip(history,flow_gate):
        check('historyReport'+a['update'],float(h['LinSolRes']),np.log10(float(a['scaledReport'])),2e-9)
        assert int(h['Linear_Solver_Iterations'])==int(a['iterations'])
        base=np.mean([float(h[k]) for k in ('rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]','rms[RhoE]')])
        first_native=max(first_native,base)
        check('nativeToleranceFactor'+a['update'],float(a['toleranceFactor']),1+max(0,1-(base-first_native)/-3),2e-9)
    for h,a in zip(history,sa_gate):
        assert int(h['LinSolIterTurb'])==int(a['iterations'])
        check('SAHistoryReport'+a['update'],float(h['LinSolResTurb']),np.log10(float(a['scaledReport'])),2e-9)
    refusal=None
    if (CASE/'linear_rejection.json').exists():
        refusal=read(CASE/'linear_rejection.json');step=refusal['attemptedUpdate']
        assert not normal and step==completed+1 and refusal['completedUpdates']==completed
        assert not refusal['flowUpdateApplied'] and not refusal['SAUpdateApplied']
        assert rejected==[dict(update=step,system='flow',report=refusal['nativeReportedResidual'],iterations=refusal['nativeIterations'])]
        tags=[p['event'] for p in phases if int(p['update'])==step]
        assert 'update_end' not in tags and 'accepted_preprocessing_begin' not in tags and step not in sa_patterns
        assert 'NK_NATIVE_LINEAR_TOLERANCE_REJECTED_BEFORE_FLOW_UPDATE' in failure_text
        refused=binary(CASE/'linear_refused_state_direction.bin',b'RAH251RJ',335726,17);refused=refused[np.argsort(refused[:,0])]
        assert np.isfinite(refused).all() and np.array_equal(refused[:,0],np.arange(335726))
        old=refused[:,1:6];vel=old[:,1:4]/old[:,0,None];energy=old[:,4]/old[:,0]-.5*np.sum(vel**2,axis=1)-300000
        T=inverse_energy(air,energy);air.evaluate(old[:,0]*air.R*T,T)
        refusal.update(savedOldStatesLegallyDecoded=len(old),savedDirectionComponents=int(refused[:,11:16].size),
          oldStateIsCompletedUpdate=completed,currentStateIsDerivativeProbe=True,physicalUpdateAbsentInStageLog=True,
          oldMinimumK=float(min(T)),oldMaximumK=float(max(T)))
        save(out/'linear_refusal_audit.json',refusal)
    linear.update(gateRecords=len(gate_rows),flowGateRecords=len(flow_gate),SAReports=len(sa_gate),rejected=rejected,
      configuredTolerance=1e-5,effectiveToleranceUpperBound=2e-5,nativeResidualScaledByToleranceFactor=True,
      SAObservedAfterUpdateOnly=True,SAFailedReports=sum(a['accepted']!='1' for a in sa_gate),independentLinearAccuracyProved=False)
    linear['SAZeroIterationEarlyExitCanReportAbsoluteNorm']=True
    result=dict(passed=True,nativeExitedNormally=normal,completedUpdates=completed,completeSavedUpdates=[x['update'] for x in summaries],
      actualFactory='CNewtonIntegration',observedMatrixFreeProducts=len(begins),completedMatrixFreeProducts=len(ends),nativeLinearHistory=linear,linearRefusal=refusal,nativeMemoryGraph=native,
      firstUpdateBitEqual=bool(first_bit),firstUpdateMaximumScaledDifference=first_error,initialBitEqual=True,
      nativeNumericalStatementsUnchanged=True,recoveryCount=len(events),recoveryRoles=roles,trialBelowLowerEnergyCount=int(trial_low),
      jointPressureTemperatureRejections=pressure_rejections,firstAuditEnergyOnlyAssumptionFailed=False,
      independentSavedStates=states,frames=summaries,terminal=terminal_result,benefitDecisionReached=bool(benefit),
      stopReason=done['stopReason'] if done else ('linear-tolerance-rejection' if refusal else 'native-domain-or-runtime-failure'),failureMessage=failure_text,originalFaces=19752,originalAreaM2=1.9646292565715553,
      newQualifiedAreaM2=0,physicalQualification=False,nonlinearConvergenceProved=False,fullJacobianAccuracyProved=False,
      supervisedTaskPassed=process['passed'],resourceContractPassed=process['stopReason'] is None and capacity['passed'],
      workingRHSIsFinalResidual=False,errorMetrics=errors)
    save(out/'report.json',result)
    print({k:result[k] for k in ('passed','nativeExitedNormally','completedUpdates','recoveryRoles','firstUpdateMaximumScaledDifference','benefitDecisionReached')})


def main():
    if CLOSED:raise SystemExit('v251 native ILU audit closed; follow NEXT_DECISION')
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/terminalCapacityV249.py',ROOT/'scripts/singleDomainFrame.py',
      ROOT/'scripts/fasRecoveryFieldsV248.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/su2PhaseAirFields.py',DOC/'PLAN.md',GEOMETRY,
      BASE/'process_rocket/report.json',BASE/'process_rocket/terminal_capacity.json',BASE/'stage01/CNewtonIntegration_original.cpp',BASE/'stage01/CNewtonIntegration_observed.cpp',
      run.PREV/'native01/state_000000.bin',run.PREV/'native01/state_000001.bin']+[p for p in CASE.iterdir() if p.is_file()]
    return supervise(run.setup_guard(),BASE,'audit',[sys.executable,str(Path(__file__).resolve()),'task'],120,1536,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v251 native ILU audit closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':audit()
    else:raise SystemExit(0 if main() else 1)
