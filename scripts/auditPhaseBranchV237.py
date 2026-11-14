"""模块职责：独立核对真实重启动初场、保存场、阶段/收益合同及原面研究出口；不执行CFD。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runPhaseBranchV237 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from singleDomainFrame import read_frame
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy
BASE=run.BASE;CASE=run.CASE;DOC=run.DOC;CLOSED=True


def rows(p):
    with p.open(encoding='utf8') as f:return list(csv.DictReader(f))


def audit(output_name='real_audit01'):
    out=BASE/output_name;out.mkdir(exist_ok=False);err={};records=[];air=PhaseAwareViscousAir()
    process=read(BASE/'process_rocket/report.json');assert process['status']=='terminal' and process['inputsStable'] and not process['liveGroupAfterCleanup']
    contract=read(BASE/'contract.json');assert sha(run.CHECKPOINT)==contract['sourceCheckpointSHA256']
    g=np.load(str(run.GEOMETRY));wall=g['wallIDs'];tri=g['triangles'];area=g['faceArea'];N=335726
    def check(name,a,b,tol=2e-9):
        a=np.asarray(a);b=np.asarray(b);assert np.isfinite(a).all() and np.isfinite(b).all(),name
        e=float(np.max(np.abs(a-b)/np.maximum(1.,np.abs(b))));err[name]=max(err.get(name,0),e);assert e<tol,(name,e)
    initial=read(CASE/'initial_check.json');assert initial['points']==N and initial['stateBitMismatches']==0 and initial['coordinateBitMismatches']==0
    assert initial['flowCFLMin']==.1 and initial['flowCFLMax']==.1
    source=np.loadtxt(str(run.CHECKPOINT),delimiter=',',skiprows=1,usecols=range(10))
    source=source[np.argsort(source[:,0])];assert np.array_equal(source[:,0],np.arange(N))
    allStages=rows(CASE/'stages.csv');stages=allStages;index={(int(r['update']),int(r['system']),r['stage']):r for r in stages}
    steps=rows(CASE/'steps.csv');n=len(steps);assert [int(r['update']) for r in steps]==list(range(1,n+1))
    expected={'after_sources','after_weak','before_strong','after_strong','before_time','after_time'}
    logs=[]
    for step in range(1,n+1):
        rms=[]
        for system in (0,2):
            assert {r['stage'] for r in stages if int(r['update'])==step and int(r['system'])==system}==expected
            before=index[(step,system,'before_time')];after=index[(step,system,'after_time')]
            for j in (range(5) if system==0 else [5]):
                value=float(after['rms%d'%j]);assert np.isfinite(value) and value>=0;rms.append(value)
                # 相消和采用明确绝对舍入界，避免近零净量的相对误判。
                delta=abs(float(before['sum%d'%j])+float(after['sum%d'%j]));bound=128*np.finfo(float).eps*max(1,N*value)
                assert delta<=bound,('IMPLICIT_SIGN',step,j,delta,bound)
        logs.append(np.log10(np.maximum(1e-300,rms)))
    logs=np.array(logs)
    for r in steps:check('acceptedWall',[float(r['maxWallSpeedComponent']),float(r['maxWallSA'])],[0,0])
    # 模块一：全场独立量热/联合域、能量坐标、实际SA与完整残差；初0单独逐位核对。
    faceframes={};last=None
    for path in sorted(CASE.glob('state_*.bin')):
        step=int(path.stem.split('_')[-1]);frame=read_frame(path,N,step==0);a=frame['values'];a=a[np.argsort(a[:,0])]
        assert np.array_equal(a[:,1:4],g['points'])
        if step==0:
            assert np.array_equal(a[:,5:11].view(np.uint64),source[:,4:10].copy().view(np.uint64))
            assert np.array_equal(a[:,1:4].copy().view(np.uint64),source[:,1:4].copy().view(np.uint64))
        u=a[:,6:9]/a[:,5,None];e=a[:,9]/a[:,5]-.5*np.sum(u*u,axis=1)-300000
        T=inverse_energy(air,e);p=a[:,5]*air.R*T;state=air.evaluate(p,T);primitive=air.evaluate(a[:,15],a[:,11])
        check('temperatureFromConserved',T,a[:,11]);check('pressureFromConserved',p,a[:,15]);check('primitiveVelocity',u,a[:,12:15])
        check('enthalpy',a[:,17]-300000,np.asarray(primitive['enthalpyJkg'])+.5*np.sum(u*u,axis=1))
        tr=air.interpolated_transport(a[:,15],a[:,11]);mu=np.asarray(tr['dynamicViscosityPas']);k=np.asarray(tr['conductivityWmK'])
        check('linearMu',a[:,19],mu);check('linearK',a[:,21],k)
        chi=a[:,10]/(mu/a[:,5]);positive=a[:,10]>=0;mut=np.zeros(N)
        mut[positive]=a[positive,5]*a[positive,10]*chi[positive]**3/(chi[positive]**3+7.1**3)
        check('updatedSAMuT',a[:,23],mut)
        if step:
            physical=a[:,27:33].copy();physical[:,4]-=300000*physical[:,0]
            actual=np.sqrt(np.mean(physical**2,axis=0));check('savedPhysicalRMS',np.log10(np.maximum(actual,1e-300)),logs[step-1])
            for j in range(6):
                r=index[(step,0 if j<5 else 2,'after_time')];delta=abs(np.sum(physical[:,j])-float(r['sum%d'%j]));bound=128*np.finfo(float).eps*max(1.,np.sum(abs(physical[:,j])))
                assert delta<=bound,('PHYSICAL_SUM',step,j,delta,bound)
        q=np.zeros(N);q[wall]=k[wall]*(a[g['normalNeighbor'],11]-300)/g['thermalDistance']
        fp=a[tri,15].mean(axis=1);fq=q[tri].mean(axis=1);ft=a[tri,11].mean(axis=1)
        faces=np.column_stack([np.arange(1,19753),area,fp,fq,ft,np.full(19752,300.)]);faceframes[step]=faces
        record=dict(update=step,frameSHA256=sha(path),pRangePa=[float(fp.min()),float(fp.max())],qRangeWm2=[float(fq.min()),float(fq.max())],
            wallFluidRangeK=[float(ft.min()),float(ft.max())],powerIntoWallW=float(area@fq),negativeHeatFaces=int(np.sum(fq<0)),
            wallTemperatureAreaRMSK=float(np.sqrt(area@((ft-300)**2)/sum(area))),
            jointDomainAccepted=True,physicalLogRMS=None if not step else logs[step-1].tolist(),physicalQualification=False)
        if last is not None:record['previousSaveChange']=[float(area@abs(faces[:,j]-last[:,j])/max(area@abs(last[:,j]),1)) for j in (2,3)]
        last=faces;records.append(record);save(out/'progress.json',records)
        del a,frame
    # 模块二：每50步的壁节点和原面变化独立回算，用同面积分母检查在线停止逻辑。
    trends=rows(CASE/'surface_trend.csv');change={};previous=None;wallrecords=[]
    for path in sorted(CASE.glob('wall_[0-9][0-9][0-9][0-9][0-9][0-9].csv')):
        step=int(path.stem.split('_')[-1]);w=np.loadtxt(str(path),delimiter=',',skiprows=1);ids=w[:,0].astype(np.int64)
        assert np.array_equal(np.sort(ids),np.sort(wall)) and len(set(ids.tolist()))==9878
        assert np.isfinite(w[:,1:4]).all() and np.isfinite(w[:,6:12]).all()
        if step:assert np.isfinite(w).all()
        else:assert np.isnan(w[:,[4,5,12,13,14]]).all()
        check('wallGradient',w[:,3],w[:,9]*(w[:,8]-300)/w[:,7])
        native=np.zeros((N,3));native[ids]=w[:,[1,3,2]];values=native[tri].mean(axis=1)
        faces=np.column_stack([np.arange(1,19753),area,values,np.full(19752,300.)])
        if step in faceframes:check('wallVsFrame',faces,faceframes[step])
        faceframes[step]=faces
        d=np.zeros(2) if previous is None else np.array([area@abs(faces[:,j]-previous[:,j])/max(area@abs(previous[:,j]),1) for j in (2,3)])
        row=next(r for r in trends if int(r['update'])==step);check('onlineAreaChanges',d,[float(row['pressureAreaL1Change']),float(row['heatAreaL1Change'])]);change[step]=d;previous=faces
        wallrecords.append(dict(update=step,yPlusNativeMax=None if not step else float(w[:,5].max()),
            nativeCacheVsCurrentGradientMaxScaled=None if not step else float(np.max(abs(w[:,4]-w[:,3])/np.maximum(1,abs(w[:,3]))))))
    decisions=rows(CASE/'benefit_decisions.csv');decisionSummary=[]
    for checkpoint in (200,300):
        if checkpoint>n:continue
        worse=False;noBenefit=0
        for end in (checkpoint-50,checkpoint):
            a=logs[end-100:end-50].mean(axis=0);b=logs[end-50:end].mean(axis=0);ratio=10**(b-a);gm=float(10**np.mean(b-a))
            old=change[end-50];now=change[end];threshold=np.array([.01,.05])
            notWorse=bool(np.all(ratio<=1.25) and np.all(now<=1.25*np.maximum(old,threshold)))
            good=notWorse and bool(gm<=.95 or np.all(b<=-8)) and bool(np.all((now<=threshold)|(now<=.95*old)))
            row=next(r for r in decisions if int(r['check'])==checkpoint and int(r['windowEnd'])==end)
            check('benefitRatios',ratio,[float(row['ratio%d'%j]) for j in range(6)]);check('benefitGeometricMean',gm,float(row['geometricMeanRatio']))
            assert int(row['noSignificantDeterioration'])==notWorse and int(row['benefit'])==good
            worse|=not notWorse;noBenefit+=not good
        stop=worse or noBenefit==2
        decisionSummary.append(dict(checkpoint=checkpoint,stop=stop,reason='significant-deterioration' if worse else 'two-windows-no-benefit' if stop else 'continue'))
        if stop:assert n==checkpoint
    # 模块三：新拒绝回调逐次核对原始试探与实际五量恢复；SA只观察、不声称联合回退。
    events=rows(CASE/'rejections.csv');eventAudit=[]
    for r in events:
        trial=np.array([float(r['trial%d'%j]) for j in range(5)]);old=np.array([float(r['old%d'%j]) for j in range(5)]);returned=np.array([float(r['returned%d'%j]) for j in range(5)])
        assert np.array_equal(old.view(np.uint64),returned.view(np.uint64)) and int(r['recovered'])==1
        assert 0<=int(r['globalID'])<N
        def valid(u):
            if not np.isfinite(u).all() or u[0]<=0:return False
            try:
                T=inverse_energy(air,np.array([u[4]/u[0]-.5*np.sum((u[1:4]/u[0])**2)-300000]));air.evaluate(u[0]*air.R*T,T);return True
            except ValueError:return False
        oldValid=valid(old);trialValid=valid(trial)
        assert oldValid and not trialValid,('REJECTION_ADMISSIBILITY',r)
        eventAudit.append(dict(update=int(r['update']),globalID=int(r['globalID']),phase=r['phase'],reason=r['reason'],oldAccepted=oldValid,trialRejected=not trialValid,fiveStateRecoveryBitEqual=True))
    save(out/'rejection_audit.json',eventAudit)
    sameStage=[r for r in rows(CASE/'wall_stages.csv') if int(r['update'])<=n];heatErrors=[];momentum=[]
    for r in sameStage:
        step=int(r['update']);before=index[(step,0,'before_strong')];after=index[(step,0,'after_strong')]
        delta=float(after['sum4'])-float(before['sum4']);q=float(r['effectiveHeatIntoWallW'])
        check('sameStageWallHeat',delta,q,2e-8);heatErrors.append(abs(delta-q))
        if r['stage']=='after_strong':check('wallMomentumReplacement',[float(r['wallRow'+d]) for d in ('X','Y','Z')],[0,0,0])
        else:
            row=np.array([float(r['wallRow'+d]) for d in ('X','Y','Z')]);force=np.array([float(r['pNout'+d])-float(r['tauNout'+d]) for d in ('X','Y','Z')])
            momentum.append(dict(update=step,row=row.tolist(),pressureMinusShear=force.tolist(),oppositeBoundaryClosureScaled=float(np.linalg.norm(row+force)/max(np.linalg.norm(force),1)),independentMomentumBalanceVerified=False))
    save(out/'same_stage_momentum.json',momentum)
    for step,faces in sorted(faceframes.items()):np.savetxt(str(out/('original_faces_%06d.csv'%step)),faces,delimiter=',',header='faceID,areaM2,pressurePa,heatFluxIntoWallWm2,wallFluidTemperatureK,givenWallTemperatureK',comments='',fmt='%.17g')
    completion=read(CASE/'native_completion.json') if (CASE/'native_completion.json').exists() else None
    if completion:
        assert completion['updates']==n and (CASE/('state_%06d.bin'%n)).exists()
        if decisionSummary and decisionSummary[-1]['stop']:assert completion['stopReason']==decisionSummary[-1]['reason']
    report=dict(passed=True,engineNormal=process['passed'],successfulUpdates=n,completeFrames=records,wallSnapshots=wallrecords,
        initialSixStateBitEqual=True,initialCoordinatesBitEqual=True,initialCFLReset=.1,exactHistoricalContinuation=False,
        errorMetrics=err,stageRows=len(stages),partialStageRows=sum(int(r['update'])>n for r in stages),newRejectionEvents=len(events),eventsIndependentlyChecked=bool(events),benefitDecisions=decisionSummary,
        nativeCompletion=completion,originalFaces=19752,originalAreaM2=float(sum(area)),physicalQualification=False,newQualifiedAreaM2=0,
        physicalMinus8All=bool(n and np.all(logs[-1]<=-8)),lastCompletePhysicalLogRMS=logs[-1].tolist(),fullSteadyBalanceVerified=False,spaceOuterDomainAndHeatValidationPassed=False,
        sameStageMomentumFinal=None if not momentum else momentum[-1],MATLABExecuted=False)
    save(out/'report.json',report);print(json.dumps(dict(passed=True,updates=n,events=len(events),physicalQualification=False)))


def main(stage='audit'):
    if CLOSED:raise SystemExit('v237 real audit closed; follow NEXT_DECISION')
    run.run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runPhaseBranchV237.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
        ROOT/'scripts/singleDomainFrame.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',
        ROOT/'data/thermo/phase_aware_dry_air_v1.json',run.GEOMETRY,run.CHECKPOINT,BASE/'contract.json',BASE/'process_rocket/report.json']
    locks += [p for p in CASE.iterdir() if p.is_file()]
    assert stage=='real_audit'
    number=1
    return run.run.guard.supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task','real_audit%02d'%number],180,1024,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v237 real audit closed; follow NEXT_DECISION')
    if len(sys.argv)>1 and sys.argv[1]=='task':audit(sys.argv[2])
    else:raise SystemExit(0 if main(sys.argv[1] if len(sys.argv)>1 else 'real_audit') else 1)
