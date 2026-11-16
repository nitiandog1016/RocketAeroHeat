"""模块职责：审核真实FAS五步控制及原面研究场，明确延拓工作数组不等于流动残差。"""
from pathlib import Path
import csv,struct,sys
import numpy as np
from scipy.optimize import brentq
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy
from singleDomainFrame import read_frame
from auditFASRocketPreflightV246 import frame
import runFASRocketControlV246 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;CASE=BASE/'control_observe01';CLOSED=False
GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz'


def rows(path):
    with Path(path).open(encoding='utf8') as f:return list(csv.DictReader(f))


def audit():
    out=BASE/'control_audit01';out.mkdir(exist_ok=False);errors={};air=PhaseAwareViscousAir()
    def check(name,a,b,tol=2e-9):
        a=np.asarray(a);b=np.asarray(b);assert np.isfinite(a).all() and np.isfinite(b).all(),name
        e=float(np.max(abs(a-b)/np.maximum(1,abs(b))));errors[name]=max(errors.get(name,0),e);assert e<tol,(name,e)
    names=['state_000000.bin','state_000001.bin','state_000005.bin','coarse_initial.bin','coarse_final.bin','history.csv','restart_flow.csv']
    for n in names:assert sha(CASE/n)==sha(BASE/'control_off01'/n),n
    for mode in ('observe','off'):
        r=read(BASE/('control_'+mode+'01/native_completion.json'));assert r['updates']==5 and r['normalReturn']
        assert read(BASE/('process_control_'+mode+'/report.json'))['passed']
    stage=rows(CASE/'stages.csv');levels=rows(CASE/'levels.csv');steps=rows(CASE/'steps.csv');assert len(stage)==80 and len(levels)==100 and len(steps)==5
    assert (CASE/'rejections.csv').stat().st_size==0 and not (CASE/'fas_rejection.txt').exists()
    for r in steps:assert all(float(r[k])==0 for k in ('lastFlowRejectionCount','totalObservedRejections','maxWallSpeedComponent','maxWallSA'))
    after={(int(r['update']),int(r['system'])):r for r in stage if r['stage']=='after_time'};assert len(after)==10
    for k,r in enumerate(stage):
        if r['stage']=='before_time':
            nextrow=stage[k+1];assert nextrow['stage']=='after_time'
            check('timeRHSSign',[float(r['sum%d'%j]) for j in range(6)],[-float(nextrow['sum%d'%j]) for j in range(6)])
    # 同一次壁处理序列配对；FAS每步有两次细层空间评估，不能用(update,label)覆盖第一轮。
    before=[r for r in stage if r['stage']=='before_strong' and r['system']=='0'];afterwall=[r for r in stage if r['stage']=='after_strong' and r['system']=='0']
    wall=rows(CASE/'wall_stages.csv');assert len(before)==len(afterwall)==10 and len(wall)==20
    for j,(b,a) in enumerate(zip(before,afterwall)):
        for w in wall[2*j:2*j+2]:
            assert w['update']==a['update']==b['update'];check('sameStageWallHeat',float(a['sum4'])-float(b['sum4']),float(w['effectiveHeatIntoWallW']),2e-8)
            if w['stage']=='after_strong':check('wallMomentumRows',[float(w['wallRow'+x]) for x in ('X','Y','Z')],[0,0,0])
    source=read(ROOT/'outputs/v222_validation/geometry01/report.json');assert sha(GEOMETRY)==source['cacheSHA256']
    g={k:v for k,v in np.load(str(GEOMETRY)).items()};tri=g['triangles'];wallids=g['wallIDs'];area=g['faceArea'];records=[];last=None
    phase_frames=[];sample_checks=0
    for step in (0,1,5):
        raw=read_frame(CASE/('state_%06d.bin'%step),335726,step==0)['values'];a=raw[np.argsort(raw[:,0])]
        assert np.array_equal(a[:,1:4],g['points']);vel=a[:,6:9]/a[:,5,None];ep=a[:,9]/a[:,5]-.5*np.sum(vel*vel,axis=1)-300000
        t=inverse_energy(air,ep);state=air.evaluate(a[:,5]*air.R*t,t);check('fineCurrentTemperature',a[:,11],t);check('fineCurrentPressure',a[:,15],state['pressurePa'])
        check('fineCurrentVelocity',a[:,12:15],vel);check('finePhysicalH',a[:,17]-300000,np.asarray(state['enthalpyJkg'])+.5*np.sum(vel*vel,axis=1))
        tr=air.interpolated_transport(a[:,15],a[:,11]);mu=np.asarray(tr['dynamicViscosityPas']);k=np.asarray(tr['conductivityWmK']);check('mu',a[:,19],mu);check('k',a[:,21],k)
        chi=a[:,10]/(mu/a[:,5]);mut=np.where(a[:,10]<0,0,a[:,5]*a[:,10]*chi**3/(chi**3+7.1**3));check('SAeddyViscosity',a[:,23],mut)
        # 抽取温度/能量极值及分布站，用Brent积分反解作独立数值路径交叉检查。
        selected=np.unique(np.r_[np.linspace(0,len(a)-1,16).astype(int),np.argmin(ep),np.argmax(ep)])
        tb=np.array([brentq(lambda x:float(air.thermal(x)[1]-air.R*x)-ep[i],60,2000,xtol=1e-10) for i in selected]);check('BrentEnergyCrosscheck',t[selected],tb);sample_checks+=len(selected)
        q=np.zeros(len(a));q[wallids]=k[wallids]*(a[g['normalNeighbor'],11]-300)/g['thermalDistance']
        faces=np.column_stack((np.arange(1,19753),area,a[tri,15].mean(axis=1),q[tri].mean(axis=1),a[tri,11].mean(axis=1),np.full(19752,300.)))
        np.savetxt(str(out/('original_faces_%06d.csv'%step)),faces,delimiter=',',header='faceID,areaM2,pressurePa,heatFluxIntoWallWm2,wallFluidTemperatureK,givenWallTemperatureK',comments='',fmt='%.17g')
        rms=None if step==0 else [float(after[(step,0 if j<5 else 2)]['rms%d'%j]) for j in range(6)]
        records.append(dict(update=step,pRangePa=[float(min(faces[:,2])),float(max(faces[:,2]))],qRangeWm2=[float(min(faces[:,3])),float(max(faces[:,3]))],
           wallFluidRangeK=[float(min(faces[:,4])),float(max(faces[:,4]))],powerIntoWallW=float(area@faces[:,3]),
           physicalSmoothingLogRMS=None if rms is None else np.log10(np.maximum(rms,1e-300)).tolist(),physicalQualification=False))
        if step==5:last=a;last_raw=raw
    # 原生末CSV与全部五流动量/SA逐位同，不使用其“Residual”列当最后流动RHS。
    csvstate=np.loadtxt(str(CASE/'restart_flow.csv'),delimiter=',',skiprows=1,usecols=range(10));csvstate=csvstate[np.argsort(csvstate[:,0])]
    assert np.array_equal(csvstate[:,:4],last[:,:4]) and np.array_equal(csvstate[:,4:10],last[:,5:11])
    coarse=frame(CASE/'coarse_final.bin');initial=frame(CASE/'coarse_initial.bin');pre=frame(BASE/'initial01/coarse_constructed.bin')
    assert np.array_equal(initial,pre) and np.array_equal(coarse[:,:5],pre[:,:5])
    vel=coarse[:,6:9]/coarse[:,5,None];ep=coarse[:,9]/coarse[:,5]-.5*np.sum(vel*vel,axis=1)-300000;t=inverse_energy(air,ep);air.evaluate(coarse[:,5]*air.R*t,t)
    with (BASE/'initial01/children.bin').open('rb') as stream:stream.read(24);pairs=np.fromfile(stream,dtype='<u8').reshape(-1,2).astype(np.int64)
    parent,child=pairs.T;weight=last_raw[child,4]/coarse[parent,4]
    check('actualSACoarseRestriction',coarse[:,10],np.bincount(parent,weights=last_raw[child,10]*weight,minlength=len(coarse)))
    mut=np.bincount(parent,weights=last_raw[child,23]*weight,minlength=len(coarse));boundary=np.loadtxt(str(BASE/'initial01/boundary.csv'),delimiter=',',skiprows=1)
    coarsewall=boundary[(boundary[:,0]==1)&(boundary[:,15]==1),3].astype(int);mut[coarsewall]=0
    check('actualCoarseEddyRestriction',coarse[:,11],mut);check('coarseNoSlip',coarse[coarsewall,6:9],np.zeros((len(coarsewall),3)))
    # 最后传回细层的修正数组：证明它是粗层修正注入，和阶段残差分开保存。
    parentof=np.empty(335726,dtype=np.int64);parentof[child]=parent
    correction=last_raw[:,27:32];first=np.zeros((len(coarse),5));seen=np.zeros(len(coarse),bool)
    for i in child:
        p=parentof[i]
        if not seen[p]:first[p]=correction[i];seen[p]=True
    check('savedFineArrayIsInjectedCorrection',correction,first[parentof])
    schema=dict(format='RAH235LE legacy binary layout',indicesZeroBased=True,
       fineFlowSlots27to31='FAS prolongation correction stored in LinSysRes, NOT flow residual',
       fineSAColumn32='last fine SA implicit RHS',initialRHSUnknown=True,
       physicalResidualSource='stages.csv after_time: fine pre-smoothing stage, before coarse correction; NOT final-cycle residual',
       coarseFinalPrimitiveCache='coarse pre-update cache; current state must be decoded from conserved variables',
       solidWallTemperatureSolved=False,physicalQualification=False)
    save(out/'field_semantics.json',schema)
    report=dict(passed=True,completeNativeFiveStepCompatibilityPassed=True,observationNeutralBitEqual=True,completeFrames=records,
       successfulUpdates=5,stageRows=len(stage),levelRows=len(levels),nativeJointDomainStateVisits=sum(int(r['nodes']) for r in levels),
       independentSavedStateChecks=3*335726+2*44484,independentBrentSamples=sample_checks,actualCoarseSAAndWallProjectionPassed=True,
       errorMetrics=errors,nativeFinalSixStateAndCoordinateBitEqual=True,savedFineFlowArrayIsCorrection=True,
       finalCycleResidualNotSaved=True,physicalMinus8All=bool(np.all(np.array(records[-1]['physicalSmoothingLogRMS'])<=-8)),
       originalFaces=19752,originalAreaM2=float(sum(area)),givenTwK=300,solidTwSolved=False,
       fullSteadyBalanceVerified=False,spaceOuterDomainHeatValidationPassed=False,physicalQualification=False,newQualifiedAreaM2=0,productionAdmitted=False)
    save(out/'report.json',report);print({k:report[k] for k in ('passed','successfulUpdates','levelRows','nativeJointDomainStateVisits','physicalMinus8All')})


def main():
    if CLOSED:raise SystemExit('v246 real control audit closed; follow NEXT_DECISION')
    g=run.prep.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda p:size(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',GEOMETRY,
      ROOT/'scripts/auditFASRocketPreflightV246.py',ROOT/'scripts/singleDomainFrame.py',ROOT/'scripts/su2PhaseAirFields.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py',
      ROOT/'data/thermo/phase_aware_dry_air_v1.json',BASE/'initial01/children.bin',BASE/'initial01/boundary.csv',BASE/'initial01/coarse_constructed.bin']
    locks += [p for mode in ('observe','off') for p in (BASE/('control_'+mode+'01')).iterdir() if p.is_file()]
    return g.supervise('control_audit',[sys.executable,str(Path(__file__).resolve()),'task'],120,1536,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v246 real control audit closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':audit()
    else:raise SystemExit(0 if main() else 1)
