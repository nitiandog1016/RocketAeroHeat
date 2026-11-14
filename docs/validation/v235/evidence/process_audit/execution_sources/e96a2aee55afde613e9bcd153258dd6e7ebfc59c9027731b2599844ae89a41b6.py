"""模块职责：核对原生5步观察中性、完整阶段/能量/保存与拒绝边界；不运行CFD。"""
from pathlib import Path
import csv,json,struct,sys
import numpy as np
from singleDomainFrame import read_frame,COLUMNS
from phaseAwareViscousAir import PhaseAwareViscousAir
import prepareSingleDomainV235 as prep
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=prep.BASE;DOC=prep.DOC;CLOSED=False


def rows(p):
    with p.open(encoding='utf8') as f:return list(csv.DictReader(f))


def audit():
    out=BASE/'audit01';out.mkdir(exist_ok=False);err={};air=PhaseAwareViscousAir()
    def check(name,a,b,tol=2e-9):
        a=np.asarray(a);b=np.asarray(b);assert np.isfinite(a).all() and np.isfinite(b).all()
        e=float(np.max(np.abs(a-b)/np.maximum(1.,np.abs(b))));err[name]=max(err.get(name,0),e);assert e<tol,(name,e)
    for n in ('state_000000.bin','state_000001.bin','history.csv'):
        assert sha(BASE/'off01'/n)==sha(BASE/'observe01'/n),n
    for n in ('state_000000.bin','state_000001.bin','state_000005.bin','history.csv','restart_flow.csv'):
        assert sha(BASE/'off02'/n)==sha(BASE/'observe02'/n),n
    for name in ('off02','observe02'):
        r=read(BASE/name/'native_completion.json');assert r['updates']==5 and r['normalReturn'] and r['initialSaved']
    stage=rows(BASE/'observe02/stages.csv');assert len(stage)==60
    indexed={(int(r['update']),int(r['system']),r['stage']):r for r in stage}
    expected={'after_sources','after_weak','before_strong','after_strong','before_time','after_time'}
    for step in range(1,6):
        for system in (0,2):assert {r['stage'] for r in stage if int(r['update'])==step and int(r['system'])==system}==expected
        for system in (0,2):
            before=indexed[(step,system,'before_time')];after=indexed[(step,system,'after_time')]
            check('implicitRHSSign',[float(after['sum%d'%j]) for j in range(6)],[-float(before['sum%d'%j]) for j in range(6)])
    fields=[]
    for step in (0,1,5):
        a=read_frame(BASE/'observe02'/('state_%06d.bin'%step),24,step==0)['values'];q=air.evaluate(a[:,15],a[:,11])
        u=a[:,6:9]/a[:,5,None];e=a[:,9]/a[:,5]-.5*np.sum(u*u,axis=1)-300000
        check('conservativeEnergy',e,np.asarray(q['internalEnergyJkg']))
        check('primitiveVelocity',a[:,12:15],u);check('density',a[:,5],q['densityKgm3'])
        check('physicalH',a[:,17]-300000,np.asarray(q['enthalpyJkg'])+.5*np.sum(u*u,axis=1))
        tr=air.interpolated_transport(a[:,15],a[:,11]);mu=np.asarray(tr['dynamicViscosityPas']);chi=a[:,10]/(mu/a[:,5])
        expected_mu=np.where(a[:,10]<0,0,a[:,5]*a[:,10]*chi**3/(chi**3+7.1**3))
        check('updatedSAMuT',a[:,23],expected_mu)
        if step:
            flow=a[:,27:32].copy();flow[:,4]-=300000*flow[:,0]
            for system,values in ((0,flow),(2,a[:,32:33])):
                r=indexed[(step,system,'after_time')];start=0 if system==0 else 5
                for j in range(values.shape[1]):
                    check('savedPhysicalRMS',np.sqrt(np.mean(values[:,j]**2)),float(r['rms%d'%(start+j)]))
                    # 对近零总和显式用分量绝对和给出舍入尺度。
                    delta=abs(np.sum(values[:,j])-float(r['sum%d'%(start+j)]));bound=64*np.finfo(float).eps*max(1.,np.sum(np.abs(values[:,j])))
                    assert delta<=bound,('EXPLICIT_SUM_BOUND',delta,bound)
        fields.append(dict(update=step,sha256=sha(BASE/'observe02'/('state_%06d.bin'%step))))
    wall=rows(BASE/'observe02/wall_stages.csv');assert len(wall)==10
    for r in wall:
        step=int(r['update']);before=indexed[(step,0,'before_strong')];after=indexed[(step,0,'after_strong')]
        check('sameStageHeatIncrement',float(after['sum4'])-float(before['sum4']),float(r['effectiveHeatIntoWallW']))
        if r['stage']=='after_strong':check('wallMomentumRows',[float(r['wallRow'+d]) for d in ('X','Y','Z')],[0,0,0])
    steps=rows(BASE/'observe02/steps.csv');assert len(steps)==5
    for r in steps:check('actualAcceptedWall',[float(r[k]) for k in ('maxWallSpeedComponent','maxWallSA')],[0,0])
    assert not rows(BASE/'observe02/rejections.csv')
    # 模块二：真实文件边界，错误身份/截断/非有限状态/错误初场RHS/重复ID必须拒绝。
    original=(BASE/'observe02/state_000005.bin').read_bytes();bad=[]
    variants=[('header',original[:7]),('magic',b'INVALID!'+original[8:]),('length',original[:-1]),
      ('columns',original[:16]+struct.pack('<Q',35)+original[24:])]
    for name,col,value in [('nan',5,float('nan')),('duplicate',34,0),('reserved',33,1)]:
        a=bytearray(original);a[24+8*col:24+8*(col+1)]=struct.pack('<d',value);variants.append((name,bytes(a)))
    for name,data in variants:
        p=out/('reject_'+name+'.bin');p.write_bytes(data)
        try:read_frame(p,24)
        except ValueError as error:bad.append(dict(case=name,error=str(error)))
        else:raise AssertionError('INVALID_FRAME_ACCEPTED_'+name)
    try:read_frame(BASE/'observe02/state_000005.bin',24,True)
    except ValueError as error:bad.append(dict(case='initial-rhs',error=str(error)))
    else:raise AssertionError('KNOWN_RHS_AS_INITIAL_ACCEPTED')
    assert len(bad)==8
    save(out/'frame_rejections.json',bad)
    # 原199文件仅提供保守态/SA，未保存自适应CFL和松弛历史，不能宣称精确迭代接续。
    old=ROOT/'outputs/v222_validation/rocket01/restart_flow_000199.csv'
    with old.open(encoding='utf8') as f:header=next(csv.reader(f))
    assert 'Nu_Tilde' in header and not any('CFL' in x or 'Relax' in x for x in header)
    save(out/'real_admission_decision.json',dict(fullRocketAdmitted=False,
      reason='native restart has five conservative fields and SA but no per-node CFL or adaptive relaxation history; no claimed exact continuation',
      next='bind checkpoint 199 physical/SA identity and native restart reset semantics, then one explicitly budgeted development branch',
      frozenFullDomainRerunRequired=False,repeatManufacturedControlRequired=False,newSourceRouteRequired=False,
      sourceHeader=header,maximumObservedSingleDomainRSSBytes=2153095168))
    save(out/'report.json',dict(passed=True,controlUpdatesEach=5,initialFirstFinalFramesBitIdentical=True,
      nativeHistoryAndRestartBitIdentical=True,stages=len(stage),wallStageRows=len(wall),errorMetrics=err,
      structuralRejections=len(bad),newRejectionCallbackDynamicallyExercised=False,
      firstInconsistentBoundaryPairPreserved=True,firstPairOneUpdateFramesBitIdentical=True,
      nativeProductionLoopExercised=True,fullRocketFlowUpdates=0,newQualifiedAreaM2=0,physicalQualification=False,
      sourceObjectsRecompiled=2,frames=fields))
    print(json.dumps(dict(passed=True,updatesEach=5,stages=len(stage),maxError=max(err.values()))))


def main():
    if CLOSED:raise SystemExit('v235 audit closed; follow NEXT_DECISION')
    g=prep.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda _:size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/singleDomainFrame.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',
      ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',ROOT/'data/thermo/phase_aware_dry_air_v1.json']
    locks += [p for name in ('off01','observe01','off02','observe02') for p in (BASE/name).iterdir() if p.is_file()]
    return g.supervise('audit',[sys.executable,str(Path(__file__).resolve()),'task'],120,128,locks,BASE)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v235 audit closed; follow NEXT_DECISION')
    if len(sys.argv)>1 and sys.argv[1]=='task':audit()
    else:raise SystemExit(0 if main() else 1)
