"""模块职责：独立核对实际HLLC工厂五步控制、相态/物理残差和强制终点，不授箭体资格。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runHLLCNativeV238 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from singleDomainFrame import read_frame
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy
BASE=run.BASE;CASE=run.CASE;DOC=run.DOC;CLOSED=False


def rows(p):
    with p.open() as f:return list(csv.DictReader(f))


def main():
    out=BASE/'native_audit01';out.mkdir(exist_ok=False);air=PhaseAwareViscousAir();errors={}
    process=read(BASE/'process_native_control/report.json');assert process['passed'] and process['inputsStable'] and not process['liveGroupAfterCleanup']
    identity=read(CASE/'native_identity.json');assert identity['actualInteriorBoundaryHLLC'] and identity['sameSAAirIdentity'] and identity['physicalGasID']==15
    summary=read(CASE/'acceptance_summary.json');assert summary['checked']==120 and summary['limited']==0 and summary['forcedFinalNativeOutput']
    completion=read(CASE/'native_completion.json');assert completion['normalReturn'] and completion['updates']==5 and completion['initialSaved']
    stages=rows(CASE/'stages.csv');assert len(stages)==60;index={(int(r['update']),int(r['system']),r['stage']):r for r in stages}
    def check(name,a,b,tol=2e-9):
        a=np.asarray(a);b=np.asarray(b);value=float(np.max(abs(a-b)/np.maximum(1,abs(b))))
        assert np.isfinite(value) and value<=tol,(name,value);errors[name]=max(errors.get(name,0),value)
    frames={}
    for step in (0,1,5):
        a=read_frame(CASE/('state_%06d.bin'%step),24,step==0)['values'];a=a[np.argsort(a[:,0])];frames[step]=a
        v=a[:,6:9]/a[:,5,None];e=a[:,9]/a[:,5]-.5*np.sum(v*v,axis=1)-300000
        t=inverse_energy(air,e);p=a[:,5]*air.R*t;physical=air.evaluate(p,t)
        check('TFromConserved',t,a[:,11]);check('pFromConserved',p,a[:,15]);check('velocity',v,a[:,12:15])
        check('physicalEnthalpy',np.asarray(physical['enthalpyJkg'])+.5*np.sum(v*v,axis=1),a[:,17]-300000)
        tr=air.interpolated_transport(p,t);check('linearMu',tr['dynamicViscosityPas'],a[:,19]);check('linearK',tr['conductivityWmK'],a[:,21])
        mu=np.asarray(tr['dynamicViscosityPas']);chi=a[:,10]*a[:,5]/mu;mut=np.where(a[:,10]>=0,a[:,5]*a[:,10]*chi**3/(chi**3+7.1**3),0)
        check('updatedSAMuT',mut,a[:,23])
        if step:
            rhs=a[:,27:33].copy();rhs[:,4]-=300000*rhs[:,0]
            for j in range(6):
                r=index[(step,0 if j<5 else 2,'after_time')]
                check('physicalRMS',np.sqrt(np.mean(rhs[:,j]**2)),float(r['rms%d'%j]))
                check('physicalSum',np.sum(rhs[:,j]),float(r['sum%d'%j]))
    for step in range(1,6):
        for system in (0,2):
            for j in (range(5) if system==0 else [5]):
                check('implicitSign',float(index[(step,system,'before_time')]['sum%d'%j]),-float(index[(step,system,'after_time')]['sum%d'%j]))
    for r in rows(CASE/'steps.csv'):assert float(r['maxWallSpeedComponent'])==float(r['maxWallSA'])==0
    assert not rows(CASE/'rejections.csv')
    walls=rows(CASE/'wall_stages.csv');assert len(walls)==10
    for r in walls:
        step=int(r['update']);delta=float(index[(step,0,'after_strong')]['sum4'])-float(index[(step,0,'before_strong')]['sum4'])
        check('physicalWallHeat',delta,float(r['effectiveHeatIntoWallW']))
    native=list(CASE.glob('restart_flow*.csv'));assert len(native)==1
    a=np.loadtxt(str(native[0]),delimiter=',',skiprows=1,usecols=range(10));a=a[np.argsort(a[:,0])]
    assert np.array_equal(a[:,0],np.arange(24))
    for left,right in ((a[:,1:4],frames[5][:,1:4]),(a[:,4:10],frames[5][:,5:11])):
        assert np.array_equal(left.copy().view(np.uint64),right.copy().view(np.uint64))
    report=dict(passed=True,schemeID=identity['schemeID'],actualInteriorAndBoundaryFactory=True,completeFlowSAUpdates=5,
        phaseGuardChecks=120,limitedUpdates=0,oldStateRecoveries=0,stageRows=60,wallStageRows=10,
        initialSaved=True,completeFrames=[0,1,5],nativeForceFinalWritten=True,nativeFinalSixStatesAndCoordinatesBitEqual=True,
        nativeFinalSHA256=sha(native[0]),errorMetrics=errors,originalRocketUpdates=0,physicalQualification=False,newQualifiedAreaM2=0,
        generalUnequalStateJacobianExactnessNotClaimed=True,stationaryContactDerivativeNotUnique=True)
    save(out/'report.json',report);print(json.dumps(report))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v238 native audit closed')
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runHLLCNativeV238.py',ROOT/'scripts/runHLLCCompatibilityV238.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
            ROOT/'scripts/singleDomainFrame.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',
            ROOT/'data/thermo/phase_aware_dry_air_v1.json',BASE/'process_native_control/report.json']+list(CASE.iterdir())
        raise SystemExit(0 if run.run.guard.supervise('native_audit',[sys.executable,str(Path(__file__).resolve()),'task'],60,1024,locks,BASE,'audit') else 1)
