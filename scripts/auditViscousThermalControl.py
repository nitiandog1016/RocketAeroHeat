"""只读热控制审核：分开进程完成、原生热流接口和独立物理解的误差门槛。"""
from pathlib import Path
import json,time,shutil,sys
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from phaseAwareViscousAir import PhaseAwareViscousAir
from openFoamFields import read_state,internal_field
from readOpenFoamWallHeat import read_wall_heat

def main(case_name='thermal02',run_name='thermal_run01',out_name='thermal_audit01'):
    start=time.monotonic();base=ROOT/'outputs/v207_validation';out=base/out_name;out.mkdir(exist_ok=False)
    case=base/case_name;run=json.loads((base/run_name/'report.json').read_text());oracle=json.loads((base/'period_repair01/oracle.json').read_text())
    if not run['passed']:raise ValueError('过程未完整，不进行正向热控制验收')
    for n in ['auditViscousThermalControl.py','readOpenFoamWallHeat.py','phaseAwareViscousAir.py']:
        shutil.copyfile(ROOT/'scripts'/n,out/n)
    gas=PhaseAwareViscousAir();T0=np.array(oracle['temperatureK']);U0=np.array(oracle['velocityMps']);dy=.001/32;area=oracle['wallAreaM2']
    initial_q={name:float(gas.transport(100000,tw)['conductivityWmK'])*(T0[i]-tw)/(dy/2) for name,tw,i in [('cold',300,0),('hot',600,-1)]}
    # 模块一：已知初态的原生一阶面质量耗散；仅计算残差，不推进任何状态。
    rho=gas.evaluate(100000,T0)['densityKgm3'];a=gas.evaluate(100000,T0)['soundSpeedMps']
    initial_mass=.5*np.maximum(a[:-1],a[1:])*(rho[:-1]-rho[1:])*area
    records=[]
    for name in sorted(run['saved'],key=float):
        f=read_state(case/name,32);wall=read_wall_heat(case,name,32)
        q={r['patch']:r['heatFluxIntoWallWm2'] for r in wall['rows']}
        qerrors={n:abs(q[n]/oracle['wallHeatIntoWallWm2'][n]-1) for n in q}
        terr=float(np.max(abs(f['T']/T0-1)));uerr=float(np.max(abs(f['U']-U0))/100)
        reconstructed={n:float(gas.interpolated_transport(100000,tw)['conductivityWmK'])*(f['T'][i]-tw)/(dy/2) for n,tw,i in [('cold',300,0),('hot',600,-1)]}
        native_error=max(abs(q[n]/reconstructed[n]-1) for n in q)
        if native_error>2e-10 or any(r['givenWallTemperatureK']!=({'cold':300,'hot':600}[r['patch']]) for r in wall['rows']):raise ValueError('实际热流符号/梯度接口不符')
        mass=internal_field(case/name/'phi',31)/area
        rows=dict(time=name,heatFluxIntoWallWm2=q,heatFluxRelativeError=qerrors,
            temperatureRelativeError=terr,velocityScaledError=uerr,rawFourierMaxRelative=native_error,
            nativeInternalNormalMassFluxMinMaxKgm2s=[float(mass.min()),float(mass.max())],
            wallPowerIntoWallW=sum(q.values())*area,referenceViscousPowerW=oracle['viscousPowerPerAreaWm2']*area,
            apparentSteadyHeatPowerDefectW=(sum(q.values())-oracle['viscousPowerPerAreaWm2'])*area,
            powerDefectIsNotFullTransientConservationResidual=True,fields=wall)
        records.append(rows)
    last=records[-1];passed=max(last['heatFluxRelativeError'].values())<=.02 and last['temperatureRelativeError']<=.005 and last['velocityScaledError']<=.005
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    unchanged=all(sha(ROOT/n.replace('\\','/'))==h for n,h in protected.items())
    if not unchanged:raise ValueError('保护文件改变')
    result=dict(passed=bool(passed),auditCompleted=True,processPassed=True,nativeSignAndGradientInterfacePassed=True,
        independentReferencePassed=bool(passed),initialDiscreteHeatFluxIntoWallWm2=initial_q,
        initialFirstOrderNormalMassFluxMinMaxKgm2s=[float(initial_mass.min()/area),float(initial_mass.max()/area)],
        exactReferenceNormalMassFlux=0,initialMassFluxCalculationIsNotUniqueCauseProof=True,
        records=records,perStepFinalStateGuardVerified=run['perStepGuardVerified'],protectedFilesUnchanged=29,
        physicalQualification=False,numericallyConverged=False,newQualifiedHeatAreaM2=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps(dict(passed=bool(passed),qError=last['heatFluxRelativeError'],Terror=last['temperatureRelativeError'],Uerror=last['velocityScaledError'],seconds=result['elapsedSeconds'])))
    return passed

if __name__=='__main__':raise SystemExit(0 if main(*sys.argv[1:]) else 1)
