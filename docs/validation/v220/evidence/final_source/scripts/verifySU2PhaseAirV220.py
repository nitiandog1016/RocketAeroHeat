"""模块职责：独立审核原生CSV、实际导数及精确构建来源，不启动求解器。"""
from pathlib import Path
import csv,json,sys
import numpy as np
from phaseAwareViscousAir import PhaseAwareViscousAir
from phaseAwareDryAir import pressure_ceiling
from prepareSU2PhaseAirV220 import ROOT,BASE,SOURCE,BUILD,OLD_BUILD,sha,read,save

def main():
    out=BASE/'audit01';out.mkdir(exist_ok=False)
    air=PhaseAwareViscousAir();data=BASE/'control03'
    native=read(data/'native_report.json');assert native['passed']
    rows=list(csv.DictReader((data/'states.csv').open(encoding='utf8')))
    temperatures=[60.,75.123,100.123,132.6312,180.123,226.50908361133003,298.15,300.123,1000.123,1600.123,2000.]
    assert len(rows)==88
    errors={};seen=set()
    for raw in rows:
        row={k:float(v) for k,v in raw.items()};case=int(row['case']);mode=int(row['mode'])
        assert 0<=case<11 and 0<=mode<8 and (case,mode) not in seen;seen.add((case,mode))
        t=temperatures[case];p=.05*float(pressure_ceiling(t));q=air.evaluate(p,t);transport=air.interpolated_transport(p,t)
        rho=float(q['densityKgm3']);cp=float(q['cpJkgK']);cv=cp-air.R
        i=min(len(air.T)-2,int(np.searchsorted(air.T,t,side='right')-1))
        dmu=(air.mu[i+1]-air.mu[i])/(air.T[i+1]-air.T[i]);dk=(air.kappa[i+1]-air.kappa[i])/(air.T[i+1]-air.T[i])
        expected=dict(p=p,T=t,rho=rho,e=float(q['internalEnergyJkg'])+300000,h=float(q['enthalpyJkg'])+300000,
            s=float(q['entropyJkgK']),cp=cp,cv=cv,a=float(q['soundSpeedMps']),mu=float(transport['dynamicViscosityPas']),
            k=float(transport['conductivityWmK']),dpdrho=air.R*t,dpde=rho*air.R/cv,dtde=1/cv,
            dhdrho=-cp*t/rho,dhdp=cp/(rho*air.R),dsdrho=-cp/rho,dsdp=cv/p,dmu=dmu,dk=dk)
        for name,value in expected.items():
            assert np.isfinite(row[name]),name
            # 小输运量也用相对误差，不能被单位量级分母掩盖。
            scale=max(abs(value),1. if name=='s' else 1e-30)
            err=abs(row[name]-value)/scale;errors[name]=max(errors.get(name,0.),err)
            assert err<=2e-10,(case,mode,name,row[name],value,err)
        assert abs(row['T']-t)<=1e-7
    fdrows=list(csv.DictReader((data/'derivatives.csv').open(encoding='utf8')))
    assert len(fdrows)==60
    derivative_errors={};keys=set()
    for row in fdrows:
        key=(float(row['T']),row['name']);assert key not in keys;keys.add(key)
        a=float(row['analytic']);f=float(row['finiteDifference']);assert np.isfinite(a) and np.isfinite(f)
        err=abs(f-a)/max(abs(a),1e-30) if a else abs(f)
        derivative_errors[row['name']]=max(derivative_errors.get(row['name'],0.),err)
        assert err<=3e-6,(row,err)
    assert native['rejections']==10 and native['restoredStates']==9
    assert native['nativeRoeCovarianceError']<=2e-10
    assert not native['actualBoundaryExecutionVerified'] and not native['physicalQualification']
    # 模块二：源依赖、旧298对象和新薄归档逐成员绑定，禁用隐式混包。
    dependencies=read(BASE/'build02/dependencies.json')
    assert all(sha(Path(n))==h for n,h in dependencies.items())
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    assert len(old)==298 and all(sha(OLD_BUILD/n)==h for n,h in old.items())
    members=read(BASE/'build02/members.json');count=0
    for group in members.values():
        for row in group:assert sha(Path(row['path']))==row['sha256'];count+=1
    assert sha(BASE/'build02/SU2_CFD.exe')==read(BASE/'build02/link_report.json')['executableSHA256']
    preserved=read(BASE/'link01/report.json')['executableSHA256'];assert sha(BUILD/'SU2_CFD.exe')==preserved
    processes=[read(p) for p in BASE.glob('process_*/report.json')]
    assert all(r['inputsStable'] and not r['liveGroupAfterCleanup'] for r in processes)
    totals={k:sum(r['elapsedSeconds'] for r in processes if r['budgetKind']==k) for k in ('offline','build','query')}
    assert totals['offline']<=120 and totals['build']<=480 and totals['query']<=60
    report=dict(passed=True,stateRows=88,finiteDifferenceRows=60,independentRelativeErrors=errors,
        native=native,derivativeRelativeErrors=derivative_errors,sourceDependencies=len(dependencies),
        oldObjectsStable=len(old),archiveMembers=count,originalUnverifiedExecutablePreserved=preserved,
        newExecutableSHA256=sha(BASE/'build02/SU2_CFD.exe'),processBudgetSecondsBeforeAudit=totals,
        completeNativeCompatibility=False,CFDAdmission=False,newQualifiedAreaM2=0)
    save(out/'report.json',report)
    print(json.dumps(dict(passed=True,maxPropertyRelativeError=max(errors.values()),
        maxDerivativeRelativeError=max(derivative_errors.values()),stateRows=88)))

if __name__=='__main__':main()
