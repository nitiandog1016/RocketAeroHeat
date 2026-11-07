"""固定原生温度反算的源审计与一次已知真解区分；不运行CFD。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from nistIntegratedAir import evaluate,AIR


def main():
    start=time.monotonic();out=ROOT/'outputs/v203_validation/source01';out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');native=Path('/opt/openfoam13')
    names=['src/thermophysicalModels/specie/thermo/thermo/thermoI.H',
        'src/thermophysicalModels/specie/thermo/thermo/thermo.H',
        'src/thermophysicalModels/specie/thermo/thermo/thermo.C',
        'src/thermophysicalModels/specie/thermo/hIcoTabulated/hIcoTabulatedThermoI.H',
        'src/thermophysicalModels/specie/thermo/sensibleInternalEnergy/sensibleInternalEnergy.H',
        'applications/modules/shockFluid/thermophysicalPredictor.C',
        'src/thermophysicalModels/specie/include/forThermo.H']
    matched=[]
    for p in (native/'src/thermophysicalModels/basic').rglob('*'):
        if 'lnInclude' in p.parts or p.suffix not in ['.C','.H']:continue
        text=p.read_text()
        if '.The(' in text or '.THE(' in text or '::calculate' in text:
            name=p.relative_to(native).as_posix();names.append(name)
            matched.extend([dict(path=name,line=i,text=s) for i,s in enumerate(text.splitlines(),1) if any(t in s for t in ['.The(','.THE(','calculate','correct()'])])
    rows=[]
    for name in sorted(set(names)):
        p=native/name;dest=out/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        rows.append(dict(path=name,sha256=sha(p)))
    boost=Path('/usr/include/boost/math/tools/toms748_solve.hpp')
    if boost.exists():
        dest=out/'toms748_solve.hpp';shutil.copyfile(boost,dest);rows.append(dict(path=str(boost),sha256=sha(boost)))
    # 模块二：只在有数据的100—2500K内给定真解；第一次试探出界即停止，不外推。
    targets=np.array([100,100.001,100.01,100.1,101,110,180,226.50908361133003,500,1000,2499.99,2500.])
    seeds=np.array([100,100.1,101,110,150,180,226.50908361133003,300,500,1000,2000,2500.])
    cases=[]
    for trueT in targets:
        target=float(evaluate(1197.0279007383137,trueT)['internalEnergyJkg'])
        for seed in seeds:
            t=float(seed);history=[];kind='iteration-limit'
            for iteration in range(102):
                g=evaluate(1197.0279007383137,t);new=t-(float(g['internalEnergyJkg'])-target)/float(g['cvJkgK'])
                history.append(dict(iteration=iteration,temperatureK=t,energyJkg=float(g['internalEnergyJkg']),trialTemperatureK=new))
                if new<100 or new>2500:kind='trial-out-of-range';break
                if abs(new-t)<=seed*1e-4:kind='returned';break
                t=new
            cases.append(dict(trueTemperatureK=float(trueT),targetEnergyJkg=target,seedK=float(seed),result=kind,history=history))
    result=dict(passed=True,files=rows,thermoCallMatches=matched,cases=cases,
        trialOutOfRangeCount=sum(r['result']=='trial-out-of-range' for r in cases),
        totalCases=len(cases),boostAvailable=boost.is_file(),exactRocketRootCauseEstablished=False,CFDRuns=0,
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result)
    print(json.dumps({k:result[k] for k in ['passed','trialOutOfRangeCount','totalCases','boostAvailable','thermoCallMatches','elapsedSeconds']}))


if __name__=='__main__':main()
