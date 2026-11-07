"""准确保存场重启核验：原生读写身份、能量重构、边界与保存时间。"""
from pathlib import Path
import json,re,time,shutil
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from openFoamFields import read_state,audit_state
from openFoamConservation import geometry,boundary_values
from openFoamAirThermo import evaluate

BASE=ROOT/'outputs/v200_validation';SOURCE=ROOT/'outputs/v199_validation/development02'


def main():
    start=time.monotonic();out=BASE/'restart_review01';out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');foam=Path('/opt/openfoam13')
    runtime=ROOT/'tools/runtime/openfoam/13-20260624/runtime.json'
    pkg=json.loads(runtime.read_text())
    for p,h in pkg['filesSHA256'].items():assert sha(p)==h,p
    names=['src/thermophysicalModels/basic/basicThermo/BasicThermo.C',
        'src/thermophysicalModels/basic/basicThermo/basicThermo.C',
        'src/thermophysicalModels/basic/psiThermo/PsiThermo.C',
        'src/OpenFOAM/db/Time/Time.C','applications/modules/shockFluid/shockFluid.C',
        'applications/modules/shockFluid/fluxPredictor.C','applications/solvers/foamRun/foamRun.C']
    texts={}
    for name in names:
        p=foam/name;q=out/'source'/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);texts[name]=p.read_text()
    basic=texts[names[0]];energy=basic[basic.index('    he_'):basic.index('    Cp_')]
    assert 'IOobject::NO_READ' in energy and 'this->p_' in energy and 'this->T_' in energy
    psi=texts[names[2]];assert 'thermoMixture.The' in psi and 'calculate();' in psi
    shock=texts[names[4]]
    assert 'IOobject::READ_IF_PRESENT' in shock[shock.index('    rho_'):shock.index('    U_')]
    assert 'IOobject::READ_IF_PRESENT' in shock[shock.index('    phi_'):shock.index('    K(')]
    assert 'fluxPredictor();' in shock and 'phi_ =' in texts[names[5]]
    clock=texts[names[3]]
    assert 'timeDict.lookup<scalar>("deltaT", userUnits())' in clock and 'timeDict.readIfPresent("index", startTimeIndex_)' in clock
    # 模块二：输入来自准确正常保存场；旧运行异常身份保持。
    run=json.loads((ROOT/'outputs/v199_validation/development_run02/report.json').read_text())
    assert not run['passed'] and run['process']['stopReason']=='supervisor-error'
    assert all(run['process']['inputStable'].values()) and not run['process']['liveGroupAfterCleanup']
    fields=run['saved']['0.002']['fieldsSHA256']
    for name,h in fields.items():assert sha(SOURCE/'0.002'/name)==h
    for path,h in run['process']['inputSHA256'].items():
        if path.startswith(str(SOURCE)+'/'):assert sha(path)==h,path
    g=geometry(SOURCE);f=read_state(SOURCE/'0.002',len(g['cells']));state=audit_state(f)
    b={name:boundary_values(SOURCE/'0.002',name,f[name],g) for name in f};boundary=audit_state(b)
    differences={}
    for label,data in [('internal',f),('boundary',b)]:
        gas=evaluate(data['p'],data['T'])
        differences[label]=dict(energyRebuildMaxScaled=float(np.max(abs(gas['internalEnergyJkg']-data['e'])/np.maximum(1,abs(data['e'])))),
            densityEOSMaxRelative=float(np.max(abs(data['rho']/gas['densityKgm3']-1))))
        assert differences[label]['energyRebuildMaxScaled']<2e-8 and differences[label]['densityEOSMaxRelative']<2e-10
    timefile=SOURCE/'0.002/uniform/time';s=timefile.read_text()
    clock_values={k:float(re.search(r'\b'+k+r'\s+([\deE.+-]+);',s).group(1)) for k in ['value','index','deltaT','deltaT0','beginTime']}
    assert clock_values['value']==.002 and clock_values['index']==7889 and clock_values['beginTime']==0
    shutil.copyfile(timefile,out/'source_time.txt')
    report=dict(passed=True,sourceFrame='outputs/v199_validation/development02/0.002',sourceFieldsSHA256=fields,
        sourceTimeSHA256=sha(timefile),sourceGeometrySHA256=sha(SOURCE/'exact_mesh.npz'),state=state,boundaryState=boundary,
        reconstruction= differences,time=clock_values,
        restartIdentity=dict(p='MUST_READ',T='MUST_READ; inverse after e reconstruction',U='MUST_READ',
            rho='READ_IF_PRESENT',e='NO_READ; rebuild sensible internal energy from saved p/T',
            phi='READ_IF_PRESENT then recomputed by constructor fluxPredictor',
            time='saved index/deltaT/deltaT0/beginTime; adjustableRunTime continues schedule'),
        continuousBitwiseEquivalent=False,oldRunPlanAcceptance=False,CFDRuns=0,engineBuilds=0,physicalQualification=False,
        sourceSHA256={str(foam/p):sha(foam/p) for p in names},runtimeSHA256=sha(runtime),elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ['passed','reconstruction','time','elapsedSeconds']}))


if __name__=='__main__':main()
