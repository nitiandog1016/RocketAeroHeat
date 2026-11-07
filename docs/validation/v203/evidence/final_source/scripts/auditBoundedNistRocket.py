"""只读核验实际压力拒绝，并计算原合同内热函数导出的参考压力需求。"""
from pathlib import Path
import json,re,shutil,time
import numpy as np
from scipy.optimize import brentq
from prepareOpenFoamControl import ROOT,sha,save
from nistIntegratedAir import evaluate,R,contract
from readOpenFoamRocketState import mapped_geometry,local_path
from auditOpenFoamDevelopment import log_progress


def main():
    start=time.monotonic();base=ROOT/'outputs/v203_validation';out=base/'rocket_audit01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');case=base/'rocket01'
    run=json.loads((base/'rocket_run01/report.json').read_text());a=json.loads((base/'rocket_run01/admission.json').read_text());process=run['process']
    assert run['profile']==a['profile']=='of13-bounded-nist-original-rocket-sync-v1' and run['inverseProfile']==a['inverseProfile']=='rah-bounded-nist-es-v1'
    assert process['nativeExitCode']==1 and process['stopReason']=='scientific-guard-after-exit'
    assert all(process['inputStable'].values()) and not process['liveGroupAfterCleanup'] and process['logExportPassed']
    assert run['actualSolverDispatchVerified'] and not run['passed'] and not run['saved']
    assert sha(base/'rocket_run01/admission.json')==run['admissionSHA256']
    for p,h in a['inputSHA256'].items():
        path=local_path(p)
        if case in path.parents:assert sha(path)==h
    assert json.loads((case/'thermo_contract.json').read_text())==contract()
    g,m,ids,owner,surface=mapped_geometry(case)
    text=(base/'rocket_run01/native/stdout.txt').read_text();progress=log_progress(text)
    regex=r'RAH_INVALID_THERMO_INPUT e=([\deE.+-]+) p=([\deE.+-]+) seed=([\deE.+-]+)'
    hits=re.findall(regex,text);assert len(hits)==1;es,p,seed=map(float,hits[0])
    assert p>100000 and 100<=seed<=2500 and np.isfinite([es,p,seed]).all()
    # 模块一：用合同内参考压力求同es的数学根；不把它当原异常场温度或生产状态。
    emin=float(evaluate(100000,100)['internalEnergyJkg']);emax=float(evaluate(100000,2500)['internalEnergyJkg'])
    assert emin<es<emax
    conditional=brentq(lambda T:float(evaluate(100000,T)['internalEnergyJkg'])-es,100,2500,xtol=1e-10)
    co=abs(run['initialMaxCo']/a['initialMaxCoExpected']-1);assert co<5e-6
    assert progress['announcedSteps']==progress['completedLogSteps']+1
    # 模块二：真实来流h0及等熵驻点的压力需求，仅为合同设计参照，不为数值场授资格。
    c=a['condition'];up=evaluate(c['pressurePa'],c['temperatureK']);U=np.array(a['velocityMps'])
    h0=float(up['enthalpyJkg'])+float(U@U)/2
    T0=brentq(lambda T:float(evaluate(100000,T)['enthalpyJkg'])-h0,100,2500,xtol=1e-10)
    p0=100000*float(np.exp((float(evaluate(100000,T0)['entropyJkgK'])-float(up['entropyJkgK']))/R))
    source_phi=dict(pressureReferencePa=100000.,source='ideal fixed-composition entropy relation at valid reference pressure',
        upstreamTotalEnthalpyJkg=h0,isentropicStagnationTemperatureK=T0,isentropicStagnationPressurePa=p0,
        pressureWithinCurrentProductionContract=bool(p0<=100000),actualShockStagnationPressure=False,
        scope='analytic pressure-coverage design reference; not an accepted out-of-contract state')
    result=dict(passed=True,candidatePassed=False,nativeExitCode=1,nativeFatalPrecededSupervisorObservation=True,
        supervisorTerminatedRunningSolver=False,progress=progress,rejectedInput=dict(internalEnergyJkg=es,pressurePa=p,seedTemperatureK=seed),
        pressureContractUpperPa=100000,pressureExcessPa=p-100000,rejectedByPressureBeforeRootSolve=True,
        inputEnergyInsideSupportedTemperatureBracket=True,energyBracketJkg=[emin,emax],
        conditionalTemperatureFromEnergyAtReferencePressureK=conditional,conditionalTemperatureIsNotSavedField=True,
        thermodynamicReference=source_phi,initialCoRelativeError=co,originalWallFaces=len(ids),originalAreaM2=float(g['area'][ids].sum()),
        originalWallExact=True,offendingCellUnknown=True,completeSavedDevelopmentFields=0,fullFieldMinimumTemperatureK=None,
        priorV202TemperatureFailureRootCauseEstablished=False,process=process,
        physicalQualification=False,numericallyConverged=False,newQualifiedHeatAreaM2=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','progress','rejectedInput','conditionalTemperatureFromEnergyAtReferencePressureK','thermodynamicReference','elapsedSeconds']}))


if __name__=='__main__':main()
