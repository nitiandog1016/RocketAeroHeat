"""模块用途。受限SST更新的小型同源对照及唯一有预算整箭候选，不自动追加。"""
import argparse
import copy
import json
from pathlib import Path
import sys
import numpy as np
from auditSU2IdealAirState import audit,load_state
from prepareSU2IdealAirRestart import prepare,validate_source
from runEmbeddedSU2 import validate_resource_budget
from su2IdealAirInterface import SST_ENERGY_FLAVOR,settings
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import run,save

BASE=ROOT/'outputs/v167_validation'


def launch(folder,seconds,parent_seconds):
    names=['scripts/runSU2SSTEnergyStudy.py','scripts/prepareSU2IdealAirRestart.py',
        'scripts/su2IdealAirInterface.py','scripts/su2IdealAirCheckpoints.py','scripts/auditSU2IdealAirState.py',
        'scripts/runEmbeddedSU2.py','scripts/su2ThermalRuntime.py','scripts/validateIdealAir.py',
        'docs/validation/v167/PLAN.md']
    for n in names:
        p=folder/'executed_sources'/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((ROOT/n).read_bytes())
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',
        maximumSeconds=seconds,threads=2,sourceFiles=names,jobToken='v167-'+folder.name)
    if seconds==900:request['resourceProfile']='v167-single-sst-energy'
    save(folder/'engine_request.json',request)
    run(folder,'parent',[sys.executable,str(ROOT/'scripts/runEmbeddedSU2.py'),'--request',str(folder/'engine_request.json')],parent_seconds)


def small():
    assert all(read(BASE/'build01/source_restore.json').values())
    folder=BASE/'small01'
    m=prepare(ROOT/'outputs/v164_validation/small01',folder,3,3,32*1024**2,checkpoint=3,runtime_flavor=SST_ENERGY_FLAVOR)
    assert not m['initialization']['samePhysicsAndNumerics'] and m['initialization']['priorSuccessfulUpdates']==4
    # 用内存反例拒绝伪称同算法、错误温区和错误来源，不修改已锁定文件。
    rejected=[]
    for name in ('unchanged-operator','wrong-temperature-range','wrong-source'):
        altered=copy.deepcopy(m)
        if name=='unchanged-operator':altered['initialization']['samePhysicsAndNumerics']=True
        elif name=='wrong-temperature-range':altered['runtimeSSTEnergy']['temperatureRangeK']=[50,2000]
        else:altered['initialization']['sourceRuntimeFlavor']='rah-ideal-air-step-v1'
        try:validate_source(folder,settings((folder/'case.cfg').read_text('utf8')),altered,altered['idealAir'])
        except AssertionError:rejected.append(name)
    assert len(rejected)==3
    save(folder/'contract_checks.json',dict(passed=True,rejected=rejected))
    launch(folder,30,120)
    r=read(folder/'engine_process.json');assert r['exitCode']==0 and all(r['inputsStable'].values())
    assert [e['localIteration'] for e in r['checkpoints']['complete']]==[0,2]
    assert not (folder/'sst_energy_updates.csv').exists() and not (folder/'admissibility.csv').exists()
    errors={};locks={};thermal={};speed=np.linalg.norm(m['expectedVelocity'])
    for step in (0,2):
        name='restart_flow_%06d.csv'%step
        oldpath=ROOT/'outputs/v166_validation/small01'/name;newpath=folder/name
        old,fields,nd=load_state(oldpath,True,captured=True);new,newfields,newnd=load_state(newpath,True,captured=True)
        assert list(fields)==list(newfields) and nd==newnd==3
        scale=np.maximum(1,abs(old))
        for j,n in enumerate(fields):
            if n in ('Velocity_x','Velocity_y','Velocity_z'):scale[:,j]=np.maximum(scale[:,j],speed)
        errors[str(step)]=float(np.max(abs(new-old)/scale));assert errors[str(step)]<=2e-13
        for p in (oldpath,newpath):locks[str(p)]=sha(p)
        thermal[str(step)]=audit(folder,folder/('audit%03d'%step),checkpoint=step if step==0 else None)['errors']
    save(BASE/'small_report.json',dict(passed=True,scaledFieldErrors=errors,thermalErrors=thermal,
        limitedUpdates=0,sourceSHA256=locks,contractRejections=rejected,physicsQualified=False,newQualifiedArea=0))
    print(json.dumps(read(BASE/'small_report.json')))


def real():
    assert read(BASE/'small_report.json')['passed']
    folder=BASE/'real01';baseline=ROOT/'outputs/v166_validation/real02'
    m=prepare(ROOT/'outputs/v164_validation/real01',folder,200,100,1536*1024**2,checkpoint=49,runtime_flavor=SST_ENERGY_FLAVOR)
    assert (m['nodes'],m['cells'])==(335726,705490) and m['initialization']['priorSuccessfulUpdates']==200
    assert all(sha(folder/n)==sha(baseline/n) for n in ('case.cfg','mesh.su2','initial_flow.csv'))
    request=dict(maximumSeconds=900,threads=2,resourceProfile='v167-single-sst-energy')
    validate_resource_budget(request,m);rejected=[]
    for key,value in (('maximumSeconds',901),('threads',3),('resourceProfile','unknown')):
        bad=dict(request);bad[key]=value
        try:validate_resource_budget(bad,m)
        except AssertionError:rejected.append(key)
    assert len(rejected)==3
    save(folder/'preflight.json',dict(passed=True,baselineDirectory=str(baseline),
        inputsIdentical=['case.cfg','mesh.su2','initial_flow.csv'],resourceRejections=rejected,
        expectedLocalCheckpoints=[0,100,199],maximumUpdates=200,maximumNativeSeconds=900,
        maximumParentSeconds=1020,maximumThreads=2,maximumMemoryBytes=4*1024**3,
        planSHA256=sha(ROOT/'docs/validation/v167/PLAN.md'),automaticAdditionalTasks=0,
        numericalOperatorChanged=True,uninterruptedEquivalent=False,physicsQualified=False))
    launch(folder,900,1020)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('case',choices=('small','real'));a=p.parse_args()
    small() if a.case=='small' else real()
