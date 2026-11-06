"""模块用途。通过有限合同和小例后，只运行一项同源整箭限步候选。"""
import argparse
import copy
import json
from pathlib import Path
import sys
import numpy as np
from auditSU2IdealAirState import audit,load_state
from prepareSU2IdealAirRestart import prepare
from su2IdealAirInterface import STEP_FLAVOR,validate
from su2ThermalRuntime import ROOT,read,resolve,sha
from validateIdealAir import run,save

BASE=ROOT/'outputs/v165_validation'
CONTROL=ROOT/'outputs/v164_validation'
INITIAL=ROOT/'outputs/embedded_workbench/24915f33-5728-439b-884c-17a471f21457/ideal_air_initial_state'


def launch(folder,seconds,parent_seconds):
    names=['scripts/runSU2ConservativeStepStudy.py','scripts/prepareSU2IdealAirRestart.py',
        'scripts/su2IdealAirInterface.py','scripts/su2IdealAirCheckpoints.py','scripts/auditSU2IdealAirState.py',
        'scripts/runEmbeddedSU2.py','scripts/su2ThermalRuntime.py','scripts/validateIdealAir.py',
        'docs/validation/v165/PLAN.md']
    for name in names:
        target=folder/'executed_sources'/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((ROOT/name).read_bytes())
    save(folder/'engine_request.json',dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',
        maximumSeconds=seconds,threads=2,sourceFiles=names,jobToken='v165-'+folder.name))
    run(folder,'parent',[sys.executable,str(ROOT/'scripts/runEmbeddedSU2.py'),'--request',str(folder/'engine_request.json')],parent_seconds)


def small():
    """小例只验证接口和合法状态保持，不能授予真实箭体物理资格。"""
    folder=BASE/'small01'
    prepare(INITIAL,folder,4,2,32*1024**2,snapshot=ROOT/'outputs/v163_validation/initial_source',runtime_flavor=STEP_FLAVOR)
    text=(folder/'case.cfg').read_text('utf8')
    _,runtime,flavor,_,_=resolve(folder,text)
    assert flavor==STEP_FLAVOR
    manifest=read(folder/'case_manifest.json');rejected=[]
    for name in ('same-numerics','wrong-parent','missing-contract','unbounded','separate-directions','clip','conservation-claim','wrong-stage'):
        m,r=copy.deepcopy(manifest),copy.deepcopy(runtime)
        if name=='same-numerics':r['samePhysicsAndNumerics']=True
        elif name=='wrong-parent':r['numericalParentFlavor']='rah-ideal-air-phase-v2'
        elif name=='missing-contract':m.pop('runtimeConservativeStep')
        elif name=='unbounded':r['conservativeStep']['maximumHalvings']=64
        elif name=='separate-directions':r['conservativeStep']['direction']='energy-only'
        elif name=='clip':r['conservativeStep']['temperatureClipping']=True
        elif name=='conservation-claim':r['conservativeStep']['globalConservationProven']=True
        else:r['conservativeStep']['currentK']='previous-output-k'
        try:validate(folder,m,text,r)
        except (AssertionError,KeyError):rejected.append(name)
    assert len(rejected)==8
    save(folder/'contract_rejections.json',dict(passed=True,rejected=rejected))
    for name in ('case.cfg','initial_flow.csv','mesh.su2'):
        assert sha(folder/name)==sha(CONTROL/'small01'/name)
    launch(folder,30,120)
    assert not (folder/'conservative_steps.csv').exists(), '合法小例意外发生限步'
    speed=float(np.linalg.norm(manifest['expectedVelocity']))
    direct,scaled,details={},{},{}
    for step in (0,2,3):
        name='restart_flow_%06d.csv'%step
        old,a,nd=load_state(CONTROL/'small01'/name,True,captured=True)
        new,b,nd_new=load_state(folder/name,True,captured=True)
        assert list(a)==list(b) and nd==nd_new==3 and old.shape==new.shape
        absolute=abs(new-old);denominator=np.maximum(1,abs(old))
        direct[str(step)]=float((absolute/denominator).max())
        for index,field in enumerate(a):
            if field in ('Velocity_x','Velocity_y','Velocity_z'):
                denominator[:,index]=np.maximum(denominator[:,index],speed)
        scaled[str(step)]=float((absolute/denominator).max())
        details[str(step)]={field:float(absolute[:,i].max()) for i,field in enumerate(a)}
        assert scaled[str(step)]<2e-13
    thermal=audit(folder)
    save(BASE/'small_report.json',dict(passed=True,rawAllColumnErrors=direct,scaledFieldErrors=scaled,
        perFieldAbsoluteErrors=details,velocityScale=speed,thermalErrors=thermal['errors'],
        contractRejections=rejected,limitingEvents=0,physicsQualified=False,newQualifiedArea=0))
    print(json.dumps(dict(smallPassed=True,raw=direct,scaled=scaled,thermal=thermal['errors'])))


def real():
    """真实候选只复用源149，原生数值配置不变，达到50次即停止。"""
    assert read(BASE/'small_report.json')['passed']
    folder=BASE/'real01'
    prepare(ROOT/'outputs/v160_validation/real01',folder,50,50,1024**3,
            checkpoint=149,linear_history=True,runtime_flavor=STEP_FLAVOR)
    same={name:sha(folder/name)==sha(CONTROL/'real01'/name) for name in ('initial_flow.csv','case.cfg','mesh.su2')}
    assert all(same.values())
    save(folder/'same_source_preflight.json',dict(passed=True,files=same,
        controlDirectory=str(CONTROL/'real01'),planSHA256=sha(ROOT/'docs/validation/v165/PLAN.md'),
        controlSHA256={str(CONTROL/'real01'/name):sha(CONTROL/'real01'/name) for name in same},
        samePhysicsAndNumerics=False,automaticContinuation=False))
    launch(folder,300,420)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',choices=('small','real'))
    args=parser.parse_args()
    small() if args.case=='small' else real()
