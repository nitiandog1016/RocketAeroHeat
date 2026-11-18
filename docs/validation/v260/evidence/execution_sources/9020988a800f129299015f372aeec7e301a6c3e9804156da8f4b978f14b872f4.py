"""模块职责：固定样本、同源薄构建和唯一量热helper查询；不运行完整CFD。"""
import sys,shutil,subprocess
from pathlib import Path
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from terminalCapacityV249 import supervise
from auditFullResidualV256 import matrix
BASE=ROOT/'outputs/v260_validation';DOC=ROOT/'docs/validation/v260';CASE=BASE/'case01';STAGE=BASE/'stage01';BUILD=BASE/'build01';CLOSED=False
PREV=ROOT/'outputs/v258_validation';CORE=ROOT/'outputs/v220_validation/source/SU2_CFD/include/fluid'
NATIVE=ROOT/'native/su2PhaseAwareAir';DRIVER=NATIVE/'caloricReplayV260.cpp';WD=ROOT/'tools/runtime/su2-thermal-build/build-v117'
OWN=[DRIVER,NATIVE/'StableCaloricInverseV260.hpp',NATIVE/'CSU2PhaseAwareAirStableV260.hpp']


def prepare():
    for p in [CASE,STAGE,BUILD]:p.mkdir(exist_ok=False)
    expected=read(PREV/'build01/dependencies.json')
    for name in ['PhaseAwareAir.hpp','CSU2PhaseAwareAir.hpp']:
        assert sha(CORE/name)==expected[str(CORE/name)] and sha(NATIVE/name)==sha(CORE/name),name
    assert sha(CORE/'PhaseAirTable.hpp')=='d0888e44ee865cb7b3a74c3135e69898e0bc21673cc2ca0cad558b60e05e102c'
    s=matrix(PREV/'rocket01/initial_state.bin',b'RAH256S1',19);r0=matrix(PREV/'rocket01/initial_full_residual.bin',b'RAH256R1',7)
    full=matrix(PREV/'rocket01/sigma_full_residual.bin',b'RAH256R1',7);half=matrix(PREV/'rocket01/half_full_residual.bin',b'RAH256R1',7)
    step=read(PREV/'rocket01/step.json')['sigma'];scale=np.genfromtxt(ROOT/'outputs/v257_validation/rocket01/scaling.csv',delimiter=',',names=True)['residualCharacteristic']
    difference=((full[:,1:]-r0[:,1:])/step-(half[:,1:]-r0[:,1:])/(step/2))/scale;score=np.linalg.norm(difference,axis=1)
    order=np.lexsort((s[:,0],-score));selected=list(order[:64]);addedWall=False
    if not np.any(s[selected,17]==1):selected.append(next(i for i in order if s[i,17]==1));addedWall=True
    rows=[];identities=[]
    for stage,label in enumerate(['initial','sigma','half']):
        state=matrix(PREV/'rocket01'/(label+'_state.bin'),b'RAH256S1',19)
        for i in selected:
            rows.append([len(rows),state[i,0],stage]+state[i,5:11].tolist())
            identities.append(dict(sample=len(rows)-1,globalID=int(state[i,0]),sourceFrame=label,sourceRow=int(i),wall=bool(state[i,17]),selectionScore=float(score[i])))
    np.savetxt(CASE/'states_input.csv',np.array(rows),delimiter=',',header='sample,globalID,stage,rho,momX,momY,momZ,physicalRhoE,nuTilde',comments='',fmt=['%d','%d','%d']+['%.17g']*6)
    save(CASE/'input_identity.json',dict(samples=identities,selectedNodes=len(selected),worstWallAdded=addedWall,selectedBeforeNewAlgorithmRun=True,
        savedPhysicalEnergyReconverted=True,exactReplayOfUnstoredOldPrimitiveTargetsClaimed=False))
    command=[]
    for part in read(PREV/'stage01/command.json'):
        if 'sundials-7.4.0' in part:continue
        if part==str(NATIVE/'directionConsistencyV258.cpp'):part=str(DRIVER)
        if part==str(PREV/'build01/preflight.exe'):part=str(BUILD/'caloric.exe')
        if part==str(PREV/'build01/preflight.d'):part=str(BUILD/'caloric.d')
        command.append(part)
    save(STAGE/'command.json',command)
    save(STAGE/'report.json',dict(passed=True,selectedNodes=len(selected),samples=len(rows),selection='top64 original scaled sigma-minus-half response, plus worst wall if absent',
        oldAndMaintainedHeadersIdentical=True,tableSHA256=sha(CORE/'PhaseAirTable.hpp'),CFDAdmitted=False,physicalQualification=False))


def build():
    subprocess.check_call(read(STAGE/'command.json'),cwd=WD)
    deps={}
    for name in (BUILD/'caloric.d').read_text(encoding='utf8').replace('\\\n',' ').split(': ',1)[1].split():
        p=Path(name);p=(p if p.is_absolute() else WD/p).resolve();deps[str(p)]=sha(p)
    save(BUILD/'dependencies.json',deps);save(BUILD/'report.json',dict(passed=True,engineSHA256=sha(BUILD/'caloric.exe'),physicsObjectsRecompiled=0,CFDExecuted=False))


def main(stage):
    if CLOSED:raise SystemExit('v260 fixed caloric helper tasks closed; follow NEXT_DECISION')
    assert stage in ['prepare','build','build2','query'];BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    measure=guard.size_bytes;guard.size_bytes=lambda p:measure(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),DOC/'PLAN.md',ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',ROOT/'scripts/auditFullResidualV256.py']+OWN
    locks += [NATIVE/'PhaseAwareAir.hpp',NATIVE/'CSU2PhaseAwareAir.hpp',CORE/'PhaseAwareAir.hpp',CORE/'CSU2PhaseAwareAir.hpp',CORE/'PhaseAirTable.hpp',
        ROOT/'data/thermo/phase_aware_dry_air_v1.json',PREV/'build01/dependencies.json',PREV/'stage01/command.json',PREV/'rocket01/step.json',ROOT/'outputs/v257_validation/rocket01/scaling.csv']
    for label in ['initial','sigma','half']:locks += [PREV/'rocket01'/(label+'_state.bin'),PREV/'rocket01'/(label+'_full_residual.bin')]
    deps=read(PREV/'build01/dependencies.json')
    for p,h in deps.items():assert sha(p)==h,p;locks.append(Path(p))
    old=ROOT/'outputs/v257_validation';index=read(ROOT/'docs/validation/v257/evidence/local_asset_index.json')
    locks += [old/'build01/members.json',ROOT/'docs/validation/v257/evidence/local_asset_index.json']
    for name,rows in read(old/'build01/members.json').items():
        archive=old/'build01'/name;assert sha(archive)==index['build01/'+name]['sha256'];locks.append(archive)
        for row in rows:assert sha(row['path'])==row['sha256'];locks.append(Path(row['path']))
    if stage=='prepare':argv=[sys.executable,str(Path(__file__).resolve()),'task_prepare'];seconds=60;kind='offline';cwd=ROOT
    else:
        assert read(STAGE/'report.json')['passed'];locks += [STAGE/'report.json',STAGE/'command.json',CASE/'states_input.csv',CASE/'input_identity.json']
        if stage in ['build','build2']:
            used=0
            if stage=='build2':
                previous=read(BASE/'process_build/report.json');assert not previous['passed'];used=previous['elapsedSeconds']
            argv=[sys.executable,str(Path(__file__).resolve()),'task_build'];seconds=90-used;kind='build';cwd=ROOT
        else:
            assert read(BUILD/'report.json')['passed'] and sha(BUILD/'caloric.exe')==read(BUILD/'report.json')['engineSHA256']
            for p,h in read(BUILD/'dependencies.json').items():assert sha(p)==h,p;locks.append(Path(p))
            locks += [BUILD/'report.json',BUILD/'dependencies.json',BUILD/'caloric.exe'];argv=[str(BUILD/'caloric.exe'),str(CASE)];seconds=30;kind='query';cwd=ROOT
    snap=BASE/('source_'+stage);snap.mkdir(exist_ok=False)
    for p in locks:
        if p in OWN+[Path(__file__).resolve(),DOC/'PLAN.md']:shutil.copyfile(p,snap/(sha(p)+p.suffix))
    return supervise(guard,BASE,stage,argv,seconds,64,list(dict.fromkeys(locks)),cwd,kind)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v260 fixed caloric helper tasks closed; follow NEXT_DECISION')
    if sys.argv[1]=='task_prepare':prepare()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
