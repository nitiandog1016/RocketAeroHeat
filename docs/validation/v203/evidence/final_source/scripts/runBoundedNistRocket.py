"""成熟有界反算的唯一真实全壁候选；物理和原几何保持。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from nistIntegratedAir import AIR,MODEL
from runOpenFoamNistRocket import Guard,prerequisites
from runOpenFoamRocketDevelopment import local_path

BASE=ROOT/'outputs/v203_validation';SOURCE=ROOT/'outputs/v202_validation/rocket01'
PROFILE='of13-bounded-nist-original-rocket-sync-v1';INVERSE='rah-bounded-nist-es-v1'


def main():
    raise ValueError('v203真实新映射候选触发既定压力上限，已关闭；不直接放宽上限续算或补场重跑。')
    from runOpenFoamProcess import run
    start=time.monotonic();old_evidence,_,_=prerequisites()
    br=json.loads((BASE/'bounded01/report.json').read_text());dispatch=json.loads((BASE/'dispatch01/report.json').read_text())
    assert br['passed'] and dispatch['passed'] and dispatch['actualSolverDispatchVerified']
    library=BASE/'bounded01/library/librahBoundedNist.so';assert sha(library)==br['registrationLibrarySHA256']
    assert (BASE/'bounded01/mixture.dict').read_text()==AIR.dictionary()
    old=json.loads((ROOT/'outputs/v202_validation/rocket_run01/admission.json').read_text())
    for p,h in old['inputSHA256'].items():
        path=local_path(p)
        if SOURCE in path.parents:assert sha(path)==h
    out=BASE/'rocket_run01';out.mkdir(exist_ok=False);case=BASE/'rocket01';case.mkdir(exist_ok=False)
    for folder in ['constant','system','0']:shutil.copytree(SOURCE/folder,case/folder)
    for name in ['exact_mesh.npz','wall_mapping.json','patches.json','thermo_contract.json']:shutil.copyfile(SOURCE/name,case/name)
    props=case/'constant/physicalProperties';text=props.read_text();assert text.count('energy sensibleInternalEnergy;')==1
    props.write_text(text.replace('energy sensibleInternalEnergy;','energy rahBoundedNistInternalEnergy;'))
    control=case/'system/controlDict';text=control.read_text();oldlib=(ROOT/'outputs/v202_validation/native01/library/librahNistAir.so').as_posix()
    assert text.count(oldlib)==1;control.write_text(text.replace(oldlib,library.as_posix()))
    for name in ['0/p','0/T','0/U','system/fvSchemes','system/fvSolution','exact_mesh.npz','wall_mapping.json']:
        assert sha(SOURCE/name)==sha(case/name)
    names=['runBoundedNistRocket.py','runOpenFoamNistRocket.py','runOpenFoamRocketDevelopment.py','nistIntegratedAir.py','openFoamFields.py',
        'openFoamSynchronousControl.py','runOpenFoamProcess.py','openFoamConservation.py']
    for name in names:shutil.copyfile(ROOT/'scripts'/name,out/name)
    locks=[*old_evidence,library,BASE/'bounded01/report.json',BASE/'dispatch01/report.json',ROOT/'data/thermo/nist_air_shomate_v1.json',
        ROOT/'docs/validation/v203/ROCKET_PLAN.md',*[ROOT/'scripts'/name for name in names],*[p for p in case.rglob('*') if p.is_file()]]
    admission=dict(profile=PROFILE,inverseProfile=INVERSE,model=MODEL,case=case.relative_to(ROOT).as_posix(),
        condition=old['condition'],velocityMps=old['velocityMps'],initialMaxCoExpected=old['initialMaxCoExpected'],
        endTimeSeconds=.00045,saveTimesSeconds=old['saveTimesSeconds'],maximumSteps=40000,maximumSeconds=900,parentSeconds=960,
        memoryMiB=2048,outputMiB=160,automaticContinuation=0,registrationLibrarySHA256=sha(library),
        inputSHA256={str(p):sha(p) for p in locks},physicalQualification=False)
    save(out/'admission.json',admission);guard=Guard(case,out/'saved_guard.json',old['initialMaxCoExpected'])
    process=run(['foamRun'],case,out/'native',900,[*locks,out/'admission.json'],160,2048,guard)
    times=sorted(float(k) for k in guard.checked);log=(out/'native/stdout.txt').read_text()
    complete=len(times)==6 and np.allclose(times,old['saveTimesSeconds'],rtol=0,atol=1e-13)
    dispatched='rahBoundedNistInternalEnergy' in log and 'hIcoTabulated' in log
    result=dict(passed=bool(process['passed'] and complete and dispatched and 'End\n' in log),
        profile=PROFILE,inverseProfile=INVERSE,model=MODEL,process=process,saved=guard.checked,
        actualSolverDispatchVerified=dispatched,initialMaxCo=guard.initial_co,admissionSHA256=sha(out/'admission.json'),
        expectedSavedTimesComplete=bool(complete),numericallyConverged=False,attachedEdgeQualified=False,
        physicalQualification=False,newQualifiedHeatAreaM2=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps(dict(passed=result['passed'],seconds=result['elapsedSeconds'],stop=process['stopReason'],saved=times)))


if __name__=='__main__':main()
