"""有限原生温变输运/热壁控制：隔离构建、32单元网格和唯一8步，各阶段只启动一次。"""
from pathlib import Path
import json,shutil,time,re,sys
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from phaseAwareViscousAir import PhaseAwareViscousAir
from openFoamFields import audit_state
from openFoamSynchronousControl import PhysicalTimeGuard

BASE=ROOT/'outputs/v207_validation';CASE=BASE/'thermal02'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def native_source_locks():
    report=read(BASE/'source01/report.json')
    if not report['passed']:raise ValueError('准确来源未过')
    paths=[]
    for row in report['sources']:
        path=Path(row['path'])
        if sha(path)!=row['sha256']:raise ValueError('原生来源改变')
        paths.append(path)
    return paths
def locks():
    return [Path(__file__).resolve(),ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',
        ROOT/'docs/validation/v207/THERMAL_CONTROL_PLAN.md',BASE/'prepare01/report.json',BASE/'prepare01/oracle.json',
        ROOT/'data/thermo/phase_aware_dry_air_v1.json',*native_source_locks()]
def checked_preparation():
    r=read(BASE/'prepare01/report.json')
    if not r['passed']:raise ValueError('准备未过')
    for name,digest in dict(r['inputSHA256'],**r['caseInputSHA256']).items():
        if sha(ROOT/name.replace('\\','/'))!=digest:raise ValueError('准备输入改变 '+name)
    if CASE.name=='thermal02':
        repaired=read(BASE/'period_repair01/report.json')
        if not repaired['passed']:raise ValueError('周期修正未过')
        for name,digest in repaired['newInputSHA256'].items():
            if sha(ROOT/name.replace('\\','/'))!=digest:raise ValueError('修正输入改变')

def build():
    from runOpenFoamProcess import run
    start=time.monotonic();checked_preparation();out=BASE/'native01';out.mkdir(exist_ok=False)
    library=out/'library';shutil.copytree(ROOT/'native/openfoamPhaseAwareViscousAir',library)
    old=ROOT/'native/openfoamPhaseAwareAir'
    for name in ['rahPhaseAwareInternalEnergy.H','rahPhaseDomain.H','BOOST_LICENSE_1_0.txt','OPENFOAM_COPYING']:
        shutil.copyfile(old/name,library/name)
    paths=locks()+[p for p in library.rglob('*') if p.is_file()]
    save(out/'admission.json',dict(seconds=100,parentSeconds=120,memoryMiB=2048,outputMiB=96,automaticRetry=0,inputSHA256={str(p):sha(p) for p in paths}))
    shutil.copyfile(__file__,out/'executed_entry.py')
    proc=run(['wmake','libso'],library,out/'build',100,paths+[out/'admission.json'],96,2048)
    result=dict(passed=proc['passed'],process=proc,elapsedSeconds=time.monotonic()-start,CFDRuns=0,
        librarySHA256=sha(library/'librahPhaseAwareViscousAir.so') if proc['passed'] else None,physicalQualification=False)
    save(out/'report.json',result);print(json.dumps(dict(passed=result['passed'],seconds=result['elapsedSeconds'])));return result['passed']

def mesh():
    from runOpenFoamProcess import run
    start=time.monotonic();checked_preparation();out=BASE/'mesh02';out.mkdir(exist_ok=False)
    if not read(BASE/'native01/report.json')['passed']:raise ValueError('隔离库未通过')
    paths=locks()+[p for p in CASE.rglob('*') if p.is_file()]
    save(out/'admission.json',dict(seconds=30,parentSeconds=45,memoryMiB=2048,outputMiB=64,automaticRetry=0,inputSHA256={str(p):sha(p) for p in paths}))
    shutil.copyfile(__file__,out/'executed_entry.py')
    block=run(['blockMesh'],CASE,out/'block',15,paths,64,2048)
    save(out/'progress.json',dict(block=block))
    if not block['passed']:return False
    mesh_files=[p for p in (CASE/'constant/polyMesh').iterdir() if p.is_file()]
    check=run(['checkMesh','-allGeometry','-allTopology'],CASE,out/'check',30-block['elapsedSeconds'],paths+mesh_files,64,2048)
    passed=check['passed'] and 'Mesh OK.' in (out/'check/stdout.txt').read_text()
    result=dict(passed=passed,block=block,check=check,meshSHA256={str(p.relative_to(ROOT)):sha(p) for p in mesh_files},
        elapsedSeconds=time.monotonic()-start,physicalQualification=False)
    save(out/'report.json',result);print(json.dumps(dict(passed=passed,seconds=result['elapsedSeconds'])));return passed

def solve():
    from runOpenFoamProcess import run
    start=time.monotonic();checked_preparation();out=BASE/'thermal_run01';out.mkdir(exist_ok=False)
    mesh_report=read(BASE/'mesh02/report.json');build_report=read(BASE/'native01/report.json')
    if not mesh_report['passed'] or not build_report['passed']:raise ValueError('网格/原生未过')
    for name,h in mesh_report['meshSHA256'].items():
        if sha(ROOT/name)!=h:raise ValueError('网格改变')
    lib=BASE/'native01/library/librahPhaseAwareViscousAir.so';guardlib=ROOT/'outputs/v206_validation/native01/library/librahPhaseAwareAir.so'
    if sha(lib)!=build_report['librarySHA256']:raise ValueError('新库改变')
    if sha(guardlib)!=read(ROOT/'outputs/v206_validation/native_audit02/report.json')['registrationLibrarySHA256']:raise ValueError('守卫库改变')
    paths=locks()+[lib,guardlib,BASE/'mesh02/report.json',BASE/'native01/report.json',BASE/'period_repair01/report.json',
        BASE/'period_repair01/oracle.json',ROOT/'docs/validation/v207/PERIOD_REPAIR_PLAN.md',*[p for p in CASE.rglob('*') if p.is_file()]]
    save(out/'admission.json',dict(seconds=60,parentSeconds=75,steps=8,deltaT=5e-9,saveEvery=2,memoryMiB=2048,outputMiB=64,
        automaticRetry=0,inputSHA256={str(p):sha(p) for p in paths}))
    for name in ['runViscousThermalControl.py','phaseAwareViscousAir.py','phaseAwareDryAir.py']:
        shutil.copyfile(ROOT/'scripts'/name,out/name)
    gas=PhaseAwareViscousAir();guard=PhysicalTimeGuard(CASE,out/'saved_guard.json',32,8,auditor=lambda f:audit_state(f,gas.evaluate))
    proc=run(['foamRun'],CASE,out/'native',60,paths+[out/'admission.json'],64,2048,guard)
    text=(out/'native/stdout.txt').read_text();times=sorted(float(k) for k in guard.checked)
    initial='RAH_PHASE_GUARD_OK index=0 time=0 checks=1' in text
    end='RAH_PHASE_GUARD_END index=8 checks=9' in text
    dispatch=bool(re.search(r'transport\s+icoTabulated\s*;',text)) and 'Fourier' in text
    complete=len(times)==4 and np.allclose(times,np.arange(1,5)*1e-8,rtol=0,atol=1e-16)
    result=dict(passed=bool(proc['passed'] and initial and end and dispatch and complete),process=proc,saved=guard.checked,
        initialGuard=initial,perStepGuardVerified=end,viscousTransportDispatch=dispatch,expectedTimesComplete=bool(complete),
        physicalQualification=False,numericallyConverged=False,newQualifiedHeatAreaM2=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps(dict(passed=result['passed'],seconds=result['elapsedSeconds'],stop=proc['stopReason'],times=times)));return result['passed']

if __name__=='__main__':
    if len(sys.argv)!=2 or sys.argv[1] not in ['build','mesh','solve']:raise ValueError('仅允许build/mesh/solve具名阶段')
    raise SystemExit(0 if globals()[sys.argv[1]]() else 1)
