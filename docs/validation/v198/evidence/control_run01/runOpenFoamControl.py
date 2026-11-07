"""唯一600步原生控制：严格准入与保存监督，不作生产箭体默认外流。"""
from pathlib import Path
import json,shutil,time
from prepareOpenFoamControl import ROOT,BASE,sha,save
from openFoamFields import SavedGuard

PROFILE='of13-nasa7-tetra-lts-600-v1'

def main():
    """模块一：新目录、一次运行、准确准入；旧几何与失败永不覆盖。"""
    from runOpenFoamProcess import run
    start=time.monotonic();evidence=[BASE/'prepared_audit01/report.json',BASE/'independent01/report.json',
        ROOT/'docs/validation/v198/CONTROL_PLAN.md']
    for p in evidence[:2]:assert json.loads(p.read_text())['passed']
    old=BASE/'control03';prep=json.loads((BASE/'preparation_repaired.json').read_text())
    for p,h in prep['inputSHA256'].items():assert sha(p)==h
    case=BASE/'development01';shutil.copytree(old,case)
    out=BASE/'control_run01';out.mkdir(exist_ok=False)
    cfg=case/'system/controlDict';text=cfg.read_text();assert text.count('endTime 0;')==1
    cfg.write_text(text.replace('endTime 0;','endTime 600;').replace('writeInterval 20;','writeInterval 100;'))
    with (case/'attempt.json').open('x') as f:json.dump(dict(profile=PROFILE,steps=600,automaticContinuation=0),f)
    scripts=['runOpenFoamControl.py','runOpenFoamProcess.py','openFoamFields.py','openFoamAirThermo.py']
    for name in scripts:shutil.copyfile(ROOT/'scripts'/name,out/name)
    locks=[*[ROOT/'scripts'/n for n in scripts],*evidence,*[p for p in case.rglob('*') if p.is_file()]]
    save(out/'admission.json',dict(profile=PROFILE,steps=600,seconds=240,addressSpaceMiB=2048,cpuCount=1,
        outputMiB=96,savedSteps=list(range(100,601,100)),inputSHA256={str(p):sha(p) for p in locks},
        physicalQualification=False,automaticContinuation=0))
    guard=SavedGuard(case,14892,out/'saved_guard.json')
    result=run(['foamRun'],case,out/'native',240,locks,96,2048,guard)
    log=(out/'native/stdout.txt').read_text();normal=result['passed'] and 'End\n' in log and set(guard.checked)==set(map(str,range(100,601,100)))
    report=dict(profile=PROFILE,process=result,completed600=normal,elapsedSeconds=time.monotonic()-start,
        outputsSHA256={str(p):sha(p) for p in case.rglob('*') if p.is_file()},physicalQualification=False)
    save(out/'report.json',report)
    print(json.dumps(dict(passed=normal,exitCode=result['nativeExitCode'],stop=result['stopReason'],elapsedSeconds=report['elapsedSeconds'])))
    return normal

if __name__=='__main__':raise SystemExit(0 if main() else 1)
