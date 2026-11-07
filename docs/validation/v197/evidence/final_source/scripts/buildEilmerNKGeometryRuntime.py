"""采用官方复数NK目标并保留已验证几何修正；不改通量、更新或共享源。"""
from pathlib import Path
import json,shutil,time
from runEilmerProcess import run,sha,PREFIX,COMMIT
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v197_validation'
PACKAGE=ROOT/'tools/runtime/eilmer-nk-geometric-v1'

def main():
    start=time.monotonic()
    old=ROOT/'outputs/v187_validation/scope03/report.json';r=json.loads(old.read_text())
    src=PREFIX/'source/src/lmr';previous=ROOT/'tools/runtime/eilmer-geometric-length-v1'
    meta=json.loads((previous/'runtime.json').read_text())
    assert all(sha(p)==h for p,h in r['originalSHA256'].items())
    assert all(sha(p)==h for p,h in meta['filesSHA256'].items())
    assert json.loads((ROOT/'outputs/v196_validation/small_validation02/report.json').read_text())['passed']
    PACKAGE.mkdir(exist_ok=False);overlay=PACKAGE/'overlay/src/lmr';overlay.mkdir(parents=True)
    shutil.copyfile(previous/'overlay/src/lmr/fluidfvcell.d',overlay/'fluidfvcell.d')
    shutil.copyfile(previous/'buildinfo.json',PACKAGE/'buildinfo.json')
    command=[]
    for a in r['compileCommand']:
        if a.startswith('-of='):a='-of='+str(PACKAGE/'lmrZ-run')
        elif a.startswith('-J=') and 'eilmer-explicit-workspace-v1' in a:a='-J='+str(PACKAGE)
        elif a.endswith('/overlay/src/lmr/fluidblock.d'):a=str(src/'fluidblock.d')
        elif a.endswith('/overlay/src/lmr/ufluidblock.d'):a=str(src/'ufluidblock.d')
        elif a==str(src/'fluidfvcell.d'):a=str(overlay/'fluidfvcell.d')
        command.append(a)
    command.insert(1,'-d-version=complex_numbers')
    assert '-d-version=newton_krylov' in command and '-boundscheck=on' in command
    assert sum(x.endswith('.d') for x in command)==347
    locks=[*r['originalSHA256'],old,previous/'runtime.json',overlay/'fluidfvcell.d',PACKAGE/'buildinfo.json',
           ROOT/'docs/validation/v197/ENTRY_DECISION.md',ROOT/'patches/eilmer/geometric-length-v1.patch',Path(__file__)]
    process=run(command,PACKAGE,BASE/'build01',240,512,locks,6144,6144)
    result=dict(passed=process['passed'],elapsedSeconds=time.monotonic()-start,process=process,physicalQualification=False)
    if process['passed']:
        files={str(p):sha(p) for p in [PACKAGE/'lmrZ-run',PACKAGE/'buildinfo.json',overlay/'fluidfvcell.d']}
        runtime=dict(schema='RAH.EilmerNKGeometryRuntime.v1',profile='rah-nk-geometric-v1',sourceCommit=COMMIT,
            binary=str(PACKAGE/'lmrZ-run'),numberType='complex',filesSHA256=files,compiledDSourceCount=347,
            buildReceipt=str(BASE/'build01/report.json'),buildReceiptSHA256=sha(BASE/'build01/report.json'),
            geometryPatchSHA256=sha(ROOT/'patches/eilmer/geometric-length-v1.patch'),
            officialFullNKStorage=True,smallControlPassed=False,physicalQualification=False)
        (PACKAGE/'runtime.json').write_text(json.dumps(runtime,indent=2)+'\n')
    (BASE/'build_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='process'}));return process['passed']

if __name__=='__main__':raise SystemExit(0 if main() else 1)
