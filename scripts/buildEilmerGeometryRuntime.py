"""修正固定Eilmer几何长度覆盖；隔离运行包、准确源与有限构建预算。"""
from pathlib import Path
import difflib
import json
import shutil
import time
import numpy as np
from runEilmerProcess import run,sha,PREFIX,COMMIT

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v196_validation'
OLD=ROOT/'tools/runtime/eilmer-explicit-workspace-v1'
PACKAGE=ROOT/'tools/runtime/eilmer-geometric-length-v1'
PROFILE='rah-geometric-length-v1'

def main():
    """模块一：先绑定实测故障、ABI/源，再登记唯一构建准入。"""
    start=time.monotonic()
    observed=BASE/'observation03/report.json';obs=json.loads(observed.read_text());assert obs['passed']
    tet,hexagon=[r['observation'] for r in obs['results']]
    assert tet['lengths']==['nan']*3 and tet['signal']==1e-99 and tet['Lmin']>0
    expected=(max(abs(np.dot(hexagon['velocity'],hexagon['normals'][i])) for i in hexagon['selectedFaces'])+
              hexagon['soundSpeed'])/min(hexagon['lengths'])
    assert abs(expected/hexagon['signal']-1)<2e-15
    prior=json.loads((ROOT/'outputs/v187_validation/scope03/report.json').read_text())
    runtime=json.loads((OLD/'runtime.json').read_text())
    assert all(sha(p)==h for p,h in {**prior['originalSHA256'],**runtime['filesSHA256']}.items())
    cell=PREFIX/'source/src/lmr/fluidfvcell.d';original=cell.read_text()
    anchor='        fvcd.lengths[id][0] = iL;\n        fvcd.lengths[id][1] = jL;\n        fvcd.lengths[id][2] = kL;\n'
    assert original.count(anchor)==1 and 'lengths= fvcd.lengths[id][0 .. 3];' in original
    patched=original.replace(anchor,'        // 各形状分支已通过属性写入共享长度；不可再用未赋值的局部量覆盖。\n')
    PACKAGE.mkdir(exist_ok=False);overlay=PACKAGE/'overlay/src/lmr';overlay.mkdir(parents=True)
    for n in ['fluidblock.d','ufluidblock.d']:
        shutil.copyfile(OLD/'overlay/src/lmr'/n,overlay/n)
    (overlay/'fluidfvcell.d').write_text(patched)
    shutil.copyfile(OLD/'buildinfo.json',PACKAGE/'buildinfo.json')
    patch=ROOT/'patches/eilmer/geometric-length-v1.patch'
    assert not patch.exists()
    patch.write_text(''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),
                         fromfile='a/src/lmr/fluidfvcell.d',tofile='b/src/lmr/fluidfvcell.d')))
    command=[]
    for arg in prior['compileCommand']:
        if arg.startswith('-of='):arg='-of='+str(PACKAGE/'lmr-run')
        elif arg=='-J='+str(OLD):arg='-J='+str(PACKAGE)
        elif arg==str(cell):arg=str(overlay/'fluidfvcell.d')
        elif str(OLD/'overlay') in arg:arg=arg.replace(str(OLD/'overlay'),str(PACKAGE/'overlay'))
        command.append(arg)
    # 旧对象是347源的单个LLVM bitcode，不能准确单模块替换；仅重编这一实数目标。
    admission=dict(profile=PROFILE,observedSHA256=sha(observed),baselineSHA256=sha(OLD/'runtime.json'),
        nativeSeconds=180,memoryMiB=4096,outputMiB=512,automaticRebuild=0,
        operatorChange='Remove post-switch overwrite of shared geometric lengths; finite CFL signal restored for non-hexahedra.',
        compileCommand=command,sourceSHA256=sha(cell),patchSHA256=sha(patch),physicalQualification=False)
    plan=ROOT/'docs/validation/v196/BUILD_ADMISSION.json';assert not plan.exists()
    plan.write_text(json.dumps(admission,indent=2)+'\n')
    locks=[*prior['originalSHA256'],*runtime['filesSHA256'],OLD/'runtime.json',observed,patch,plan,
           Path(__file__),PACKAGE/'buildinfo.json',*overlay.glob('*.d')]
    result=run(command,PACKAGE,BASE/'build01',180,512,locks,4096,4096)
    report=dict(passed=result['passed'],process=result,elapsedSeconds=time.monotonic()-start,profile=PROFILE,
                sameNumericsAsBaseline=False,physicalQualification=False)
    if result['passed']:
        package=dict(schema='RAH.EilmerGeometryRuntime.v1',profile=PROFILE,sourceCommit=COMMIT,
            binary=str(PACKAGE/'lmr-run'),baseRuntime=str(OLD/'runtime.json'),baseRuntimeSHA256=sha(OLD/'runtime.json'),
            filesSHA256={str(p):sha(p) for p in [PACKAGE/'lmr-run',PACKAGE/'buildinfo.json',*overlay.glob('*.d')]},
            buildReceipt=str(BASE/'build01/report.json'),buildReceiptSHA256=sha(BASE/'build01/report.json'),
            patch=str(patch),patchSHA256=sha(patch),sameNumericsAsBaseline=False,
            change=admission['operatorChange'],physicalQualification=False)
        (PACKAGE/'runtime.json').write_text(json.dumps(package,indent=2)+'\n')
    (BASE/'build_result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='process'}))
    return result['passed']

if __name__=='__main__':raise SystemExit(0 if main() else 1)
