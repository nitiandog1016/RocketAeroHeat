"""隔离构建显式工作区包；容量/调用范围先核对，旧安装和对象只读。"""
from pathlib import Path
import argparse
import difflib
import json
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from runEilmerProcess import PREFIX, COMMIT, run, sha

BASE = ROOT / 'outputs/v187_validation'
PACKAGE = ROOT / 'tools/runtime/eilmer-explicit-workspace-v1'
PROFILE = 'rah-explicit-workspace-v1'


def audit():
    """模块一：精确分派/空访问与编译器类型布局，避免把载荷当总内存。"""
    started = time.monotonic()
    out = BASE / 'scope03'; out.mkdir(parents=True, exist_ok=False)
    prior = json.loads((ROOT / 'outputs/v185_validation/scope02/report.json').read_text())
    runtime_path = ROOT / 'tools/runtime/eilmer-explicit-dense-v1/runtime.json'
    runtime = json.loads(runtime_path.read_text())
    assert sha(runtime_path) == '8362d3615c9334039919e56042ec77a2dd429166e0f6ad20ebfc6dd504f82631'
    assert all(sha(p) == h for p, h in prior['originalSHA256'].items())
    assert all(sha(p) == h for p, h in runtime['filesSHA256'].items())
    source = PREFIX / 'source/src'; lmr = source / 'lmr'
    # 保存实际调用和定义；NK/check-jacobian路径不属于该显式运行包的开放范围。
    direct = []
    pattern = r'\b(?:saved_gradients|saved_lsqgradients|saved_source_terms|cell_jacobians|grad_save|gradients_save|Q_save|dRdU|invBlockDiag|doNotPerturb|evaluate_jacobian|evalRHS|SolverMode.transient)\b'
    for path in sorted(lmr.rglob('*.d')):
        if path.is_file():
            for n, line in enumerate(path.read_text().splitlines(), 1):
                if re.search(pattern, line): direct.append(dict(path=str(path.relative_to(source)), line=n, text=line))
    (out / 'accesses.json').write_text(json.dumps(direct, indent=2) + '\n')
    spans = {'commands/runsim.d': [(145,195),(340,402)], 'fluidfvcell.d': [(179,253),(1113,1129)],
        'fluidblock.d': [(237,355),(1481,1543),(1717,1824)], 'ufluidblock.d': [(300,355)],
        'flowgradients.d': [(35,112)], 'lsqinterp.d': [(250,344)], 'coredata.d': [(20,54)],
        'simcore.d': [(1,80)]}
    for name, intervals in spans.items():
        lines = (lmr / name).read_text().splitlines()
        (out / (name.replace('/','_') + '.txt')).write_text('\n'.join(
            '%d: %s' % (i,lines[i-1]) for lo,hi in intervals for i in range(lo,min(hi,len(lines))+1))+'\n')
    cell = (lmr / 'fluidfvcell.d').read_text()
    assert all(x in cell for x in ['if (fvcd.saved_gradients)', 'if (fvcd.saved_lsqgradients)',
        'if (fvcd.saved_source_terms)', 'if (fvcd.cell_jacobians)', 'invBlockDiag = new Matrix!double(ncq,ncq);'])
    PACKAGE.mkdir(parents=True, exist_ok=True)
    assert not (PACKAGE/'overlay').exists() and not (PACKAGE/'lmr-run').exists()
    # 只编译类型查询，不链接执行文件；保持所有影响布局的版本开关。
    probe = out / 'layout.d'
    types = {'FlowGradients':'FlowGradients.sizeof', 'LSQInterpGradients':'LSQInterpGradients.sizeof',
        'WLSQGradWorkspace':'WLSQGradWorkspace.sizeof', 'LSQInterpWorkspace':'LSQInterpWorkspace.sizeof',
        'FlowState':'FlowState.sizeof', 'FluidFVCell':'__traits(classInstanceSize, FluidFVCell)',
        'FVInterface':'__traits(classInstanceSize, FVInterface)', 'Matrix':'__traits(classInstanceSize, Matrix!double)',
        'Vector3':'Vector3.sizeof', 'Slice':'(double[]).sizeof', 'Scalar':'number.sizeof'}
    probe.write_text('module rah_layout;\nimport lmr.flowgradients; import lmr.lsqinterp; import lmr.flowstate;\n'
        'import lmr.fluidfvcell; import lmr.fvinterface; import nm.bbla; import nm.number; import geom;\n'+
        '\n'.join('pragma(msg, "RAH-LAYOUT %s=", %s);' % (k,v) for k,v in types.items())+'\n')
    flags = [x for x in prior['compileCommand'][1:] if x.startswith(('-d-version=','--d-version=','-I')) and x!='-d-version=run_main']
    # 官方命令直接列出的外部模块，在只查询类型时须补其模块搜索路径。
    flags += ['-I'+str(source/'extern'/p) for p in ['eqc/source','D-YAML/source','gzip','tinyendian/source']]
    flags += [x for x in prior['compileCommand'] if x.startswith('-J=')]
    command = [prior['compiler'], '-c', '-o-', '-boundscheck=on', *flags, str(probe)]
    layout = run(command, out, BASE / 'layout03', 30, 64, [probe, Path(__file__), *prior['originalSHA256']])
    assert layout['passed'], '布局查询失败'
    sizes = {k:int(v) for k,v in re.findall(r'RAH-LAYOUT (\w+)=(\d+)', (BASE/'layout03/stdout.txt').read_text())}
    assert set(sizes)==set(types) and sizes['Scalar']==8 and sizes['Slice']==16
    geometry_path = ROOT/'outputs/v178_validation/real01_geometry.json'
    geom = json.loads(geometry_path.read_text()); n=geom['cells']+sum(geom['boundaryCounts'].values())
    nf=geom['internalFaces']+sum(geom['boundaryCounts'].values()); neq=7
    # 精确nsp=1/nmodes=0/nturb=2；动态梯度载荷按构造器的原数组长度计。
    removable = dict(denseJacobian=n*neq*neq*8, savedSource=n*neq*8, savedPerturb=n,
        savedFlowGradients=n*(sizes['FlowGradients']+24), savedLSQGradients=n*(sizes['LSQInterpGradients']+48),
        cellJacobianRowViews=n*neq*16, inverseMatrixObjects=n*sizes['Matrix'], inverseMatrixData=n*neq*neq*8)
    necessary = dict(stageArrays=n*neq*8*8, sourceTerms=n*neq*8,
        flowStateInline=n*sizes['FlowState'], flowGradients=n*(sizes['FlowGradients']+24),
        flowGradientWorkspace=n*sizes['WLSQGradWorkspace'], interpolationWorkspace=n*sizes['LSQInterpWorkspace'],
        interpolationGradients=n*(sizes['LSQInterpGradients']+48), cellObjects=n*sizes['FluidFVCell'],
        faceObjects=nf*sizes['FVInterface'])
    capacity=dict(sizes=sizes, cellsIncludingGhostCapacity=n, faces=nf, removablePayloadBytes=removable,
        removablePayloadSumBytes=sum(removable.values()), necessaryPartialPayloadBytes=necessary,
        necessaryPartialSumBytes=sum(necessary.values()), budgetBytes=8*1024**3,
        fullProcessEstimate=False, excluded='Additional face/vertex/connectivity/state data, array spare capacity, allocator/GC/runtime/temporary buffers and virtual mappings; sums are payload only.',
        necessaryLowerBoundExceedsBudget=sum(necessary.values())>8*1024**3)
    assert not capacity['necessaryLowerBoundExceedsBudget'], '已知必要载荷已超预算，停止构建'
    # 模块二：在旧包准确覆盖层上应用完整工作区改动，不改共享源。
    overlay = PACKAGE/'overlay/src/lmr'; overlay.mkdir(parents=True)
    shutil.copyfile(ROOT/'tools/runtime/eilmer-explicit-dense-v1/overlay/src/lmr/fluidblock.d', overlay/'fluidblock.d')
    shutil.copyfile(lmr/'ufluidblock.d', overlay/'ufluidblock.d')
    patch = ROOT/'patches/eilmer/explicit-workspace-v1.patch'
    assert sha(patch)=='fda235021f0c9ecf287c55b9a72ee65511cc51d53baa398007d2446cf358fd2c'
    applied=subprocess.check_output(['patch','--batch','--forward','-p1','-i',str(patch)],cwd=PACKAGE/'overlay',text=True,timeout=10)
    (out/'patch_application.txt').write_text(applied)
    # 记录真实已分配长度及构造后单元对象，不能只记录预期值。
    ufile=overlay/'ufluidblock.d'; text=ufile.read_text()
    anchor='        facedata.dL.length = nfaces;'
    assert text.count(anchor)==1
    observation='''        // 观察完整工作区分配，正常梯度和工作区保留。
        if (myConfig.solverMode == SolverMode.transient) {
            writeln("RAH-WORKSPACE jacobian=", celldata.cell_jacobians.length,
                " sourceBackup=", celldata.saved_source_terms.length,
                " viscousBackup=", celldata.saved_gradients.length,
                " interpBackup=", celldata.saved_lsqgradients.length,
                " perturbBackup=", celldata.doNotPerturb.length,
                " viscous=", celldata.gradients.length, " interp=", celldata.lsqgradients.length,
                " interpWorkspace=", celldata.lsqws.length);
        }
'''
    text=text.replace(anchor,observation+anchor)
    anchor='        // Bind the interfaces, vertices and cells together,'
    assert text.count(anchor)==1
    observation='''        // 观察已构造的内部单元，不触发隐式计算、不修改任何状态。
        if (myConfig.solverMode == SolverMode.transient) {
            size_t implicitObjects = 0;
            foreach (cell; cells) {
                if (cell.dRdU.length || cell.Q_save.length || cell.grad_save !is null ||
                    cell.gradients_save !is null || cell.invBlockDiag !is null) implicitObjects++;
            }
            writeln("RAH-WORKSPACE-CELLS internal=", cells.length, " implicitObjects=", implicitObjects);
        }
'''
    ufile.write_text(text.replace(anchor,observation+anchor))
    built_patch=ROOT/'patches/eilmer/explicit-workspace-checked-v1.patch'
    assert not built_patch.exists()
    built_patch.write_text(''.join(''.join(difflib.unified_diff((lmr/name).read_text().splitlines(True),
        (overlay/name).read_text().splitlines(True),fromfile='a/src/lmr/'+name,tofile='b/src/lmr/'+name)) for name in ['fluidblock.d','ufluidblock.d']))
    for name in ['fluidblock.d','ufluidblock.d']:shutil.copyfile(overlay/name,out/('executed_'+name))
    command=[]
    for arg in prior['compileCommand']:
        if arg.startswith('-of='):command.append('-of='+str(PACKAGE/'lmr-run'))
        elif arg=='-J='+str(ROOT/'tools/runtime/eilmer-explicit-dense-v1'):command.append('-J='+str(PACKAGE))
        elif arg==str(ROOT/'tools/runtime/eilmer-explicit-dense-v1/overlay/src/lmr/fluidblock.d'):command.append(str(overlay/'fluidblock.d'))
        elif arg==str(lmr/'ufluidblock.d'):command.append(str(overlay/'ufluidblock.d'))
        else:command.append(arg)
    assert '-d-version=newton_krylov' in command and '-boundscheck=on' in command
    shutil.copyfile(ROOT/'tools/runtime/eilmer-explicit-dense-v1/buildinfo.json',PACKAGE/'buildinfo.json')
    locks={**prior['originalSHA256'],**runtime['filesSHA256']}
    locks[str(runtime_path)]=sha(runtime_path)
    locks[str(geometry_path)]=sha(geometry_path)
    report=dict(passed=True,profile=PROFILE,sourceCommit=COMMIT,capacity=capacity,
        originalSHA256=locks,compileCommand=command,compiledDSourceCount=prior['compiledDSourceCount'],
        layoutReceiptSHA256=sha(BASE/'layout03/report.json'),patch=str(built_patch),patchSHA256=sha(built_patch),
        patchedSHA256={str(overlay/name):sha(overlay/name) for name in ['fluidblock.d','ufluidblock.d']},
        buildInfoSHA256=sha(PACKAGE/'buildinfo.json'),elapsedSeconds=time.monotonic()-started,physicalQualification=False)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=True,capacity=capacity,elapsedSeconds=report['elapsedSeconds'])))


def build():
    """模块三：固定预算重编必要流动目标，准确回执随运行包绑定。"""
    started=time.monotonic();scope_path=BASE/'scope03/report.json';scope=json.loads(scope_path.read_text())
    assert scope['passed'] and not scope['capacity']['necessaryLowerBoundExceedsBudget']
    assert all(sha(p)==h for p,h in scope['originalSHA256'].items())
    assert all(sha(p)==h for p,h in scope['patchedSHA256'].items())
    assert sha(PACKAGE/'buildinfo.json')==scope['buildInfoSHA256']
    locks=[*scope['originalSHA256'],*scope['patchedSHA256'],scope_path,scope['patch'],PACKAGE/'buildinfo.json',Path(__file__)]
    process=run(scope['compileCommand'],PACKAGE,BASE/'build01',300-(time.monotonic()-started),512,locks)
    result=dict(passed=process['passed'],process=process,profile=PROFILE,elapsedSeconds=time.monotonic()-started,
        compiledDSourceCount=scope['compiledDSourceCount'],sharedFilesStable=process['inputsStable'],physicalQualification=False)
    if process['passed']:
        files={str(PACKAGE/n):sha(PACKAGE/n) for n in ['lmr-run','buildinfo.json','overlay/src/lmr/fluidblock.d','overlay/src/lmr/ufluidblock.d']}
        package=dict(schema='RAH.EilmerWorkspaceRuntime.v1',profile=PROFILE,sourceCommit=COMMIT,binary=str(PACKAGE/'lmr-run'),
            filesSHA256=files,patch=scope['patch'],patchSHA256=scope['patchSHA256'],boundsCheck='on',
            buildReceipt=str(BASE/'build01/report.json'),buildReceiptSHA256=sha(BASE/'build01/report.json'),
            sameUpdateFormulas=True,smallCaseCompatibility='pending',physicalQualification=False)
        (PACKAGE/'runtime.json').write_text(json.dumps(package,indent=2)+'\n')
        result.update(runtimePackage=str(PACKAGE/'runtime.json'),runtimePackageSHA256=sha(PACKAGE/'runtime.json'))
    (BASE/'build_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['passed','elapsedSeconds','compiledDSourceCount','sharedFilesStable']}))
    return result['passed']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['audit','build']);a=p.parse_args()
    if a.stage=='audit':audit()
    else:raise SystemExit(0 if build() else 1)
