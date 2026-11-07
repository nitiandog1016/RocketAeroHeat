"""构建隔离显式存储包：复用官方实数目标和固定依赖，不修改共享安装或对象。"""
from pathlib import Path
import argparse
import difflib
import json
import re
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from runEilmerProcess import PREFIX, COMMIT, run, sha

PROFILE = 'rah-explicit-dense-v1'
PACKAGE = ROOT / 'tools/runtime/eilmer-explicit-dense-v1'
BASE = ROOT / 'outputs/v185_validation'


def audit():
    """模块一：准确源码的直接访问、别名、ghost范围及官方编译命令只读核对。"""
    start = time.monotonic()
    out = BASE / 'scope01'
    out.mkdir(parents=True, exist_ok=False)
    source = PREFIX / 'source'
    src = source / 'src'
    lmr = src / 'lmr'
    assert subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip() == COMMIT
    subprocess.run(['git', '-C', str(source), 'diff', '--exit-code'], check=True)
    direct = []
    for path in sorted(src.rglob('*.d')):
        for line, text in enumerate(path.read_text().splitlines(), 1):
            if re.search(r'\.(?:U|dUdt)[0-4]\b', text):
                direct.append(dict(path=str(path.relative_to(source)), line=line, text=text))
    relevant = [r for r in direct if r['path'].startswith('src/lmr/')]
    assert len(relevant) == 37 and {r['path'] for r in relevant} == {'src/lmr/fluidblock.d', 'src/lmr/fluidfvcell.d'}
    assert all(r['path'].startswith('src/slf/') for r in direct if r not in relevant)
    spans = {'coredata.d': [(20, 94)], 'conservedquantities.d': [(24, 54)],
             'fluidfvcell.d': [(180, 235)], 'ufluidblock.d': [(317, 371), (419, 487)],
             'fluidblock.d': [(237, 304), (2012, 2138)],
             'globalconfig.d': [(2420, 2433)],
             'simcore_gasdynamic_step.d': [(1103, 1228), (1394, 1407)]}
    for name, ranges in spans.items():
        lines = (lmr / name).read_text().splitlines()
        selected = [f'{i}: {lines[i-1]}' for lo, hi in ranges for i in range(lo, hi + 1)]
        (out / (name + '.txt')).write_text('\n'.join(selected) + '\n')
    cell = (lmr / 'fluidfvcell.d').read_text()
    block = (lmr / 'ufluidblock.d').read_text()
    core = (lmr / 'coredata.d').read_text()
    assert 'immutable alias ConservedQuantities = number[];' in (lmr / 'conservedquantities.d').read_text()
    for family in ('U', 'dUdt'):
        for stage in range(5):
            assert f'fvcd.{family}{stage}[id*ncq + 0 .. id*ncq + ncq]' in cell
    assert 'foreach (bndry; grid.boundaries) nghost += bndry.face_id_list.length;' in block
    assert 'int ghost_cell_id = to!int(cells.length);' in block and 'ghost_cell_id++;' in block
    assert 'if (bc[i].ghost_cell_data_available)' in block
    assert 'void allocate_minimal' in core  # CHT独立分配器未修改，不能声称整个引擎都已压缩。
    dry = subprocess.check_output(['make', '-n', '-W', 'fluidblock.d', 'FLAVOUR=fast', 'DMD=ldc2',
                                   'OPT_OUT_MPI=1', 'lmr-run'], cwd=lmr, text=True, timeout=20)
    (out / 'official_make_dry_run.txt').write_text(dry)
    commands = [shlex.split(line) for line in dry.splitlines() if line.startswith('ldc2 ')]
    assert len(commands) == 1
    original = commands[0]
    assert '-boundscheck=off' in original and '-flto=full' in original and '-d-version=run_main' in original
    assert not any('complex_numbers' in x or 'mpi_parallel' in x for x in original)
    PACKAGE.mkdir(parents=True, exist_ok=False)
    overlay = PACKAGE / 'overlay/src/lmr'
    overlay.mkdir(parents=True)
    shutil.copyfile(lmr / 'fluidblock.d', overlay / 'fluidblock.d')
    patch = ROOT / 'patches/eilmer/explicit-dense-storage-v1.patch'
    assert sha(patch) == '0dd224f7470df1e125c9caf703522c5fea9910ea74e811f1cccda14aa44b069f'
    applied = subprocess.check_output(['patch', '--batch', '--forward', '-p1', '-i', str(patch)],
                                     cwd=PACKAGE / 'overlay', text=True, timeout=10)
    (out / 'patch_application.txt').write_text(applied)
    patched = overlay / 'fluidblock.d'
    text = patched.read_text()
    anchor = '        celldata.source_terms.length = (ncells + nghost)*neq;'
    assert text.count(anchor) == 1
    observation = '''        // 记录实际分配容量；仅观察，不参与任何数值更新。
        if (myConfig.solverMode == SolverMode.transient) {
            writeln("RAH-STORAGE cells=", ncells, " ghosts=", nghost, " neq=", neq,
                    " nftl=", nftl, " scalarBytes=", number.sizeof,
                    " stageLength=", stage_storage, " payloadBytes=",
                    (celldata.U0.length + celldata.U1.length + celldata.U2.length + celldata.U3.length + celldata.U4.length
                     + celldata.dUdt0.length + celldata.dUdt1.length + celldata.dUdt2.length + celldata.dUdt3.length + celldata.dUdt4.length)*number.sizeof);
        }
'''
    patched.write_text(text.replace(anchor, observation + anchor))
    built_patch = ROOT / 'patches/eilmer/explicit-dense-checked-v1.patch'
    assert not built_patch.exists()
    built_patch.write_text(''.join(difflib.unified_diff((lmr / 'fluidblock.d').read_text().splitlines(True),
        patched.read_text().splitlines(True), fromfile='a/src/lmr/fluidblock.d', tofile='b/src/lmr/fluidblock.d')))
    shutil.copyfile(patched, out / 'executed_fluidblock.d')
    compiler = PREFIX / 'ldc2-1.41.0-linux-x86_64/bin/ldc2'
    command = [str(compiler)]
    dependencies = []
    for arg in original[1:]:
        if arg == '-boundscheck=off': command.append('-boundscheck=on')
        elif arg.startswith('-of='): command.append('-of=' + str(PACKAGE / 'lmr-run'))
        elif arg == '-I..': command.append('-I' + str(src))
        elif arg == '-J=.': command.extend(['-J=' + str(PACKAGE), '-J=' + str(lmr)])
        elif not arg.startswith('-'):
            path = (lmr / arg).resolve()
            assert path.is_file()
            dependencies.append(path)
            command.append(str(patched if path == lmr / 'fluidblock.d' else path))
        else: command.append(arg)
    buildinfo = dict(buildDate=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()), revisionId=COMMIT[:7],
                     fullRevisionId=COMMIT, revisionDate='Mon Dec 22 17:46:51 2025 +1000', parallelType='OpenMPI')
    (PACKAGE / 'buildinfo.json').write_text(json.dumps(buildinfo, indent=2) + '\n')
    original_locks = dependencies + [lmr / 'lmr-run.o', lmr / 'lmr-run', compiler]
    original_locks += [p for p in (PREFIX / 'install/bin').glob('lmr*') if p.is_file()]
    report = dict(passed=True, profile=PROFILE, sourceCommit=COMMIT, compiler=str(compiler), compilerSHA256=sha(compiler),
        directAccesses=direct, relevantDirectAccesses=len(relevant),
        aliasContract='ConservedQuantities is number[]; cells have neq slices; residuals index cidx*neq+j; ghost IDs start at ncells and increase at most once per boundary face.',
        untouchedMinimalCHTAllocator=True, allGeometryDynamicCoverage=False,
        originalSHA256={str(p): sha(p) for p in original_locks},
        package=str(PACKAGE), patch=str(built_patch), patchSHA256=sha(built_patch),
        patchedSource=str(patched), patchedSourceSHA256=sha(patched), compileCommand=command,
        buildInfoSHA256=sha(PACKAGE / 'buildinfo.json'), boundsCheck='on',
        compiledDSourceCount=sum(p.suffix == '.d' for p in dependencies),
        elapsedSeconds=time.monotonic() - start, physicalQualification=False)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['passed', 'relevantDirectAccesses', 'compiledDSourceCount', 'elapsedSeconds']}))


def build():
    """模块二：完整重编必要流动目标，固定源/编译器/库并保留真实退出。"""
    started = time.monotonic()
    report_path = BASE / 'scope01/report.json'
    scope = json.loads(report_path.read_text())
    assert scope['passed'] and scope['sourceCommit'] == COMMIT
    assert all(sha(p) == h for p, h in scope['originalSHA256'].items())
    assert sha(scope['patchedSource']) == scope['patchedSourceSHA256']
    assert sha(PACKAGE / 'buildinfo.json') == scope['buildInfoSHA256']
    locks = list(scope['originalSHA256']) + [report_path, scope['patchedSource'], scope['patch'],
                                            PACKAGE / 'buildinfo.json', Path(__file__)]
    r = run(scope['compileCommand'], PACKAGE, BASE / 'build01', 300-(time.monotonic()-started), 512, locks)
    result = dict(passed=r['passed'], sourceCommit=COMMIT, profile=PROFILE,
                  process=r, elapsedSeconds=time.monotonic()-started, productionKernelSourceCompiles=1,
                  compiledDSourceCount=scope['compiledDSourceCount'], sharedFilesStable=r['inputsStable'],
                  sourceAuditSHA256=sha(report_path), physicalQualification=False)
    if r['passed']:
        files = {str(PACKAGE / n): sha(PACKAGE / n) for n in ['lmr-run', 'buildinfo.json', 'overlay/src/lmr/fluidblock.d']}
        package = dict(schema='RAH.EilmerDenseRuntime.v1', profile=PROFILE, sourceCommit=COMMIT,
                       binary=str(PACKAGE / 'lmr-run'), filesSHA256=files,
                       patch=scope['patch'], patchSHA256=scope['patchSHA256'],
                       boundsCheck='on', buildReceipt=str(BASE / 'build01/report.json'),
                       buildReceiptSHA256=sha(BASE / 'build01/report.json'),
                       sameUpdateFormulas=True, smallCaseCompatibility='pending', physicalQualification=False)
        (PACKAGE / 'runtime.json').write_text(json.dumps(package, indent=2)+'\n')
        result['runtimePackage']=str(PACKAGE / 'runtime.json')
        result['runtimePackageSHA256']=sha(PACKAGE / 'runtime.json')
    (BASE / 'build_result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['passed', 'elapsedSeconds', 'sharedFilesStable', 'compiledDSourceCount']}))
    return result['passed']


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('stage', choices=['audit', 'build']); a = p.parse_args()
    if a.stage == 'audit': audit()
    else: raise SystemExit(0 if build() else 1)
