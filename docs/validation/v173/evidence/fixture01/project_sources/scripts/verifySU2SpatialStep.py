"""模块用途：有限编译真实层流CNS夹具，独立审核限步并恢复全部共享源；不启动CFD。"""
import argparse
import csv
import hashlib
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from prepareSU2SpatialStepPatch import ROOT, SOURCE, BUILD, FILE, prepared_sources
from su2ThermalRuntime import read, sha
from validateIdealAir import BIN, run, save

DEFAULT_OUTPUT = ROOT/'outputs/v173_validation/fixture01'
CONFIG = ROOT/'outputs/v171_validation/spatial_profile01/case.cfg'


def verify_events(folder):
    """模块一：真实类事件以独立物性重新核对；合成上下文不充当真实网格或失败定位。"""
    from idealAirProperties import IdealAir
    with (folder/'conservative_steps.csv').open(encoding='utf8', newline='') as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 25 and [int(row['event']) for row in rows] == list(range(1,26))
    assert sum(row['old_valid'] == '0' for row in rows) == 2
    assert sum(row['old_valid'] == '1' and float(row['theta']) == 0 for row in rows) == 2
    core = IdealAir(); maximum = 0.; positive = []; effective = 0
    for row in rows:
        assert int(row['dimension']) == 3 and float(row['k']) == 0
        assert row['context_known'] == '0', '此实际类夹具没有真实网格上下文，不得冒充节点定位'
        assert 'RAH_PROPERTY_QUERY mode=rhoe' in row['reason']
        if row['old_valid'] == '0':
            assert 'phase=node.old' in row['reason']
            continue
        old, trial, accepted = [np.array([float(row[prefix+'_'+name]) for name in ('rho','rhou','rhov','rhow','rhoE')])
            for prefix in ('old','trial','accepted')]
        theta = float(row['theta']); attempts = int(row['attempts'])
        assert np.isfinite(accepted).all() and accepted[0] > 0 and 1 <= attempts <= 32
        assert theta == 0 or theta == 2.**(-attempts)
        error = float(np.max(abs(accepted-(old+theta*(trial-old)))/np.maximum(1.,abs(accepted))))
        assert error < 2e-13; maximum = max(maximum,error)
        energy = accepted[4]/accepted[0]-.5*np.sum((accepted[1:4]/accepted[0])**2)
        state = core.state('rhoe',accepted[0],energy)
        assert 60 <= state['T'][0] <= 2000
        if theta == 0:
            assert attempts == 32 and np.array_equal(accepted,old)
        else:
            assert not np.array_equal(accepted,old), '正theta没有实际状态改变，不能称正进展'
            positive.append(theta); effective += 1
    return dict(events=len(rows), invalidOldEvents=2, zeroStepEvents=2, positiveEffectiveEvents=effective,
        maximumConservativeSegmentError=maximum, minimumPositiveTheta=min(positive),
        laminarKExactlyZero=True, independentEOSChecked=True, syntheticContextOnly=True,
        sourceSSTEventPairsUsed=0, physicsQualified=False, newQualifiedArea=0)


def main(output):
    """模块二：准确源临时覆盖、独立fixture.exe、finally恢复；共享对象和引擎exe不得写入。"""
    folder = Path(output).resolve(); folder.relative_to(ROOT/'outputs'); folder.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); compile_seconds = test_seconds = 0.; snapshots = {}; originals = {}; objects = {}
    report = dict(passed=False, CFDLaunches=0, sourceSSTEventPairsUsed=0, physicsQualified=False,
        maximumCompileSeconds=40, maximumTestSeconds=10, maximumMiscSeconds=10, maximumStageSeconds=60)
    project_files = [ROOT/name for name in ('scripts/verifySU2SpatialStep.py','scripts/prepareSU2SpatialStepPatch.py',
        'scripts/prepareSU2ConservativeStepPatch.py','tests/testSU2SpatialStep.cpp',
        'third_party/ideal_air/ConservativeStep.inc','third_party/ideal_air/ConservativeStepTrace.hpp',
        'third_party/ideal_air/PropertyQueryContext.hpp','scripts/validateIdealAir.py','docs/validation/v173/PLAN.md')]
    locks = {str(path):sha(path) for path in project_files+[CONFIG]}
    save(folder/'input_locks.json',locks)
    for path in project_files:
        target = folder/'project_sources'/path.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True); target.write_bytes(path.read_bytes())
    try:
        previous, tree, inherited, snapshots = prepared_sources()
        assert len(snapshots) == 21, '准确待覆盖源数不是已登记21'
        save(folder/'parent_source_tree.json',tree)
        save(folder/'runtime_parent.json',previous)
        originals = {path:path.read_bytes() if path.exists() else None for path in snapshots}
        before_sources = {str(path):hashlib.sha256(raw).hexdigest() if raw is not None else None for path,raw in originals.items()}
        save(folder/'shared_sources_before.json',before_sources)
        objects = {str(path.relative_to(BUILD)):sha(path) for path in BUILD.rglob('*.obj')}
        assert len(objects) == 298; save(folder/'objects_before.json',objects)
        binaries = [BUILD/'SU2_CFD/src/libSU2core.a', BUILD/'Common/src/libSU2Common.a', BUILD/'SU2_CFD/src/SU2_CFD.exe']
        binary_locks = {str(path):sha(path) for path in binaries}; save(folder/'shared_binary_locks.json',binary_locks)
        for path, raw in snapshots.items():
            path.resolve().relative_to(SOURCE.resolve())
            target = folder/'executed_sources'/path.relative_to(SOURCE)
            target.parent.mkdir(parents=True,exist_ok=True); target.write_bytes(raw)
            path.write_bytes(raw)
        save(folder/'executed_source_tree.json',{str(path):sha(path) for path in snapshots})
        command = [str(BIN/'g++.exe'),'-std=c++17','-O3','-fopenmp','-ffast-math','-fno-finite-math-only',
            '-I'+str(SOURCE/'SU2_CFD/include'),'-I'+str(BUILD/'Common/src'),
            '-I'+str(SOURCE/'externals/eigen'),'-I'+str(SOURCE/'externals/mel'),
            str(ROOT/'tests/testSU2SpatialStep.cpp'),str(SOURCE/FILE),'-Wl,--start-group',
            str(BUILD/'SU2_CFD/src/libSU2core.a'),str(BUILD/'Common/src/libSU2Common.a'),
            '-Wl,--end-group','-o',str(folder/'fixture.exe')]
        save(folder/'command.json',dict(compile=command,configPath=str(CONFIG),configSHA256=sha(CONFIG)))
        assert time.monotonic()-started <= 10, '编译前辅助工作超过独立10秒预算'
        run(folder,'compile',command,40)
        compile_seconds = read(folder/'compile_process.json')['elapsedSeconds']
        assert time.monotonic()-started-compile_seconds <= 10, '测试前辅助工作超预算'
        run(folder,'test',[sys.executable,str(Path(__file__).resolve()),'--worker','--output',str(folder)],10)
        test_seconds = read(folder/'test_process.json')['elapsedSeconds']
        assert all(sha(path)==hashlib.sha256(raw).hexdigest() for path,raw in snapshots.items()), '夹具改变了待恢复执行源'
        report['actualClass'] = read(folder/'results.json'); assert report['actualClass']['passed']
        report['events'] = verify_events(folder)
        report['sharedBinariesStable'] = {path:sha(Path(path))==digest for path,digest in binary_locks.items()}
        assert all(report['sharedBinariesStable'].values())
        report['passed'] = True
    except Exception as error:
        report['error'] = type(error).__name__+': '+str(error)
        raise
    finally:
        # 单文件恢复范围逐项限制在原SOURCE；不递归移动/删除，不碰OpenCode或用户文件。
        restored = {}
        for path, raw in originals.items():
            path.resolve().relative_to(SOURCE.resolve())
            if raw is None:
                if path.exists(): path.unlink()
                restored[str(path)] = not path.exists()
            else:
                path.write_bytes(raw); restored[str(path)] = path.read_bytes()==raw
        save(folder/'source_restore.json',restored)
        stable_objects = {name:sha(BUILD/name)==digest for name,digest in objects.items()}
        save(folder/'objects_unchanged.json',stable_objects)
        stable_inputs = {path:sha(Path(path))==digest for path,digest in locks.items()}
        total = time.monotonic()-started
        # 编译/测试失败也使用真实回执，不将其执行时间错误归为辅助检查。
        if (folder/'compile_process.json').exists(): compile_seconds=read(folder/'compile_process.json').get('elapsedSeconds',0)
        if (folder/'test_process.json').exists(): test_seconds=read(folder/'test_process.json').get('elapsedSeconds',0)
        report.update(inputSHA256=locks,inputsStable=stable_inputs,sourcesRestored=restored,
            engineObjectsUnchanged=len(stable_objects) if all(stable_objects.values()) else False,
            elapsedSeconds=total,compileSeconds=compile_seconds,testSeconds=test_seconds,
            miscSeconds=total-compile_seconds-test_seconds)
        safe = all(restored.values()) and all(stable_objects.values()) and all(stable_inputs.values())
        budget = compile_seconds<=40 and test_seconds<=10 and report['miscSeconds']<=10 and total<=60
        report['passed'] = bool(report['passed'] and safe and budget)
        report['budgetPassed']=budget; save(folder/'report.json',report)
        assert safe, '共享源/对象或输入未保持；准确记录已保存'
        assert budget, '阶段A实际执行超过原预算；不得追加或改写'
    assert report['passed']
    print(dict(passed=True,legalCases=report['actualClass']['legalCases'],limitedCases=report['actualClass']['limitedCases'],
        zeroStepCases=report['actualClass']['zeroStepCases'],elapsedSeconds=report['elapsedSeconds']))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--output',default=str(DEFAULT_OUTPUT))
    parser.add_argument('--worker',action='store_true'); args=parser.parse_args()
    if args.worker:
        folder=Path(args.output).resolve(); folder.relative_to(ROOT/'outputs')
        raise SystemExit(subprocess.call([str(folder/'fixture.exe'),str(CONFIG),str(folder/'results.json')],cwd=str(folder)))
    main(args.output)
