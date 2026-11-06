"""独占离线夹具验证正常/异常保存帧身份；不修改原任务、不重算 CFD。"""
import argparse
import copy
import csv
import os
from pathlib import Path
import re
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from nacaSavedFrameStudy import select_saved_frames
from su2ThermalRuntime import read, sha
from validateIdealAir import save


def check(output):
    """模块一：保存准确源；模块二：只读复用完整小场；模块三：正向与拒绝测试。"""
    started = time.monotonic()
    out = Path(output).resolve()
    out.relative_to(ROOT / 'outputs')
    out.mkdir(parents=True, exist_ok=False)
    source = ROOT / 'outputs/v171_validation/spatial_scalar01'
    receipt = read(source / 'engine_process.json')
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0
    assert receipt['checkpoints']['normalScheduleComplete'] and all(receipt['inputsStable'].values())
    frames = [source / row[name] for row in receipt['checkpoints']['complete'] for name in ('volume', 'surface')]
    assert [row['localIteration'] for row in receipt['checkpoints']['complete']] == [0, 2, 3]
    originals = [source / name for name in ('engine_process.json', 'case_manifest.json', 'mesh.su2', 'history.csv')] + frames
    production = [ROOT / 'scripts' / name for name in (
        'nacaSavedFrameStudy.py', 'su2NACASpatialInterface.py', 'su2NACALaminarInterface.py',
        'su2ThermalRuntime.py', 'auditSU2SSTState.py', 'validateIdealAir.py')]
    locks = {str(path): sha(path) for path in originals + production + [Path(__file__).resolve()]}
    save(out / 'input_manifest.json', locks)
    for path in production + [Path(__file__).resolve()]:
        (out / ('executed_' + path.name)).write_bytes(path.read_bytes())
    with (source / 'history.csv').open(encoding='utf8') as stream:
        raw = list(csv.reader(stream, skipinitialspace=True))
    header, rows = raw[0], raw[1:]
    step_column = [name.strip() for name in header].index('Inner_Iter')
    results = []

    def fixture(name, steps=(0, 2, 3), abnormal=False, half_last=False, changed_rows=None):
        folder = out / name
        folder.mkdir(exist_ok=False)
        for file in ('case_manifest.json', 'mesh.su2'):
            os.link(source / file, folder / file)
        local = copy.deepcopy(receipt)
        for file in frames:
            step = int(re.search(r'_(\d+)\.', file.name)[1])
            if step in steps and not (half_last and step == 3 and file.suffix == '.vtk'):
                os.link(file, folder / file.name)
            else:
                local['outputSHA256'].pop(file.name, None)
        # 可变历史单独创建；所有硬链接仅供读取，绝不改写或设置其属性。
        with (folder / 'history.csv').open('w', encoding='utf8', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            writer.writerows(rows if changed_rows is None else changed_rows)
        local['outputSHA256']['history.csv'] = sha(folder / 'history.csv')
        if abnormal:
            local.update(exitCode=1, timedOut=True, cancelled=False)
        return folder, local

    def exercise(name, expected, opt_in=False, mutate=None, verify=None, **kwargs):
        folder, local = fixture(name, **kwargs)
        if mutate:
            mutate(local)
        save(folder / 'fixture_receipt.json', local)
        before = copy.deepcopy(local)
        file_locks = {str(path): sha(path) for path in folder.iterdir() if path.is_file()}
        try:
            checkpoints, identity = select_saved_frames(folder, local, saved_frames_only=opt_in)
        except (AssertionError, KeyError, OSError, ValueError, RuntimeError) as exc:
            assert not expected, '应该接受的夹具被拒绝: ' + name + ': ' + str(exc)
            result = dict(name=name, passed=True, accepted=False, reason=type(exc).__name__ + ': ' + str(exc))
        else:
            assert expected, '应拒绝的夹具被接受: ' + name
            assert identity['savedFrameSelectionPassed'] and not identity['thermalInterfaceAudited']
            assert not identity['physicalSolutionQualified'] and not identity['experimentalQualification']
            assert not identity['abnormalInstantStateAvailable'] and not identity['unsavedTerminalStateKnown']
            assert identity['newQualifiedArea_m2'] == 0
            if verify:
                verify(checkpoints, identity)
            result = dict(name=name, passed=True, accepted=True, checkpoints=checkpoints, identity=identity)
        assert local == before, '选择函数修改了传入回执'
        assert all(sha(Path(path)) == digest for path, digest in file_locks.items()), '选择函数修改了夹具文件'
        save(folder / 'result.json', result)
        results.append(result)

    def require(condition, message):
        assert condition, message

    exercise('normal-complete', True, verify=lambda c, i: require(
        i['normalRunCompleted'] and c['normalScheduleComplete'] and i['auditScope'] == 'normal-complete-saved-frames',
        '正常身份丢失'))
    exercise('normal-complete-explicit', True, opt_in=True, verify=lambda c, i: require(
        i['normalRunCompleted'] and c['normalScheduleComplete'], '开关不能把正常完整任务标成异常'))
    exercise('abnormal-default-rejected', False, abnormal=True, steps=(0, 2))
    exercise('abnormal-saved-only', True, abnormal=True, steps=(0, 2), opt_in=True,
        mutate=lambda r: r.update(checkpoints={'complete': [], 'normalScheduleComplete': True}),
        verify=lambda c, i: require([p['localIteration'] for p in c['complete']] == [0, 2]
            and i['lastCompleteLocalIteration'] == 2 and i['lastFiniteHistoryLocalIteration'] == 3
            and not i['normalRunCompleted'] and not i['normalScheduleComplete'] and i['actualExitCode'] == 1
            and i['terminalFlags']['timedOut'] and i['auditScope'] == 'complete-saved-frames-from-abnormal-task',
            '没有独立重清点或混淆了历史和末场'))
    exercise('abnormal-half-last', True, abnormal=True, half_last=True, opt_in=True,
        verify=lambda c, i: require([p['localIteration'] for p in c['complete']] == [0, 2]
            and [p['localIteration'] for p in i['rejectedSavedPairs']] == [3], '半写帧必须拒绝并保留旧完整帧'))
    exercise('abnormal-missing-periodic', True, abnormal=True, steps=(0,), opt_in=True,
        verify=lambda c, i: require(i['missingScheduledLocalIterations'] == [2]
            and i['lastCompleteLocalIteration'] == 0 and not i['normalScheduleComplete'], '缺周期帧必须明示'))
    for opt_in in (False, True):
        exercise('normal-missing-final-' + str(opt_in), False, steps=(0, 2), opt_in=opt_in)
        exercise('normal-missing-periodic-' + str(opt_in), False, steps=(0, 3), opt_in=opt_in)
    exercise('abnormal-history-gap', False, abnormal=True, opt_in=True,
        changed_rows=[row for row in rows if int(float(row[step_column])) != 1])
    exercise('abnormal-no-complete', False, abnormal=True, steps=(), opt_in=True)
    exercise('not-terminal', False, abnormal=True, opt_in=True, mutate=lambda r: r.update(status='running'))
    exercise('input-unstable', False, abnormal=True, opt_in=True,
        mutate=lambda r: r['inputsStable'].update({next(iter(r['inputsStable'])): False}))
    exercise('input-empty', False, abnormal=True, opt_in=True, mutate=lambda r: r.update(inputsStable={}))
    exercise('input-lock-missing', False, abnormal=True, opt_in=True,
        mutate=lambda r: r['inputsStable'].pop(next(iter(r['inputsStable']))))
    exercise('supervisor-error', False, abnormal=True, opt_in=True,
        mutate=lambda r: r.update(supervisorError='offline integrity fixture'))
    exercise('history-hash-mismatch', False, abnormal=True, opt_in=True,
        mutate=lambda r: r['outputSHA256'].update({'history.csv': '0' * 64}))
    exercise('mesh-input-hash-mismatch', False, abnormal=True, opt_in=True,
        mutate=lambda r: r['inputSHA256'].update({next(p for p in r['inputSHA256'] if Path(p).name == 'mesh.su2'): '0' * 64}))
    exercise('only-frame-hash-mismatch', False, abnormal=True, steps=(0,), opt_in=True,
        mutate=lambda r: r['outputSHA256'].update({'restart_flow_000000.csv': '0' * 64}))
    stable = {path: sha(Path(path)) == digest for path, digest in locks.items()}
    assert all(stable.values())
    report = dict(passed=True, cases=results, acceptedCases=sum(row['accepted'] for row in results),
        rejectedCases=sum(not row['accepted'] for row in results), inputSHA256=locks, inputsStable=stable,
        originalTasksModified=False, newCFDTasks=0, newPhysicalQualification=False,
        elapsedSeconds=time.monotonic() - started)
    save(out / 'report.json', report)
    print(dict(passed=True, cases=len(results), accepted=report['acceptedCases'], rejected=report['rejectedCases'],
               elapsedSeconds=report['elapsedSeconds']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    check(parser.parse_args().output)
