"""模块用途：独立比较新输出包与原缺首帧任务；严格字节/全部数值列相等，温域另判。"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf8')


def numeric_comparison(previous, current, expected_columns=None):
    """模块一：原字节独立比较；数值差只报告，不能取代或放宽逐字相等门槛。"""
    raw_before, raw_after = previous.read_bytes(), current.read_bytes()
    parsed = []
    for path in (previous, current):
        with path.open(encoding='utf8', newline='') as stream:
            reader = csv.reader(stream, skipinitialspace=True)
            columns = [value.strip() for value in next(reader)]
            data = np.asarray([[float(value) for value in row] for row in reader], dtype=float)
        assert data.ndim == 2 and data.shape[1] == len(columns) and np.isfinite(data).all()
        parsed.append((columns, data))
    columns, before = parsed[0]
    other_columns, after = parsed[1]
    assert columns == other_columns and before.shape == after.shape, '两任务数值表结构不同'
    if expected_columns is not None:
        assert len(columns) == expected_columns
    delta = after - before
    scale = np.maximum(1., np.maximum(abs(before), abs(after)))
    per_column = []
    for index, column in enumerate(columns):
        raw_index = int(np.argmax(abs(delta[:, index])))
        scaled_index = int(np.argmax(abs(delta[:, index]) / scale[:, index]))
        per_column.append(dict(column=column, exactNumericEqual=bool(np.array_equal(before[:, index], after[:, index])),
            differentNumericValues=int(np.count_nonzero(delta[:, index])),
            maximumAbsoluteDifference=float(abs(delta[raw_index, index])),
            maximumAbsoluteDifferenceRow=raw_index,
            valuesAtMaximumAbsoluteDifference=[float(before[raw_index, index]), float(after[raw_index, index])],
            maximumScaledDifference=float(abs(delta[scaled_index, index]) / scale[scaled_index, index]),
            maximumScaledDifferenceRow=scaled_index))
    first = next((index for index, (a, b) in enumerate(zip(raw_before, raw_after)) if a != b), None)
    if first is None and len(raw_before) != len(raw_after):
        first = min(len(raw_before), len(raw_after))
    return dict(previousFile=str(previous), currentFile=str(current), previousSHA256=sha(previous),
        currentSHA256=sha(current), previousBytes=len(raw_before), currentBytes=len(raw_after),
        exactByteEqual=raw_before == raw_after, firstDifferingByteOffset=first,
        exactAllColumnsNumericEqual=bool(np.array_equal(before, after)), rows=len(before), columns=len(columns),
        maximumAbsoluteDifference=float(abs(delta).max()), maximumScaledDifference=float((abs(delta) / scale).max()),
        scaledDifferenceDefinition='abs(current-previous)/max(1,abs(previous),abs(current)); diagnostic only',
        loosenedToleranceApplied=False, perColumn=per_column)


def compare(folder, previous, audit_path, out):
    """模块二：绑定真实终态、旧失败及严格热审核，再完成共享2/3和完整历史比较。"""
    started = time.monotonic()
    folder, previous, audit_path, out = [Path(path).resolve() for path in (folder, previous, audit_path, out)]
    out.relative_to(ROOT / 'outputs')
    out.mkdir(parents=True, exist_ok=False)
    locks, rebinding = {}, []
    report = dict(comparisonCompleted=False, comparisonPassed=False, exactOutputEqualityPassed=False,
        physicalQualified=False, experimentalQualification=False, newQualifiedArea_m2=0,
        oldMissingFirstFramePreserved=True, oldSupervisorErrorPreserved=True,
        newCFDTasks=0, builds=0, originalReceiptsModified=False)

    def bind(path, expected=None):
        path = Path(path).resolve()
        digest = sha(path)
        if expected is not None:
            assert digest == expected, '输入身份不符: ' + str(path)
        if str(path) in locks:
            assert locks[str(path)] == digest
        locks[str(path)] = digest
        return digest

    def bind_receipt_input(receipt, path):
        path = Path(path).resolve()
        matches = [value for key, value in receipt['inputSHA256'].items() if Path(key).resolve() == path]
        assert len(matches) == 1, '原运行回执缺输入绑定: ' + str(path)
        return bind(path, matches[0])

    try:
        source = Path(__file__).resolve()
        bind(source)
        (out / ('executed_' + source.name)).write_bytes(source.read_bytes())
        for path in (folder / 'engine_process.json', previous / 'engine_process.json', audit_path):
            bind(path)
        current, old, audit = read(folder / 'engine_process.json'), read(previous / 'engine_process.json'), read(audit_path)
        assert current['status'] == old['status'] == 'terminal' and current['exitCode'] == old['exitCode'] == 0
        assert not any(current.get(key, False) for key in ('timedOut', 'cancelled', 'outputBudgetExceeded', 'supervisorError'))
        assert not any(old.get(key, False) for key in ('timedOut', 'cancelled', 'outputBudgetExceeded'))
        assert old['supervisorError'].startswith('Checkpoint audit failed:')
        assert current['checkpoints']['normalScheduleComplete'] and current['checkpoints']['normalFinalPairPresent']
        assert [item['localIteration'] for item in current['checkpoints']['complete']] == [0, 2, 3]
        assert current['checkpoints']['expectedLocalIterations'] == [0, 2, 3]
        assert current['checkpoints']['missingScheduledLocalIterations'] == []
        assert [item['localIteration'] for item in old['checkpoints']['complete']] == [2, 3]
        assert old['checkpoints']['missingScheduledLocalIterations'] == [0]
        assert old['checkpoints']['normalScheduleComplete'] is False
        assert not (previous / 'restart_flow_000000.csv').exists()
        assert not (previous / 'surface_flow_000000.vtk').exists()
        for receipt in (current, old):
            assert set(receipt['inputsStable']) == set(receipt['inputSHA256']) and all(receipt['inputsStable'].values())

        # 已发生旧任务的准确源使用其事先保存快照；不以当前分派源追改旧回执。
        snapshot = ROOT / 'outputs/v174_validation/small02_source_snapshot'
        bind(snapshot / 'engine_process.json', sha(previous / 'engine_process.json'))
        for name, digest in old['inputSHA256'].items():
            path = Path(name).resolve()
            if path.is_file() and sha(path) == digest:
                bind(path, digest)
            else:
                relative = path.relative_to(ROOT)
                archived = snapshot / 'project_sources' / relative
                bind(archived, digest)
                rebinding.append(dict(originalPath=str(path), exactArchivedPath=str(archived), SHA256=digest))
        for name, digest in current['inputSHA256'].items():
            bind(name, digest)

        # 绑定根代理刚完成的严格审核；本比较不重复热计算或以研究读取绕过生产入口。
        assert audit['reviewCompleted'] and audit['normalRunCompleted'] and audit['interfacePassed']
        assert audit['savedFrameInterfacePassed'] and all(audit['inputsStable'].values())
        assert [frame['localIteration'] for frame in audit['frames']] == [0, 2, 3]
        for name, digest in audit['inputSHA256'].items():
            bind(name, digest)
        assert audit['inputSHA256'][str(folder / 'engine_process.json')] == sha(folder / 'engine_process.json')
        for frame in audit['frames']:
            assert frame['strictThermalConsistencyPassed']
            assert max(frame['energyMaximumRelativeError'], frame['eosMaximumRelativeError'],
                frame['viscosityMaximumRelativeError'], frame['constantCpAlgebraMaximumRelativeError']) < 2e-11
            assert frame['nativeHeatMaximumScaledError'] < 2e-10 and frame['powerScaledDifference'] < 2e-9
            for filename, digest_key in [('wallCSV', 'wallCSV_SHA256'),
                    ('measurementJSON', 'measurementSHA256'), ('segmentCSV', 'segmentCSV_SHA256')]:
                bind(audit_path.parent / frame[filename], frame[digest_key])
        old_review_path = ROOT / 'outputs/v174_validation/review_small02/report.json'
        bind(old_review_path)
        old_review = read(old_review_path)
        assert old_review['researchReviewPassed'] and old_review['normalTaskPassed'] is False
        assert old_review['actualSupervisorError'] == old['supervisorError']
        assert old_review['inputSHA256'][str(previous / 'engine_process.json')] == sha(previous / 'engine_process.json')

        identities = []
        for name in ('case.cfg', 'mesh.su2', 'geometry.json', 'mesh_arrays.npz', 'wall_temperature.csv', 'initial_flow.csv'):
            left = bind_receipt_input(old, previous / name)
            right = bind_receipt_input(current, folder / name)
            identities.append(dict(file=name, previousSHA256=left, currentSHA256=right, exactByteEqual=left == right))
        old_manifest, manifest = read(previous / 'case_manifest.json'), read(folder / 'case_manifest.json')
        assert manifest['publicBenchmark'] == old_manifest['publicBenchmark']
        assert manifest['initialization']['priorSuccessfulUpdates'] == old_manifest['initialization']['priorSuccessfulUpdates'] == 0
        assert manifest['runtimeIdentity'] != old_manifest['runtimeIdentity']
        assert manifest['nativeIdealGasOutputContract']['nativeSharedFrameEqualityVerified'] is False
        assert manifest['nativeIdealGasOutputContract']['numericalUpdateSourceUnchanged'] is True
        assert manifest['nativeIdealGasOutputContract']['outputScheduleChanged'] is True
        comparisons = []
        surfaces = []
        for step in (2, 3):
            name = 'restart_flow_%06d.csv' % step
            bind(previous / name, old['outputSHA256'][name])
            bind(folder / name, current['outputSHA256'][name])
            comparisons.append(dict(localIteration=step, **numeric_comparison(previous / name, folder / name, 27)))
            surface = 'surface_flow_%06d.vtk' % step
            before = bind(previous / surface, old['outputSHA256'][surface])
            after = bind(folder / surface, current['outputSHA256'][surface])
            surfaces.append(dict(localIteration=step, previousSHA256=before, currentSHA256=after, exactByteEqual=before == after))
        for case, receipt in ((previous, old), (folder, current)):
            bind(case / 'history.csv', receipt['outputSHA256']['history.csv'])
            bind(case / 'rah_wall_temperature_applied.csv', receipt['outputSHA256']['rah_wall_temperature_applied.csv'])
        history = numeric_comparison(previous / 'history.csv', folder / 'history.csv')
        assert history['rows'] == 4
        passed = all(item['exactByteEqual'] for item in identities)
        passed = passed and all(item['exactByteEqual'] and item['exactAllColumnsNumericEqual'] for item in comparisons)
        passed = passed and history['exactByteEqual'] and history['exactAllColumnsNumericEqual']
        report.update(comparisonCompleted=True, comparisonPassed=passed, exactOutputEqualityPassed=passed,
            currentNativeExitCode=0, currentNormalTaskPassed=True, currentCompleteSavedLocalIterations=[0, 2, 3],
            previousNativeExitCode=0, previousNormalTaskPassed=False, previousSupervisorError=old['supervisorError'],
            previousMissingScheduledLocalIterations=[0], physicalInputIdentities=identities,
            previousRuntimeIdentity=old_manifest['runtimeIdentity'], currentRuntimeIdentity=manifest['runtimeIdentity'],
            sharedSavedFields=comparisons, sharedSurfaceByteComparisons=surfaces, fullHistory=history,
            appliedTemperatureReceiptByteEqual=sha(previous / 'rah_wall_temperature_applied.csv') == sha(folder / 'rah_wall_temperature_applied.csv'),
            strictThermalAuditPassed=True, boundThermalAudit=str(audit_path), boundThermalAuditSHA256=sha(audit_path),
            allSavedFramesWithinDevelopmentModelDomain=audit['allSavedFramesWithinDevelopmentModelDomain'],
            developmentTemperatureRangeK=[150., 360.], modelDomainSeparateFromOutputAndThermalInterface=True,
            temperatureFrames=[{key: frame[key] for key in ('localIteration', 'temperatureRangeK',
                'developmentModelDomainPassed', 'belowDevelopmentRangeNodes', 'aboveDevelopmentRangeNodes',
                'originalResidualGatePassed', 'givenTwMaximumFluidDifferenceK')} for frame in audit['frames']],
            exactHistoricalInputRebinding=rebinding, noToleranceRelaxation=True,
            onlySharedSavedStatesCompared=True, arbitraryFutureNumericalEqualityProven=False)
        assert passed, '逐字/全列严格相等未过；原始差异已保留，禁止改成宽容差通过'
    except Exception as error:
        report.update(comparisonPassed=False, error=type(error).__name__ + ': ' + str(error), traceback=traceback.format_exc())
        raise
    finally:
        stable = {name: Path(name).is_file() and sha(Path(name)) == digest for name, digest in locks.items()}
        if not all(stable.values()):
            report.update(comparisonPassed=False, inputMutationDetected=True)
        report.update(inputSHA256=locks, inputsStable=stable, elapsedSeconds=time.monotonic() - started)
        save(out / 'input_manifest.json', locks)
        save(out / 'report.json', report)
    assert report['comparisonPassed'] and all(report['inputsStable'].values())
    print(dict(comparisonPassed=True, strictThermalAuditPassed=True, sharedFrames=[2, 3],
        allSavedFramesWithinDevelopmentModelDomain=report['allSavedFramesWithinDevelopmentModelDomain'],
        elapsedSeconds=report['elapsedSeconds']))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--previous', required=True)
    parser.add_argument('--audit', required=True)
    parser.add_argument('--output', required=True)
    arguments = parser.parse_args()
    compare(arguments.folder, arguments.previous, arguments.audit, arguments.output)
