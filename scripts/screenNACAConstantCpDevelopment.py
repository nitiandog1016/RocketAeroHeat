"""模块用途：按v174事前计划判定唯一实际定比热候选；开发收益不等于物理资格。"""
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
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf8')


def history_table(path):
    """模块一：只收完整有限历史；缺行或损坏单列，不能缩短序列冒充连续发展。"""
    accepted, rejected = [], []
    with Path(path).open(encoding='utf8', newline='') as stream:
        reader = csv.reader(stream, skipinitialspace=True)
        names = [value.strip() for value in next(reader, [])]
        if not all(name in names for name in ('Inner_Iter', 'rms[Rho]', 'rms[RhoE]')):
            return [], [dict(line=1, reason='缺少完整有限历史表头')], False
        for line, row in enumerate(reader, 2):
            try:
                values = np.asarray(row, dtype=float)
                assert len(values) == len(names) and np.isfinite(values).all()
                item = dict(zip(names, values.tolist()))
                assert item['Inner_Iter'] == int(item['Inner_Iter'])
                accepted.append(item)
            except (AssertionError, ValueError, KeyError) as error:
                rejected.append(dict(line=line, reason=type(error).__name__ + ': ' + str(error)))
    ids = [int(row['Inner_Iter']) for row in accepted]
    return accepted, rejected, ids == list(range(len(ids)))


def relative_change(before, after, area):
    """模块二：按原生投影面积比较最后两保存壁场；零基准明确未知，不造量纲阈值。"""
    numerator = float(np.sum(area * abs(after - before)))
    denominator = float(np.sum(area * abs(before)))
    return dict(weightedAbsoluteDifference=numerator, weightedReferenceMagnitude=denominator,
        relativeL1=numerator / denominator if denominator > 0 else None,
        defined=denominator > 0, acceptanceThresholdAssigned=False)


def plot_history_temperature(rows, frames, target):
    """绘图模块：有限历史可连线；稀疏保存温区只画点，未保存的中间温度未知。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with plt.rc_context({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False}):
        figure, axes = plt.subplots(2, 1, figsize=(8, 5.7), constrained_layout=True)
        indices = [row['Inner_Iter'] for row in rows]
        for key, label, color in [('rms[Rho]', 'Density residual', '#255B8A'), ('rms[RhoE]', 'Energy residual', '#B15A22')]:
            axes[0].plot(indices, [row[key] for row in rows], color=color, linewidth=1.2, label=label)
        axes[0].axhline(-8, color='#777777', linewidth=.8, linestyle='--', label='Original -8 gate')
        axes[0].set(xlabel='Native local iteration (0 = first update)', ylabel='Residual log10')
        axes[0].legend(frameon=False, fontsize=8)
        axes[0].grid(alpha=.18)
        if frames:
            positions = [frame['localIteration'] for frame in frames]
            axes[1].plot(positions, [frame['temperatureRangeK'][0] for frame in frames], 'o', color='#255B8A', linestyle='None', label='Saved minimum')
            axes[1].plot(positions, [frame['temperatureRangeK'][1] for frame in frames], 's', color='#B15A22', linestyle='None', label='Saved maximum')
        for boundary in (150, 360):
            axes[1].axhline(boundary, color='#777777', linewidth=.8, linestyle='--')
        axes[1].set(xlabel='Complete saved local iteration', ylabel='Temperature [K]')
        axes[1].text(.5, 1.025, 'Unsaved intermediate temperatures are unknown; no interpolation.',
            transform=axes[1].transAxes, ha='center', va='bottom', fontsize=8, color='#555555')
        if frames:
            axes[1].legend(frameon=False, fontsize=8)
        axes[1].grid(alpha=.18)
        figure.suptitle('NACA native ideal-gas development | physical qualification: not established', fontsize=10)
        figure.savefig(target, dpi=160)
        plt.close(figure)


def screen(folder, audit_path, out):
    """模块三：冻结回执及严格热证据；保留异常身份，对每个事前条件分别作判定。"""
    started = time.monotonic()
    folder, audit_path, out = [Path(path).resolve() for path in (folder, audit_path, out)]
    out.relative_to(ROOT / 'outputs')
    out.mkdir(parents=True, exist_ok=False)
    locks, gates, reasons = {}, {}, []
    report = dict(reviewCompleted=False, candidateAccepted=False, physicalQualified=False,
        experimentalQualification=False, newQualifiedArea_m2=0, newCFDTasks=0, builds=0,
        automaticContinuationAuthorized=False, candidate='v174-unique-fresh-native-ideal-gas-validation-profile')

    def bind(path, expected=None):
        path = Path(path).resolve()
        value = sha(path)
        assert expected is None or value == expected, '原始身份改变: ' + str(path)
        assert str(path) not in locks or locks[str(path)] == value
        locks[str(path)] = value
        return value

    def gate(name, value, explanation):
        gates[name] = dict(passed=value, reason=explanation)
        if value is not True:
            reasons.append(dict(gate=name, status='failed' if value is False else 'pending', reason=explanation))

    try:
        source = Path(__file__).resolve()
        bind(source)
        (out / ('executed_' + source.name)).write_bytes(source.read_bytes())
        receipt_path = folder / 'engine_process.json'
        bind(receipt_path)
        receipt = read(receipt_path)
        assert receipt['status'] == 'terminal'
        manifest_path = folder / 'case_manifest.json'
        bind(manifest_path, receipt['inputSHA256'][str(manifest_path)])
        manifest = read(manifest_path)
        assert manifest['caseKind'] == 'naca-tn3344-constant-cp-spatial-wall-v1'
        assert manifest['specification']['variant'] == 'validation-profile'
        assert manifest['specification']['actualValidationMesh'] is True
        assert manifest['nodes'] == 190334 and manifest['elements'] == 179568
        assert manifest['initialization']['priorSuccessfulUpdates'] == 0
        assert manifest['initialization']['reusesPreviousConservativeState'] is False
        contract = manifest['idealGas']
        assert contract['energyReference'] == 'e=cv*T; h=cp*T; no-OTH-offset'
        assert contract['developmentTemperatureRangeK'] == [150., 360.]
        for name, digest in receipt['inputSHA256'].items():
            path = Path(name).resolve()
            if path == (ROOT / 'docs/validation/v174/PLAN.md').resolve() and sha(path) != digest:
                bind(folder / 'plan_at_prepare.md', digest)
            else:
                bind(path, digest)
        bind(folder / 'plan_at_prepare.md', manifest['planSHA256'])
        (out / 'plan_at_run.md').write_bytes((folder / 'plan_at_prepare.md').read_bytes())
        spec = manifest['specification']
        gate('registeredBudget', [spec[key] for key in ('updates', 'every', 'maximumSeconds', 'threads', 'maximumOutputBytes')]
            == [500, 100, 600, 2, 1536 * 1024**2], '固定500/100/600秒/2线程/1536MiB，不增加预算')
        stable = receipt.get('inputsStable', {})
        gate('nativeInputsStable', bool(stable) and set(stable) == set(receipt['inputSHA256']) and all(value is True for value in stable.values()),
            '原运行回执及本次读取逐项绑定')
        flags = {key: receipt.get(key, False) for key in ('timedOut', 'cancelled', 'outputBudgetExceeded', 'modelDomainExceeded')}
        assert all(isinstance(value, bool) for value in flags.values())
        normal = receipt['exitCode'] == 0 and not any(flags.values()) and not receipt.get('supervisorError')
        gate('normalNativeTerminal', normal, '必须正常原生终态；越域停止、异常和监督错误不能接续')
        checkpoints = receipt.get('checkpoints', {})
        complete = checkpoints.get('complete', [])
        gate('completeSavedSchedule', checkpoints.get('normalScheduleComplete') is True,
            '完整500计划或真正−8提前收敛的首帧/中间计划帧/末帧')

        # 守卫实际运行身份与温域审核分别判；缺守卫或未检查某完整帧不能当通过。
        guard = receipt.get('savedFieldDomainGuard', {})
        expected_guard = contract.get('savedFieldDomainGuard')
        guard_enabled = bool(expected_guard and guard.get('contract') == expected_guard)
        guard_enabled = guard_enabled and expected_guard.get('temperatureRangeK') == [150., 360.]
        gate('guardActuallyEnabled', guard_enabled and bool(guard.get('helperAttempts')),
            '要求当前回执具有准确守卫合同和真实helper尝试，不能仅靠配置开关')
        gate('guardCompletedWithoutError', bool(guard) and guard.get('error') is None and not guard.get('active', True)
            and not guard.get('cleanupErrors'), '助手错误/预算耗尽/未清理均与物理越域分开')
        checked = {int(item['step']): item for item in guard.get('checkedFrames', [])}
        covered = bool(complete)
        for pair in complete:
            step = pair['localIteration']
            for key, digest_key in [('volume', 'volumeSHA256'), ('surface', 'surfaceSHA256')]:
                bind(folder / pair[key], pair[digest_key])
                assert receipt['outputSHA256'][pair[key]] == pair[digest_key]
            item = checked.get(step)
            covered = covered and bool(item and item.get('status') == 'complete'
                and item.get('volumeSHA256') == pair['volumeSHA256'] and item.get('surfaceSHA256') == pair['surfaceSHA256'])
        gate('guardCoversEveryCompleteSavedFrame', bool(covered), '守卫已检查结果须逐帧绑定相同原始CSV/VTK')
        guard_domain = bool(checked) and all(item.get('exceeded') is False for item in checked.values())
        gate('guardModelDomain', guard_domain and not flags['modelDomainExceeded'], '150–360K开发域固定；仅限完整保存场')

        audit = read(audit_path) if audit_path.is_file() else {}
        frames = audit.get('frames', [])
        if audit:
            bind(audit_path)
            assert audit['inputSHA256'][str(receipt_path)] == sha(receipt_path)
            for name, digest in audit['inputSHA256'].items():
                bind(name, digest)
        strict = bool(audit.get('reviewCompleted') and audit.get('savedFrameInterfacePassed'))
        strict = strict and [frame['localIteration'] for frame in frames] == [pair['localIteration'] for pair in complete]
        for frame in frames:
            strict = strict and frame.get('strictThermalConsistencyPassed') is True
            strict = strict and max(frame['energyMaximumRelativeError'], frame['eosMaximumRelativeError'],
                frame['viscosityMaximumRelativeError'], frame['constantCpAlgebraMaximumRelativeError']) < 2e-11
            strict = strict and frame['nativeHeatMaximumScaledError'] < 2e-10 and frame['powerScaledDifference'] < 2e-9
        gate('strictThermalAudit', bool(strict), '能量/EOS/黏度2e−11、热梯度2e−10、功率2e−9原门槛保持')
        domain = bool(frames) and all(frame['developmentModelDomainPassed'] is True for frame in frames)
        gate('allSavedFieldsWithinDevelopmentDomain', domain, '热一致性不授予物性开发温域资格')

        history_path = folder / 'history.csv'
        rows, rejected, contiguous = [], [], False
        if history_path.is_file():
            bind(history_path, receipt['outputSHA256']['history.csv'])
            rows, rejected, contiguous = history_table(history_path)
        gate('finiteContiguousHistory', bool(rows) and contiguous and not rejected, '0起连续有限历史；拒绝缺失/重复/损坏行')
        final_residuals = {key: rows[-1][key] for key in ('rms[Rho]', 'rms[RhoE]')} if rows else {}
        last = int(rows[-1]['Inner_Iter']) if rows else None
        converged = bool(rows and last >= 10 and all(value <= -8 for value in final_residuals.values()))
        full = len(rows) == 500 and last == 499 and contiguous and not rejected
        trend = dict(available=len(rows) >= 100, windowRows=50, qualifiedSteadySolution=False)
        if len(rows) >= 100:
            trend['fields'] = {}
            for key in ('rms[Rho]', 'rms[RhoE]'):
                values = np.asarray([row[key] for row in rows])
                first, previous, final = [float(np.mean(part)) for part in (values[:50], values[-100:-50], values[-50:])]
                trend['fields'][key] = dict(first50MeanLog=first, preceding50MeanLog=previous,
                    last50MeanLog=final, firstToLastLogImprovement=first-final,
                    lastMinusPrecedingMeanLog=final-previous, improvementAtLeastHalfLog=first-final >= .5,
                    noLastWindowRebound=final <= previous)
        benefit = bool(full and trend['available'] and all(item['improvementAtLeastHalfLog'] and item['noLastWindowRebound']
            for item in trend.get('fields', {}).values()))
        gate('registeredDevelopmentBenefit', bool(converged or benefit) if len(rows) >= 100 or converged else None,
            '原−8已达，或完整500且rho/E首末50均改善≥0.5log并无末窗反弹；不足窗口保持待验')

        # 热/压力稳定指标没有新增门槛；保持其物理意义与相对变化分母。
        changes = dict(available=False, acceptanceThresholdAssigned=False)
        if len(frames) >= 2:
            pair_frames = frames[-2:]
            tables = []
            for frame in pair_frames:
                path = audit_path.parent / frame['wallCSV']
                bind(path, frame['wallCSV_SHA256'])
                with path.open(encoding='utf8') as stream:
                    tables.append(list(csv.DictReader(stream)))
            before, after = tables
            assert [item['PointID'] for item in before] == [item['PointID'] for item in after]
            area = np.asarray([float(item['nativeProjectedArea_m2']) for item in before])
            other_area = np.asarray([float(item['nativeProjectedArea_m2']) for item in after])
            assert np.array_equal(area, other_area) and np.all(area > 0)
            changes = dict(available=True, localIterations=[frame['localIteration'] for frame in pair_frames],
                weight='nativeProjectedArea_m2', summedProjectedAreaM2=float(area.sum()), acceptanceThresholdAssigned=False)
            for column, label in [('pressure_Pa', 'pressure'), ('qIntoWall_Wm2', 'heatFlux')]:
                left = np.asarray([float(item[column]) for item in before])
                right = np.asarray([float(item[column]) for item in after])
                changes[label] = relative_change(left, right, area)
            left, right = [frame['nativeNodalPowerIntoWallW'] for frame in pair_frames]
            changes['power'] = dict(previousIntoWallW=left, currentIntoWallW=right, absoluteDifferenceW=abs(right-left),
                relativeAbsoluteChange=abs(right-left)/abs(left) if left != 0 else None, defined=left != 0)
        last_frame = frames[-1] if frames else {}
        yplus = last_frame.get('independentFirstNeighborYPlus99')
        physical = dict(originalResidualMinus8Passed=converged, independentYPlus99=yplus,
            independentYPlus99AtMost1Passed=yplus <= 1 if yplus is not None else None,
            lastTwoSavedWallFieldChanges=changes, heatStabilityQualified=False, gridAndOuterDomainQualified=False,
            localExternalEdgeQualified=False, experimentalComparisonQualified=False,
            givenTwIsSolvedSolidWallTemperature=False)
        report.update(reviewCompleted=True, candidateAccepted=all(item['passed'] is True for item in gates.values()),
            gates=gates, decisionReasons=reasons, terminalExitCode=receipt['exitCode'], terminalFlags=flags,
            originalSupervisorError=receipt.get('supervisorError'), originalStopReason=receipt.get('stopReason'),
            domainGuardTerminationRequested=receipt.get('domainGuardTerminationRequested'),
            terminalOnlyDomainDetection=guard.get('terminalOnlyDetection'), guardActualStatus=guard,
            completeSavedLocalIterations=[pair['localIteration'] for pair in complete],
            historyRows=len(rows), lastFiniteHistoryLocalIteration=last, rejectedHistoryRows=rejected,
            finalResiduals=final_residuals, originalResidualConvergencePassed=converged,
            full500UpdatesObserved=full, developmentWindowTrend=trend,
            originalPhysicalGates=physical, strictThermalAuditPath=str(audit_path),
            temperatureFrames=[{key: frame[key] for key in ('localIteration', 'temperatureRangeK',
                'developmentModelDomainPassed', 'aboveDevelopmentRangeNodes', 'belowDevelopmentRangeNodes')} for frame in frames],
            decision='retain-single-development-evidence-no-automatic-continuation' if all(item['passed'] is True for item in gates.values())
                else 'candidate-not-accepted-no-continuation-or-parameter-scan')

        # 模块四：克制的研究图，只画实际有限历史/保存温区，不插值未知瞬态或理论真值。
        if rows:
            plot_history_temperature(rows, frames, out / 'development_history_temperature.png')
            report['figure'] = 'development_history_temperature.png'
            report['figureSHA256'] = sha(out / report['figure'])
    except Exception as error:
        report.update(reviewCompleted=False, candidateAccepted=False,
            error=type(error).__name__ + ': ' + str(error), traceback=traceback.format_exc())
        raise
    finally:
        stable = {path: Path(path).is_file() and sha(Path(path)) == value for path, value in locks.items()}
        if not all(stable.values()):
            report.update(reviewCompleted=False, candidateAccepted=False, inputMutationDetected=True)
        report.update(inputSHA256=locks, inputsStable=stable, elapsedSeconds=time.monotonic() - started)
        save(out / 'input_manifest.json', locks)
        save(out / 'report.json', report)
    assert report['reviewCompleted'] and all(report['inputsStable'].values())
    print(dict(reviewCompleted=True, candidateAccepted=report['candidateAccepted'], historyRows=report['historyRows'],
        lastFiniteHistoryLocalIteration=report['lastFiniteHistoryLocalIteration'], elapsedSeconds=report['elapsedSeconds']))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--audit', required=True)
    parser.add_argument('--output', required=True)
    arguments = parser.parse_args()
    screen(arguments.folder, arguments.audit, arguments.output)
