"""模块用途：只读汇总v173唯一实际候选的事前发展门槛；不启动CFD、不读取大场。"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time


ROOT = Path(__file__).resolve().parents[1]
FLAVOR = 'rah-ideal-air-spatial-step-v1'
MAXIMUM_OUTPUT_BYTES = 1024 * 1024 * 1024
THERMAL_LIMITS = {
    'energyMaximumRelativeError': 2e-11,
    'eosMaximumRelativeError': 2e-11,
    'paperCpMaximumRelativeError': 2e-8,
    'viscosityMaximumRelativeError': 2e-11,
    'nativeHeatMaximumScaledError': 2e-10,
    'powerScaledDifference': 2e-9,
}


def sha(path):
    """模块一：只散列本脚本实际读取的小报告，不重读报告引用的大场。"""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def stable_record(report):
    hashes = report.get('inputSHA256')
    stable = report.get('inputsStable')
    return (isinstance(hashes, dict) and bool(hashes) and isinstance(stable, dict)
            and set(stable) == set(hashes) and all(value is True for value in stable.values()))


def gate(status, value=None, threshold=None, evidence=None, reason=None):
    assert status in ('passed', 'failed', 'pending')
    return dict(status=status, value=value, threshold=threshold, evidence=evidence, reason=reason)


def screen(folder, audit, development, verification, output):
    """模块二：绑定终态回执与三份既存报告；缺报告记待验，身份冲突记失败。"""
    started = time.monotonic()
    folder = Path(folder).resolve()
    output = Path(output).resolve()
    output.relative_to(ROOT / 'outputs')
    assert not output.exists(), '必须使用新的筛选输出目录，不能覆盖旧证据'
    locks, reports, binding_errors, missing = {}, {}, [], []

    def read_small(path, required=False):
        path = Path(path).resolve()
        if not path.is_file():
            if required:
                raise FileNotFoundError('缺必要终态回执: ' + str(path))
            missing.append(str(path))
            return None
        data = path.read_bytes()
        assert len(data) <= 32 * 1024 * 1024, '只允许读取小型JSON报告: ' + str(path)
        locks[str(path)] = hashlib.sha256(data).hexdigest()
        return json.loads(data.decode('utf-8-sig'))

    def report_path(path):
        path = Path(path).resolve()
        return path / 'report.json' if path.is_dir() or path.suffix.lower() != '.json' else path

    def bound(report, path):
        target = Path(path).resolve()
        values = [value for key, value in report.get('inputSHA256', {}).items()
                  if Path(key).resolve() == target]
        return bool(values) and len(set(values)) == 1 and values[0] == locks[str(target)]

    receipt_path = folder / 'engine_process.json'
    receipt = read_small(receipt_path, required=True)
    assert isinstance(receipt, dict) and receipt.get('status') == 'terminal', '禁止筛选运行中的任务'
    assert type(receipt.get('exitCode')) is int, '终态回执必须包含明确原生退出码'
    assert receipt.get('runtimeFlavor') == FLAVOR, '本筛选仅适用于v173空间壁温限步包'
    manifest_path = folder / 'case_manifest.json'
    manifest = read_small(manifest_path, required=True)
    spec = manifest.get('specification', {})
    assert spec.get('actualValidationMesh') is True, '禁止把小型组合当作实际候选'
    assert spec.get('updates') == 1000 and spec.get('every') == 250
    assert spec.get('maximumOutputBytes') == MAXIMUM_OUTPUT_BYTES
    assert bound(receipt, manifest_path), '实际回执未绑定准确任务清单'
    if manifest.get('requiredRuntimeFlavor') != FLAVOR:
        binding_errors.append('case_manifest.requiredRuntimeFlavor与实际包不同')
    paths = {'audit': report_path(audit), 'development': report_path(development),
             'verification': report_path(verification)}
    for name, path in paths.items():
        item = read_small(path)
        reports[name] = item
        if item is None:
            continue
        if not isinstance(item, dict):
            binding_errors.append(name + ': 报告不是对象')
            reports[name] = None
            continue
        if not stable_record(item):
            binding_errors.append(name + ': 原审核未记录完整稳定输入')
        if not bound(item, receipt_path):
            binding_errors.append(name + ': 未绑定本次准确原生回执')
    a, d, v = (reports[name] or {} for name in ('audit', 'development', 'verification'))
    for name in ('development', 'verification'):
        if reports[name] is not None and a and not bound(reports[name], paths['audit']):
            binding_errors.append(name + ': 未绑定选定热审核报告')
    if a and not bound(a, manifest_path):
        binding_errors.append('audit: 未绑定本次准确任务清单')

    # 模块三：正常终态、容量与完整保存只认原生回执，不能用离线清点冒充事中监督。
    flags = {key: receipt.get(key, False) for key in (
        'timedOut', 'cancelled', 'outputBudgetExceeded', 'memoryBudgetExceeded', 'supervisorError')}
    normal = receipt['exitCode'] == 0 and not any(flags.values())
    cp = receipt.get('checkpoints') or {}
    pairs = cp.get('complete', [])
    steps = [item['localIteration'] for item in pairs]
    assert steps == sorted(set(steps)), '完整帧身份重复或乱序'
    last = cp.get('lastFiniteHistoryLocalIteration')
    last_complete = steps[-1] if steps else None
    gates = {}
    gates['normalTermination'] = gate('passed' if normal else 'failed',
        dict(exitCode=receipt['exitCode'], flags=flags), 'exitCode=0且无终止标志', 'engine_process.json',
        None if normal else '实际任务未正常完成；退出码不自动解释为某种具体异常')
    gates['stableNativeInputs'] = gate('passed' if stable_record(receipt) else 'failed',
        dict(count=len(receipt.get('inputSHA256', {})), allStable=stable_record(receipt)),
        '非空逐项输入稳定回执', 'engine_process.json.inputsStable')
    gates['reportEvidenceBinding'] = gate('failed' if binding_errors else 'pending' if missing else 'passed',
        dict(errors=binding_errors, missing=missing), '三份报告绑定同一回执/热审核且来源稳定',
        {name: str(path) for name, path in paths.items()})
    verified = bool(v.get('reviewCompleted') is True and v.get('passed') is True)
    gates['sourceAndRunVerification'] = gate('passed' if verified else 'failed' if v else 'pending',
        dict(reviewCompleted=v.get('reviewCompleted'), passed=v.get('passed')),
        '只读来源/运行合同复核通过；不等于求解正常', 'verification.report')
    budget_fields = ('generatedOutputBytes' in receipt and 'outputBudgetExceeded' in receipt)
    byte_count = receipt.get('generatedOutputBytes')
    budget_values_ok = (type(byte_count) is int and 0 <= byte_count <= MAXIMUM_OUTPUT_BYTES
                        and receipt.get('outputBudgetExceeded') is False)
    capacity_verified = (verified and v.get('nativeOutputBudgetRecorded') is True
        and v.get('nativeOutputCapacityEnforcementRecorded') is True
        and v.get('nativeGeneratedOutputBytes') == byte_count
        and v.get('maximumOutputBytes') == MAXIMUM_OUTPUT_BYTES)
    capacity_status = ('failed' if not budget_fields or not budget_values_ok else
                       'pending' if not v else 'passed' if capacity_verified else 'failed')
    gates['nativeOutputCapacity'] = gate(capacity_status,
        dict(nativeMetadataPresent=budget_fields, generatedOutputBytes=byte_count,
             enforcementPathRecorded=v.get('nativeOutputCapacityEnforcementRecorded'),
             overCapacityStopDemonstrated=v.get('nativeOutputCapacityEnforcementDemonstrated')),
        dict(maximumOutputBytes=MAXIMUM_OUTPUT_BYTES, outputBudgetExceeded=False),
        'engine_process.json + verification.report',
        '检查事中监督路径及终态容量；不要求本次故意超容量，也不宣称超容量停止已实测')
    cp_recorded = 'checkpoints' in receipt
    cp_verified = verified and v.get('nativeCheckpointsRecorded') is True
    gates['nativeCheckpointRecord'] = gate('failed' if not cp_recorded else 'pending' if not v
        else 'passed' if cp_verified else 'failed',
        dict(nativePresent=cp_recorded, completeLocalIterations=steps), '原生完整帧清点与独立复核相同',
        'engine_process.json.checkpoints + verification.report')
    schedule_ok = (normal and cp.get('normalScheduleComplete') is True
                   and cp.get('normalFinalPairPresent') is True and bool(pairs))
    gates['normalCompleteSchedule'] = gate('passed' if schedule_ok else 'failed',
        dict(normalFinalPairPresent=cp.get('normalFinalPairPresent'),
             normalScheduleComplete=cp.get('normalScheduleComplete'),
             expected=cp.get('expectedLocalIterations'), missing=cp.get('missingScheduledLocalIterations')),
        '正常终态及应有保存完整', 'engine_process.json.checkpoints')

    # 模块四：严格热检查逐保存帧保持原阈值，异常任务也可以有合格的已存首帧。
    thermal_frames = a.get('frames', [])
    thermal_details = []
    thermal_ok = bool(a and pairs and a.get('savedFrameInterfacePassed') is True
                      and [item.get('localIteration') for item in thermal_frames] == steps)
    for frame in thermal_frames:
        checks = {key: bool(finite(frame.get(key)) and 0 <= frame[key] < limit)
                  for key, limit in THERMAL_LIMITS.items()}
        frame_ok = frame.get('interfacePassed') is True and all(checks.values())
        thermal_ok = thermal_ok and frame_ok
        thermal_details.append(dict(localIteration=frame.get('localIteration'),
            passed=frame_ok, checks=checks, values={key: frame.get(key) for key in THERMAL_LIMITS}))
    gates['allSavedStrictThermalAudit'] = gate('pending' if not a or not pairs else
        'passed' if thermal_ok else 'failed', thermal_details, THERMAL_LIMITS, 'audit.report.frames')

    # 模块五：发展量只使用明确身份的实际窗口和保存帧，不把首帧移名为末帧。
    history = d.get('history', {})
    windows = history.get('windows', {})
    first_window, last_window = windows.get('first20'), windows.get('last20')
    windows_ok = bool(first_window and last_window and
        first_window.get('firstLocalIteration') == 0 and first_window.get('lastLocalIteration') == 19
        and type(last) is int and last >= 39
        and last_window.get('firstLocalIteration') == last - 19
        and last_window.get('lastLocalIteration') == last
        and history.get('lastFiniteLocalIteration') == last)
    for key, name in (('rms[Rho]', 'densityResidualDevelopment'), ('rms[RhoE]', 'energyResidualDevelopment')):
        first_value = first_window.get('meanLogResidual', {}).get(key) if first_window else None
        last_value = last_window.get('meanLogResidual', {}).get(key) if last_window else None
        delta = last_value - first_value if windows_ok and finite(first_value) and finite(last_value) else None
        gates[name] = gate('pending' if delta is None else 'passed' if delta <= -.25 else 'failed',
            dict(first20MeanLog=first_value, last20MeanLog=last_value, lastMinusFirst=delta),
            dict(lastMinusFirstMaximum=-.25), 'development.report.history.windows',
            None if delta is not None else '真实首末20窗不足或身份不完整，不能补造发展量')
    terminal_available = bool(normal and cp.get('normalFinalPairPresent') is True
                              and last_complete is not None and last_complete == last)
    terminal = next((item for item in d.get('frames', []) if
                     item.get('localIteration') == last), None) if terminal_available else None
    for name, key, limit, inclusive in (
        ('terminalWallTemperature', 'wallTemperatureMaximumDifferenceK', 19.9253, False),
        ('terminalIndependentYPlus', 'independentFirstNeighborYPlus99', 1., True)):
        value = terminal.get(key) if terminal else None
        known = finite(value) and value >= 0
        passed = known and (value <= limit if inclusive else value < limit)
        gates[name] = gate('pending' if not known else 'passed' if passed else 'failed',
            dict(localIteration=last if terminal else None, value=value),
            dict(operator='<=' if inclusive else '<', value=limit), 'development.report.frames',
            None if known else '缺正常终态完整场；保存首帧不代表最终壁流体状态')
    change = None
    if len(steps) >= 2:
        matching = [item for item in d.get('consecutiveSurfaceChanges', [])
                    if item.get('fromLocalIteration') == steps[-2] and item.get('toLocalIteration') == steps[-1]]
        if len(matching) == 1:
            change = matching[0]
    for name, field, limit in (('lastTwoPressureAreaL1', 'p', .05), ('lastTwoHeatAreaL1', 'q', .10)):
        item = change.get('fields', {}).get(field, {}) if change else {}
        value = item.get('relativeL1')
        known = finite(value) and value >= 0 and item.get('zeroReference') is False
        gates[name] = gate('pending' if not known else 'passed' if value <= limit else 'failed',
            dict(fromLocalIteration=steps[-2] if len(steps) >= 2 else None,
                 toLocalIteration=last_complete, relativeL1=value), dict(maximum=limit),
            'development.report.consecutiveSurfaceChanges',
            None if known else '缺最近两份完整场或非零参照积分；不能用同一首帧自比较')

    # 模块六：事件只看已独立重算的摘要；阶段未知或截断不能当作末20无零步。
    events = v.get('conservativeStepEvents', {})
    invalid_old = events.get('invalidOldEvents')
    no_trace_verified = verified and events.get('traceFilePresent') is False and events.get('eventsObserved') == 0
    event_trace_verified = (verified and events.get('actualStepEventsVerified') is True
                            and events.get('traceComplete') is True)
    old_safe = no_trace_verified or (event_trace_verified and invalid_old == 0)
    gates['noInvalidOldState'] = gate('passed' if old_safe else 'failed' if finite(invalid_old)
        and invalid_old > 0 else 'pending', dict(invalidOldEvents=0 if no_trace_verified else invalid_old),
        '独立核验事件中非法旧态为0', 'verification.report.conservativeStepEvents')
    event_path = paths['verification'].parent / 'step_event_summary.json'
    event_rows = read_small(event_path) if events.get('eventsObserved', 0) else []
    event_window = dict(firstLocalIteration=last_window.get('firstLocalIteration') if last_window else None,
                        lastLocalIteration=last, zeroStepEvents=None, unknownEventTimes=None,
                        eventsAfterLastFiniteUpdate=None)
    zero_status, zero_reason = 'pending', '缺已核验完整事件或末20窗口'
    if windows_ok and no_trace_verified:
        zero_status, zero_reason = 'passed', '已核验运行合同中无事件文件/记录，不表示限步机制被实际激活'
        event_window.update(zeroStepEvents=0, unknownEventTimes=0, eventsAfterLastFiniteUpdate=0)
    elif windows_ok and event_trace_verified and isinstance(event_rows, list) and len(event_rows) == events.get('eventsObserved'):
        unknown = [item for item in event_rows if not item.get('geometryContextKnown') or
                   type(item.get('innerIteration')) is not int]
        after = [item for item in event_rows if type(item.get('innerIteration')) is int
                 and item['innerIteration'] > last]
        zero = [item for item in event_rows if item.get('oldValid') is True and item.get('theta') == 0
                and type(item.get('innerIteration')) is int and last - 19 <= item['innerIteration'] <= last]
        event_window.update(zeroStepEvents=len(zero), unknownEventTimes=len(unknown),
                            eventsAfterLastFiniteUpdate=len(after), zeroStepEventIDs=[item.get('event') for item in zero])
        zero_status = 'failed' if zero else 'pending' if unknown or after else 'passed'
        zero_reason = ('末20发现零步' if zero else '存在未知阶段或末有限历史之后事件，不能宣称终态稳定'
                       if unknown or after else '完整事件按实际innerIteration核对，末20无零步')
    gates['noZeroStepInLast20'] = gate(zero_status, event_window, '真实末20更新无零步',
        'verification.report + step_event_summary.json', zero_reason)

    # 模块七：通过发展筛选不授予实验/全面积资格；失败和待验均阻止采纳与自动续算。
    locks[str(Path(__file__).resolve())] = sha(Path(__file__).resolve())
    stable = {path: Path(path).is_file() and sha(path) == digest for path, digest in locks.items()}
    gates['screeningInputsStable'] = gate('passed' if all(stable.values()) else 'failed', stable,
        '本次实际读取的小报告逐项SHA保持', 'screen.inputSHA256')
    failed = [name for name, item in gates.items() if item['status'] == 'failed']
    pending = [name for name, item in gates.items() if item['status'] == 'pending']
    adopted = not failed and not pending
    result = dict(screeningCompleted=True, developmentScreeningPassed=adopted, candidateAdopted=adopted,
        adoptionScope='仅本次事前发展筛选；默认配置不由本脚本修改',
        candidateDecision='passed-development-screening' if adopted else 'not-adopted',
        failedCriteria=failed, pendingCriteria=pending, criteria=gates,
        nativeExitCode=receipt['exitCode'], normalRunCompleted=normal,
        completeSavedLocalIterations=steps, lastCompleteLocalIteration=last_complete,
        lastCompleteSeriesSuccessfulUpdates=pairs[-1].get('seriesSuccessfulUpdates') if pairs else None,
        lastFiniteHistoryLocalIteration=last,
        lastFiniteHistorySeriesSuccessfulUpdates=cp.get('seriesUpdatesAtLastFiniteHistory'),
        terminalCompleteFieldAvailable=terminal_available, abnormalInstantStateAvailable=False,
        finalTerminationCause='未由本筛选推断；退出码、恢复事件均不单独证明最后终止原因',
        originalDensityEnergyMinus8Passed=history.get('originalDensityEnergyMinus8Passed'),
        temperatureIdentities=dict(givenTw='已给定空间壁温输入',
            wallFluidTemperature='仅实际完整保存场中的流体温度；最终场缺失则未知', solvedSolidWallTemperature=False),
        fieldMappingNotes=[
            'audit.frames.givenTwMaximumFluidDifferenceK与development.frames.wallTemperatureMaximumDifferenceK均指给定Tw和壁流体T之差；终态仅取准确末场。',
            'independentFirstNeighborYPlus99为独立首邻诊断，不使用层流原生Y_Plus=0作为通过证据。',
            'history meanLogResidual为已取log的残差均值；末减首<=-0.25，不再取log。',
            'p/q面积L1取最近两份不同完整场，以较早场绝对面积积分为分母；分母为0则待验。',
            'checkpoint阶段标签与热审核阶段标签不同，不能按字符串相等合并。',
            '事件innerIteration属于原生求解阶段；未知或末有限历史之后事件不补入已完成窗口。',
            'nativeOutputCapacityEnforcementDemonstrated=false表示未实测超容量停止，不等于事中路径未启用。'],
        physicalQualified=False, experimentalQualification=False, newQualifiedArea_m2=0,
        fullRocket20DegreeQualified=False, defaultsChanged=False, automaticContinuation=False,
        newCFDTasks=0, largeFieldReads=0, inputSHA256=locks, inputsStable=stable,
        elapsedSeconds=time.monotonic() - started)
    output.mkdir(parents=True, exist_ok=False)
    (output / 'executed_screen.py').write_bytes(Path(__file__).read_bytes())
    (output / 'report.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf8')
    print(json.dumps(dict(screeningCompleted=True, candidateAdopted=adopted,
        failedCriteria=failed, pendingCriteria=pending, elapsedSeconds=result['elapsedSeconds']), ensure_ascii=False))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('folder', 'audit', 'development', 'verification', 'output'):
        parser.add_argument('--' + argument, required=True)
    args = parser.parse_args()
    screen(args.folder, args.audit, args.development, args.verification, args.output)
