"""只读选择 NACA 已保存完整帧；异常任务须显式选择，永不补造终止瞬间场。"""
from pathlib import Path

from su2NACASpatialInterface import inventory
from su2ThermalRuntime import sha


def select_saved_frames(folder, receipt, saved_frames_only=False):
    """模块一：核验终态与来源；模块二：重清点；模块三：返回准确研究身份。

    返回 (checkpoints, identity)，不修改传入回执或任何文件。
    saved_frames_only 只增加异常终态的离线研究入口；正常任务仍须完整。
    此处只选择可审核的帧，热量、物性和实验资格仍由后续审核分别判断。
    """
    folder = Path(folder).resolve()
    assert isinstance(saved_frames_only, bool), '保存帧开关必须为布尔量'
    assert isinstance(receipt, dict) and receipt.get('status') == 'terminal', '任务尚无真实终态回执'
    code = receipt.get('exitCode')
    assert isinstance(code, int) and not isinstance(code, bool), '终态缺少明确整数退出码'
    assert not receipt.get('supervisorError'), '监督完整性错误不能作为保存帧研究入口'
    stable = receipt.get('inputsStable')
    assert isinstance(stable, dict) and stable and all(value is True for value in stable.values()), '输入稳定记录缺失或失败'
    input_hashes = receipt.get('inputSHA256')
    assert isinstance(input_hashes, dict) and input_hashes and set(stable) == set(input_hashes), '输入稳定记录与锁定来源不一致'
    output_hashes = receipt.get('outputSHA256')
    assert isinstance(output_hashes, dict) and output_hashes, '保存输出哈希缺失'
    assert sha(folder / 'history.csv') == output_hashes.get('history.csv'), '有限历史与终态回执哈希不符'

    # 允许离线夹具更换目录；物理网格和清点合同必须仍与本次锁定输入一致。
    for name in ('mesh.su2', 'case_manifest.json'):
        target = folder / name
        exact = [value for path, value in input_hashes.items() if Path(path).resolve() == target]
        # 实际任务可同时锁定本目录和上游网格；必须先按完整路径确认本任务身份。
        # 仅离线搬迁夹具没有精确路径时允许同名回退，且所有候选哈希必须一致。
        expected = exact or [value for path, value in input_hashes.items() if Path(path).name == name]
        assert expected and len(set(expected)) == 1 and sha(target) == expected[0], '网格或清点合同偏离运行来源: ' + name
    flags = {}
    for key in ('timedOut', 'cancelled', 'outputBudgetExceeded'):
        value = receipt.get(key, False)
        assert isinstance(value, bool), '终态标记不是布尔量: ' + key
        flags[key] = value
    normal = code == 0 and not any(flags.values())
    assert normal or saved_frames_only, '默认审核拒绝异常任务；仅显式保存帧研究可读取'

    # 不信任旧 checkpoints.complete。重新检查完整历史、双文件哈希和场身份。
    checkpoints = inventory(folder, receipt)
    assert checkpoints['normalTerminal'] == normal, '终态和重新清点身份冲突'
    assert checkpoints['historyContiguous'], '历史不连续；不能确认保存帧发展阶段'
    assert checkpoints['complete'], '没有可独立重验的完整保存帧'
    if normal:
        assert checkpoints['normalScheduleComplete'], '正常任务缺失或拒绝计划保存帧'

    last = checkpoints['lastComplete']
    reasons = [name for name, value in flags.items() if value]
    if code != 0:
        reasons.append('nonzero-exit-code')
    identity = dict(
        auditScope='normal-complete-saved-frames' if normal else 'complete-saved-frames-from-abnormal-task',
        savedFramesOnlyRequested=saved_frames_only,
        normalRunCompleted=normal,
        taskCompletedNormally=normal,
        actualTerminalStatus=receipt['status'],
        actualExitCode=code,
        terminalFlags=flags,
        terminalReasons=reasons,
        savedFrameSelectionPassed=True,
        thermalInterfaceAudited=False,
        lastCompleteLocalIteration=last['localIteration'],
        lastFiniteHistoryLocalIteration=checkpoints['lastFiniteHistoryLocalIteration'],
        missingScheduledLocalIterations=list(checkpoints['missingScheduledLocalIterations']),
        rejectedSavedPairs=checkpoints['rejected'],
        normalScheduleComplete=checkpoints['normalScheduleComplete'],
        abnormalInstantStateAvailable=False,
        unsavedTerminalStateKnown=False,
        physicalSolutionQualified=False,
        gridAndDomainQualified=False,
        experimentalQualification=False,
        newQualifiedArea_m2=0,
    )
    return checkpoints, identity
