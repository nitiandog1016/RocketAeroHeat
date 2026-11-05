"""监督器专用进程回归。模拟输出写入中断，核对取消、超时、失败和容量终止。"""
from pathlib import Path
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from validateIdealAir import run, save, sha
from su2ThermalRuntime import read
from su2IdealAirCheckpoints import inventory

OUT = ROOT/'outputs/v156_validation/terminal01'
SOURCE = ROOT/'outputs/v156_validation/cases01/sparse'


def fixture_supervisor(request):
    """仅在专用测试进程替换引擎解析器，其余监督和终态路径执行生产代码。"""
    import su2ThermalRuntime
    import runEmbeddedSU2
    def fixture_runtime(folder, config):
        locked = {str(folder/name): sha(folder/name) for name in ('case_manifest.json', 'ideal_air_contract.json')}
        return Path(sys.executable), dict(version='TEST-FIXTURE-NOT-A-FLOW-SOLVER', archiveSHA256='test-only'), \
            'rah-ideal-air-phase-v2', None, locked
    su2ThermalRuntime.resolve = fixture_runtime
    sys.argv = [str(ROOT/'scripts/runEmbeddedSU2.py'), '--request', request]
    return runEmbeddedSU2.main()


def run_case(name):
    """在真实子进程中写一对完整样本及一份截断样本，随后触发指定终态。"""
    folder = OUT/name
    folder.mkdir()
    for filename in ('mesh.su2', 'case_manifest.json', 'ideal_air_contract.json'):
        shutil.copyfile(str(SOURCE/filename), str(folder/filename))
    if name == 'capacity':
        contract = read(folder/'ideal_air_contract.json')
        contract['checkpointContract']['maximumOutputBytes'] = 1024
        save(folder/'ideal_air_contract.json', contract)
    program = '''from pathlib import Path
import json, shutil, time
source = Path(%r)
folder = Path.cwd()
for filename in ('restart_flow_000000.csv', 'surface_flow_000000.vtk'):
    shutil.copyfile(str(source/filename), str(folder/filename))
history = (source/'history.csv').read_text().splitlines(True)
(folder/'history.csv').write_text(''.join(history[:4]))
(folder/'restart_flow_000002.csv').write_bytes((source/'restart_flow_000002.csv').read_bytes()[:90])
if %r == 'cancel':
    (folder/'engine_cancel.json').write_text(json.dumps(dict(jobToken='fixture-token')))
if %r == 'failure':
    raise SystemExit(7)
time.sleep(8)
''' % (str(SOURCE), name, name)
    (folder/'case.cfg').write_text(program)
    save(folder/'engine_request.json', dict(workDirectory=str(folder), configFile='case.cfg', meshFile='mesh.su2',
        threads=2, maximumSeconds=1 if name == 'timeout' else 5, jobToken='fixture-token',
        sourceFiles=['tests/testSU2IdealAirCheckpointTerminals.py', 'scripts/su2IdealAirCheckpoints.py']))
    try:
        run(OUT, name, [sys.executable, '-X', 'utf8', str(Path(__file__).resolve()),
                       '--fixture-supervisor', str(folder/'engine_request.json')], 10)
    except AssertionError:
        pass
    receipt = read(folder/'engine_process.json')
    process = read(OUT/(name+'_process.json'))
    expected = dict(cancel=125, timeout=124, failure=7, capacity=126)[name]
    assert process['status'] == 'terminal' and process['exitCode'] == expected and process['stoppedByBudget'] is None
    assert receipt['status'] == 'terminal' and receipt['executableVersion'] == 'TEST-FIXTURE-NOT-A-FLOW-SOLVER'
    assert all(receipt['inputsStable'].values())
    saved = receipt['checkpoints']
    assert saved['lastComplete']['localIteration'] == 0 and saved['lastFiniteHistoryLocalIteration'] == 2
    assert [row['localIteration'] for row in saved['rejected']] == [2]
    assert saved['missingScheduledLocalIterations'] == [2]
    assert not saved['normalTerminal'] and not saved['abnormalInstantStateAvailable']
    assert receipt['cancelled'] == (name == 'cancel') and receipt['timedOut'] == (name == 'timeout')
    assert receipt['outputBudgetExceeded'] == (name == 'capacity')
    return dict(name=name, supervisorExitCode=process['exitCode'], nativeFixtureExitCode=receipt['exitCode'],
        elapsedSeconds=receipt['elapsedSeconds'], lastComplete=0, lastFiniteHistory=2,
        rejectedPartial=2, flowSolverRun=False, passed=True)


def main():
    OUT.mkdir(exist_ok=False)
    sources = [Path(__file__).resolve(), ROOT/'scripts/runEmbeddedSU2.py', ROOT/'scripts/su2IdealAirCheckpoints.py']
    locks = {str(path): sha(path) for path in sources}
    for path in sources:
        target = OUT/'executed_sources'/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
    results = [run_case(name) for name in ('cancel', 'timeout', 'failure', 'capacity')]
    # 对已经完成的真实稀疏任务只读复核完整保存计划，不增加原生求解。
    native = inventory(SOURCE, read(SOURCE/'engine_process.json'))
    assert native['normalScheduleComplete'] and native['expectedLocalIterations'] == [0, 2, 4, 5]
    stable = {name: sha(Path(name)) == digest for name, digest in locks.items()}
    assert all(stable.values())
    save(OUT/'report.json', dict(passed=True, results=results, sourceSHA256=locks, inputsStable=stable,
        nativeScheduleReadOnly=native, additionalNativeFlowRuns=0, physicsQualified=False))
    print(json.dumps(dict(passed=True, terminalChecks=4, nativeFlowRuns=0)))


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--fixture-supervisor':
        raise SystemExit(fixture_supervisor(sys.argv[2]))
    main()
