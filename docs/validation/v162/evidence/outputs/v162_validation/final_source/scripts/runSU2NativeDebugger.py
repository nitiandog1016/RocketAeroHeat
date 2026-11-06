"""固定内嵌引擎的有限异常观察，分别记录调试器状态和原生保存场。"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time

import psutil

from runEmbeddedSU2 import atomic_write_json, sha
from su2ThermalRuntime import ROOT, resolve, read

GDB = ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/gdb.exe'


def observe(command, folder, locks, seconds, maximum_bytes, cancel_file=None):
    """模块一。只管理当前创建及观察到的子进程，终态后回收剩余自有进程。"""
    receipt = folder/'debugger_process.json'
    with receipt.open('x', encoding='utf8') as stream:
        json.dump(dict(status='starting', command=command), stream)
    started = time.monotonic()
    record = dict(status='starting', command=command, inputSHA256=locks,
        maximumSeconds=seconds, maximumRSSBytes=4*1024**3, maximumOutputBytes=maximum_bytes,
        supervisorProcessID=os.getpid(), physicsQualified=False, nativeExitCode=None)
    family = {}
    peak_rss, peak_bytes, reason, error = 0, 0, None, None
    env = dict(os.environ, OMP_NUM_THREADS='2')
    for key in ('RAH_WALL_TEMPERATURE_FILE', 'RAH_FROZEN_ENERGY_MODE', 'RAH_SST_UPDATE_TRACE'):
        env.pop(key, None)
    with (folder/'debugger.log').open('wb') as log:
        process = subprocess.Popen(command, cwd=str(folder), env=env, stdout=log, stderr=log,
            creationflags=subprocess.CREATE_NO_WINDOW)
        root = psutil.Process(process.pid)
        family[root.pid] = root
        record.update(status='running', processID=process.pid)
        try:
            atomic_write_json(receipt, record)
            while process.poll() is None:
                try:
                    for child in root.children(recursive=True):
                        family[child.pid] = child
                except psutil.NoSuchProcess:
                    pass
                rss = 0
                for child in family.values():
                    try:
                        if child.is_running():
                            rss += child.memory_info().rss
                    except psutil.NoSuchProcess:
                        pass
                peak_rss = max(peak_rss, rss)
                length = sum(p.stat().st_size for p in folder.iterdir()
                    if p.is_file() and str(p) not in locks and p.name != receipt.name)
                peak_bytes = max(peak_bytes, length)
                if cancel_file and cancel_file.exists():
                    reason = 'explicit-cancel'
                elif time.monotonic()-started >= seconds:
                    reason = 'time-budget'
                elif rss > record['maximumRSSBytes']:
                    reason = 'memory-budget'
                elif length > maximum_bytes:
                    reason = 'output-budget'
                if reason:
                    break
                time.sleep(.05)
        except Exception as exception:
            error = type(exception).__name__+': '+str(exception)
        finally:
            # 先回收仍存活的根进程树，之后按带创建时间的Process对象核对已观察子进程。
            if process.poll() is None:
                subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'],
                    stdout=log, stderr=log, creationflags=subprocess.CREATE_NO_WINDOW)
            for child in reversed(list(family.values())):
                try:
                    if child.is_running():
                        child.kill()
                except psutil.NoSuchProcess:
                    pass
            gone, alive = psutil.wait_procs(list(family.values()), timeout=5)
            code = process.wait(timeout=5)
    record.update(status='terminal', debuggerExitCode=code, elapsedSeconds=time.monotonic()-started,
        stopReason=reason, supervisorError=error, peakFamilyRSSBytes=peak_rss,
        peakOutputBytes=peak_bytes, observedProcessIDs=list(family), survivingProcessIDs=[p.pid for p in alive],
        inputsStable={name:sha(Path(name)) == digest for name, digest in locks.items()})
    atomic_write_json(receipt, record)
    assert not alive and error is None and all(record['inputsStable'].values())
    return record


def run(folder, plan):
    """模块二。固定真实输入和观察命令，独立保存异常终态身份。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    assert not (folder/'engine_process.json').exists()
    plan = Path(plan).resolve()
    exe, runtime, flavor, _, locks = resolve(folder, (folder/'case.cfg').read_text('utf8'))
    assert flavor == 'rah-ideal-air-phase-v2'
    manifest = read(folder/'case_manifest.json')
    assert manifest['initialization']['sourceLocalIteration'] == 149
    assert manifest['flight']['iterations'] == 50
    assert manifest['idealAir']['checkpointContract']['maximumOutputBytes'] == 1024**3
    commands = folder/'observe.gdb'
    assert not commands.exists()
    commands.write_text('''set pagination off
set confirm off
set print demangle on
break __cxa_call_terminate
break _ZN9__gnu_cxx24__trap_terminate_handlerEv
run
echo RAH_DEBUG_STOP\\n
thread apply all bt 30
info registers rcx rdx r8 r9 rip
echo RAH_EXCEPTION_STORAGE\\n
x/16gx $rcx-96
x/s *(char**)($rcx+40)
echo RAH_OBSERVATION_COMPLETE\\n
quit
''', encoding='ascii')
    for path in (exe, GDB, commands, plan, Path(__file__).resolve(), folder/'case.cfg',
                 folder/'mesh.su2', folder/'initial_flow.csv', ROOT/'scripts/runEmbeddedSU2.py'):
        locks[str(path)] = sha(path)
    snapshot = folder/'executed_sources'
    # 输入表和二进制按SHA绑定，准确执行源码单独归档。
    for path in (Path(__file__).resolve(), ROOT/'scripts/prepareSU2IdealAirRestart.py',
                 ROOT/'scripts/su2ThermalRuntime.py', ROOT/'scripts/su2IdealAirInterface.py',
                 ROOT/'scripts/su2IdealAirCheckpoints.py', plan):
        target = snapshot/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
        locks[str(path)] = sha(path)
    command = [str(GDB), '--batch', '-nx', '-q', '-x', str(commands), '--args', str(exe), 'case.cfg', '--threads', '2']
    record = observe(command, folder, locks, 300, 1024**3, folder/'debug_cancel.json')
    text = (folder/'debugger.log').read_text('utf8', errors='replace')
    record.update(runtimeFlavor=flavor, runtimeSHA256=sha(exe), runtimeVersion=runtime['version'],
        terminationStopCaptured=bool(re.search(r'hit Breakpoint \d+[^\n]*(?:__cxa_call_terminate|__trap_terminate_handler)', text)),
        observationComplete='RAH_OBSERVATION_COMPLETE' in text,
        thermodynamicModel=manifest['idealAir'], normalNativeCompletion=False)
    record['outputSHA256'] = {p.name:sha(p) for p in folder.iterdir()
        if p.is_file() and str(p) not in locks and p.name != 'debugger_process.json'}
    from su2IdealAirCheckpoints import inventory
    record['checkpoints'] = inventory(folder, dict(record, exitCode=None))
    atomic_write_json(folder/'debugger_process.json', record)
    print(json.dumps({key:record[key] for key in ('debuggerExitCode', 'elapsedSeconds', 'terminationStopCaptured',
        'observationComplete', 'stopReason', 'peakFamilyRSSBytes')}))
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--plan', required=True)
    args = parser.parse_args()
    result = run(args.folder, args.plan)
    success = result['terminationStopCaptured'] and result['observationComplete'] and result['stopReason'] is None
    raise SystemExit(0 if success else 1)
