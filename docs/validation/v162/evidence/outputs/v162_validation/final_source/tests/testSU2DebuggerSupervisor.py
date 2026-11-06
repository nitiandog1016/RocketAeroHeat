"""验证调试监督器的时间/取消终态和自有子进程回收，零流场任务。"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import threading
import time

import psutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from runEmbeddedSU2 import atomic_write_json, sha
from runSU2NativeDebugger import observe


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--child', action='store_true')
    parser.add_argument('--folder', required=True)
    args = parser.parse_args()
    folder = Path(args.folder).resolve()
    folder.relative_to(ROOT/'outputs')
    if args.child:
        child = subprocess.Popen([sys.executable, '-c', 'import time;time.sleep(60)'],
            creationflags=subprocess.CREATE_NO_WINDOW)
        atomic_write_json(folder/'child.json', dict(pid=child.pid))
        child.wait()
        return
    folder.mkdir(parents=True, exist_ok=False)
    records = []
    for mode in ('timeout', 'cancel'):
        case = folder/mode
        case.mkdir()
        cancel = case/'cancel.json'
        def request_cancel():
            deadline = time.monotonic()+10
            while not (case/'child.json').exists() and time.monotonic() < deadline:
                time.sleep(.02)
            atomic_write_json(cancel, dict(requested=True))
        worker = None
        if mode == 'cancel':
            worker = threading.Thread(target=request_cancel)
            worker.start()
        command = [sys.executable, str(Path(__file__).resolve()), '--child', '--folder', str(case)]
        locks = {str(p):sha(p) for p in (Path(__file__).resolve(), ROOT/'scripts/runSU2NativeDebugger.py')}
        record = observe(command, case, locks, 2 if mode == 'timeout' else 10, 1024**2, cancel)
        if worker:
            worker.join(timeout=12)
            assert not worker.is_alive()
        child = json.loads((case/'child.json').read_text())['pid']
        assert child in record['observedProcessIDs'] and not psutil.pid_exists(child)
        assert record['stopReason'] == ('time-budget' if mode == 'timeout' else 'explicit-cancel')
        assert not record['survivingProcessIDs'] and record['nativeExitCode'] is None
        records.append(record)
    atomic_write_json(folder/'report.json', dict(passed=True, cases=records, nativeFlowRuns=0))
    print('Process-tree timeout and cancellation verified')


if __name__ == '__main__':
    main()
