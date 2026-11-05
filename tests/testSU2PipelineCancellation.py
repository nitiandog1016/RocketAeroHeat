"""独立验证后台准备/审核取消时回收本次子孙进程，不启动物性或SU2。"""
from pathlib import Path
import json
import os
import subprocess
import sys
import threading
import time

import psutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from runSU2WorkbenchPipeline import Pipeline, PipelineStopped
from runEmbeddedSU2 import atomic_write_json


def main():
    folder = ROOT/'outputs/v157_validation/process_tree'
    folder.mkdir(parents=True, exist_ok=False)
    request = dict(workDirectory=str(folder), jobToken='v157-owned-process-tree',
                   idealAirInitialization=dict(model='ideal-air-SST', mode='wallDistance'))
    atomic_write_json(folder/'engine_request.json', request)
    pipeline = Pipeline(folder/'engine_request.json')
    child = folder/'child.py'
    child.write_text(
        'import subprocess,sys,json,os\n'
        'from pathlib import Path\n'
        'p=subprocess.Popen([sys.executable,"-c","import time; time.sleep(20)"],creationflags=subprocess.CREATE_NO_WINDOW)\n'
        'Path("owned_pids.json").write_text(json.dumps([os.getpid(),p.pid]))\n'
        'p.wait()\n', encoding='utf8')
    seen = []

    def cancel_when_started():
        started = time.monotonic()
        while time.monotonic()-started < 4:
            path = folder/'owned_pids.json'
            if path.exists():
                try:
                    pids = json.loads(path.read_text())
                except ValueError:
                    time.sleep(.01)
                    continue
                seen.extend((pid, psutil.Process(pid).create_time()) for pid in pids)
                atomic_write_json(folder/'engine_cancel.json', dict(jobToken=request['jobToken']))
                return
            time.sleep(.01)

    watcher = threading.Thread(target=cancel_when_started)
    watcher.start()
    stopped = False
    try:
        pipeline.phase('preparation', 5, lambda: pipeline.run_property(pipeline.out, 'fixture',
            [sys.executable, str(child)], folder, dict(os.environ), 5))
    except PipelineStopped:
        stopped = True
    watcher.join()
    assert stopped and len(seen) == 2
    for pid, created in seen:
        assert not psutil.pid_exists(pid) or psutil.Process(pid).create_time() != created
    receipt = json.loads((pipeline.out/'fixture_process.json').read_text())
    assert receipt['status'] == 'terminal' and receipt['cancelled'] and receipt['exitCode'] != 0
    atomic_write_json(folder/'report.json', dict(passed=True, ownedProcesses=2,
        allOwnedProcessesTerminated=True, childExitCode=receipt['exitCode'],
        elapsedSeconds=receipt['elapsedSeconds'], newSU2Runs=0))
    print('Owned child and grandchild terminated; native solver runs 0')


if __name__ == '__main__':
    main()
