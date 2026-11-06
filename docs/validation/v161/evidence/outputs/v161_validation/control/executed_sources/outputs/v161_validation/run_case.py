"""本轮实际求解父监督，限制RSS并保留取消来源。"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import psutil

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from runEmbeddedSU2 import atomic_write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('case')
    parser.add_argument('--seconds', type=int, default=60)
    parser.add_argument('--cancel', action='store_true')
    args = parser.parse_args()
    folder = (ROOT/'outputs/v161_validation'/args.case).resolve()
    folder.relative_to(ROOT/'outputs/v161_validation')
    receipt = folder/'parent_process.json'
    assert not receipt.exists() and not (folder/'engine_process.json').exists()
    source_files = ['scripts/prepareSU2IdealAirRestart.py', 'scripts/su2IdealAirInterface.py',
        'scripts/su2IdealAirCheckpoints.py', 'scripts/auditSU2IdealAirState.py',
        'scripts/runEmbeddedSU2.py', 'scripts/su2ThermalRuntime.py', 'docs/validation/v161/PLAN.md',
        'outputs/v161_validation/run_case.py']
    snapshot = folder/'executed_sources'
    for name in source_files:
        dest = snapshot/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes((ROOT/name).read_bytes())
    request = dict(workDirectory=str(folder), configFile='case.cfg', meshFile='mesh.su2',
        maximumSeconds=args.seconds, threads=2, sourceFiles=source_files,
        jobToken='v161-'+args.case)
    atomic_write_json(folder/'engine_request.json', request)
    started = time.monotonic()
    record = dict(status='starting', maximumRSSBytes=4*1024**3, maximumSeconds=args.seconds+120)
    atomic_write_json(receipt, record)
    with (folder/'parent.log').open('wb') as stream:
        proc = subprocess.Popen([sys.executable, str(ROOT/'scripts/runEmbeddedSU2.py'),
            '--request', str(folder/'engine_request.json')], cwd=str(ROOT), stdout=stream, stderr=stream,
            creationflags=subprocess.CREATE_NO_WINDOW)
        record.update(status='running', processID=proc.pid)
        atomic_write_json(receipt, record)
        peak, cancelled, reason = 0, False, None
        while proc.poll() is None:
            try:
                parent = psutil.Process(proc.pid)
                rss = 0
                for child in [parent]+parent.children(recursive=True):
                    try:
                        rss += child.memory_info().rss
                    except psutil.NoSuchProcess:
                        pass
                peak = max(peak, rss)
            except psutil.NoSuchProcess:
                pass
            if peak > record['maximumRSSBytes']:
                reason = 'memory'
            if time.monotonic()-started > record['maximumSeconds']:
                reason = 'supervisor-time'
            if reason:
                subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'], stdout=stream, stderr=stream)
                break
            if args.cancel and not cancelled and (folder/'surface_flow_000000.vtk').exists():
                try:
                    from auditSU2IdealAirState import load_state
                    from auditSU2SSTState import vtk_tokens
                    data, fields, _ = load_state(folder/'restart_flow_000000.csv', True, True, True)
                    xyz, surface = vtk_tokens(folder/'surface_flow_000000.vtk')
                    native = json.loads((folder/'engine_process.json').read_text('utf8'))
                    assert len(xyz) and len(data) and psutil.pid_exists(native['processID'])
                    atomic_write_json(folder/'engine_cancel.json', dict(jobToken=request['jobToken']))
                    cancelled = True
                    record['cancelRequestedAfterCompleteLocal0'] = True
                except (OSError, ValueError, AssertionError, KeyError, StopIteration):
                    pass
            time.sleep(.025)
        code = proc.wait()
    record.update(status='terminal', exitCode=code, elapsedSeconds=time.monotonic()-started,
        sampledPeakRSSBytes=peak, stoppedByBudget=reason, cancellationRequested=cancelled)
    atomic_write_json(receipt, record)
    print(json.dumps(record))
    raise SystemExit(0 if reason is None and code==(125 if args.cancel else 0) else 1)


if __name__ == '__main__':
    main()
