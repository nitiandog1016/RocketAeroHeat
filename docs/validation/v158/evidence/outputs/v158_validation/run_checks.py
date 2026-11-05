"""本轮有限准备监督，记录准确执行源、时间、内存及终态。"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import psutil

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'outputs/v158_validation'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(name, command, seconds, memory):
    receipt = OUT/(name+'_process.json')
    assert not receipt.exists()
    source_paths = ['scripts/prepareSU2PublicValidation.py', 'scripts/idealAirFreestream.py',
        'scripts/idealAirProperties.py', 'scripts/prepareSU2AxisymmetricBenchmark.py',
        'scripts/su2AxisymmetricAnalysis.py', 'scripts/su2SSTInterface.py', 'scripts/su2IdealAirInterface.py',
        'tests/testSU2PublicPreparation.py', 'tests/testSU2PublicFreestream.m',
        'tests/data/tmr_mach7/verify_contract.py', 'scripts/validateIdealAir.py',
        'idealAirThermodynamicState.m', 'third_party/ideal_air/air_ideal.json']
    source_dir = OUT/(name+'_source')
    source_dir.mkdir()
    locks = {}
    for name_in_root in source_paths:
        p = ROOT/name_in_root
        locks[name_in_root] = digest(p)
        destination = source_dir/name_in_root
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(p.read_bytes())
    state = dict(command=command, maximumSeconds=seconds, maximumMemoryBytes=memory,
                 maximumOutputBytes=64*1024**2, sourceSHA256=locks, status='starting')

    def save():
        receipt.write_text(json.dumps(state, indent=2)+'\n', encoding='utf8')

    save()
    began = time.monotonic()
    with (OUT/(name+'.log')).open('wb') as log:
        proc = subprocess.Popen(command, cwd=str(ROOT), stdout=log, stderr=log,
            creationflags=subprocess.CREATE_NO_WINDOW,
            env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
        state.update(status='running', processID=proc.pid)
        save()
        peak = 0
        reason = None
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
            size = sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
            if time.monotonic()-began > seconds:
                reason = 'time'
            if peak > memory:
                reason = 'memory'
            if size > state['maximumOutputBytes']:
                reason = 'storage'
            if reason:
                subprocess.run(['taskkill', '/PID', str(proc.pid), '/T', '/F'], stdout=log, stderr=log)
                break
            time.sleep(.1)
        code = proc.wait()
    state.update(status='terminal', exitCode=code, elapsedSeconds=time.monotonic()-began,
                 sampledPeakMemoryBytes=peak, stoppedByBudget=reason,
                 sourceStable=all(digest(ROOT/p)==expected for p, expected in locks.items()))
    save()
    print(json.dumps({key: state[key] for key in ('status', 'exitCode', 'elapsedSeconds', 'stoppedByBudget', 'sourceStable')}), flush=True)
    if code or reason or not state['sourceStable']:
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('step', choices=['prepare', 'review', 'matlab'])
    parser.add_argument('--name')
    args = parser.parse_args()
    commands = {
        'prepare': [sys.executable, 'scripts/prepareSU2PublicValidation.py', '--out', 'outputs/v158_validation/preparation01'],
        'review': [sys.executable, 'tests/testSU2PublicPreparation.py', '--folder', 'outputs/v158_validation/preparation01'],
        'matlab': [sys.executable, 'scripts/runMatlabResearch.py', 'freestream01',
                   "addpath('tests');testSU2PublicFreestream('outputs/v158_validation/preparation01');", 
                   '--input', 'outputs/v158_validation/preparation01/preparation.json',
                   '--folder', 'outputs/v158_validation/matlab']}
    run(args.name or args.step+'01', commands[args.step], 180 if args.step=='matlab' else 60,
        (4 if args.step=='matlab' else 2)*1024**3)
