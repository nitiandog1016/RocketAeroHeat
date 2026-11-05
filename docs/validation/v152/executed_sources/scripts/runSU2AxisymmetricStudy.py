"""从已审核公开算例检查点进行一次受控收敛计算，保存独立局部步号。"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
import psutil

from prepareSU2AxisymmetricBenchmark import ROOT, sha, save_json
from prepareSU2CheckpointContinuation import config
from su2SSTInterface import state


def prepare(source, folder, plan):
    """核对源终态及相同模型，逐字复用网格和守恒场。"""
    source = Path(source).resolve()
    folder = Path(folder).resolve()
    plan = Path(plan).resolve()
    source.relative_to(ROOT/'outputs')
    folder.relative_to(ROOT/'outputs')
    assert not folder.exists() and plan.is_file()
    receipt = json.loads((source/'engine_process.json').read_text(encoding='utf-8'))
    benchmark = json.loads((source/'benchmark_manifest.json').read_text(encoding='utf-8'))
    audit = json.loads((source/'audit02.json').read_text(encoding='utf-8'))
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0
    assert not receipt['timedOut'] and all(receipt['inputsStable'].values())
    assert benchmark['kind'] == 'NASA-ACC-SBLI-M7-axisymmetric-startup-v1'
    assert audit['thermalReadbackPassed'] and audit['checkpointLocalIteration'] == 1
    for path, h in receipt['sourceSHA256'].items():
        assert sha(path) == h, path
    for name, h in receipt['outputSHA256'].items():
        assert sha(source/name) == h, name
    a, keys, nd, rho, velocity, thermal = state(source/'restart_flow.csv', True, True)
    assert nd == 2 and len(a) == 16261 and np.min(thermal) > 0
    old = config(source/'case.cfg')
    changes = dict(RESTART_SOL='YES', SOLUTION_FILENAME='initial_flow', READ_BINARY_RESTART='NO',
                   ITER='12000', OUTPUT_WRT_FREQ='1000, 3000', WRT_RESTART_OVERWRITE='NO',
                   WRT_SURFACE_OVERWRITE='NO', SCREEN_WRT_FREQ_INNER='500',
                   VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL, TIMESTEP)')
    current = dict(old)
    current.update(changes)
    assert {k: v for k, v in current.items() if k not in changes} == {
        k: v for k, v in old.items() if k not in changes}
    folder.mkdir(parents=True)
    os.link(str(source/'mesh.su2'), str(folder/'mesh.su2'))
    os.link(str(source/'restart_flow.csv'), str(folder/'initial_flow.csv'))
    (folder/'case.cfg').write_bytes(('\n'.join(k+'= '+v for k, v in current.items())+'\n').encode('ascii'))
    locked = [source/n for n in ('engine_process.json', 'benchmark_manifest.json', 'audit02.json',
                                  'mesh.su2', 'restart_flow.csv', 'case.cfg')]
    manifest = dict(kind='NASA-ACC-SBLI-M7-axisymmetric-convergence-v1',
                    benchmark=benchmark, sourceDirectory=str(source), sourceLocalIteration=1,
                    sourceOriginalStartupPassed=audit['interfacePassed'],
                    sourceSavedWallFluidTemperatureRangeK=audit['savedWallFluidTemperatureRangeK'],
                    controlsChanged=changes, samePhysicsAndNumerics=True,
                    maximumSeconds=720, maximumRSSBytes=4*1024**3, maximumOutputBytes=120*1024**2,
                    sourceSHA256={str(p): sha(p) for p in locked},
                    plan=str(plan), planSHA256=sha(plan),
                    physicsQualified=False, newQualifiedArea=0)
    save_json(folder/'study_manifest.json', manifest)
    return folder, manifest


def run(folder, manifest):
    """监督原生进程及其子进程，只结束本任务，任何终态均保留回执。"""
    target = folder/'engine_process.json'
    assert not target.exists()
    engine = ROOT/'tools/runtime/su2/8.6.0/bin/SU2_CFD.exe'
    files = [Path(__file__).resolve(), engine, folder/'case.cfg', folder/'mesh.su2',
             folder/'initial_flow.csv', folder/'study_manifest.json', Path(manifest['plan']),
             ROOT/'scripts/su2SSTInterface.py', ROOT/'scripts/prepareSU2AxisymmetricBenchmark.py',
             ROOT/'scripts/prepareSU2CheckpointContinuation.py']
    locked = {str(p): sha(p) for p in files}
    command = [str(engine), 'case.cfg', '--threads', '2']
    started = time.monotonic()
    peak_rss, peak_bytes = 0, 0
    reason = None
    checked_states = {}
    with (folder/'engine.log').open('wb') as log:
        process = subprocess.Popen(command, cwd=str(folder), stdout=log, stderr=subprocess.STDOUT,
                                   env=dict(os.environ, OMP_NUM_THREADS='2'),
                                   creationflags=subprocess.CREATE_NO_WINDOW)
        save_json(target, dict(status='running', processID=process.pid, command=command,
                               sourceSHA256=locked, maximumSeconds=manifest['maximumSeconds']))
        while process.poll() is None:
            try:
                family = [psutil.Process(process.pid)] + psutil.Process(process.pid).children(recursive=True)
                rss = sum(p.memory_info().rss for p in family if p.is_running())
                peak_rss = max(peak_rss, rss)
            except psutil.NoSuchProcess:
                pass
            size = sum(p.stat().st_size for p in folder.iterdir() if p.is_file())
            peak_bytes = max(peak_bytes, size)
            if time.monotonic()-started > manifest['maximumSeconds']:
                reason = 'wall-clock budget'
            if peak_rss > manifest['maximumRSSBytes']:
                reason = 'family memory budget'
            if size > manifest['maximumOutputBytes']:
                reason = 'output file length budget'
            if manifest.get('validateSavedState'):
                # 下一百步表面已落盘时，较早体场的writer已经结束，避免读取半写文件。
                surfaces = [int(p.stem.split('_')[-1]) for p in folder.glob('surface_flow_*.vtk')]
                latest_surface = max(surfaces) if surfaces else -1
                for saved in folder.glob('restart_flow_*.csv'):
                    step = int(saved.stem.split('_')[-1])
                    if step >= latest_surface or saved.name in checked_states:
                        continue
                    try:
                        state(saved, True, True, True)
                        checked_states[saved.name] = True
                    except (AssertionError, ValueError) as error:
                        checked_states[saved.name] = str(error)
                        reason = 'invalid saved conservative state'
            if reason and process.poll() is None:
                subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'],
                               stdout=log, stderr=log, creationflags=subprocess.CREATE_NO_WINDOW)
                break
            time.sleep(.25)
        code = process.wait()
    if manifest.get('validateSavedState'):
        # 进程终态后检查末帧和所有尚未检查的场，原生退出码与审核原因分别保留。
        for saved in folder.glob('restart_flow_*.csv'):
            if saved.name in checked_states:
                continue
            try:
                state(saved, True, True, True)
                checked_states[saved.name] = True
            except (AssertionError, ValueError) as error:
                checked_states[saved.name] = str(error)
                reason = reason or 'invalid saved conservative state'
    outputs = {p.name: sha(p) for p in folder.iterdir()
               if p.name in ('engine.log', 'history.csv') or p.name.startswith(('restart_flow', 'surface_flow'))}
    result = dict(status='terminal', exitCode=code, stopReason=reason,
                  elapsedSeconds=time.monotonic()-started, peakFamilyRSSBytes=peak_rss,
                  peakOutputFileBytes=peak_bytes, command=command, sourceSHA256=locked,
                  inputsStable={p: sha(p) == h for p, h in locked.items()}, outputSHA256=outputs,
                  savedStateChecks=checked_states)
    save_json(target, result)
    print(json.dumps({k: result[k] for k in ('status', 'exitCode', 'stopReason', 'elapsedSeconds',
                                            'peakFamilyRSSBytes', 'peakOutputFileBytes')}), flush=True)
    return 0 if code == 0 and reason is None and all(result['inputsStable'].values()) else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', required=True)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--plan', required=True)
    args = parser.parse_args()
    folder, manifest = prepare(args.source, args.folder, args.plan)
    raise SystemExit(run(folder, manifest))
