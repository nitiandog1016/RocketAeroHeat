"""在明确预算内运行公开轴对称基准，保存终态和输入/输出身份。"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from prepareSU2AxisymmetricBenchmark import ROOT, sha, save_json


def run(folder):
    """专用两更新入口独立记录 SST-Vm，拒绝覆盖已有原生任务。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    receipt = folder/'engine_process.json'
    assert not receipt.exists()
    manifest = json.loads((folder/'benchmark_manifest.json').read_text(encoding='utf-8'))
    assert manifest['kind'] == 'NASA-ACC-SBLI-M7-axisymmetric-startup-v1'
    for name, expected in manifest['inputSHA256'].items():
        assert sha(folder/name) == expected
    config = (folder/'case.cfg').read_text()
    assert '\nITER= 2\n' in config and '\nAXISYMMETRIC= YES\n' in config
    assert '\nSST_OPTIONS= (V1994m, VORTICITY)\n' in config
    engine = ROOT/'tools/runtime/su2/8.6.0/bin/SU2_CFD.exe'
    locked = {str(p): sha(p) for p in (Path(__file__), engine, folder/'case.cfg',
                                     folder/'mesh.su2', folder/'benchmark_manifest.json')}
    command = [str(engine), 'case.cfg', '--threads', '2']
    start = time.monotonic()
    timed_out = False
    with (folder/'engine.log').open('wb') as log:
        process = subprocess.Popen(command, cwd=str(folder), stdout=log, stderr=subprocess.STDOUT,
                                   env=dict(os.environ, OMP_NUM_THREADS='2'),
                                   creationflags=subprocess.CREATE_NO_WINDOW)
        save_json(receipt, dict(status='running', processID=process.pid, command=command,
                                sourceSHA256=locked))
        try:
            code = process.wait(timeout=60)
        except subprocess.TimeoutExpired:
            timed_out = True
            subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'],
                           stdout=log, stderr=log, creationflags=subprocess.CREATE_NO_WINDOW)
            code = process.wait()
    result = dict(status='terminal', exitCode=code, timedOut=timed_out,
                  elapsedSeconds=time.monotonic()-start, command=command, sourceSHA256=locked,
                  inputsStable={p: sha(p) == h for p, h in locked.items()},
                  outputSHA256={p.name: sha(p) for p in folder.iterdir()
                                if p.name in ('engine.log', 'history.csv', 'restart_flow.csv', 'surface_flow.vtk')})
    save_json(receipt, result)
    print(json.dumps(result))
    return 0 if code == 0 and not timed_out and all(result['inputsStable'].values()) else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    raise SystemExit(run(parser.parse_args().folder))
