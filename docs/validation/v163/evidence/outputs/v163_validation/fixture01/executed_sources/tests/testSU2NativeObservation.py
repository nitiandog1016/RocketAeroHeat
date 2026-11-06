"""模块用途。限时构建原生夹具，核对异常终止和C++文本精度，保存准确来源。"""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from validateIdealAir import BIN, run, save, sha


def main():
    folder = ROOT/'outputs/v163_validation/fixture01'
    folder.mkdir(parents=True, exist_ok=False)
    paths = [Path(__file__).resolve(), ROOT/'tests/testSU2NativeObservation.cpp',
             ROOT/'third_party/su2_thermal/NativeTermination.hpp', ROOT/'docs/validation/v163/PLAN.md']
    locks = {str(path): sha(path) for path in paths}
    for path in paths:
        target = folder/'executed_sources'/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
    executable = folder/'fixture.exe'
    run(folder, 'compile', [str(BIN/'g++.exe'), '-std=c++17', '-O3', '-fopenmp',
        str(ROOT/'tests/testSU2NativeObservation.cpp'), '-o', str(executable)], 30)
    results = {}
    env = dict(os.environ, PATH=str(BIN)+os.pathsep+os.environ['PATH'])
    for mode in ('precision', 'std', 'omp', 'nonstd'):
        started = time.monotonic()
        process = subprocess.run([str(executable), mode], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            env=env, timeout=10, creationflags=subprocess.CREATE_NO_WINDOW)
        (folder/(mode+'.log')).write_bytes(process.stdout)
        output = process.stdout.decode('utf8', errors='replace')
        results[mode] = dict(exitCode=process.returncode, elapsedSeconds=time.monotonic()-started)
        save(folder/'runs.json', results)
        if mode == 'precision':
            assert process.returncode == 0
            results[mode].update(json.loads(output))
            assert results[mode]['newExact'] == 9 and results[mode]['maximumScalarBytes'] <= 32
        else:
            assert process.returncode != 0
            assert output.count('RAH_NATIVE_TERMINATION_BEGIN') == 1
            assert output.count('RAH_NATIVE_TERMINATION_END') == 1
            assert 'exit_action=delegate-original-terminate' in output and 'frame=0 address=' in output
            expected = {'std': 'what=rah-fixture-domain-error', 'omp': 'what=rah-fixture-openmp-error',
                        'nonstd': 'exception_type=non-std'}[mode]
            assert expected in output
    stable = {name: sha(Path(name)) == digest for name, digest in locks.items()}
    assert all(stable.values())
    save(folder/'report.json', dict(passed=True, sourceSHA256=locks, inputsStable=stable,
                                  executableSHA256=sha(executable), cases=results))
    print(json.dumps(results))


if __name__ == '__main__':
    main()
