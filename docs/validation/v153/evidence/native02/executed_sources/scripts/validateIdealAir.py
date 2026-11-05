"""v153有限原生热物性检查，绑定来源并保存准确执行源和真实终态。"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import time

import numpy as np
import psutil
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/v153_validation'
SOURCE = ROOT / 'tools/runtime/su2-thermal-build/source/subprojects/CoolProp'
BIN = ROOT / 'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin'
LIBRARY = ROOT / 'tools/runtime/su2-coolprop-build/cfd-build/subprojects/CoolProp/libCoolProp.a'
R = 8.31451 / .02896546


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf8')


def run(folder, name, command, seconds):
    """只管理本次创建的进程树，检查墙钟和总工作集，保留失败回执。"""
    receipt = folder / (name + '_process.json')
    assert not receipt.exists(), 'Existing process receipt'
    state = dict(command=list(map(str, command)), maximumSeconds=seconds,
                 maximumMemoryBytes=4 * 1024**3, status='starting')
    save(receipt, state)
    env = dict(os.environ, PATH=str(BIN) + os.pathsep + os.environ['PATH'],
               OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    started = time.monotonic()
    with (folder / (name + '.log')).open('wb') as log:
        proc = subprocess.Popen(command, cwd=str(ROOT), env=env, stdout=log,
                                stderr=log, creationflags=subprocess.CREATE_NO_WINDOW)
        state.update(status='running', processID=proc.pid)
        save(receipt, state)
        peak = 0
        stopped = None
        while proc.poll() is None:
            try:
                parent = psutil.Process(proc.pid)
                family = [parent] + parent.children(recursive=True)
                memory = 0
                for child in family:
                    try:
                        memory += child.memory_info().rss
                    except psutil.NoSuchProcess:
                        pass
                peak = max(peak, memory)
            except psutil.NoSuchProcess:
                family = []
            if time.monotonic() - started > seconds or peak > 4 * 1024**3:
                stopped = 'time' if time.monotonic() - started > seconds else 'memory'
                subprocess.run(['taskkill', '/PID', str(proc.pid), '/T', '/F'], stdout=log, stderr=log)
                break
            time.sleep(.05)
        code = proc.wait()
    state.update(status='terminal', exitCode=code, elapsedSeconds=time.monotonic()-started,
                 sampledPeakMemoryBytes=peak, stoppedByBudget=stopped)
    save(receipt, state)
    print(json.dumps(dict(name=name, exitCode=code, elapsedSeconds=state['elapsedSeconds'])), flush=True)
    assert code == 0 and stopped is None, 'Process failed; retain receipt'


def write_rows(path, rows):
    with path.open('w', newline='', encoding='ascii') as stream:
        writer = csv.writer(stream)
        writer.writerow(['mode', 'a', 'b'])
        writer.writerows(rows)


def read_rows(path):
    return np.genfromtxt(str(path), delimiter=',', names=True)


def relative(actual, expected):
    return float(np.max(np.abs(actual-expected) / np.maximum(1., np.abs(expected))))


def paper_cp(temperature):
    """独立原论文式18/Table12，以温度表达式复核Helmholtz导数。"""
    t = temperature
    u, v, w = 3364.011/t, 2242.45/t, 11580.4/t
    return R*(3.490888032 + 2.395525583e-6*t + 7.172111248e-9*t*t
              - 3.115413101e-13*t**3 + .223806688/t**1.5
              + .791309509*u*u*np.exp(-u)/(-np.expm1(-u))**2
              + .212236768*v*v*np.exp(-v)/(-np.expm1(-v))**2
              + (2./3.)*.197938904*w*w*np.exp(-w)/(1+(2./3.)*np.exp(-w))**2)


def validate(folder, executable):
    """先独立原生参照，再双向恢复、差分、积分与拒绝测试。"""
    temperatures = np.unique(np.r_[np.linspace(60, 2000, 101), 60+1e-8, 80, 150, 200, 300, 1000, 2000-1e-8])
    rows = [('rhoT', float(rho), float(t)) for rho in (1e-5, .02, 1., 10.) for t in temperatures]
    write_rows(folder / 'states.csv', rows)
    for mode in ('core', 'oracle'):
        run(folder, mode, [str(executable), '--'+mode, str(folder/'states.csv'), str(folder/(mode+'.csv'))], 60)
    core, oracle = read_rows(folder/'core.csv'), read_rows(folder/'oracle.csv')
    assert np.all(core['accepted'] == 1) and np.all(oracle['accepted'] == 1)
    fields = ('T', 'rho', 'p', 'e', 'h', 's', 'cp', 'cv', 'sound')
    native_errors = {name: relative(core[name], oracle[name]) for name in fields}
    paper_error = relative(core['cp'], paper_cp(core['T']))
    identities = {'cp_cv': relative(core['cp']-core['cv'], np.full(len(core), R)),
                  'enthalpy': relative(core['h']-core['e'], R*core['T']),
                  'eos': relative(core['p'], core['rho']*R*core['T'])}
    inverse_rows = []
    for item in core:
        inverse_rows.extend([('PT', item['p'], item['T']), ('Prho', item['p'], item['rho']),
                             ('rhoe', item['rho'], item['e']), ('hs', item['h'], item['s']),
                             ('Ps', item['p'], item['s'])])
    write_rows(folder/'inverse.csv', inverse_rows)
    run(folder, 'inverse', [str(executable), '--core', str(folder/'inverse.csv'), str(folder/'recovered.csv')], 60)
    recovered = read_rows(folder/'recovered.csv')
    accepted = recovered['accepted'] == 1
    inverse_errors = {name: relative(recovered[name][accepted], np.repeat(core[name], 5)[accepted]) for name in fields}

    # 每个偏导均在其指定自变量保持条件下作独立中心差分。
    difference_rows, checks = [], []
    selected = core[(core['T'] > 61) & (core['T'] < 1999)][::13]
    for item in selected:
        rho, t, p, e, s = [float(item[name]) for name in ('rho', 'T', 'p', 'e', 's')]
        definitions = [
            ('rhoT', rho, t, 1, 'e', item['cv']),
            ('rhoT', rho, t, 1, 'h', item['cp']),
            ('rhoT', rho, t, 1, 'cp', item['dcp_dt']),
            ('rhoT', rho, t, 1, 'mu', item['dmu_dt']),
            ('rhoT', rho, t, 1, 'lambda', item['dlambda_dt']),
            ('rhoe', rho, e, 0, 'p', item['dp_drho_e']),
            ('rhoe', rho, e, 1, 'p', item['dp_de_rho']),
            ('rhoe', rho, e, 0, 'T', item['dt_drho_e']),
            ('rhoe', rho, e, 1, 'T', item['dt_de_rho']),
            ('Prho', p, rho, 1, 'h', item['dh_drho_p']),
            ('Prho', p, rho, 0, 'h', item['dh_dp_rho']),
            ('Prho', p, rho, 1, 's', item['ds_drho_p']),
            ('Prho', p, rho, 0, 's', item['ds_dp_rho']),
            ('Ps', p, s, 0, 'rho', 1/item['sound']**2)]
        for mode, a, b, variable, output, expected in definitions:
            values = [a, b]
            step = abs(values[variable])*1e-5
            minus, plus = values[:], values[:]
            minus[variable] -= step
            plus[variable] += step
            difference_rows.extend([(mode, *minus), (mode, *plus)])
            checks.append((output, step, expected))
    write_rows(folder/'difference.csv', difference_rows)
    run(folder, 'difference', [str(executable), '--core', str(folder/'difference.csv'), str(folder/'difference_output.csv')], 60)
    differences = read_rows(folder/'difference_output.csv')
    assert np.all(differences['accepted'] == 1)
    derivative_errors = {}
    for index, (field, step, expected) in enumerate(checks):
        actual = (differences[field][2*index+1]-differences[field][2*index])/(2*step)
        error = abs(actual-expected)/max(1., abs(expected))
        derivative_errors[field] = max(derivative_errors.get(field, 0.), float(error))

    invalid = [('rhoT', 0, 300), ('rhoT', -1, 300), ('rhoT', 1, 59.999),
               ('rhoT', 1, 2000.001), ('PT', 0, 300), ('PT', 1, 0), ('Prho', -1, 1),
               ('Prho', 1, 0), ('rhoe', 1, core['e'].min()-1),
               ('rhoe', 1, core['e'].max()+1), ('hs', core['h'].min()-1, 0),
               ('hs', core['h'].max()+1, 0), ('hs', 300000, 1e300),
               ('Ps', 1, -1e10), ('Ps', 1, 1e10), ('rhoT', 'nan', 300),
               ('PT', 1, 'inf'), ('invalid', 1, 1), ('rhoT', 1e308, 300)]
    write_rows(folder/'invalid.csv', invalid)
    run(folder, 'invalid', [str(executable), '--core', str(folder/'invalid.csv'), str(folder/'invalid_output.csv')], 60)
    rejected = read_rows(folder/'invalid_output.csv')

    integration_errors = dict(e=0., h=0., s=0.)
    reference = core[0]
    for item in core[:len(temperatures):7]:
        t = float(item['T'])
        expected = [quad(lambda v: paper_cp(v)-R, 60, t, epsabs=1e-6)[0],
                    quad(paper_cp, 60, t, epsabs=1e-6)[0],
                    quad(lambda v: (paper_cp(v)-R)/v, 60, t, epsabs=1e-8)[0]]
        for field, integral in zip(('e', 'h', 's'), expected):
            error = abs((item[field]-reference[field])-integral)/max(1., abs(integral))
            integration_errors[field] = max(integration_errors[field], float(error))

    run(folder, 'benchmark', [str(executable), '--benchmark', str(folder/'states.csv'), str(folder/'benchmark.csv')], 60)
    benchmark = read_rows(folder/'benchmark.csv')[1:]
    timings = {name: float(np.median(benchmark[name])) for name in ('analytic_rhoT_s', 'analytic_rhoe_s', 'HEOS_rhoT_s')}
    ratios = {name: timings['HEOS_rhoT_s']/timings[name] for name in ('analytic_rhoT_s', 'analytic_rhoe_s')}
    result = dict(states=len(core), inverseStates=len(recovered), rejectedInverseRows=np.flatnonzero(~accepted).tolist(),
                  nativeErrors=native_errors, paperCpError=paper_error, identityErrors=identities,
                  inverseErrors=inverse_errors, derivativeCases=len(checks), derivativeErrors=derivative_errors,
                  integralErrors=integration_errors, invalidCount=len(invalid),
                  invalidRejected=bool(np.all(rejected['accepted'] == 0)),
                  medianSecondsPer20000=timings, speedupVsHEOS=ratios,
                  performanceGoalPassed=all(value >= 5 for value in ratios.values()))
    result['numericalPassed'] = bool(max(native_errors.values()) <= 2e-11 and paper_error <= 2e-8
        and max(identities.values()) <= 2e-12 and np.all(accepted) and max(inverse_errors.values()) <= 2e-10
        and max(derivative_errors.values()) <= 3e-6 and max(integration_errors.values()) <= 3e-8
        and result['invalidRejected'])
    save(folder/'summary.json', result)
    print(json.dumps(result, indent=2), flush=True)
    assert result['numericalPassed'], 'Numerical gates failed; preserve exact sources'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--attempt', required=True)
    args = parser.parse_args()
    assert args.attempt.isdigit()
    folder = OUT / ('native' + args.attempt)
    folder.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__).resolve(), ROOT/'third_party/ideal_air/IdealAir.hpp', ROOT/'tests/idealAirNativeProbe.cpp',
             ROOT/'third_party/ideal_air/air_ideal.json', ROOT/'third_party/ideal_air/LICENSE.CoolProp',
             ROOT/'docs/validation/v153/PLAN.md', SOURCE/'dev/fluids/Air.json', SOURCE/'include/AbstractState.h',
             LIBRARY, BIN/'g++.exe']
    lock = {p.relative_to(ROOT).as_posix(): sha(p) for p in files}
    save(folder/'source_lock.json', lock)
    for path in files[:6]:
        destination = folder/'executed_sources'/path.relative_to(ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(str(path), str(destination))
    assert sha(SOURCE/'dev/fluids/Air.json') == 'c5e0d064d83fc91ed35d7354cea1fa38fea832a09d2e85faf4336e0017c33616'
    executable = folder/'idealAirProbe.exe'
    try:
        command = [str(BIN/'g++.exe'), '-std=c++17', '-O2', '-I'+str(SOURCE/'include'),
                   '-I'+str(SOURCE/'externals/fmtlib/include'), str(ROOT/'tests/idealAirNativeProbe.cpp'),
                   str(LIBRARY), '-o', str(executable)]
        run(folder, 'compile', command, 120)
        save(folder/'runtime.json', dict(path=executable.relative_to(ROOT).as_posix(),
                                       sha256=sha(executable), bytes=executable.stat().st_size))
        validate(folder, executable)
    finally:
        stable = {path: sha(ROOT/path) == value for path, value in lock.items()}
        save(folder/'source_stability.json', stable)
        assert all(stable.values()), 'Source changed during validation'


if __name__ == '__main__':
    main()
