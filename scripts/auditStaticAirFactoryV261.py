"""核对静止空气配置、节点回读和异常恢复，也可从已归档的小文件复核。"""
import argparse
import csv
import itertools
import json
import sys
from decimal import Decimal as D, localcontext
from pathlib import Path, PureWindowsPath

import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT, read, save, sha
from terminalCapacityV249 import supervise

BASE = ROOT / 'outputs/v261_validation'
DOC = ROOT / 'docs/validation/v261'
TABLE = ROOT / 'data/thermo/phase_aware_dry_air_v1.json'
STAGES = ('prepare', 'build', 'default', 'stable', 'unknown', 'incompatible')
METHODS = {'default': 'TOMS748', 'stable': 'STABLE_PIECEWISE'}
ERRORS = {'unknown': ('SetConfig_Parsing', 'RAH_CALORIC_INVERSE expects STABLE_PIECEWISE'),
          'incompatible': ('SetPostprocessing', 'RAH_CALORIC_INVERSE requires FLUID_MODEL=RAH_PHASE_AIR')}
T, P = D('300.123'), D(100000)
LIMIT = D('2e-10')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def reference():
    """从原始 Cp 节点做梯形积分，使用配置的十进制温度和独立的零焓参考。"""
    require(sha(TABLE) == '7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985', 'Cp input version')
    data = json.loads(TABLE.read_text(encoding='utf8'), parse_float=D)
    grid = data['temperatureK']
    cp = data['cpJkgK']
    start = D('298.15')
    require(start in grid and grid[-1] >= T > start, 'reference temperature grid')
    h = D(0)
    for i in range(grid.index(start), len(grid) - 1):
        left, right = grid[i:i + 2]
        end = min(right, T)
        end_cp = cp[i] + (cp[i + 1] - cp[i]) * (end - left) / (right - left)
        h += (cp[i] + end_cp) * (end - left) / 2
        if end == T:
            break
    rt = data['gasConstantJkgK'] * T
    return P / rt, h - rt + D(300000)


def csv_rows(path):
    with path.open(encoding='utf8', newline='') as stream:
        return list(csv.DictReader(stream))


def inspect_frames(folder, method, rho_ref, energy_ref):
    """逐节点比对两个快照，同时检查编号、空间覆盖和密度加权能量的单位。"""
    maxima = {key: D(0) for key in ('temperatureAbsK', 'pressureRel', 'densityRel', 'energyRel', 'energyAbsJkg')}
    coords = {}
    fields = ('rho', 'rhoUx', 'rhoUy', 'rhoUz', 'rhoE_internal')
    first = folder / 'initial.csv'
    require(first.read_bytes() == (folder / 'terminal.csv').read_bytes(), f'{method}: frame changed')
    for row in csv_rows(first):
        point = int(row['id'])
        require(point not in coords and row['method'] == method, f'{method}: node identity {point}')
        coords[point] = tuple(D(row[axis]) for axis in ('x', 'y', 'z'))
        rho, mx, my, mz, density_energy = [D.from_float(float(row[key])) for key in fields]
        require(rho.is_finite() and rho > 0 and mx == my == mz == 0, f'{method}: static state {point}')
        t = D.from_float(float(row['T_K']))
        p = D.from_float(float(row['p_Pa']))
        de = abs(density_energy / rho - energy_ref)
        errors = (abs(t - T), abs(p - P) / P, abs(rho - rho_ref) / rho_ref, de / abs(energy_ref), de)
        for key, error in zip(maxima, errors):
            require(error.is_finite(), f'{method}: nonfinite {key}')
            maxima[key] = max(maxima[key], error)
    axis = [D(i) / 4 for i in range(5)]
    require(set(coords) == set(range(125)) and set(coords.values()) == set(itertools.product(axis, repeat=3)),
            f'{method}: missing nodes or coordinates')
    for key in ('temperatureAbsK', 'pressureRel', 'densityRel', 'energyRel'):
        require(maxima[key] <= LIMIT, f'{method}: {key} {maxima[key]}')
    return coords, {k: float(v) for k, v in maxima.items()}


def audit(base):
    prepared = read(base / 'prepared.json')
    binaries = read(base / 'build/report.json')['executables']
    rho, energy = reference()
    cases = {}
    processes = {}
    default_coords = None
    for stage in STAGES:
        folder = base / stage
        proc = base / ('process_' + stage)
        receipt = read(proc / 'report.json')
        cap = read(proc / 'terminal_capacity.json')
        expected_exit = 1 if stage in ERRORS else 0
        require(receipt['status'] == 'terminal' and receipt['exitCode'] == expected_exit
                and receipt['inputsStable'] and receipt['stopReason'] is None
                and receipt['liveGroupAfterCleanup'] == [] and cap['passed'], f'{stage}: process receipt')
        processes[stage] = {k: receipt[k] for k in ('exitCode', 'elapsedSeconds', 'peakObservedGroupRSSBytes')}
        if stage in ('prepare', 'build'):
            continue
        exe = PureWindowsPath(receipt['argv'][0])
        expected_exe = 'SU2_CFD.exe' if stage in ERRORS else 'inspect.exe'
        require(exe.name == expected_exe and receipt['inputSHA256'][str(exe)] == binaries[expected_exe]
                and receipt['argv'][1:] == ['case.cfg'], f'{stage}: executable')
        cfg = folder / 'case.cfg'
        run_cfg = exe.parent.parent / stage / 'case.cfg'
        require(sha(cfg) == prepared['configs'][stage] == receipt['inputSHA256'][str(run_cfg)],
                f'{stage}: executed input')
        log = (proc / 'stdout.txt').read_text(encoding='utf8')
        if stage in ERRORS:
            where, why = ERRORS[stage]
            require(where in log and why in log and 'Error Exit' in log and 'RAH caloric initialization:' not in log,
                    f'{stage}: unexpected rejection')
            cases[stage] = {'rejectionLocation': where, 'reason': why}
            continue
        method = METHODS[stage]
        require(log.splitlines().count('RAH caloric initialization: ' + method) == 1
                and 'Exit Success (SU2_CFD)' in log, f'{stage}: initialization or finalization')
        obs = read(folder / 'initialization.json')
        require(obs['workerMethod'] == method and obs['threadsObserved'] == 1 and obs['nodes'] == 125
                and obs['nonlinearUpdates'] == obs['timeUpdates'] == 0, f'{stage}: scope')
        stable = stage == 'stable'
        require(obs['stableRecoveryCallsInFrames'] == (250 if stable else 0)
                and obs['inversePath'] == ('stable-piecewise-quadratic' if stable else 'legacy-TOMS748'),
                f'{stage}: method dispatch')
        coords, errors = inspect_frames(folder, method, rho, energy)
        if default_coords is None:
            default_coords = coords
        require(coords == default_coords, f'{stage}: geometry mismatch')
        events = csv_rows(folder / 'admissibility.csv')
        require(len(events) == 1 and events[0]['reason'] == 'RAH_PHASE_ENERGY_ENTROPY_DOMAIN'
                and events[0]['recovered'] == '1' and D(events[0]['trial_e']) == -1000000
                and obs['invalidCandidateAccepted'] is False and obs['oldConservedStateRestored'] is True
                and abs(D.from_float(obs['restoredTemperatureK']) - T) <= LIMIT, f'{stage}: exception recovery')
        cases[stage] = {'observation': obs, 'maxErrors': errors, 'rejectionReason': events[0]['reason']}
    return {'passed': True, 'reference': {'temperatureK': str(T), 'pressurePa': str(P), 'densityKgm3': str(rho),
            'internalEnergyJkg': str(energy), 'tableSHA256': sha(TABLE)}, 'cases': cases, 'processes': processes,
            'newQualifiedSurfaceAreaM2': 0, 'flowAdvanced': False}


def run_supervised():
    """锁定配置、原始帧和进程回执，30 秒内完成离线核对。"""
    guard.BASE = guard.CASE = BASE
    guard.DOC = DOC
    count_bytes = guard.size_bytes
    guard.size_bytes = lambda p: count_bytes(BASE) if p.name.startswith('process_') else 0
    script = Path(__file__).resolve()
    locks = [script, Path(guard.__file__), ROOT / 'scripts/terminalCapacityV249.py', DOC / 'PLAN.md',
             TABLE, BASE / 'prepared.json', BASE / 'build/report.json']
    for stage in STAGES:
        proc = BASE / ('process_' + stage)
        locks.extend(proc / name for name in ('report.json', 'terminal_capacity.json', 'stdout.txt'))
        if stage in METHODS or stage in ERRORS:
            locks.append(BASE / stage / 'case.cfg')
        if stage in METHODS:
            locks.extend(BASE / stage / name for name in ('initialization.json', 'initial.csv', 'terminal.csv',
                                                         'admissibility.csv'))
    return supervise(guard, BASE, 'audit', [sys.executable, str(script), '--worker'], 30, 160, locks, ROOT, 'audit')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--worker', action='store_true')
    mode.add_argument('--replay', type=Path, help='只读复核保留原目录结构的证据')
    args = parser.parse_args()
    if args.worker or args.replay:
        with localcontext() as ctx:
            ctx.prec = 60
            result = audit(args.replay or BASE)
        if args.worker:
            save(BASE / 'audit/report.json', result)
        print(json.dumps(result, indent=2))
    else:
        raise SystemExit(0 if run_supervised() else 1)
