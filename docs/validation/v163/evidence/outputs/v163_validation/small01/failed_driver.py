"""模块用途。复用已验证输入，完成具名观察包的小型门槛及一次有限真实完整域任务。"""
import argparse
import copy
import json
from pathlib import Path
import sys

import numpy as np

from auditSU2IdealAirState import audit, load_state
from prepareSU2IdealAirRestart import prepare
from su2IdealAirInterface import OBSERVED_FLAVOR, validate
from su2ThermalRuntime import ROOT, read, resolve, sha
from validateIdealAir import run, save

BASE = ROOT/'outputs/v163_validation'
REFERENCE = ROOT/'outputs/embedded_workbench/24915f33-5728-439b-884c-17a471f21457'


def launch(folder, seconds, parent_seconds):
    """原监督器保存真实退出与完整帧，外层复用进程树RSS/墙钟限制。"""
    names = ['scripts/runSU2ObservedStudy.py', 'scripts/prepareSU2IdealAirRestart.py',
        'scripts/su2IdealAirInterface.py', 'scripts/su2IdealAirCheckpoints.py',
        'scripts/auditSU2IdealAirState.py', 'scripts/runEmbeddedSU2.py',
        'scripts/su2ThermalRuntime.py', 'scripts/validateIdealAir.py', 'docs/validation/v163/PLAN.md']
    for name in names:
        target = folder/'executed_sources'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((ROOT/name).read_bytes())
    save(folder/'engine_request.json', dict(workDirectory=str(folder), configFile='case.cfg', meshFile='mesh.su2',
        maximumSeconds=seconds, threads=2, sourceFiles=names, jobToken='v163-'+folder.name))
    run(folder, 'parent', [sys.executable, str(ROOT/'scripts/runEmbeddedSU2.py'),
        '--request', str(folder/'engine_request.json')], parent_seconds)


def small():
    """小型数值保持、保存态和运行包身份检查，失败立即阻断真实任务。"""
    folder = BASE/'small01'
    prepare(REFERENCE/'ideal_air_initial_state', folder, 4, 2, 32*1024**2,
        snapshot=ROOT/'outputs/v157_validation/executed_sources03', runtime_flavor=OBSERVED_FLAVOR)
    text = (folder/'case.cfg').read_text('utf8')
    _, runtime, flavor, _, _ = resolve(folder, text)
    assert flavor == OBSERVED_FLAVOR
    manifest = read(folder/'case_manifest.json')
    rejected = []
    for name in ('precision', 'baseline', 'missing-observation', 'wrong-runtime'):
        m, r = copy.deepcopy(manifest), copy.deepcopy(runtime)
        if name == 'precision':
            r['asciiSignificantDigits'] = 16
        elif name == 'baseline':
            r['numericalBaselineFlavor'] = 'rah-ideal-air-axisymmetric-v1'
        elif name == 'missing-observation':
            m.pop('runtimeObservation')
        else:
            r['flavor'] = 'rah-ideal-air-phase-v2'
        try:
            validate(folder, m, text, r)
        except (AssertionError, KeyError):
            rejected.append(name)
    assert len(rejected) == 4
    save(folder/'contract_rejections.json', dict(passed=True, rejected=rejected))
    launch(folder, 30, 120)
    errors = {}
    locks = {}
    for step in (0, 2, 3):
        name = 'restart_flow_%06d.csv' % step
        old, old_fields, nd = load_state(REFERENCE/name, True, captured=True)
        new, new_fields, new_nd = load_state(folder/name, True, captured=True)
        assert list(old_fields) == list(new_fields) and new_nd == nd == 3 and old.shape == new.shape
        difference = float(np.max(abs(new-old)/np.maximum(1, abs(old))))
        errors[str(step)] = difference
        assert difference < 2e-13, '同源小型字段变化超出输出舍入门槛'
        locks[str(REFERENCE/name)] = sha(REFERENCE/name)
        locks[str(folder/name)] = sha(folder/name)
    thermal = audit(folder)
    save(BASE/'small_report.json', dict(passed=True, fieldRelativeErrors=errors,
        thermalErrors=thermal['errors'], contractRejections=rejected, sourceSHA256=locks,
        physicsQualified=False, newQualifiedArea=0))
    print(json.dumps(dict(smallPassed=True, fieldRelativeErrors=errors, thermalErrors=thermal['errors'])))


def real():
    """单次真实域从源149冷重启，结果与异常原样保留，禁止自动续算。"""
    assert read(BASE/'small_report.json')['passed']
    folder = BASE/'real01'
    prepare(ROOT/'outputs/v160_validation/real01', folder, 50, 49, 1024**3,
            checkpoint=149, linear_history=True, runtime_flavor=OBSERVED_FLAVOR)
    old = ROOT/'outputs/v161_validation/control'
    assert sha(folder/'initial_flow.csv') == sha(old/'initial_flow.csv')
    assert sha(folder/'case.cfg') == sha(old/'case.cfg')
    assert sha(folder/'mesh.su2') == sha(old/'mesh.su2')
    save(folder/'same_source_preflight.json', dict(passed=True,
        sourceDirectory=str(old), fields=['initial_flow.csv', 'case.cfg', 'mesh.su2'],
        sourceSHA256={str(old/name): sha(old/name) for name in ('initial_flow.csv', 'case.cfg', 'mesh.su2')}))
    launch(folder, 300, 420)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=('small', 'real'))
    args = parser.parse_args()
    small() if args.case == 'small' else real()
