"""模块用途。新数值回退合同的小型门槛与唯一同源整箭有限计算。"""
import argparse
import copy
import json
from pathlib import Path
import sys

import numpy as np

from auditSU2IdealAirState import audit, load_state
from prepareSU2IdealAirRestart import prepare
from su2IdealAirInterface import ADMISSIBLE_FLAVOR, validate
from su2ThermalRuntime import ROOT, read, resolve, sha
from validateIdealAir import run, save

BASE = ROOT/'outputs/v164_validation'
OLD = ROOT/'outputs/v163_validation'
INITIAL = ROOT/'outputs/embedded_workbench/24915f33-5728-439b-884c-17a471f21457/ideal_air_initial_state'


def launch(folder, seconds, parent_seconds):
    """模块一。绑定准确脚本/计划；原生任务和父进程均有墙钟/内存预算。"""
    names = ['scripts/runSU2AdmissibleStudy.py', 'scripts/prepareSU2IdealAirRestart.py',
        'scripts/su2IdealAirInterface.py', 'scripts/su2IdealAirCheckpoints.py',
        'scripts/auditSU2IdealAirState.py', 'scripts/runEmbeddedSU2.py',
        'scripts/su2ThermalRuntime.py', 'scripts/validateIdealAir.py', 'docs/validation/v164/PLAN.md']
    for name in names:
        target = folder/'executed_sources'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((ROOT/name).read_bytes())
    save(folder/'engine_request.json', dict(workDirectory=str(folder), configFile='case.cfg', meshFile='mesh.su2',
        maximumSeconds=seconds, threads=2, sourceFiles=names, jobToken='v164-'+folder.name))
    run(folder, 'parent', [sys.executable, str(ROOT/'scripts/runEmbeddedSU2.py'),
        '--request', str(folder/'engine_request.json')], parent_seconds)


def small():
    """模块二。合法路径保持、独立热审核及数值身份反例，阻止未验证包进入整箭。"""
    folder = BASE/'small01'
    prepare(INITIAL, folder, 4, 2, 32*1024**2,
        snapshot=OLD/'initial_source', runtime_flavor=ADMISSIBLE_FLAVOR)
    text = (folder/'case.cfg').read_text('utf8')
    _, runtime, flavor, _, _ = resolve(folder, text)
    assert flavor == ADMISSIBLE_FLAVOR
    manifest = read(folder/'case_manifest.json')
    rejected = []
    for name in ('precision', 'wrong-baseline', 'no-numerical-change', 'unbounded-trace',
                 'temperature-clipping', 'missing-contract', 'same-numerics'):
        m, r = copy.deepcopy(manifest), copy.deepcopy(runtime)
        if name == 'precision':
            r['asciiSignificantDigits'] = 16
        elif name == 'wrong-baseline':
            r['numericalBaselineFlavor'] = 'rah-ideal-air-observed-v1'
        elif name == 'no-numerical-change':
            r['samePhysicsAndNumerics'] = True
        elif name == 'unbounded-trace':
            r['admissibility']['maximumTraceEvents'] = 0
        elif name == 'temperature-clipping':
            r['admissibility']['temperatureClipping'] = True
        elif name == 'missing-contract':
            m.pop('runtimeAdmissibility')
        else:
            m['initialization']['samePhysicsAndNumerics'] = True
        try:
            validate(folder, m, text, r)
        except (AssertionError, KeyError):
            rejected.append(name)
    assert len(rejected) == 7
    save(folder/'contract_rejections.json', dict(passed=True, rejected=rejected))
    launch(folder, 30, 120)
    errors, locks = {}, {}
    for step in (0, 2, 3):
        name = 'restart_flow_%06d.csv' % step
        old, old_fields, nd = load_state(OLD/'small02'/name, True, captured=True)
        new, new_fields, new_nd = load_state(folder/name, True, captured=True)
        assert list(old_fields) == list(new_fields) and new_nd == nd == 3 and old.shape == new.shape
        errors[str(step)] = float(np.max(abs(new-old)/np.maximum(1, abs(old))))
        assert errors[str(step)] < 2e-13
        for path in (OLD/'small02'/name, folder/name):
            locks[str(path)] = sha(path)
    assert not (folder/'admissibility.csv').exists(), '小型合法路径发生未预期回退'
    thermal = audit(folder)
    save(BASE/'small_report.json', dict(passed=True, fieldRelativeErrors=errors,
        thermalErrors=thermal['errors'], contractRejections=rejected, sourceSHA256=locks,
        newFallbackEvents=0, physicsQualified=False, newQualifiedArea=0))
    print(json.dumps(dict(smallPassed=True, fieldRelativeErrors=errors, thermalErrors=thermal['errors'])))


def real():
    """模块三。只使用同源149、原几何与原设置，有限50更新后停止，不自动续算。"""
    assert read(BASE/'small_report.json')['passed']
    folder = BASE/'real01'
    prepare(ROOT/'outputs/v160_validation/real01', folder, 50, 50, 1024**3,
            checkpoint=149, linear_history=True, runtime_flavor=ADMISSIBLE_FLAVOR)
    old = OLD/'real02'
    names = ('initial_flow.csv', 'case.cfg', 'mesh.su2')
    for name in names:
        assert sha(folder/name) == sha(old/name)
    save(folder/'same_source_preflight.json', dict(passed=True, samePhysicsAndNumerics=False,
        sourceDirectory=str(old), fields=list(names),
        sourceSHA256={str(old/name): sha(old/name) for name in names}))
    launch(folder, 300, 420)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=('small', 'real'))
    args = parser.parse_args()
    small() if args.case == 'small' else real()
