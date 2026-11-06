"""模块用途：空间Tw限步组合的严格正向/拒绝检查；不以新增保护声明替代实际能力。"""
import argparse
import copy
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from prepareNACASpatialStep import FLAVOR, source_files
from prepareSU2NACALaminar import STAGE
from runEmbeddedSU2 import validate_resource_budget
from su2NACASpatialStepInterface import validate
from su2ThermalRuntime import read, sha, resolve
from validateIdealAir import save


def check(output):
    """模块一：固定准确源与已准备小例；模块二：真实resolve；模块三：独占拒绝夹具。"""
    started = time.monotonic(); out = Path(output).resolve(); out.relative_to(ROOT/'outputs')
    out.mkdir(parents=True, exist_ok=False)
    source = ROOT/'outputs/v173_validation/small01'
    names = ('step_preparation.json', 'step_initialization.json', 'case_manifest.json', 'ideal_air_contract.json',
        'case.cfg', 'mesh.su2', 'mesh_arrays.npz', 'geometry.json', 'wall_temperature.csv', 'initial_flow.csv', 'request.json')
    originals = [source/name for name in names]
    manifest = read(source/'case_manifest.json'); request = read(source/'request.json')
    text = (source/'case.cfg').read_text(encoding='ascii')
    runtime = read(ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)/'runtime.json')
    paths = source_files() + [Path(__file__).resolve(), ROOT/'scripts/validateIdealAir.py']
    locks = dict(manifest['initialization']['sourceSHA256'])
    locks.update({str(path): sha(path) for path in originals + paths})
    save(out/'input_manifest.json', locks)
    for path in paths:
        target = out/'executed_source'/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(path.read_bytes())
    exe, _, flavor, profile, resolved = resolve(source, text)
    budget = validate_resource_budget(request, manifest)
    assert flavor == FLAVOR and profile == source/'wall_temperature.csv' and budget == (30., 2)
    assert sha(exe) != manifest['initialization']['sourceExecutableSHA256']
    assert runtime['numericalParentFlavor'] == manifest['initialization']['sourceRuntimeFlavor']
    source_audit = read(Path(manifest['initialization']['sourceAuditPath']))
    # 报告用途身份与原生保存列阶段是两个合同；不能再误作同一个字符串。
    assert source_audit['stage'] == 'same-flow-update-laminar-no-SST'
    assert manifest['initialization']['sourceStage'] == STAGE
    assert source_audit['stage'] != STAGE
    results = [dict(name='small01-real-resolve-and-budget', passed=True, accepted=True,
        locksValidated=len(resolved), auditStage=source_audit['stage'], savedColumnStage=STAGE,
        maximumSeconds=budget[0], threads=budget[1], mocked=False)]

    def fixture(name):
        folder = out/name; folder.mkdir(exist_ok=False)
        # 夹具可变输入均为独立字节副本，不使用硬链接，避免误写历史原件。
        for filename in names:
            (folder/filename).write_bytes((source/filename).read_bytes())
        return folder, read(folder/'case_manifest.json'), text, copy.deepcopy(runtime)

    def bind_initialization(folder, local, change):
        initialization = read(folder/'step_initialization.json')
        change(initialization)
        save(folder/'step_initialization.json', initialization)
        contract = read(folder/'ideal_air_contract.json')
        contract['initializationSHA256'] = sha(folder/'step_initialization.json')
        save(folder/'ideal_air_contract.json', contract)
        prepared = read(folder/'step_preparation.json'); prepared['initializationSHA256'] = contract['initializationSHA256']
        save(folder/'step_preparation.json', prepared)
        local['initialization'] = initialization; local['idealAir'] = contract
        local['idealAirContractSHA256'] = sha(folder/'ideal_air_contract.json')
        save(folder/'case_manifest.json', local)

    def reject(name, mutation):
        folder, local, config, run = fixture(name)
        changed = mutation(folder, local, config, run)
        if changed is not None:
            config = changed
        save(folder/'tested_manifest.json', local); save(folder/'tested_runtime.json', run)
        (folder/'tested_config.cfg').write_text(config, encoding='ascii')
        try:
            validate(folder, local, config, run)
        except (AssertionError, ValueError, KeyError, OSError, RuntimeError) as exc:
            record = dict(name=name, passed=True, accepted=False, reason=type(exc).__name__+': '+str(exc))
        else:
            raise AssertionError('未拒绝非法入口: '+name)
        save(folder/'rejection.json', record); results.append(record)

    reject('extra-CFL-parameter', lambda f, m, c, r: c.replace(
        'CFL_ADAPT_PARAM= (0.1, 1.2, 0.4, 100, 0.001, 0)', 'CFL_ADAPT_PARAM= (0.1, 1.2, 0.4, 100, 0.001, 0, 1)'))
    reject('relaxed-original-residual', lambda f, m, c, r: c.replace('CONV_RESIDUAL_MINVAL= -8', 'CONV_RESIDUAL_MINVAL= -7'))
    reject('changed-Mach', lambda f, m, c, r: c.replace('MACH_NUMBER= 1.97', 'MACH_NUMBER= 2.5'))
    reject('unregistered-SST', lambda f, m, c, r: m.update(sst={'model': 'SST-2003m'}))
    reject('changed-runtime-flavor', lambda f, m, c, r: r.update(flavor='rah-ideal-air-observed-v1'))
    reject('missing-step-contract', lambda f, m, c, r: r.update(conservativeStep={}))
    reject('changed-half-step-limit', lambda f, m, c, r: r['conservativeStep'].update(maximumHalvings=64))
    reject('false-CFL-coupled', lambda f, m, c, r: m['numericalChange'].update(effectiveBacktrackingFedToCFL=True))
    reject('false-boundary-recovery', lambda f, m, c, r: m['numericalChange'].update(riemannBoundaryRecovery=True))
    reject('false-whole-domain-protection', lambda f, m, c, r: m['numericalChange'].update(wholeQueryPathAdmissibilityGuaranteed=True))
    reject('changed-query-contract', lambda f, m, c, r: r['queryContext'].update(domainErrorTypePreserved=False))
    reject('false-numerical-equivalence', lambda f, m, c, r: r.update(samePhysicsAndNumerics=True))


    def changed_wall(folder, local, config, run):
        lines = (folder/'wall_temperature.csv').read_text(encoding='ascii').splitlines()
        x, tw = lines[1].split(','); lines[1] = x+','+format(float(tw)+1., '.17g')
        (folder/'wall_temperature.csv').write_text('\n'.join(lines)+'\n', encoding='ascii')
    reject('changed-original-Tw', changed_wall)
    for name, change in (
        ('false-source-step', lambda init: init.update(sourceLocalIteration=2)),
        ('false-cumulative-offset', lambda init: init.update(priorSuccessfulUpdates=3)),
        ('cancelled-source-claim', lambda init: init['sourceTermination'].update(cancelled=True)),
        ('unreviewed-source-claim', lambda init: init.update(sourceSavedFrameInterfacePassed=False)),
        ('false-source-audit-hash', lambda init: init.update(sourceAuditSHA256='0'*64)),
        ('audit-stage-as-saved-stage', lambda init: init.update(sourceStage='same-flow-update-laminar-no-SST')),
        ('false-uninterrupted-equivalence', lambda init: init.update(uninterruptedTrajectoryEquivalent=True)),
    ):
        reject(name, lambda f, m, c, r, change=change: bind_initialization(f, m, change))

    def changed_initial(folder, local, config, run):
        path = folder/'initial_flow.csv'; header = path.read_text(encoding='ascii').splitlines()[0]
        data = np.loadtxt(str(path), delimiter=',', skiprows=1); data[0, 4] *= 1.01
        np.savetxt(str(path), data, delimiter=',', fmt='%.17g', comments='', header=header)
        bind_initialization(folder, local, lambda init: init.update(initialStateSHA256=sha(path)))
    reject('changed-conservative-initial-with-rebound-hash', changed_initial)

    # 模块四：完整源proof不改；这里只证明旧固定配置拒绝新增自适应，不伪称旧resolve全来源重验。
    legacy = ROOT/'outputs/v171_validation/spatial_profile01'
    legacy_config = (legacy/'case.cfg').read_text(encoding='ascii')
    altered = legacy_config.replace('CFL_ADAPT= NO', 'CFL_ADAPT= YES')
    assert altered != legacy_config
    try:
        resolve(legacy, altered)
    except (AssertionError, ValueError, KeyError, OSError, RuntimeError) as exc:
        results.append(dict(name='old-v171-rejects-new-adaptive-config', passed=True, accepted=False,
            reason=type(exc).__name__+': '+str(exc), fullHistoricalResolveRequalified=False, mocked=False))
    else:
        raise AssertionError('旧v171合同放行自适应配置')
    for key, value in (('maximumSeconds', 31), ('threads', 3), ('resourceProfile', 'v171-naca-spatial-validation')):
        bad = copy.deepcopy(request); bad[key] = value
        try:
            validate_resource_budget(bad, manifest)
        except (AssertionError, ValueError, KeyError, RuntimeError) as exc:
            results.append(dict(name='reject-resource-'+key, passed=True, accepted=False, reason=type(exc).__name__+': '+str(exc)))
        else:
            raise AssertionError('预算入口未拒绝: '+key)
    stable = {path: sha(Path(path)) == digest for path, digest in locks.items()}; assert all(stable.values())
    report = dict(passed=True, cases=results, inputSHA256=locks, inputsStable=stable,
        acceptedCases=sum(row['accepted'] for row in results), rejectedCases=sum(not row['accepted'] for row in results),
        newCFDTasks=0, historicalTasksModified=False, fullOldResolveRequalified=False,
        physicalQualification=False, elapsedSeconds=time.monotonic()-started)
    save(out/'report.json', report)
    print(dict(passed=True, accepted=report['acceptedCases'], rejected=report['rejectedCases'], elapsedSeconds=report['elapsedSeconds']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', required=True)
    check(parser.parse_args().output)
