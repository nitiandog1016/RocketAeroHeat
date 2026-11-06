"""模块用途：只在独占小副本上检验新定比热入口及真实拒绝，不启动CFD或修改原任务。"""
import argparse
import copy
import csv
import io
from pathlib import Path
import sys
import time
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from nacaConstantCpModel import conditions
from prepareNACAConstantCpCase import KIND, FLAVOR
from su2NACAConstantCpInterface import validate, inlet_state_check, resource_budget
from su2IdealAirInterface import PHASE_FLAVORS
from su2ThermalRuntime import read, sha, resolve
from validateIdealAir import save


def verify(folder, out):
    """模块一：冻结已准备小任务，复制到独占目录；测试不修改真实准备物或父运行包。"""
    started = time.monotonic()
    folder, out = Path(folder).resolve(), Path(out).resolve()
    out.relative_to(ROOT / 'outputs')
    assert not out.exists(), '拒绝覆盖已有合同证据'
    assert not (folder / 'engine_process.json').exists(), '本项只在唯一原生小任务启动前执行'
    out.mkdir(parents=True, exist_ok=False)
    work = out / 'fixture'
    work.mkdir()
    locks, rejections = {}, []
    report = dict(passed=False, productionContractRefusalsTested=0, physicalQualified=False,
                  experimentalQualification=False, newQualifiedArea_m2=0, newCFDTasks=0, builds=0)

    def bind(path):
        path = Path(path).resolve()
        value = sha(path)
        if str(path) in locks:
            assert locks[str(path)] == value, '检查期间来源改变'
        locks[str(path)] = value
        return value

    try:
        source_names = ['tests/verifyNACAConstantCpContractV174.py',
            'scripts/prepareNACAConstantCpCase.py', 'scripts/su2NACAConstantCpInterface.py',
            'scripts/nacaConstantCpModel.py', 'scripts/su2ThermalRuntime.py',
            'scripts/su2IdealAirInterface.py', 'scripts/su2IdealAirCheckpoints.py',
            'scripts/runEmbeddedSU2.py']
        for name in source_names:
            path = ROOT / name
            bind(path)
            (out / ('executed_' + path.name)).write_bytes(path.read_bytes())
        names = ('case_manifest.json', 'ideal_gas_contract.json', 'case.cfg', 'mesh.su2',
                 'geometry.json', 'mesh_arrays.npz', 'wall_temperature.csv', 'initial_flow.csv',
                 'request.json', 'plan_at_prepare.md')
        original = {}
        for name in names:
            bind(folder / name)
            original[name] = (folder / name).read_bytes()
            (work / name).write_bytes(original[name])
        manifest = read(work / 'case_manifest.json')
        assert manifest['caseKind'] == KIND and manifest['specification']['variant'] == 'small-profile'
        assert manifest['nodes'] == 549 and manifest['initialization']['priorSuccessfulUpdates'] == 0
        runtime_path = ROOT / 'tools/runtime/su2' / ('8.6.0-' + FLAVOR) / 'runtime.json'
        bind(runtime_path)
        runtime = read(runtime_path)
        bind(Path(runtime['executable']))
        # 与真实runner一样使用通用换行读取；原文件bytes/SHA不作任何转换。
        text = (work / 'case.cfg').read_text(encoding='ascii')
        valid_locks = validate(work, copy.deepcopy(manifest), text, copy.deepcopy(runtime))
        assert valid_locks and all(sha(Path(name)) == digest for name, digest in valid_locks.items())
        executable, actual_runtime, actual_flavor, profile, route_locks = resolve(work, text)
        assert actual_flavor == FLAVOR and actual_runtime == runtime and profile == work / 'wall_temperature.csv'
        assert executable == Path(runtime['executable']).resolve() and route_locks
        assert FLAVOR in PHASE_FLAVORS
        request = read(work / 'request.json')
        assert resource_budget(request, manifest) == (30., 2)

        # 模块二：读取先前50位独立真值；不重新调用生产来流函数自造预期。
        oracle_path = ROOT / 'outputs/v174_validation/independent_anchor01/report.json'
        bind(oracle_path)
        oracle = read(oracle_path)
        assert oracle['passed'] and oracle['scope'] == 'independent-50-digit-anchored-algebra-no-production-validator'
        for name, digest in oracle['inputSHA256'].items():
            assert sha(Path(name)) == digest
            bind(name)
        c = conditions()
        field_map = {'staticTemperature_K': 'staticTemperatureK', 'pressure_Pa': 'pressurePa',
            'density_kg_m3': 'densityKgM3', 'speed_m_s': 'speedMS', 'soundSpeed_m_s': 'soundSpeedMS',
            'viscosity_Pas': 'viscosityPas', 'totalEnthalpy_Jkg': 'totalEnthalpyCpT0Jkg',
            'cp_JkgK': 'cp', 'cv_JkgK': 'cv', 'gasConstant_JkgK': 'R'}
        anchor_errors = {key: abs(c[key] / oracle['values'][truth] - 1) for key, truth in field_map.items()}
        assert max(anchor_errors.values()) < 2e-13

        def restore():
            for name, data in original.items():
                (work / name).write_bytes(data)
            assert all((work / name).read_bytes() == data for name, data in original.items())

        def contract_change(item, key, value):
            item['idealGas'][key] = value
            save(work / 'ideal_gas_contract.json', item['idealGas'])
            item['idealGasContractSHA256'] = sha(work / 'ideal_gas_contract.json')

        def initial_energy_change(item):
            data = np.loadtxt(str(work / 'initial_flow.csv'), delimiter=',', skiprows=1)
            # 错用旧OTH总焓零点的rhoE；同步登记SHA，确保不是仅被文件hash拒绝。
            old = read(ROOT / 'tests/data/naca_tn3344/contract.json')['freeStream']
            difference = old['totalEnthalpy_OTH_Jkg'] - c['totalEnthalpy_Jkg']
            data[:, 8] += data[:, 4] * difference
            header = original['initial_flow.csv'].decode('ascii').splitlines()[0]
            np.savetxt(str(work / 'initial_flow.csv'), data, delimiter=',', fmt='%.17g', comments='', header=header)
            item['initialization']['initialStateSHA256'] = sha(work / 'initial_flow.csv')

        def initial_columns_change(item):
            data = original['initial_flow.csv'].decode('ascii')
            assert '"Energy"' in data.splitlines()[0]
            (work / 'initial_flow.csv').write_text(data.replace('"Energy"', '"SpecificEnergy"', 1), encoding='ascii')
            item['initialization']['initialStateSHA256'] = sha(work / 'initial_flow.csv')

        def table_change(item):
            table_text = original['wall_temperature.csv'].decode('ascii')
            (work / 'wall_temperature.csv').write_text(table_text.replace('X_m,T_K', 'X_mm,T_C', 1), encoding='ascii')
            item['thermalProfileSHA256'] = sha(work / 'wall_temperature.csv')

        def reject(name, mutator, runtime_mutator=None, config_mutator=None):
            restore()
            item, engine = copy.deepcopy(manifest), copy.deepcopy(runtime)
            mutator(item)
            if runtime_mutator:
                runtime_mutator(engine)
            case_text = config_mutator(text) if config_mutator else text
            save(work / 'case_manifest.json', item)
            try:
                validate(work, item, case_text, engine)
            except (AssertionError, ValueError) as error:
                rejections.append(dict(name=name, rejected=True, errorType=type(error).__name__, reason=str(error)))
            else:
                raise AssertionError('真实validator未拒绝: ' + name)

        # 模块三：每例改变一个独立合同约束，运行的都是实际新validator。
        reject('old-case-kind', lambda item: item.update(caseKind='naca-tn3344-spatial-wall-laminar-v1'))
        reject('old-selected-model', lambda item: item.update(selectedFluidModel='IDEAL_AIR_VARIABLE_CP'))
        reject('wrong-runtime-flavor', lambda item: None, lambda engine: engine.update(flavor='rah-ideal-air-spatial-step-v1'))
        reject('low-precision-runtime', lambda item: None, lambda engine: engine.update(asciiSignificantDigits=16))
        reject('OTH-energy-reference', lambda item: contract_change(item, 'energyReference', 'OTH'))
        reject('wrong-gamma', lambda item: contract_change(item, 'gamma', 1.39))
        reject('implicit-standard-air-R', lambda item: contract_change(item, 'gasConstant', 287.058))
        reject('old-thermodynamic-model', lambda item: contract_change(item, 'thermodynamicModel', 'Lemmon2000-OTH'))
        reject('runtime-metadata-masquerades-as-selected-model',
               lambda item: contract_change(item, 'reusedRuntimeMetadataIsSelectedModel', True))
        # 新具名输出包只声明源更新不变；还未运行的相等性不能由元数据提前授予。
        reject('inherited-parent-model-masquerades-as-selected-model',
               lambda item: contract_change(item, 'inheritedPackageMetadataIsSelectedModel', True))
        reject('selected-runtime-model-identity-missing',
               lambda item: contract_change(item, 'selectedRuntimeMetadataIsSelectedModel', False))
        reject('wrong-native-output-contract-first-frame',
               lambda item: item['nativeIdealGasOutputContract'].update(steadyFirstCompletedUpdateLocalIteration=1))
        reject('manifest-equality-verified-before-run',
               lambda item: item['nativeIdealGasOutputContract'].update(nativeSharedFrameEqualityVerified=True))
        reject('runtime-equality-verified-before-run', lambda item: None,
               lambda engine: engine.update(nativeSharedFrameEqualityVerified=True))
        reject('nested-runtime-equality-verified-before-run', lambda item: None,
               lambda engine: engine['nativeIdealGasOutputContract'].update(nativeSharedFrameEqualityVerified=True))
        reject('inherited-OTH-selected-runtime-model', lambda item: None,
               lambda engine: engine.update(thermodynamicModel=engine['inheritedPackageThermodynamicModel']))
        reject('inherited-temperature-range-promoted-to-selected-range', lambda item: None,
               lambda engine: engine.update(temperatureRangeK=[60, 2000]))
        reject('missing-output-schedule-change-declaration', lambda item: None,
               lambda engine: engine.update(outputScheduleChanged=False))
        reject('same-physics-and-numerics-premature-claim', lambda item: None,
               lambda engine: engine.update(samePhysicsAndNumerics=True))
        reject('old-parent-runtime-use-contract',
               lambda item: item.update(runtimeUse='native-ideal-gas-v171-parent-without-output-extension'))
        reject('wrong-native-output-contract-model', lambda item: None,
               lambda engine: engine['nativeIdealGasOutputContract'].update(supportedModel='IDEAL_AIR_VARIABLE_CP'))
        reject('old-series-offset', lambda item: item['initialization'].update(priorSuccessfulUpdates=1501))
        reject('old-conservative-state-reuse', lambda item: item['initialization'].update(reusesPreviousConservativeState=True))
        reject('required-source-lock-omitted', lambda item: item['initialization']['sourceSHA256'].pop(
            str((ROOT / 'scripts/nacaConstantCpModel.py').resolve())))
        reject('old-OTH-rhoE-with-consistent-hash', initial_energy_change)
        reject('wrong-initial-energy-column-with-consistent-hash', initial_columns_change)
        reject('wrong-Tw-units-with-consistent-hash', table_change)
        reject('SST-not-laminar', lambda item: None, config_mutator=lambda value: value.replace('KIND_TURB_MODEL= NONE', 'KIND_TURB_MODEL= SST'))
        reject('axisymmetric-source-enabled', lambda item: None, config_mutator=lambda value: value.replace('AXISYMMETRIC= NO', 'AXISYMMETRIC= YES'))
        reject('old-free-stream-derived-pressure', lambda item: item['publicBenchmark']['conditions'].update(pressure_Pa=26135.37495648134))
        reject('checkpoint-stage-mismatch', lambda item: contract_change(item, 'checkpointContract',
            dict(item['idealGas']['checkpointContract'], stage='old-OTH-stage')))
        reject('undeclared-output-cap', lambda item: contract_change(item, 'checkpointContract',
            dict(item['idealGas']['checkpointContract'], maximumOutputBytes=64 * 1024**2)))
        restore()

        # 模块四：Mach向量语义直接检查真实入口函数，避免只被配置逐字比较挡住。
        lines = text.splitlines()
        index = next(i for i, line in enumerate(lines) if line.startswith('MARKER_RIEMANN='))
        tokens = [value.strip() for value in lines[index].split('=', 1)[1].strip('() ').split(',')]
        tokens[4] = '-1'
        lines[index] = 'MARKER_RIEMANN= (' + ', '.join(tokens) + ')'
        try:
            inlet_state_check('\n'.join(lines) + '\n', c, manifest['expectedVelocity'])
        except AssertionError as error:
            rejections.append(dict(name='unit-direction-not-Mach-vector', rejected=True,
                                   errorType=type(error).__name__, reason=str(error)))
        else:
            raise AssertionError('真实入口未拒绝单位Mach方向')
        for name, changes in [('undeclared-time', dict(maximumSeconds=31)),
                              ('undeclared-threads', dict(threads=3))]:
            changed = dict(request, **changes)
            try:
                resource_budget(changed, manifest)
            except AssertionError as error:
                rejections.append(dict(name=name, rejected=True, errorType=type(error).__name__, reason=str(error)))
            else:
                raise AssertionError('真实资源合同未拒绝: ' + name)
        assert len(rejections) == 35
        report.update(passed=True, positiveProductionValidationCases=1,
            actualResolverDispatchPassed=True, registeredPhaseFlavor=True,
            nativeOutputEnforcementExecuted=False, nativeCheckpointsExecuted=False,
            productionContractRefusalsTested=len(rejections), rejections=rejections,
            independentAnchorRelativeErrors=anchor_errors, independentAnchorMaximumRelativeError=max(anchor_errors.values()),
            parentRuntimeMetadataPreserved=True, originalPreparedCaseModified=False)
    except Exception as error:
        report.update(passed=False, error=type(error).__name__ + ': ' + str(error), traceback=traceback.format_exc())
        raise
    finally:
        stable = {path: Path(path).is_file() and sha(Path(path)) == digest for path, digest in locks.items()}
        if not all(stable.values()):
            report.update(passed=False, inputMutationDetected=True)
        report.update(inputSHA256=locks, inputsStable=stable, elapsedSeconds=time.monotonic() - started)
        save(out / 'input_manifest.json', locks)
        save(out / 'report.json', report)
    assert report['passed'] and all(report['inputsStable'].values())
    print(dict(passed=True, positiveProductionValidationCases=1, productionRefusals=len(rejections),
               elapsedSeconds=report['elapsedSeconds']))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--output', required=True)
    arguments = parser.parse_args()
    verify(arguments.folder, arguments.output)
