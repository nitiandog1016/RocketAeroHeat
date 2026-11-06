"""模块用途：独立小副本验证实际候选资源、守卫与保存终态；不启动实际500步或CFD。"""
import argparse
import copy
import csv
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from prepareNACAConstantCpCase import specification, checkpoint_plan, PLAN, ACTUAL_RESOURCE_PROFILE, GUARD_REPORT
from su2NACAConstantCpInterface import validate_development_contract, resource_budget, inventory
from nacaConstantCpDomainGuard import domain_guard_contract
from nacaSavedFrameStudy import select_saved_frames
from su2ThermalRuntime import read, sha
from validateIdealAir import save


def verify(folder, out):
    """模块一：冻结小例真实场及当前准确源码；所有实际规格/历史组合均明确标为合成。"""
    started = time.monotonic()
    folder, out = Path(folder).resolve(), Path(out).resolve()
    out.relative_to(ROOT / 'outputs')
    out.mkdir(parents=True, exist_ok=False)
    locks, refusals, positives = {}, [], []
    report = dict(passed=False, syntheticOnly=True, actual500NativeExecuted=False,
        actual190334MeshFullyValidated=False, originalFilesModified=False, newCFDTasks=0,
        builds=0, physicalQualified=False, experimentalQualification=False, newQualifiedArea_m2=0)

    def bind(path, expected=None):
        path = Path(path).resolve()
        digest = sha(path)
        assert expected is None or digest == expected, '输入身份变化: ' + str(path)
        assert str(path) not in locks or locks[str(path)] == digest
        locks[str(path)] = digest
        return digest

    def refuse(label, operation):
        try:
            operation()
        except (AssertionError, KeyError, ValueError) as error:
            refusals.append(dict(name=label, rejected=True, errorType=type(error).__name__, reason=str(error)))
        else:
            raise AssertionError('实际生产接口未拒绝: ' + label)

    try:
        for path in [Path(__file__).resolve()] + [ROOT / 'scripts' / name for name in
            ('prepareNACAConstantCpCase.py', 'su2NACAConstantCpInterface.py', 'nacaConstantCpDomainGuard.py',
             'nacaSavedFrameStudy.py', 'runEmbeddedSU2.py', 'su2IdealAirCheckpoints.py',
             'su2NACALaminarInterface.py', 'su2ThermalRuntime.py', 'auditSU2SSTState.py', 'validateIdealAir.py')]:
            bind(path)
            (out / ('executed_' + path.name)).write_bytes(path.read_bytes())
        bind(PLAN)
        (out / 'plan_at_check.md').write_bytes(PLAN.read_bytes())
        bind(GUARD_REPORT)
        guard_proof = read(GUARD_REPORT)
        assert guard_proof['passed'] and all(guard_proof['inputsStable'].values())
        original_receipt = read(folder / 'engine_process.json')
        bind(folder / 'engine_process.json')
        assert original_receipt['exitCode'] == 0 and not original_receipt.get('supervisorError')
        for name in ('case_manifest.json', 'mesh.su2'):
            expected = original_receipt['inputSHA256'][str(folder / name)]
            bind(folder / name, expected)
        for name in ('restart_flow_000000.csv', 'surface_flow_000000.vtk'):
            bind(folder / name, original_receipt['outputSHA256'][name])

        # 模块二：准确资源及低成本实际开发合同；不以小网格假称通过完整几何validator。
        manifest = copy.deepcopy(read(folder / 'case_manifest.json'))
        spec = specification('validation-profile')
        assert [spec[key] for key in ('updates', 'every', 'maximumSeconds', 'threads', 'maximumOutputBytes')] == [500, 100, 600, 2, 1536 * 1024**2]
        assert spec['actualValidationMesh'] is True
        manifest['specification'] = spec
        manifest['nodes'], manifest['elements'] = 190334, 179568
        manifest['publicBenchmark']['fixture'] = False
        manifest['developmentPrerequisites'] = dict(coarseFixtureDevelopmentDomainPassed=False)
        manifest['idealGas']['savedFieldDomainGuard'] = domain_guard_contract()
        manifest['idealGas']['checkpointContract'] = checkpoint_plan(spec, 190334, 3226, 3096)
        manifest['planSHA256'] = sha(PLAN)
        request = dict(maximumSeconds=600, threads=2, resourceProfile=ACTUAL_RESOURCE_PROFILE)
        assert validate_development_contract(manifest) == spec
        assert resource_budget(request, manifest) == (600., 2)
        assert manifest['idealGas']['checkpointContract']['estimatedMaximumBytes'] == 1176607488
        save(out / 'synthetic_actual_spec_manifest.json', manifest)
        save(out / 'synthetic_actual_resource_request.json', request)
        positives.append(dict(name='actual-spec-and-resource-structure', productionFunctions=['validate_development_contract', 'resource_budget']))

        def changed_manifest(mutator):
            value = copy.deepcopy(manifest)
            mutator(value)
            return value

        mutations = [
            ('missing-guard', lambda value: value['idealGas'].pop('savedFieldDomainGuard')),
            ('widened-guard-temperature-range', lambda value: value['idealGas']['savedFieldDomainGuard'].update(temperatureRangeK=[150., 400.])),
            ('widened-model-temperature-range', lambda value: value['idealGas'].update(developmentTemperatureRangeK=[60., 2000.])),
            ('disabled-live-history-protection', lambda value: value['idealGas']['savedFieldDomainGuard'].update(historyStrictlyAfterForLive=False)),
            ('increased-helper-budget', lambda value: value['idealGas']['savedFieldDomainGuard'].update(helperMaximumSeconds=16.)),
            ('nonfresh-prior', lambda value: value['initialization'].update(priorSuccessfulUpdates=4)),
            ('reused-conservative-state', lambda value: value['initialization'].update(reusesPreviousConservativeState=True)),
            ('old-checkpoint-metadata', lambda value: value['initialization'].update(sourceCheckpointIdentity={'old': True})),
            ('coarse-domain-failure-erased', lambda value: value['developmentPrerequisites'].update(coarseFixtureDevelopmentDomainPassed=True)),
            ('small-mesh-masquerades-as-actual', lambda value: value.update(nodes=549)),
            ('unregistered-extra-update', lambda value: value['specification'].update(updates=501)),
        ]
        for name, mutator in mutations:
            changed = changed_manifest(mutator)
            refuse(name, lambda changed=changed: validate_development_contract(changed))
        for name, updates in [('missing-resource-profile', {'resourceProfile': None}),
                ('extra-native-second', {'maximumSeconds': 601}), ('extra-thread', {'threads': 3})]:
            value = dict(request, **updates)
            refuse(name, lambda value=value: resource_budget(value, manifest))
        value = changed_manifest(lambda item: item['idealGas']['checkpointContract'].update(maximumOutputBytes=1537 * 1024**2))
        refuse('extra-output-capacity', lambda: resource_budget(request, value))

        # 模块三：准确0帧只是文件夹具；改名/合成history从不代表真实演化或收敛。
        def fixture(label, last, saved, rho=-3., energy=2., flags=None, missing_history=None, broken_surface=False):
            work = out / label
            work.mkdir()
            (work / 'mesh.su2').write_bytes((folder / 'mesh.su2').read_bytes())
            save(work / 'case_manifest.json', manifest)
            with (work / 'history.csv').open('w', encoding='ascii', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(['Inner_Iter', 'rms[Rho]', 'rms[RhoE]'])
                for index in range(last + 1):
                    if index != missing_history:
                        writer.writerow([index, rho if index == last else -3., energy if index == last else 2.])
            for step in saved:
                for prefix, extension in [('restart_flow', 'csv'), ('surface_flow', 'vtk')]:
                    data = (folder / ('%s_000000.%s' % (prefix, extension))).read_bytes()
                    if broken_surface and step == saved[-1] and prefix == 'surface_flow':
                        data = data[:len(data) // 2]
                    (work / ('%s_%06d.%s' % (prefix, step, extension))).write_bytes(data)
            inputs = {str(work / name): sha(work / name) for name in ('mesh.su2', 'case_manifest.json')}
            outputs = {path.name: sha(path) for path in work.iterdir() if path.name == 'history.csv' or path.name.startswith(('restart_flow_', 'surface_flow_'))}
            receipt = dict(status='terminal', exitCode=0, timedOut=False, cancelled=False, outputBudgetExceeded=False,
                modelDomainExceeded=False, inputSHA256=inputs, inputsStable={key: True for key in inputs}, outputSHA256=outputs,
                syntheticOnly=True, actualNativeTask=False, originalStep0Relabelled=True, physicalHistoryFabricatedForContractTest=True)
            if flags:
                receipt.update(flags)
            save(work / 'synthetic_receipt.json', receipt)
            return work, receipt

        for label, last, saved, density, energy in [('synthetic_full500', 499, [0, 100, 200, 300, 400, 499], -3., 2.),
                                                   ('synthetic_early10', 10, [0, 10], -8., -8.)]:
            work, receipt = fixture(label, last, saved, density, energy)
            checked, identity = select_saved_frames(work, receipt)
            assert checked['normalScheduleComplete'] and identity['normalRunCompleted']
            assert checked['earlyConvergenceAccepted'] == (last == 10)
            save(work / 'inventory.json', checked)
            positives.append(dict(name=label, syntheticOnly=True, normalScheduleAccepted=True))
        cases = [
            ('synthetic_false-early-flag', 3, [0, 3], -3., 2., {'earlyConvergenceAccepted': True}, None, False),
            ('synthetic_density-only-converged', 10, [0, 10], -8., -7.9, None, None, False),
            ('synthetic_energy-only-converged', 10, [0, 10], -7.9, -8., None, None, False),
            ('synthetic_before-startiter', 9, [0, 9], -8., -8., None, None, False),
            ('synthetic_missing-final', 10, [0], -8., -8., None, None, False),
            ('synthetic_missing-first', 10, [10], -8., -8., None, None, False),
            ('synthetic_incomplete-final-surface', 10, [0, 10], -8., -8., None, None, True),
            ('synthetic_history-gap', 10, [0, 10], -8., -8., None, 5, False),
        ]
        for arguments in cases:
            work, receipt = fixture(*arguments)
            checked = inventory(work, receipt)
            assert not checked['normalScheduleComplete']
            save(work / 'inventory.json', checked)
            refuse(arguments[0], lambda work=work, receipt=receipt: select_saved_frames(work, receipt))
        work, receipt = fixture('synthetic_terminal-only-domain-discovery', 10, [0, 10], -8., -8.,
            {'modelDomainExceeded': True, 'savedFieldDomainGuard': {'exceeded': True, 'terminalOnlyDetection': True}})
        refuse('default-rejects-domain-stop-even-exit0', lambda: select_saved_frames(work, receipt))
        checked, identity = select_saved_frames(work, receipt, True)
        assert identity['normalRunCompleted'] is False and identity['terminalFlags']['modelDomainExceeded'] is True
        assert identity['terminalFlags']['cancelled'] is False and identity['actualExitCode'] == 0
        assert not identity['unsavedTerminalStateKnown'] and not checked['normalScheduleComplete']
        save(work / 'selection_identity.json', identity)
        positives.append(dict(name='explicit-domain-stop-saved-only', syntheticOnly=True, userCancellation=False))
        broken = dict(receipt, supervisorError='Saved field domain guard failed: synthetic parser failure')
        refuse('guard-error-not-valid-saved-only-entry', lambda: select_saved_frames(work, broken, True))
        broken = dict(receipt, modelDomainExceeded='true')
        refuse('nonboolean-domain-stop', lambda: select_saved_frames(work, broken, True))
        work, receipt = fixture('synthetic_no-complete-frame', 10, [], -8., -8., {'modelDomainExceeded': True})
        refuse('abnormal-without-complete-frame', lambda: select_saved_frames(work, receipt, True))
        report.update(passed=True, actualSpec=spec, conservativeMaximumOutputBytes=1176607488,
            positiveCases=positives, productionRefusals=len(refusals), refusals=refusals,
            actualNormalEarlyConvergenceNotClaimed=True, actualGuardTerminationNotClaimed=True,
            syntheticGeometryNote='549-node copied geometry with structural actual-spec manifest; full prepare validator intentionally not invoked')
    except Exception as error:
        report.update(passed=False, error=type(error).__name__ + ': ' + str(error), traceback=traceback.format_exc())
        raise
    finally:
        stable = {path: Path(path).is_file() and sha(Path(path)) == value for path, value in locks.items()}
        if not all(stable.values()):
            report.update(passed=False, inputMutationDetected=True)
        report.update(inputSHA256=locks, inputsStable=stable, elapsedSeconds=time.monotonic() - started)
        save(out / 'input_manifest.json', locks)
        save(out / 'report.json', report)
    assert report['passed'] and all(report['inputsStable'].values())
    print(dict(passed=True, positiveCases=len(positives), productionRefusals=len(refusals), elapsedSeconds=report['elapsedSeconds']))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--output', required=True)
    arguments = parser.parse_args()
    verify(arguments.folder, arguments.output)
