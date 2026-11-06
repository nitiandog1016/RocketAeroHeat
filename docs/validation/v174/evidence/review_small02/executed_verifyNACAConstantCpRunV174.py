"""模块用途：研究复核v174缺首帧任务的真实2/3场；不使原任务通过、不启动CFD。"""
import argparse
import copy
import csv
from pathlib import Path
import sys
import time
import traceback

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from nacaConstantCpModel import GAMMA, R, CP, CV, MODEL
from nacaSpatialWallFlux import geometry, evaluate, read_temperature_table, temperature_at, verify_applied
from su2NACAConstantCpInterface import load_state, inventory, inlet_state_check
from su2NACALaminarInterface import mesh
from auditSU2SSTState import vtk_tokens
from su2ThermalRuntime import read, sha
from validateIdealAir import save


def verify(folder, out):
    """模块一：仅接受已发生的缺0任务，绑定原回执、准确源码与真实文件。"""
    started = time.monotonic()
    folder, out = Path(folder).resolve(), Path(out).resolve()
    out.relative_to(ROOT / 'outputs')
    out.mkdir(parents=True, exist_ok=False)
    locks, rebinding, frames = {}, [], []
    report = dict(reviewCompleted=False, researchReviewPassed=False, normalTaskPassed=False,
        originalSupervisorErrorPreserved=True, originalReceiptModified=False,
        productionAuditBypassed=False, selectorUsedToAuthorizeResearch=False,
        physicalQualified=False, experimentalQualification=False, newQualifiedArea_m2=0,
        newCFDTasks=0, builds=0, frames=frames)

    def bind(path, expected=None):
        path = Path(path).resolve()
        digest = sha(path)
        if expected is not None:
            assert digest == expected, '文件身份不符: ' + str(path)
        if str(path) in locks:
            assert locks[str(path)] == digest
        locks[str(path)] = digest
        return digest

    try:
        for path in [Path(__file__).resolve()] + [ROOT / 'scripts' / name for name in (
            'nacaConstantCpModel.py', 'su2NACAConstantCpInterface.py', 'su2NACALaminarInterface.py',
            'nacaSpatialWallFlux.py', 'nacaWallFlux.py', 'auditSU2SSTState.py',
            'su2ThermalRuntime.py', 'validateIdealAir.py')]:
            bind(path)
            (out / ('executed_' + path.name)).write_bytes(path.read_bytes())
        receipt_path = folder / 'engine_process.json'
        original_receipt_hash = bind(receipt_path)
        receipt = read(receipt_path)
        assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0
        assert receipt.get('supervisorError', '').startswith('Checkpoint audit failed:')
        assert receipt.get('timedOut') is False and receipt.get('cancelled') is False
        assert receipt.get('outputBudgetExceeded') is False
        assert set(receipt['inputsStable']) == set(receipt['inputSHA256']) and all(receipt['inputsStable'].values())
        for name, digest in receipt['inputSHA256'].items():
            path = Path(name).resolve()
            if path == (ROOT / 'docs/validation/v174/PLAN.md').resolve() and sha(path) != digest:
                snapshot = folder / 'plan_at_prepare.md'
                bind(snapshot, digest)
                rebinding.append(dict(originalPath=str(path), historicalSHA256=digest,
                    exactSnapshotPath=str(snapshot), reason='PLAN追加后保留准确执行快照，原回执不改'))
            else:
                bind(path, digest)
        manifest = read(folder / 'case_manifest.json')
        assert manifest['caseKind'] == 'naca-tn3344-constant-cp-spatial-wall-v1'
        assert manifest['specification']['updates'] == 4 and manifest['initialization']['priorSuccessfulUpdates'] == 0
        contract = read(folder / 'ideal_gas_contract.json')
        assert manifest['idealGas'] == contract and contract['thermodynamicModel'] == MODEL
        assert contract['energyReference'] == 'e=cv*T; h=cp*T; no-OTH-offset'
        cp = receipt['checkpoints']
        assert [pair['localIteration'] for pair in cp['complete']] == [2, 3]
        assert cp['missingScheduledLocalIterations'] == [0] and cp['normalScheduleComplete'] is False
        assert cp['lastFiniteHistoryLocalIteration'] == 3 and cp['historyRows'] == 4
        assert not (folder / 'restart_flow_000000.csv').exists() and not (folder / 'surface_flow_000000.vtk').exists()
        report.update(actualNativeExitCode=0, actualSupervisorError=receipt['supervisorError'],
            nativeProcessExitedZero=True, actualNativeRecordedNormalScheduleComplete=False,
            missingScheduledLocalIterations=[0], completeSavedLocalIterations=[2, 3],
            originalSavedFirstFieldAvailable=False, lastFiniteLocalIteration=3,
            originalRecordedCheckpointNormalTerminal=cp['normalTerminal'],
            recordedCheckpointWasCreatedBeforeSupervisorError=True,
            actualGeneratedOutputBytes=receipt['generatedOutputBytes'])
        bind(folder / 'history.csv', receipt['outputSHA256']['history.csv'])
        with (folder / 'history.csv').open(encoding='utf8') as stream:
            history = {int(float(item['Inner_Iter'])): item for item in
                [{key.strip(): value for key, value in row.items()}
                 for row in csv.DictReader(stream, skipinitialspace=True)]}
        assert sorted(history) == [0, 1, 2, 3]
        points, cells, markers = mesh(folder / 'mesh.su2')
        meta = read(folder / 'geometry.json')
        with np.load(str(folder / 'mesh_arrays.npz')) as arrays:
            assert np.array_equal(arrays['points'], points)
            labels = arrays['nodeIndices']
        geom = geometry(points, cells, markers, np.asarray(meta['s'])[labels[:, 1]])
        table = read_temperature_table(folder / 'wall_temperature.csv')
        tw = temperature_at(points[geom['ids'], 0], table)
        applied_path = folder / 'rah_wall_temperature_applied.csv'
        bind(applied_path, receipt['outputSHA256'][applied_path.name])
        applied = verify_applied(applied_path, points, geom, tw)
        assert applied['nodes'] == 61
        conditions = manifest['publicBenchmark']['conditions']
        inlet = inlet_state_check((folder / 'case.cfg').read_text(), conditions, manifest['expectedVelocity'])

        # 模块二：已保存2/3场的研究审核，不能追补原未保存0或改变任务失败身份。
        for pair in cp['complete']:
            step = pair['localIteration']
            volume, surface_path = folder / pair['volume'], folder / pair['surface']
            bind(volume, pair['volumeSHA256'])
            bind(surface_path, pair['surfaceSHA256'])
            assert pair['volumeSHA256'] == receipt['outputSHA256'][volume.name]
            assert pair['surfaceSHA256'] == receipt['outputSHA256'][surface_path.name]
            with volume.open(encoding='utf8') as stream:
                columns = next(csv.reader(stream))
            assert len(columns) == 27
            data, values = load_state(volume)
            assert len(data) == 549 and np.max(abs(data[:, 1:4] - points)) < 2e-13
            rho = values['Density']
            velocity = np.column_stack([values['Momentum_' + axis] / rho for axis in 'xyz'])
            temperature = values['Temperature']
            internal = values['Energy'] / rho - .5 * np.sum(velocity * velocity, axis=1)
            energy = float(np.max(abs(internal - CV * temperature) / np.maximum(1., CV * temperature)))
            eos = float(np.max(abs(rho * R * temperature / values['Pressure'] - 1)))
            viscosity = 1.716e-5 * (temperature / 273.15) ** 1.5 * 383.55 / (temperature + 110.4)
            mu = float(np.max(abs(values['Laminar_Viscosity'] / viscosity - 1)))
            flux = evaluate(geom, temperature, CP * values['Laminar_Viscosity'] / .72, tw)
            ids, q = geom['ids'], flux['qIntoWall']
            heat = float(np.max(abs(q - values['Heat_Flux'][ids]) / np.maximum(1., abs(q))))
            power = abs(flux['powerIntoWall'] - float(history[step]['HF'])) / max(1., abs(flux['powerIntoWall']))
            xyz, surface = vtk_tokens(surface_path)
            distance, matched = cKDTree(points[ids]).query(xyz)
            assert max(distance) < 1e-6 and len(np.unique(matched)) == len(ids) == len(xyz)
            assert all(surface[key].shape == (len(ids), 1) for key in ('Pressure', 'Temperature', 'Heat_Flux'))
            surface_differences = {key: float(np.max(abs(surface[key][:, 0] - values[key][ids[matched]])
                / np.maximum(1., abs(values[key][ids[matched]])))) for key in ('Pressure', 'Temperature', 'Heat_Flux')}
            cf = np.column_stack([values['Skin_Friction_Coefficient_' + axis][ids] for axis in 'xyz'])
            shear = np.linalg.norm(cf, axis=1) * .5 * conditions['density_kg_m3'] * conditions['speed_m_s'] ** 2
            yplus = geom['euclideanDistance'] * np.sqrt(rho[ids] * shear) / values['Laminar_Viscosity'][ids]
            passed = bool(max(energy, eos, mu) < 2e-11 and heat < 2e-10 and power < 2e-9)
            frames.append(dict(localIteration=step, seriesSuccessfulUpdates=step + 1, nativeColumns=columns,
                availableFieldThermalConsistencyPassed=passed, normalTaskPassed=False,
                energyMaximumRelativeError=energy, eosMaximumRelativeError=eos,
                viscosityMaximumRelativeError=mu, nativeHeatMaximumScaledError=heat,
                powerScaledDifference=power, nativeNodalPowerIntoWallW=flux['powerIntoWall'],
                nativeHistoryPowerW=float(history[step]['HF']), surfaceToVolumeScaledDifferences=surface_differences,
                surfaceMappingMaximumDistanceM=float(max(distance)),
                temperatureRangeK=[float(temperature.min()), float(temperature.max())],
                developmentModelDomainPassed=bool(np.all((temperature >= 150.) & (temperature <= 360.))),
                givenTwMaximumFluidDifferenceK=float(np.max(abs(tw - temperature[ids]))),
                independentFirstNeighborYPlus99=float(np.percentile(yplus, 99)),
                originalResidualGatePassed=bool(max(float(history[step][key]) for key in ('rms[Rho]', 'rms[RhoE]')) <= -8),
                physicalQualified=False, experimentalQualification=False, newQualifiedArea_m2=0))
        assert all(frame['availableFieldThermalConsistencyPassed'] for frame in frames)

        # 模块三：仅在独占副本中模拟少算/少保存；把实际2帧重命名为0的动作明确标合成。
        synthetic = []
        for label, saved_steps, history_rows in [('synthetic_only0', [0], 1), ('synthetic_without_final3', [0, 2], 3)]:
            fixture = out / label
            fixture.mkdir()
            for name in ('case_manifest.json', 'mesh.su2'):
                (fixture / name).write_bytes((folder / name).read_bytes())
            raw_history = (folder / 'history.csv').read_bytes().splitlines(keepends=True)
            (fixture / 'history.csv').write_bytes(b''.join(raw_history[:history_rows + 1]))
            synthetic_hashes = {'history.csv': sha(fixture / 'history.csv')}
            for saved_step in saved_steps:
                for prefix, extension in [('restart_flow', 'csv'), ('surface_flow', 'vtk')]:
                    destination = fixture / ('%s_%06d.%s' % (prefix, saved_step, extension))
                    destination.write_bytes((folder / ('%s_000002.%s' % (prefix, extension))).read_bytes())
                    synthetic_hashes[destination.name] = sha(destination)
            synthetic_receipt = dict(exitCode=0, timedOut=False, cancelled=False,
                outputBudgetExceeded=False, outputSHA256=synthetic_hashes)
            save(fixture / 'synthetic_receipt.json', synthetic_receipt)
            observed = inventory(fixture, synthetic_receipt)
            assert observed['normalTerminal'] is True and observed['normalScheduleComplete'] is False
            assert observed['normalFinalPairPresent'] is False and 3 in observed['missingScheduledLocalIterations']
            assert observed['expectedLocalIterations'] == [0, 2, 3]
            save(fixture / 'inventory.json', observed)
            synthetic.append(dict(name=label, syntheticOnly=True, realNativeTask=False,
                relabelledActualStep2AsFixtureStep0=True, noPhysicalTrajectoryClaim=True,
                normalFinalPairRejected=True, normalCompleteScheduleRejected=True,
                missingScheduledLocalIterations=observed['missingScheduledLocalIterations']))
        assert sha(receipt_path) == original_receipt_hash
        report.update(reviewCompleted=True, researchReviewPassed=True,
            availableSavedFieldsConsistencyPassed=True, normalTaskPassed=False,
            originalMissingFirstFrameStillUnresolved=True, all27ColumnsFinite=True,
            appliedTemperatureCheck=applied, inletStateCheck=inlet, syntheticInventoryRefusals=synthetic,
            historicalInputRebinding=rebinding, selectedThermodynamicModel=MODEL,
            scope='read-only research of authentic saved2/3; not production audit acceptance',
            developmentDomainDoesNotGrantQualification=True, noNewNativeTaskRequiredForThisReview=True)
    except Exception as error:
        report.update(researchReviewPassed=False, error=type(error).__name__ + ': ' + str(error), traceback=traceback.format_exc())
        raise
    finally:
        stable = {name: Path(name).is_file() and sha(Path(name)) == digest for name, digest in locks.items()}
        if not all(stable.values()):
            report.update(researchReviewPassed=False, inputMutationDetected=True)
        report.update(inputSHA256=locks, inputsStable=stable, elapsedSeconds=time.monotonic() - started)
        save(out / 'input_manifest.json', locks)
        save(out / 'report.json', report)
    assert report['researchReviewPassed'] and all(report['inputsStable'].values())
    print(dict(researchReviewPassed=True, normalTaskPassed=False, missingScheduledLocalIterations=[0],
               frames=[2, 3], elapsedSeconds=report['elapsedSeconds']))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--output', required=True)
    arguments = parser.parse_args()
    verify(arguments.folder, arguments.output)
