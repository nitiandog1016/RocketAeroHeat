"""模块用途：独立审核定比热NACA保存态的能量、弱壁热流、原段功率与开发域。"""
import argparse
import csv
from pathlib import Path
import time
import traceback

import numpy as np

from nacaConstantCpModel import MODEL, GAMMA, R, CP, CV, PR, state
from nacaMeasurementSegments import segments, area
from nacaSavedFrameStudy import select_saved_frames
from nacaSpatialWallFlux import (geometry, evaluate, read_temperature_table,
    temperature_at, segment_integrals, verify_applied)
from prepareNACAConstantCpCase import KIND, FLAVOR, PLAN
from su2NACAConstantCpInterface import load_state, inlet_state_check
from su2NACALaminarInterface import mesh
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save


def column_identity(meta, points, arrays_path):
    """模块一：核对原全局节点与法向列；这里只取原s标签，不求模型外缘。"""
    shape = (len(meta['offsets']), len(meta['s']), len(meta['azimuth']))
    with np.load(str(arrays_path)) as arrays:
        assert np.array_equal(arrays['points'], points)
        labels = arrays['nodeIndices'].copy()
        indices = arrays['columnNodeIDs'].copy()
    assert labels.shape == (len(points), 3) and indices.shape == shape
    assert np.issubdtype(labels.dtype, np.integer) and np.issubdtype(indices.dtype, np.integer)
    assert indices.min() >= 0 and indices.max() < len(points) and len(np.unique(indices)) == len(points)
    for node, (j, i, k) in enumerate(labels):
        assert 0 <= j < shape[0] and 0 <= i < shape[1] and 0 <= k < shape[2]
        assert indices[j, i, k] == node
    return labels, indices


def history_rows(path):
    """模块二：保留有限历史及坏行身份；只从同一保存步取得原生有符号功率。"""
    accepted, rejected = {}, []
    with Path(path).open(encoding='utf8') as stream:
        reader = csv.reader(stream, skipinitialspace=True)
        keys = [key.strip() for key in next(reader)]
        assert len(keys) == len(set(keys)) and {'Inner_Iter', 'HF', 'rms[Rho]', 'rms[RhoE]'}.issubset(keys)
        for number, row in enumerate(reader, 2):
            try:
                values = np.asarray(row, float)
                assert len(values) == len(keys) and np.isfinite(values).all()
                item = dict(zip(keys, values))
                step = item['Inner_Iter']
                assert step == int(step) and int(step) >= 0 and int(step) not in accepted
                accepted[int(step)] = item
            except (ValueError, AssertionError, OverflowError) as error:
                rejected.append(dict(line=number, reason=type(error).__name__ + ': ' + str(error)))
    return accepted, rejected


def original_segment_powers(meta, geom, flux, rows, summary):
    """模块三：20个重叠原测量段只作真实有符号功率映射，不偷用旧Lemmon归约。"""
    bounds = np.asarray([[float(row[key]) for key in ('s_left_m', 's_right_m')] for row in rows])
    base = dict(selectedThermodynamicModel=MODEL, experimentalQualification=False,
        newQualifiedArea_m2=0, comparisonReady=False, modelAwareEdgeAnalysisImplemented=False,
        edgeAndReductionStatus='not-implemented-for-selected-constant-cp-model',
        freestreamFallbackUsed=False, originalTransportTableEquivalent=False,
        all20CenterIdentitiesClosed=summary['strict_all20_local_abscissa_identity_closed'],
        overlappingSegmentPowersAreNotTotal=True, measuredParameterIsDimensionalHeatFlux=False)
    if bounds.min() < np.min(geom['patchS']) - 1e-14 or bounds.max() > np.max(geom['patchS']) + 1e-14:
        return dict(base, status='wall-does-not-cover-original-measurement-segments', segments=[],
                    nativeNodalPowerDecomposition=False)
    actual = segment_integrals(geom, flux, bounds, meta['azimuth'][-1] - meta['azimuth'][0])
    for row, item in zip(rows, actual):
        left, right = item['boundsM']
        nominal = float(area(left, right))
        assert abs(nominal / float(row['area_full_circumference_m2']) - 1) < 2e-12
        item.update(row=int(row['heat_row_1based']), nominalFullAreaM2=nominal,
            meanQByNominalAreaWm2=item['symmetryScaledPowerIntoWallFullW'] / nominal,
            discreteVsNominalAreaRelativeDifference=item['discreteGeometricAreaSectorM2']
                * item['fullRevolutionScale'] / nominal - 1,
            strictLocalAbscissaIdentityClosed=row['strict_local_abscissa_identity_closed'] == 'True',
            measuredParameter=float(row['original_Nu_s_over_sqrt_Re_s']),
            selectedThermodynamicModel=MODEL, edgeAndReductionStatus=base['edgeAndReductionStatus'],
            experimentalQualification=False, comparisonReady=False)
    return dict(base, status='discrete-segment-power-only-constant-cp', segments=actual,
                nativeNodalPowerDecomposition=True)


def write_segment_csv(path, measurement):
    """模块四：CSV同时给出扇区功率、对称放大、两类面积与明确未完成的归约状态。"""
    fields = ['row', 's_left_m', 's_right_m', 'powerIntoWallSectorW',
        'symmetryScaledPowerIntoWallFullW', 'discreteGeometricAreaSectorM2',
        'nativeProjectedAreaSectorM2', 'nominalFullAreaM2', 'meanQByNominalAreaWm2',
        'meanQByDiscreteGeometricAreaWm2', 'meanQByNativeProjectedAreaWm2',
        'measuredDimensionlessNuOverSqrtRe', 'strictLocalAbscissaIdentityClosed',
        'method', 'thermodynamicModel', 'edgeAndReductionStatus', 'experimentalQualification']
    with Path(path).open('w', encoding='utf8', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        for item in measurement['segments']:
            writer.writerow([item['row'], *item['boundsM'], item['powerIntoWallSectorW'],
                item['symmetryScaledPowerIntoWallFullW'], item['discreteGeometricAreaSectorM2'],
                item['nativeProjectedAreaSectorM2'], item['nominalFullAreaM2'],
                item['meanQByNominalAreaWm2'], item['meanQByDiscreteGeometricAreaWm2'],
                item['meanQByNativeProjectedAreaWm2'], item['measuredParameter'],
                item['strictLocalAbscissaIdentityClosed'], item['method'], MODEL,
                item['edgeAndReductionStatus'], False])


def audit(folder, out, saved_frames_only=False):
    """模块五：严格同阶段热审核与模型开发域分别判断；全部失败保存准确报告。"""
    started = time.monotonic()
    folder, out = Path(folder).resolve(), Path(out).resolve()
    out.relative_to(ROOT / 'outputs')
    out.mkdir(parents=True, exist_ok=False)
    locks, frames = {}, []
    report = dict(reviewCompleted=False, interfacePassed=False, savedFrameInterfacePassed=False,
        selectedThermodynamicModel=MODEL, selectedFluidModel='IDEAL_GAS',
        energyReference='e=cv*T; h=cp*T; no-OTH-offset',
        physicalQualified=False, experimentalQualification=False, newQualifiedArea_m2=0,
        normalRunCompleted=False, newCFDTasks=0, frames=frames)

    def bind(path, expected=None):
        path = Path(path).resolve()
        digest = sha(path)
        if expected is not None:
            assert digest == expected, '输入或输出身份改变: ' + str(path)
        if str(path) in locks:
            assert locks[str(path)] == digest, '审核期间源变化: ' + str(path)
        locks[str(path)] = digest
        return digest

    def bind_recorded(path, receipt):
        path = Path(path).resolve()
        found = [digest for name, digest in receipt['inputSHA256'].items() if Path(name).resolve() == path]
        assert found and len(set(found)) == 1, '缺本任务准确输入路径绑定: ' + str(path)
        return bind(path, found[0])

    try:
        names = ('auditNACAConstantCp.py', 'nacaConstantCpModel.py', 'prepareNACAConstantCpCase.py',
            'su2NACAConstantCpInterface.py', 'nacaSavedFrameStudy.py', 'nacaSpatialWallFlux.py',
            'nacaWallFlux.py', 'nacaMeasurementSegments.py', 'naca3344Model.py',
            'su2NACALaminarInterface.py', 'su2SSTInterface.py', 'su2IdealAirCheckpoints.py',
            'su2IdealAirInterface.py', 'su2ThermalRuntime.py', 'auditSU2SSTState.py',
            'analyzeSU2DualVolumes.py', 'validateIdealAir.py')
        for name in names:
            path = ROOT / 'scripts' / name
            bind(path)
            (out / ('executed_' + name)).write_bytes(path.read_bytes())
        receipt_path = folder / 'engine_process.json'
        bind(receipt_path)
        receipt = read(receipt_path)
        manifest_path = folder / 'case_manifest.json'
        bind_recorded(manifest_path, receipt)
        manifest = read(manifest_path)
        assert manifest['caseKind'] == KIND and manifest['requiredRuntimeFlavor'] == FLAVOR
        assert manifest['selectedFluidModel'] == 'IDEAL_GAS' and not manifest['physicalQualified']
        assert 'idealAir' not in manifest and 'sst' not in manifest
        assert not (folder / 'ideal_air_contract.json').exists() and not (folder / 'sst_contract.json').exists()
        contract_path = folder / 'ideal_gas_contract.json'
        bind_recorded(contract_path, receipt)
        bind(contract_path, manifest['idealGasContractSHA256'])
        contract = read(contract_path)
        assert contract == manifest['idealGas'] and contract['schema'] == 'RAH.NACA.NativeIdealGas.v1'
        assert contract['selectedFluidModel'] == 'IDEAL_GAS' and contract['thermodynamicModel'] == MODEL
        assert contract['gamma'] == GAMMA and contract['gasConstant'] == R
        assert contract['energyReference'] == report['energyReference'] and not contract['turbulent']
        assert contract['developmentTemperatureRangeK'] == [150., 360.]
        assert not contract['physicalQualified'] and contract['newQualifiedArea_m2'] == 0
        init = manifest['initialization']
        assert init['kind'] == KIND and init['priorSuccessfulUpdates'] == 0
        assert init['reusesPreviousConservativeState'] is False
        assert init['energyReference'] == report['energyReference']
        bind_recorded(folder / 'initial_flow.csv', receipt)
        bind(folder / 'initial_flow.csv', init['initialStateSHA256'])
        snapshot = folder / 'plan_at_prepare.md'
        if snapshot.exists():
            bind(snapshot, manifest['planSHA256'])
        else:
            bind(PLAN, manifest['planSHA256'])

        # 真实终态选择器自行重验完整双文件；异常只在显式模式下允许研究保存帧。
        checkpoints, identity = select_saved_frames(folder, receipt, saved_frames_only)
        assert receipt.get('runtimeFlavor') == FLAVOR and 'checkpoints' in receipt
        assert checkpoints == receipt['checkpoints'], '原生与当前独立清点不一致'
        assert 'generatedOutputBytes' in receipt and 'outputBudgetExceeded' in receipt
        for name, key in [('mesh.su2', 'meshSHA256'), ('geometry.json', 'geometrySHA256'),
                          ('mesh_arrays.npz', 'meshArraysSHA256'), ('case.cfg', 'configurationSHA256')]:
            bind_recorded(folder / name, receipt)
            bind(folder / name, contract[key])
        bind(folder / 'history.csv', receipt['outputSHA256']['history.csv'])
        points, cells, markers = mesh(folder / 'mesh.su2')
        meta = read(folder / 'geometry.json')
        labels, indices = column_identity(meta, points, folder / 'mesh_arrays.npz')
        node_s = np.asarray(meta['s'])[labels[:, 1]]
        geom = geometry(points, cells, markers, node_s)
        del indices
        from su2IdealAirInterface import settings
        text = (folder / 'case.cfg').read_text(encoding='utf8')
        cfg = settings(text)
        assert cfg['FLUID_MODEL'] == 'IDEAL_GAS' and float(cfg['GAMMA_VALUE']) == GAMMA
        assert float(cfg['GAS_CONSTANT']) == R and float(cfg['PRANDTL_LAM']) == PR
        assert cfg['SOLVER'] == 'NAVIER_STOKES' and cfg['KIND_TURB_MODEL'] == 'NONE' and cfg['AXISYMMETRIC'] == 'NO'
        assert cfg['VISCOSITY_MODEL'] == 'SUTHERLAND' and cfg['CONDUCTIVITY_MODEL'] == 'CONSTANT_PRANDTL'
        assert float(cfg['MU_REF']) == 1.716e-5 and float(cfg['MU_T_REF']) == 273.15
        assert float(cfg['SUTHERLAND_CONSTANT']) == 110.4
        conditions = manifest['publicBenchmark']['conditions']
        assert 'totalEnthalpy_OTH_Jkg' not in conditions
        inlet = inlet_state_check(text, conditions, manifest['expectedVelocity'])
        for name, digest in conditions['sourceSHA256'].items():
            bind(name, digest)
        table_path = folder / 'wall_temperature.csv'
        bind_recorded(table_path, receipt)
        bind(table_path, manifest['thermalProfileSHA256'])
        table = read_temperature_table(table_path)
        given_tw = temperature_at(points[geom['ids'], 0], table)
        assert np.all((given_tw >= 150.) & (given_tw <= 360.))
        applied_path = folder / 'rah_wall_temperature_applied.csv'
        assert applied_path.name in receipt['outputSHA256'], '原生实际Tw回执缺输出身份'
        bind(applied_path, receipt['outputSHA256'][applied_path.name])
        applied = verify_applied(applied_path, points, geom, given_tw)
        history, rejected_history = history_rows(folder / 'history.csv')
        if identity['normalRunCompleted']:
            assert not rejected_history, '正常终态存在不完整历史行'
        rows, summary = segments()
        data_folder = ROOT / 'tests/data/naca_tn3344'
        bind(data_folder / 'segments_v170/summary.json')
        for name, digest in summary['generated_csv_sha256'].items():
            bind(data_folder / 'segments_v170' / name, digest)
        for name, digest in summary['old_inputs_sha256_before_and_after_identical'].items():
            if name.startswith('tests/'):
                bind(ROOT / name, digest)

        for pair in checkpoints['complete']:
            path = folder / pair['volume']
            digest = bind(path, pair['volumeSHA256'])
            bind(folder / pair['surface'], pair['surfaceSHA256'])
            data, values = load_state(path)
            assert len(data) == len(points) and np.max(abs(data[:, 1:4] - points)) < 2e-13
            rho = values['Density']
            velocity = np.column_stack([values['Momentum_' + axis] / rho for axis in 'xyz'])
            temperature = values['Temperature']
            thermo = state('PT', values['Pressure'], temperature)
            thermal_energy = values['Energy'] / rho - .5 * np.sum(velocity * velocity, axis=1)
            # 与保存守恒量独立比较固定cv*T，绝不使用OTH零点或人为修正能量。
            energy = float(np.max(abs(thermal_energy - CV * temperature) / np.maximum(1., abs(CV * temperature))))
            eos = float(np.max(abs(rho * R * temperature / values['Pressure'] - 1)))
            algebra = float(max(np.max(abs(thermo['e'] / (CV * temperature) - 1)),
                np.max(abs(thermo['h'] / (CP * temperature) - 1)), np.max(abs(thermo['cp'] / CP - 1))))
            viscosity = 1.716e-5 * (temperature / 273.15) ** 1.5 * (273.15 + 110.4) / (temperature + 110.4)
            mu_error = float(np.max(abs(values['Laminar_Viscosity'] / viscosity - 1)))
            conductivity = CP * values['Laminar_Viscosity'] / PR
            flux = evaluate(geom, temperature, conductivity, given_tw)
            ids = geom['ids']
            q = flux['qIntoWall']
            q_error = float(np.max(abs(q - values['Heat_Flux'][ids]) / np.maximum(1., abs(q))))
            step = pair['localIteration']
            assert pair['successfulUpdates'] == pair['seriesSuccessfulUpdates'] == step + 1
            native_power = float(history[step]['HF'])
            power_error = abs(flux['powerIntoWall'] - native_power) / max(1., abs(flux['powerIntoWall']))
            cf = np.column_stack([values['Skin_Friction_Coefficient_' + axis][ids] for axis in 'xyz'])
            shear = np.linalg.norm(cf, axis=1) * .5 * conditions['density_kg_m3'] * conditions['speed_m_s'] ** 2
            yplus = geom['euclideanDistance'] * np.sqrt(rho[ids] * shear) / values['Laminar_Viscosity'][ids]
            domain_passed = bool(np.all((temperature >= 150.) & (temperature <= 360.)))
            strict_passed = bool(max(energy, eos, mu_error, algebra) < 2e-11 and q_error < 2e-10 and power_error < 2e-9)
            filename = 'wall_constant_cp_%06d.csv' % step
            with (out / filename).open('w', encoding='utf8', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(['PointID', 'X_m', 'Y_m', 'Z_m', 'normalNeighbor', 'pressure_Pa',
                    'givenTw_K', 'wallFluidT_K', 'conductivity_WmK', 'normalDistance_m',
                    'euclideanFirstDistance_m', 'normalDistanceFloorTriggered', 'nativeProjectedArea_m2',
                    'qIntoWall_Wm2', 'nodalPowerIntoWall_W', 'independentFirstNeighborYPlus',
                    'thermodynamicModel', 'method', 'solidWallTemperatureSolved', 'physicalQualified'])
                for i, node in enumerate(ids):
                    writer.writerow([int(node), *points[node], int(geom['normalNeighbor'][i]),
                        float(values['Pressure'][node]), float(given_tw[i]), float(temperature[node]),
                        float(conductivity[node]), float(geom['distance'][i]), float(geom['euclideanDistance'][i]),
                        bool(geom['distanceFloorTriggered'][i]), float(geom['area'][i]), float(q[i]),
                        float(flux['nodalPowerIntoWall'][i]), float(yplus[i]), MODEL, flux['method'], False, False])
            measurement = original_segment_powers(meta, geom, flux, rows, summary)
            measurement.update(fieldSource=str(path), fieldSHA256=digest, localIteration=step,
                               seriesSuccessfulUpdates=step + 1)
            measurement_name = 'measurement_%06d.json' % step
            segment_name = 'segments_%06d.csv' % step
            save(out / measurement_name, measurement)
            write_segment_csv(out / segment_name, measurement)
            residuals = {key: float(history[step][key]) for key in ('rms[Rho]', 'rms[RhoE]')}
            frame = dict(localIteration=step, successfulUpdates=step + 1, seriesSuccessfulUpdates=step + 1,
                interfacePassed=strict_passed, strictThermalConsistencyPassed=strict_passed,
                energyMaximumRelativeError=energy, eosMaximumRelativeError=eos,
                constantCpAlgebraMaximumRelativeError=algebra, viscosityMaximumRelativeError=mu_error,
                nativeHeatMaximumScaledError=q_error, powerScaledDifference=power_error,
                nativeNodalPowerIntoWallW=flux['powerIntoWall'], nativeHistoryPowerW=native_power,
                geometricAreaSectorM2=geom['discreteGeometricArea'], nativeProjectedAreaSectorM2=geom['nativeProjectedArea'],
                temperatureRangeK=[float(temperature.min()), float(temperature.max())],
                developmentTemperatureRangeK=[150., 360.], developmentModelDomainPassed=domain_passed,
                belowDevelopmentRangeNodes=int(np.sum(temperature < 150.)),
                aboveDevelopmentRangeNodes=int(np.sum(temperature > 360.)), temperatureClipped=False,
                givenTwRangeK=[float(given_tw.min()), float(given_tw.max())],
                wallFluidTemperatureRangeK=[float(temperature[ids].min()), float(temperature[ids].max())],
                givenTwMaximumFluidDifferenceK=float(np.max(abs(given_tw - temperature[ids]))),
                independentFirstNeighborYPlus99=float(np.percentile(yplus, 99)),
                independentFirstNeighborYPlusMaximum=float(yplus.max()),
                independentFirstNeighborYPlusDefinition='Euclidean original normal-neighbor distance * sqrt(rho_wall*abs(tau_native))/mu_wall',
                normalDistanceFloorTriggered=int(geom['distanceFloorTriggered'].sum()), residuals=residuals,
                originalResidualGatePassed=bool(max(residuals.values()) <= -8),
                strictThermalAuditDoesNotGrantModelDomain=True, heatStabilityQualified=False,
                physicalQualified=False, experimentalQualification=False, newQualifiedArea_m2=0,
                wallCSV=filename, wallCSV_SHA256=sha(out / filename),
                measurementJSON=measurement_name, measurementSHA256=sha(out / measurement_name),
                segmentCSV=segment_name, segmentCSV_SHA256=sha(out / segment_name),
                modelAwareEdgeAnalysisImplemented=False, edgeAndReductionStatus=measurement['edgeAndReductionStatus'])
            frames.append(frame)
        saved_passed = bool(frames and all(frame['interfacePassed'] for frame in frames))
        identity.update(thermalInterfaceAudited=True, savedFrameInterfacePassed=saved_passed)
        report.update(reviewCompleted=True, interfacePassed=bool(identity['normalRunCompleted'] and saved_passed),
            savedFrameInterfacePassed=saved_passed, normalRunCompleted=identity['normalRunCompleted'],
            auditScope=identity['auditScope'], savedFrameIdentity=identity, frames=frames,
            allSavedFramesWithinDevelopmentModelDomain=bool(frames and all(frame['developmentModelDomainPassed'] for frame in frames)),
            developmentTemperatureRangeK=[150., 360.], modelDomainIsSeparateFromThermalInterface=True,
            appliedTemperatureCheck=applied, inletStateCheck=inlet, rejectedHistoryRows=rejected_history,
            iterationOrigin={key: init.get(key) for key in ('kind', 'priorSuccessfulUpdates',
                'reusesPreviousConservativeState', 'energyReference')},
            stage='same-flow-update-laminar-no-SST', nativeLaminarYPlusIsQualification=False,
            givenTwIsSolvedSolidTemperature=False, solvedSolidWallTemperature=False,
            modelAwareEdgeAnalysisImplemented=False, experimentalComparisonReady=False,
            gridAndDomainQualified=False, originalTransportTableEquivalent=False)
        assert saved_passed, '严格定比热空间壁热接口未过；保留失败，不重算CFD'
    except Exception as error:
        report.update(interfacePassed=False, error=type(error).__name__ + ': ' + str(error), traceback=traceback.format_exc())
        raise
    finally:
        stable = {name: Path(name).is_file() and sha(Path(name)) == digest for name, digest in locks.items()}
        if not all(stable.values()):
            report.update(interfacePassed=False, savedFrameInterfacePassed=False, inputMutationDetected=True)
        report.update(inputSHA256=locks, inputsStable=stable, elapsedSeconds=time.monotonic() - started)
        save(out / 'input_manifest.json', locks)
        save(out / 'report.json', report)
    assert all(report['inputsStable'].values()), '审核期间输入改变'
    print(dict(interfacePassed=report['interfacePassed'], savedFrameInterfacePassed=report['savedFrameInterfacePassed'],
        normalRunCompleted=report['normalRunCompleted'], frames=len(frames),
        allSavedFramesWithinDevelopmentModelDomain=report['allSavedFramesWithinDevelopmentModelDomain']))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--saved-frames-only', action='store_true',
                        help='仅显式允许审核已终态异常任务的完整帧；不补造终止瞬间场')
    arguments = parser.parse_args()
    audit(arguments.folder, arguments.output, arguments.saved_frames_only)
