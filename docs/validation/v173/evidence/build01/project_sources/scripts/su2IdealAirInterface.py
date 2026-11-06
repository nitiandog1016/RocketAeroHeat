"""变比热理想空气的显式配置合同，独立于旧定比热SST和HEOS入口。"""
from pathlib import Path
import math
import re

from su2ThermalRuntime import read, sha
from su2SSTInterface import PHASE_OUTPUT

FLAVOR = 'rah-ideal-air-v1'
PHASE_FLAVOR = 'rah-ideal-air-phase-v2'
OBSERVED_FLAVOR = 'rah-ideal-air-observed-v1'
ADMISSIBLE_FLAVOR = 'rah-ideal-air-admissible-v1'
STEP_FLAVOR = 'rah-ideal-air-step-v1'
SST_ENERGY_FLAVOR = 'rah-sst-energy-v1'
SST_ENERGY = dict(model='IDEAL_AIR_VARIABLE_CP', stage='after-native-SST-shared-relaxation',
    update='clip((Solution_Old*rho_current+alpha*LinSysSol)/rho_current)',
    variables=['k','omega'], flowConservativeUnchanged=True, temperatureRangeK=[60,2000],
    roundingInset='32*eps*(abs(E)+velocity2/2+emax+abs(old_k)); anchored at legal zero step',
    maximumBisections=64, tinyAlphaThreshold=1e-10, invalidZeroStep='record-and-rethrow',
    traceFile='sst_energy_updates.csv', maximumTraceEvents=65536, pointNumbering='local-RCM',
    changesNumerics=True, temperatureClipping=False, kineticEnergyReset=False,
    globalConservationProven=False)
SST_ENERGY_CHANGE = 'SST shared k/omega relaxation constrained by current conserved thermal energy; admissible-v1 flow fallback retained'
CONSERVATIVE_STEP = dict(model='IDEAL_AIR_VARIABLE_CP',
    candidateRejection='domain_error-before-primitive-read',
    direction='Solution_Old-to-trial-all-flow-conservative-components',
    currentK='actual-SetPrimVar-turb_ke', stepRule='first-EOS-valid-dyadic-step',
    maximumHalvings=32, noPositiveStep='restore-old-and-record-theta-zero',
    invalidOldState='record-and-rethrow', traceFile='conservative_steps.csv',
    maximumTraceEvents=4096, pointNumbering='native-global-plus-local-RCM',
    changesNumerics=True, temperatureClipping=False, kineticEnergyReset=False,
    globalConservationProven=False)
STEP_CHANGE = 'ideal-air domain rejection halves the old-to-trial conservative segment with current k; invalid old state terminates'
ADMISSIBILITY = dict(model='IDEAL_AIR_VARIABLE_CP',
    candidateRejection='domain_error-before-primitive-read',
    fallback='Solution_Old-with-current-turb_ke', invalidOldState='record-and-rethrow',
    traceFile='admissibility.csv', maximumTraceEvents=4096, pointNumbering='local-RCM',
    changesNumerics=True, temperatureClipping=False, kineticEnergyReset=False)
OBSERVATION = dict(numericalBaselineFlavor=PHASE_FLAVOR, asciiSignificantDigits=17,
                   terminationObservation='cpp-current-exception-windows-stack-delegate-v1')
AXISYMMETRIC_FLAVOR = 'rah-ideal-air-axisymmetric-v1'
PHASE_FLAVORS = ('rah-ideal-air-phase-v1', PHASE_FLAVOR, OBSERVED_FLAVOR, ADMISSIBLE_FLAVOR, STEP_FLAVOR, SST_ENERGY_FLAVOR, AXISYMMETRIC_FLAVOR)
PHASE_FLAVORS = PHASE_FLAVORS + ('rah-ideal-air-spatial-wall-v1',)
MODEL = 'Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1'
R = 8.31451 / .02896546


def settings(text):
    pairs = re.findall(r'^\s*([A-Z0-9_]+)\s*=\s*([^%\r\n]+)', text, re.M)
    assert len(pairs) == len(dict(pairs)), '重复配置项'
    return {key: value.strip() for key, value in pairs}


def validate(folder, manifest, text, runtime):
    """创建进程前核对模型、参考态、边界、SST阶段和逐字输入身份。"""
    c = settings(text)
    if manifest.get('caseKind') == 'naca-tn3344-spatial-step-cold-restart-v1':
        from su2NACASpatialStepInterface import validate as validate_spatial_step_naca
        return validate_spatial_step_naca(folder, manifest, text, runtime)
    if manifest.get('caseKind') == 'naca-tn3344-adaptive-cold-restart-v1':
        from su2NACAAdaptiveInterface import validate as validate_adaptive_naca
        return validate_adaptive_naca(folder, manifest, text, runtime)
    if manifest.get('caseKind') == 'naca-tn3344-spatial-wall-laminar-v1':
        from su2NACASpatialInterface import validate as validate_spatial_naca
        return validate_spatial_naca(folder, manifest, text, runtime)
    if manifest.get('caseKind') == 'naca-tn3344-quarter-3d-laminar-v1':
        # 具名层流小例单独严格核验；旧正半径SST合同和运行包身份保持。
        from su2NACALaminarInterface import validate as validate_naca
        return validate_naca(folder, manifest, text, runtime)
    path = folder / 'ideal_air_contract.json'
    contract = read(path)
    flavor = manifest['requiredRuntimeFlavor']
    assert runtime['flavor'] == flavor and flavor in (FLAVOR, PHASE_FLAVOR, OBSERVED_FLAVOR, ADMISSIBLE_FLAVOR, STEP_FLAVOR, SST_ENERGY_FLAVOR, AXISYMMETRIC_FLAVOR)
    if flavor == OBSERVED_FLAVOR:
        # 模块用途。诊断和文本精度单独具名，阻止数值基线或运行包身份被静默替换。
        assert manifest['runtimeObservation'] == OBSERVATION
        assert all(runtime[key] == value for key, value in OBSERVATION.items())
    else:
        assert 'runtimeObservation' not in manifest
    if flavor in (ADMISSIBLE_FLAVOR, SST_ENERGY_FLAVOR):
        # 新数值身份仅允许已登记的原态回退，禁止伪称数值相同或静默放宽温区。
        assert manifest['runtimeAdmissibility'] == runtime['admissibility'] == ADMISSIBILITY
        assert runtime['samePhysicsAndNumerics'] is False
        assert runtime['numericalBaselineFlavor'] == flavor
        assert runtime['numericalParentFlavor'] == (OBSERVED_FLAVOR if flavor == ADMISSIBLE_FLAVOR else ADMISSIBLE_FLAVOR)
        for key in ('asciiSignificantDigits', 'terminationObservation'):
            assert runtime[key] == OBSERVATION[key]
    else:
        assert 'runtimeAdmissibility' not in manifest
    if flavor == STEP_FLAVOR:
        # 模块用途。限步独立具名，不把共同方向、有限事件记录误称为全局守恒证明。
        assert manifest['runtimeConservativeStep'] == runtime['conservativeStep'] == CONSERVATIVE_STEP
        assert runtime['samePhysicsAndNumerics'] is False
        assert runtime['numericalBaselineFlavor'] == STEP_FLAVOR
        assert runtime['numericalParentFlavor'] == ADMISSIBLE_FLAVOR
        for key in ('asciiSignificantDigits', 'terminationObservation'):
            assert runtime[key] == OBSERVATION[key]
    else:
        assert 'runtimeConservativeStep' not in manifest
    if flavor == SST_ENERGY_FLAVOR:
        assert manifest['runtimeSSTEnergy'] == runtime['sstEnergy'] == SST_ENERGY
    else:
        assert 'runtimeSSTEnergy' not in manifest
    assert runtime['thermodynamicModel'] == contract['thermodynamicModel'] == MODEL
    assert contract == manifest['idealAir'] and sha(path) == manifest['idealAirContractSHA256']
    assert contract['temperatureRangeK'] == [60, 2000] and contract['gasConstant'] == R
    assert contract['energyEntropyReference'] == 'OTH'
    general = contract['initialization'] == 'general-mesh-ideal-air-SST-wall-distance-v1'
    axisymmetric = flavor == AXISYMMETRIC_FLAVOR
    sparse = flavor in (PHASE_FLAVOR, OBSERVED_FLAVOR, ADMISSIBLE_FLAVOR, STEP_FLAVOR, SST_ENERGY_FLAVOR, AXISYMMETRIC_FLAVOR)
    assert general or contract['initialization'] == 'native-uniform-freestream' or sparse
    if axisymmetric:
        assert contract['initialization'] == 'axisymmetric-ideal-air-SST-Vm-wall-distance-v1'
        assert contract['schemaVersion'] == 4 and contract['turbulent']
    elif sparse:
        assert contract['initialization'] == 'same-mesh-ideal-air-state-v1'
        assert contract['schemaVersion'] == 3 and contract['turbulent']
    if general:
        assert contract['schemaVersion'] == 2 and contract['turbulent']
    assert contract['physicsQualified'] is False and contract['newQualifiedArea'] == 0
    assert sha(folder/'case.cfg') == contract['configurationSHA256']
    assert sha(folder/'mesh.su2') == contract['meshSHA256']
    required = dict(FLUID_MODEL='IDEAL_AIR_VARIABLE_CP', MATH_PROBLEM='DIRECT',
        SYSTEM_MEASUREMENTS='SI', REF_DIMENSIONALIZATION='DIMENSIONAL',
        INIT_OPTION='TD_CONDITIONS', RESTART_SOL='YES' if general or sparse else 'NO', MGLEVEL='0',
        VISCOSITY_MODEL='SUTHERLAND', CONDUCTIVITY_MODEL='CONSTANT_PRANDTL',
        CONV_NUM_METHOD_FLOW='ROE', WRT_RESTART_COMPACT='NO',
        WRT_RESTART_OVERWRITE='NO', WRT_SURFACE_OVERWRITE='NO',
        VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE, FLOW_PHASE, RESIDUAL)' if axisymmetric else
            '(COORDINATES, SOLUTION, PRIMITIVE, FLOW_PHASE)' if sparse else '(COORDINATES, SOLUTION, PRIMITIVE)',
        OUTPUT_WRT_FREQ=str(contract['checkpointContract']['every']) if sparse else '1')
    for key, value in required.items():
        assert c.get(key) == value, key
    if general:
        assert c.get('SOLUTION_FILENAME') == 'initial_flow' and c.get('READ_BINARY_RESTART') == 'NO'
        assert 'RESTART_ITER' not in c and 'continuation' not in manifest
    assert c.get('TIME_DOMAIN', 'NO') == 'NO' and c.get('LOW_MACH_CORR', 'NO') == 'NO'
    assert float(c['GAS_CONSTANT']) == R and float(c['CONV_RESIDUAL_MINVAL']) == -8
    for key in ('MARKER_WALL_FUNCTIONS', 'MARKER_PYTHON_CUSTOM', 'KIND_TRANS_MODEL', 'CONFIG_LIST',
                'NZONES', 'MARKER_GILES', 'GRID_MOVEMENT', 'FLUID_NAME', 'MARKER_HEATFLUX',
                'MARKER_FAR', 'MARKER_INLET', 'MARKER_OUTLET'):
        assert key not in c, key
    assert 60 <= float(c['FREESTREAM_TEMPERATURE']) <= 2000
    for key in ('FREESTREAM_PRESSURE', 'MU_REF', 'MU_T_REF', 'PRANDTL_LAM'):
        assert math.isfinite(float(c[key])) and float(c[key]) > 0, key
    wall = re.fullmatch(r'\(\s*wall\s*,\s*([-+\deE.]+)\s*\)', c['MARKER_ISOTHERMAL'])
    assert wall and 60 <= float(wall[1]) <= 2000
    boundary = [] if axisymmetric else [s.strip() for s in c['MARKER_RIEMANN'].strip('() ').split(',')]
    assert len(boundary) % 7 == 0
    for index in range(0, len(boundary), 7):
        kind = boundary[index+1]
        values = list(map(float, boundary[index+2:index+7]))
        assert all(map(math.isfinite, values)) and values[0] > 0
        assert kind in ('STATIC_SUPERSONIC_INFLOW_PT', 'STATIC_PRESSURE')
        if kind == 'STATIC_SUPERSONIC_INFLOW_PT':
            assert 60 <= values[1] <= 2000
    turbulent = c['SOLVER'] == 'RANS'
    assert c['SOLVER'] in ('RANS', 'NAVIER_STOKES')
    assert contract['turbulent'] == turbulent
    locks = {str(p): sha(p) for p in (path, folder/'case_manifest.json', Path(__file__).resolve())}
    if turbulent:
        sst_path = folder/'sst_contract.json'
        sst = read(sst_path)
        assert sst == manifest['sst'] and sha(sst_path) == manifest['sstContractSHA256']
        assert sst['model'] == ('SST-Vm' if axisymmetric else 'SST-2003m') and sst['thermodynamics'] == MODEL
        assert sst['outputStateContract'] == PHASE_OUTPUT and sst['provenance'].strip()
        assert sst['physicsQualified'] is False and sst['newQualifiedArea'] == 0
        for key, value in dict(KIND_TURB_MODEL='SST', SST_OPTIONS='(V1994m, VORTICITY)' if axisymmetric else 'V2003m',
            TURBULENT_CONDUCTIVITY_MODEL='CONSTANT_PRANDTL_TURB', CONV_NUM_METHOD_TURB='SCALAR_UPWIND',
            MUSCL_TURB='NO', TIME_DISCRE_TURB='EULER_IMPLICIT').items():
            assert c.get(key) == value, key
        for key, name in (('FREESTREAM_TURBULENCEINTENSITY', 'intensity'),
                          ('FREESTREAM_TURB2LAMVISCRATIO', 'viscosityRatio'), ('PRANDTL_TURB', 'Prt')):
            assert 0 < sst[name] == float(c[key]) and math.isfinite(sst[name])
        assert sst['intensity'] < 1
        assert set(re.findall(r'RMS_\w+', c['CONV_FIELD'])) == {'RMS_DENSITY', 'RMS_ENERGY', 'RMS_TKE', 'RMS_DISSIPATION'}
        locks[str(sst_path)] = sha(sst_path)
    else:
        assert c['KIND_TURB_MODEL'] == 'NONE' and 'sst' not in manifest
        assert not (folder/'sst_contract.json').exists()
    if general:
        from prepareSU2IdealAirWallDistance import validate_source
        locks.update(validate_source(folder, c, manifest, contract))
    if axisymmetric:
        from prepareSU2AxisymmetricAir import validate_source
        locks.update(validate_source(folder, c, manifest, contract))
    elif sparse:
        from prepareSU2IdealAirRestart import validate_source
        locks.update(validate_source(folder, c, manifest, contract))
    return locks
