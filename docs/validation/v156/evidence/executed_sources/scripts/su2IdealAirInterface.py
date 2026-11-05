"""变比热理想空气的显式配置合同，独立于旧定比热SST和HEOS入口。"""
from pathlib import Path
import math
import re

from su2ThermalRuntime import read, sha
from su2SSTInterface import PHASE_OUTPUT

FLAVOR = 'rah-ideal-air-v1'
PHASE_FLAVOR = 'rah-ideal-air-phase-v1'
MODEL = 'Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1'
R = 8.31451 / .02896546


def settings(text):
    pairs = re.findall(r'^\s*([A-Z0-9_]+)\s*=\s*([^%\r\n]+)', text, re.M)
    assert len(pairs) == len(dict(pairs)), '重复配置项'
    return {key: value.strip() for key, value in pairs}


def validate(folder, manifest, text, runtime):
    """创建进程前核对模型、参考态、边界、SST阶段和逐字输入身份。"""
    c = settings(text)
    path = folder / 'ideal_air_contract.json'
    contract = read(path)
    flavor = manifest['requiredRuntimeFlavor']
    assert runtime['flavor'] == flavor and flavor in (FLAVOR, PHASE_FLAVOR)
    assert runtime['thermodynamicModel'] == contract['thermodynamicModel'] == MODEL
    assert contract == manifest['idealAir'] and sha(path) == manifest['idealAirContractSHA256']
    assert contract['temperatureRangeK'] == [60, 2000] and contract['gasConstant'] == R
    assert contract['energyEntropyReference'] == 'OTH'
    general = contract['initialization'] == 'general-mesh-ideal-air-SST-wall-distance-v1'
    sparse = flavor == PHASE_FLAVOR
    assert general or contract['initialization'] == 'native-uniform-freestream' or sparse
    if sparse:
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
        VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE, FLOW_PHASE)' if sparse else '(COORDINATES, SOLUTION, PRIMITIVE)',
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
    boundary = [s.strip() for s in c['MARKER_RIEMANN'].strip('() ').split(',')]
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
        assert sst['model'] == 'SST-2003m' and sst['thermodynamics'] == MODEL
        assert sst['outputStateContract'] == PHASE_OUTPUT and sst['provenance'].strip()
        assert sst['physicsQualified'] is False and sst['newQualifiedArea'] == 0
        for key, value in dict(KIND_TURB_MODEL='SST', SST_OPTIONS='V2003m',
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
    if sparse:
        from prepareSU2IdealAirRestart import validate_source
        locks.update(validate_source(folder, c, manifest, contract))
    return locks
