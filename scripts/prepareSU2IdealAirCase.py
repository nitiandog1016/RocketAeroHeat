"""显式准备变比热空气算例，支持原生均匀来流和一般三维壁距初场。"""
from pathlib import Path
import argparse
import json
import math
import re

from prepareSU2SSTCase import set_value
from su2ThermalRuntime import read, sha
from su2IdealAirInterface import FLAVOR, MODEL, R, settings, validate
from su2SSTInterface import PHASE_OUTPUT
from validateIdealAir import ROOT, paper_cp
from runEmbeddedSU2 import atomic_write_json


def configuration(text, turbulent, provenance, freestream=None, turbulence=None):
    """模块一。组装共用配置与方法清单，计算阶段不改写输入文件。"""
    original = settings(text)
    assert original['RESTART_SOL'] == 'NO' and original['SOLVER'] == 'NAVIER_STOKES'
    assert original['KIND_TURB_MODEL'] == 'NONE' and provenance.strip()
    text = re.sub(r'^\s*(FLUID_NAME|GAMMA_VALUE)\s*=.*\n?', '', text, flags=re.M)
    changes = dict(FLUID_MODEL='IDEAL_AIR_VARIABLE_CP', GAS_CONSTANT=format(R, '.17g'),
        CONV_NUM_METHOD_FLOW='ROE', WRT_RESTART_COMPACT='NO', WRT_RESTART_OVERWRITE='NO',
        WRT_SURFACE_OVERWRITE='NO', OUTPUT_WRT_FREQ='1',
        VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE)')
    temperature = float(original['FREESTREAM_TEMPERATURE'])
    if freestream is None:
        cp = float(paper_cp(temperature))
        sound = (cp/(cp-R)*R*temperature)**.5
    else:
        cp, sound = float(freestream['cp'][0]), float(freestream['sound'][0])
    speed = float(original['MACH_NUMBER'])*sound
    contract = dict(schemaVersion=1, thermodynamicModel=MODEL, gasConstant=R,
        temperatureRangeK=[60, 2000], energyEntropyReference='OTH',
        initialization='native-uniform-freestream', turbulent=bool(turbulent),
        flowFlux='general-EOS Roe', conductivity='molecular mu*cp(T)/Pr; native eddy term once',
        freestreamSpeed=speed, freestreamSoundSpeed=sound, freestreamCp=cp,
        provenance=provenance, physicsQualified=False, newQualifiedArea=0)
    sst = None
    if turbulent:
        parameters = turbulence or dict(intensity=.01, viscosityRatio=10., Prt=.9)
        intensity, ratio, prt = [float(parameters[key]) for key in ('intensity', 'viscosityRatio', 'Prt')]
        assert all(math.isfinite(value) and value > 0 for value in (intensity, ratio, prt)) and intensity < 1
        mu = float(original['MU_REF'])*(temperature/float(original['MU_T_REF']))**1.5
        mu *= (float(original['MU_T_REF'])+float(original['SUTHERLAND_CONSTANT']))/(temperature+float(original['SUTHERLAND_CONSTANT']))
        rho = float(original['FREESTREAM_PRESSURE'])/(R*temperature)
        kinetic = 1.5*(speed*intensity)**2
        sst = dict(schemaVersion=1, model='SST-2003m', thermodynamics=MODEL,
            outputStateContract=PHASE_OUTPUT, intensity=intensity, viscosityRatio=ratio, Prt=prt,
            kInfinity=kinetic, omegaInfinity=rho*kinetic/(mu*ratio),
            assumption='fully turbulent; transition unknown', provenance=provenance,
            physicsQualified=False, newQualifiedArea=0)
        changes.update(SOLVER='RANS', KIND_TURB_MODEL='SST', SST_OPTIONS='V2003m',
            FREESTREAM_TURBULENCEINTENSITY=str(intensity), FREESTREAM_TURB2LAMVISCRATIO=str(ratio),
            TURBULENT_CONDUCTIVITY_MODEL='CONSTANT_PRANDTL_TURB', PRANDTL_TURB=str(prt),
            CONV_NUM_METHOD_TURB='SCALAR_UPWIND', MUSCL_TURB='NO', TIME_DISCRE_TURB='EULER_IMPLICIT',
            CFL_REDUCTION_TURB='1.0', CONV_FIELD='(RMS_DENSITY, RMS_ENERGY, RMS_TKE, RMS_DISSIPATION)')
    for key, value in changes.items():
        text = set_value(text, key, value)
    return text, contract, sst


def publish(folder, text, manifest, contract, sst):
    """模块二。先发布数据和合同，最后原子发布输入清单，部分写出由验证器拒绝。"""
    pending = folder/'case.cfg.pending'
    pending.write_text(text, encoding='utf8')
    pending.replace(folder/'case.cfg')
    contract.update(configurationSHA256=sha(folder/'case.cfg'), meshSHA256=sha(folder/'mesh.su2'))
    if sst is not None:
        atomic_write_json(folder/'sst_contract.json', sst)
        manifest.update(sst=sst, sstContractSHA256=sha(folder/'sst_contract.json'))
    atomic_write_json(folder/'ideal_air_contract.json', contract)
    manifest.update(requiredRuntimeFlavor=FLAVOR, idealAir=contract,
                    idealAirContractSHA256=sha(folder/'ideal_air_contract.json'))
    atomic_write_json(folder/'case_manifest.json', manifest)


def prepare(folder, turbulent=False, provenance='', wall_distance=None):
    """模块三。保留既有三参数入口，新初场采用显式选项和独立来源审核。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    assert not (folder/'engine_process.json').exists() and not (folder/'ideal_air_contract.json').exists()
    assert not (folder/'ideal_air_preparation.json').exists(), '已有初场准备记录，使用新目录保留失败来源'
    if wall_distance is not None:
        assert turbulent, '一般壁距初场当前显式支持SST-2003m'
        from prepareSU2IdealAirWallDistance import prepare_wall_distance
        return prepare_wall_distance(folder, wall_distance, provenance)
    text = (folder/'case.cfg').read_text(encoding='utf-8-sig')
    manifest = read(folder/'case_manifest.json')
    assert not any(key in manifest for key in ('sst', 'initialization', 'continuation', 'thermalProfileSHA256'))
    text, contract, sst = configuration(text, turbulent, provenance)
    publish(folder, text, manifest, contract, sst)
    runtime = read(ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)/'runtime.json')
    validate(folder, manifest, text, runtime)
    return contract


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    parser.add_argument('--sst', action='store_true')
    parser.add_argument('--provenance', required=True)
    parser.add_argument('--wall-distance-m', type=float)
    parser.add_argument('--intensity', type=float, default=.01)
    parser.add_argument('--viscosity-ratio', type=float, default=10.)
    parser.add_argument('--prt', type=float, default=.9)
    args = parser.parse_args()
    spec = None
    if args.wall_distance_m is not None:
        spec = dict(thicknessM=args.wall_distance_m, intensity=args.intensity,
                    viscosityRatio=args.viscosity_ratio, Prt=args.prt)
    print(json.dumps(prepare(args.folder, args.sst, args.provenance, spec), indent=2))
