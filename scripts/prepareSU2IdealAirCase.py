"""将已有单区小算例显式转为变比热空气；原能量参考的重启动需另行适配。"""
from pathlib import Path
import argparse
import json
import re

from prepareSU2SSTCase import set_value
from su2ThermalRuntime import read, sha
from su2IdealAirInterface import FLAVOR, MODEL, R, settings, validate
from su2SSTInterface import PHASE_OUTPUT
from validateIdealAir import ROOT, paper_cp


def prepare(folder, turbulent=False, provenance=''):
    """保留几何与数值设置，建立相容新初态及完整方法清单。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    assert not (folder/'engine_process.json').exists() and not (folder/'ideal_air_contract.json').exists()
    text = (folder/'case.cfg').read_text(encoding='utf-8-sig')
    original = settings(text)
    assert original['RESTART_SOL'] == 'NO' and original['SOLVER'] == 'NAVIER_STOKES'
    assert original['KIND_TURB_MODEL'] == 'NONE' and provenance.strip()
    manifest = read(folder/'case_manifest.json')
    assert not any(key in manifest for key in ('sst', 'initialization', 'continuation', 'thermalProfileSHA256'))
    text = re.sub(r'^\s*(FLUID_NAME|GAMMA_VALUE)\s*=.*\n?', '', text, flags=re.M)
    changes = dict(FLUID_MODEL='IDEAL_AIR_VARIABLE_CP', GAS_CONSTANT=format(R, '.17g'),
        WRT_RESTART_COMPACT='NO', WRT_RESTART_OVERWRITE='NO', WRT_SURFACE_OVERWRITE='NO',
        OUTPUT_WRT_FREQ='1', VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE)')
    temperature = float(original['FREESTREAM_TEMPERATURE'])
    cp = float(paper_cp(temperature))
    speed = float(original['MACH_NUMBER'])*(cp/(cp-R)*R*temperature)**.5
    contract = dict(schemaVersion=1, thermodynamicModel=MODEL, gasConstant=R,
        temperatureRangeK=[60, 2000], energyEntropyReference='OTH',
        initialization='native-uniform-freestream', turbulent=bool(turbulent),
        flowFlux='general-EOS Roe', conductivity='molecular mu*cp(T)/Pr; native eddy term once',
        freestreamSpeed=speed, freestreamSoundSpeed=speed/float(original['MACH_NUMBER']),
        freestreamCp=cp, provenance=provenance, physicsQualified=False, newQualifiedArea=0)
    if turbulent:
        intensity, ratio, prt = .01, 10., .9
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
        (folder/'sst_contract.json').write_text(json.dumps(sst, indent=2)+'\n', encoding='utf8')
        manifest.update(sst=sst, sstContractSHA256=sha(folder/'sst_contract.json'))
    for key, value in changes.items():
        text = set_value(text, key, value)
    (folder/'case.cfg').write_text(text, encoding='utf8')
    contract.update(configurationSHA256=sha(folder/'case.cfg'), meshSHA256=sha(folder/'mesh.su2'))
    (folder/'ideal_air_contract.json').write_text(json.dumps(contract, indent=2)+'\n', encoding='utf8')
    manifest.update(requiredRuntimeFlavor=FLAVOR, idealAir=contract,
                    idealAirContractSHA256=sha(folder/'ideal_air_contract.json'))
    (folder/'case_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf8')
    runtime = read(ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)/'runtime.json')
    validate(folder, manifest, text, runtime)
    return contract


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    parser.add_argument('--sst', action='store_true')
    parser.add_argument('--provenance', required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.folder, args.sst, args.provenance), indent=2))
