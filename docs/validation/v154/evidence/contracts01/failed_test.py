"""真实运行合同的预启动负例，检查新模型不会放宽旧SST规则。"""
from pathlib import Path
import copy
import json
import os
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from su2ThermalRuntime import resolve, validate_sst_manifest, read, sha
from su2IdealAirInterface import validate
from prepareSU2SSTCase import set_value
from validateIdealAir import save


def main():
    source = ROOT/'outputs/v154_validation/cases01/sst'
    out = ROOT/'outputs/v154_validation/contracts01'
    out.mkdir(exist_ok=False)
    text = (source/'case.cfg').read_text()
    runtime = read(ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-v1/runtime.json')
    resolve(source, text)
    rejected = []
    changes = [('restart', 'RESTART_SOL', 'YES'), ('temperature', 'FREESTREAM_TEMPERATURE', '2001'),
               ('flux', 'CONV_NUM_METHOD_FLOW', 'HLLC'), ('units', 'SYSTEM_MEASUREMENTS', 'US'),
               ('compact', 'WRT_RESTART_COMPACT', 'YES'), ('turbulence', 'KIND_TURB_MODEL', 'SA'),
               ('unsteady', 'TIME_DOMAIN', 'YES'), ('wall_temperature', 'MARKER_ISOTHERMAL', '(wall, 59)')]
    for name, key, value in changes:
        folder = out/name
        folder.mkdir()
        os.link(str(source/'mesh.su2'), str(folder/'mesh.su2'))
        cfg = set_value(text, key, value)
        (folder/'case.cfg').write_text(cfg)
        manifest = copy.deepcopy(read(source/'case_manifest.json'))
        contract = manifest['idealAir']
        contract['configurationSHA256'] = sha(folder/'case.cfg')
        save(folder/'ideal_air_contract.json', contract)
        save(folder/'sst_contract.json', manifest['sst'])
        manifest['idealAirContractSHA256'] = sha(folder/'ideal_air_contract.json')
        save(folder/'case_manifest.json', manifest)
        try:
            validate(folder, manifest, cfg, runtime)
        except AssertionError:
            rejected.append(name)
    assert len(rejected) == len(changes)
    try:
        validate_sst_manifest(source, read(source/'case_manifest.json'), text)
    except AssertionError:
        rejected.append('old_calorically_perfect_contract_rejects_variable_cp')
    assert len(rejected) == len(changes)+1
    save(out/'report.json', dict(passed=True, rejected=rejected, nativeRuns=0, physicsQualified=False))
    print(json.dumps(rejected))


if __name__ == '__main__':
    main()
