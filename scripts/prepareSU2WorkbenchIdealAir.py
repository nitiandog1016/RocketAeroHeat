"""主工作台专用准备。独立保存一般初场来源，再原子发布主目录阶段合同。"""
from pathlib import Path
import argparse
import math
import os
import shutil

from prepareSU2IdealAirCase import prepare as prepare_initial
from prepareSU2IdealAirRestart import prepare as prepare_stage
from su2ThermalRuntime import ROOT, read


def prepare(folder, spec):
    """模块用途。继承正常导出的三维工况、原网格和明确湍流参数，拒绝重复任务。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    assert spec['model'] == 'ideal-air-SST' and spec['mode'] == 'wallDistance'
    # MATLAB JSON可把整数容量写为指数形式，先严格检查再统一为Python整数。
    maximum_bytes = spec['maximumOutputBytes']
    assert isinstance(maximum_bytes, (int, float)) and not isinstance(maximum_bytes, bool)
    assert math.isfinite(maximum_bytes) and maximum_bytes > 0 and maximum_bytes == int(maximum_bytes)
    maximum_bytes = int(maximum_bytes)
    assert not (folder/'engine_process.json').exists() and not (folder/'ideal_air_contract.json').exists()
    manifest = read(folder/'case_manifest.json')
    initial = folder/'ideal_air_initial_state'
    initial.mkdir(exist_ok=False)
    os.link(str(folder/'mesh.su2'), str(initial/'mesh.su2'))
    for name in ('case.cfg', 'case_manifest.json'):
        shutil.copyfile(str(folder/name), str(initial/name))
    parameters = {key: spec[key] for key in ('thicknessM', 'intensity', 'viscosityRatio', 'Prt')}
    prepare_initial(initial, True, spec['provenance'], parameters)
    flight = manifest['flight']
    result = prepare_stage(initial, folder, int(flight['iterations']),
        min(int(flight['outputEvery']), int(flight['iterations'])), maximum_bytes, exported_folder=True)
    assert result['flight']['Mach'] == flight['Mach']
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--request', required=True)
    args = parser.parse_args()
    request = read(Path(args.request).resolve())
    prepare(request['workDirectory'], request['idealAirInitialization'])
