"""实际启动前的轴对称合同反例、旧入口兼容和截断场拒绝。"""
import copy
import os
from pathlib import Path
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from runEmbeddedSU2 import atomic_write_json
from su2IdealAirCheckpoints import inventory, output_plan
from su2ThermalRuntime import read, resolve, sha


def main():
    base = ROOT/'outputs/v159_validation/contracts01'
    base.mkdir(exist_ok=False)
    source = ROOT/'outputs/v159_validation/public01'
    resolve(source, (source/'case.cfg').read_text('utf8'))
    rejected = []
    for name in ('planar', 'wrongModel', 'wrongInlet', 'wrongVelocity', 'missingPhase'):
        folder = base/name
        folder.mkdir()
        for file in source.iterdir():
            if file.is_file():
                if file.name in ('mesh.su2', 'initial_flow.csv'):
                    os.link(str(file), str(folder/file.name))
                else:
                    shutil.copyfile(str(file), str(folder/file.name))
        manifest = read(folder/'case_manifest.json')
        text = (folder/'case.cfg').read_text('utf8')
        if name == 'planar':
            text = text.replace('AXISYMMETRIC= YES', 'AXISYMMETRIC= NO')
        elif name == 'wrongModel':
            text = text.replace('SST_OPTIONS= (V1994m, VORTICITY)', 'SST_OPTIONS= V2003m')
        elif name == 'wrongInlet':
            text = text.replace('MARKER_SUPERSONIC_OUTLET= (outlet)', 'MARKER_SUPERSONIC_OUTLET= (inlet)')
        elif name == 'wrongVelocity':
            manifest['expectedVelocity'][0] *= .5
        else:
            text = text.replace(', FLOW_PHASE', '')
        (folder/'case.cfg').write_text(text, encoding='utf8')
        manifest['idealAir']['configurationSHA256'] = sha(folder/'case.cfg')
        atomic_write_json(folder/'ideal_air_contract.json', manifest['idealAir'])
        manifest['idealAirContractSHA256'] = sha(folder/'ideal_air_contract.json')
        atomic_write_json(folder/'case_manifest.json', manifest)
        try:
            resolve(folder, text)
        except AssertionError as error:
            rejected.append(dict(name=name, reason=str(error)))
        assert not (folder/'engine_process.json').exists()
    assert len(rejected) == 5
    # 正常旧合同只读解析，保持其已验证模型入口。
    for relative in ('outputs/v154_validation/cases01/laminar', 'outputs/v154_validation/cases01/sst',
                     'outputs/v145_validation/fixture02/case'):
        old = ROOT/relative
        resolve(old, (old/'case.cfg').read_text('utf8'))
    old = ROOT/'outputs/v159_validation/small01'
    broken = base/'truncated'
    broken.mkdir()
    for name in ('case_manifest.json', 'mesh.su2', 'history.csv', 'surface_flow_000000.vtk'):
        shutil.copyfile(str(old/name), str(broken/name))
    rows = (old/'restart_flow_000000.csv').read_bytes().splitlines(keepends=True)
    (broken/'restart_flow_000000.csv').write_bytes(b''.join(rows[:-1]))
    receipt = copy.deepcopy(read(old/'engine_process.json'))
    receipt['outputSHA256'] = {p.name: sha(p) for p in broken.iterdir() if p.is_file()}
    inspected = inventory(broken, receipt)
    assert not inspected['complete'] and len(inspected['rejected']) == 1
    atomic_write_json(broken/'inventory.json', inspected)
    capacity = output_plan(64521, 321, 320, 1200, 100, 1024**3, extra_columns=7)
    assert capacity['estimatedMaximumBytes'] == 994382592
    # 小壳全局单元方向和边配对独立检查。
    from su2AxisymmetricAnalysis import mesh_data, median_volumes
    from collections import Counter
    p, cells, markers = mesh_data(old/'mesh.su2')
    a, v = median_volumes(p, cells)
    edges = Counter((int(a), int(b)) for cell in cells for a, b in zip(cell, np.roll(cell, -1)))
    assert max(edges.values()) == 1
    exposed = set(edge for edge in edges if edge[::-1] not in edges)
    assert exposed == set(tuple(edge) for values in markers.values() for edge in values)
    expected_volume = np.pi*(max(p[:, 1])**2-min(p[:, 1])**2)*(max(p[:, 0])-min(p[:, 0]))
    assert abs(sum(v)/expected_volume-1) < 2e-12 and min(a) > 0
    report = dict(passed=True, rejected=rejected, truncatedFrameRejected=True, legacyContracts=3,
                  fixtureVolumeRelativeError=abs(sum(v)/expected_volume-1), capacity=capacity)
    atomic_write_json(base/'report.json', report)
    print(report)


if __name__ == '__main__':
    main()
