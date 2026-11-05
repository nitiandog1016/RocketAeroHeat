"""公开域续算与表面采样的定向回归，复用保存场避免数值重跑。"""
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from prepareSU2AxisymmetricBenchmark import ROOT
from prepareSU2AxisymmetricContinuation import series_budget
from su2AxisymmetricAnalysis import mesh_data
from su2AxisymmetricSurface import surface_frame, compare_volume
from su2BenchmarkAcceptance import stability


def verify():
    source = ROOT/'outputs/v151_validation/acc_sbli01'
    p, _, markers = mesh_data(source/'mesh.su2')
    contract = json.loads((ROOT/'tests/data/tmr_mach7/contract.json').read_text(encoding='utf8'))
    surface = source/'surface_flow_001199.vtk'
    frame = surface_frame(surface, p, markers, contract)
    errors = compare_volume(frame, source/'restart_flow_001199.csv', p, markers)
    assert series_budget(1200, 10800)['maximumSeriesUpdates'] == 12000
    rejected = 0
    for completed, requested in ((1200, 10801), (12000, 1), (-1, 1), (1200, 0)):
        try:
            series_budget(completed, requested)
        except AssertionError:
            rejected += 1
    assert rejected == 4
    # 对同一小表面改坐标和连接，必须拒绝，避免错误映射生成看似可信的稳定性。
    raw = surface.read_text()
    corrupt = [raw.replace('7.838999e-01', '8.838999e-01', 1),
               raw.replace('2\t1\t0\t', '2\t0\t1\t', 1)]
    with tempfile.TemporaryDirectory() as folder:
        target = Path(folder)/surface.name
        for text in corrupt:
            assert text != raw
            target.write_text(text)
            try:
                surface_frame(target, p, markers, contract)
            except AssertionError:
                rejected += 1
    assert rejected == 6
    frames = [dict(iteration=i, p=frame['p'], q=frame['q'], power=frame['power'])
              for i in list(range(9700, 10800, 100))+[10799]]
    assert stability(frames, 10799, 1.)['spanUpdates'] == 1099
    assert not stability(frames[:-1], 10799, 1.)['available']
    assert not stability(frames[:3]+frames[4:], 10799, 1.)['available']
    assert not stability([frames[0], frames[-1]], 10799, 1.)['available']
    return dict(passed=True, pairedErrors=errors, coordinateErrorM=frame['coordinateErrorM'],
                rejectionChecks=9, scope='saved surface pairing, exact connectivity, lineage budget and sampling')


if __name__ == '__main__':
    print(json.dumps(verify()))
