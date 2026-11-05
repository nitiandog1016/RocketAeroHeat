"""独立核对公开准备，包含总态极限、拒绝输入、原始网格与21点坐标。"""
import argparse
import csv
import json
from pathlib import Path
import sys
import time
import zipfile

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
sys.path.insert(0, str(ROOT/'tests/data/tmr_mach7'))
from idealAirFreestream import total_to_static, static_to_total
from idealAirProperties import IdealAir
from prepareSU2AxisymmetricBenchmark import sha
from su2AxisymmetricAnalysis import mesh_data, median_volumes
from su2IdealAirInterface import R
from validateIdealAir import paper_cp
from verify_contract import read_grid


def verify(folder):
    """模块一。论文式cp独立积分求根，覆盖正常总态、零速和无域内根。"""
    started = time.monotonic()
    folder = Path(folder)
    report_path = folder/'independent_review.json'
    assert not report_path.exists()
    model = IdealAir()
    maximum_error = 0
    cases = [(7.11, 1600*5/9, 364.8*6894.757293168), (0, 60, 1000),
             (0, 2000, 1e5), (5, 700, 1e5), (1, 300, 2e5), (2, 600, 3e4), (7, 1600, 1e6)]
    for mach, t0, p0 in cases:
        value = total_to_static(model, mach, t0, p0)

        def balance(t):
            cp = paper_cp(t)
            return quad(paper_cp, t, t0, epsabs=1e-7)[0] - .5*mach**2*cp/(cp-R)*R*t

        expected_t = t0 if mach == 0 else brentq(balance, 60, t0, xtol=1e-11)
        expected_p = p0*np.exp(-quad(lambda t: paper_cp(t)/t, expected_t, t0, epsabs=1e-8)[0]/R)
        recovered = static_to_total(model, value['p'], value['T'], value['speed'])
        errors = [abs(value['T']/expected_t-1), abs(value['p']/expected_p-1),
                  abs(recovered['totalTemperature']/t0-1), abs(recovered['totalPressure']/p0-1)]
        maximum_error = max(maximum_error, max(errors))
    assert maximum_error < 2e-10
    refused = 0
    bad = [(-1, 900, 1e5), (float('nan'), 900, 1e5), (7, 59, 1e5),
           (7, 2001, 1e5), (7, 900, 0), (100, 900, 1e5)]
    for values in bad:
        try:
            total_to_static(model, *values)
        except ValueError:
            refused += 1
    assert refused == len(bad)

    # 模块二。直接读取原PLOT3D的Fortran排列，逐点比对选定网格和全部连接。
    preparation = json.loads((folder/'preparation.json').read_text('utf8'))
    grid = preparation['grid']
    with zipfile.ZipFile(str(ROOT/grid['archive'])) as archive:
        raw_blocks = read_grid(archive.read(grid['member']))
    blocks = [np.stack((-b[0].T, b[1].T), axis=2) for b in raw_blocks]
    points, cells, markers = mesh_data(folder/'mesh.su2')
    expected = np.concatenate((blocks[0].reshape(-1, 2), blocks[1][:, 1:].reshape(-1, 2)))
    assert np.array_equal(points, expected)
    nj, ni = blocks[0].shape[:2]
    index0 = np.arange(ni*nj).reshape(nj, ni)
    index1 = np.column_stack((index0[:, -1], np.arange(ni*nj, len(points)).reshape(nj, ni-1)))
    expected_cells = np.concatenate([np.stack((ix[:-1, :-1], ix[1:, :-1], ix[1:, 1:], ix[:-1, 1:]),
                                              axis=2).reshape(-1, 4) for ix in (index0, index1)])
    assert np.array_equal(cells, expected_cells)
    planar, volume = median_volumes(points, cells)
    volume_error = abs(sum(volume)/preparation['geometry']['revolvedVolumeM3']-1)
    assert min(planar) > 0 and min(volume) > 0 and volume_error < 2e-10
    boundary_ids = np.unique(np.concatenate(list(markers.values())))
    assert len(boundary_ids) == sum(len(v) for v in markers.values())

    # 模块三。独立回代每个插值组合，核对测点身份和固定旧摩擦尺度估计。
    with (folder/'profile_locations.csv').open(encoding='utf8') as stream:
        profile = list(csv.DictReader(stream))
    maximum_location_error = 0
    contract = json.loads((ROOT/'tests/data/tmr_mach7/contract.json').read_text('utf8'))
    assert len(profile) == len(contract['measurements']['upstream_profile']) == 21
    for row, experiment in zip(profile, contract['measurements']['upstream_profile']):
        ids = [int(row['node%d' % i]) for i in range(4)]
        weights = np.array([float(row['weight%d' % i]) for i in range(4)])
        assert min(weights) >= 0 and abs(sum(weights)-1) < 1e-14
        xy = weights @ points[ids]
        actual = np.array([.06, blocks[0][0, 0, 1]+experiment['y_from_wall_m']])
        maximum_location_error = max(maximum_location_error, max(abs(xy-actual)))
        assert int(row['sourceRow']) == experiment['source_row']
        assert float(row['distanceFromWall']) == experiment['y_from_wall_m']
    assert maximum_location_error < 2e-14
    estimates = np.genfromtxt(str(folder/'near_wall_estimates.csv'), delimiter=',', names=True)
    assert len(estimates) == 321 and np.isfinite(estimates.tolist()).all()
    quantile = float(np.quantile(estimates['fixedOldFrictionYPlus'], .99))
    assert quantile <= 1 and abs(quantile-preparation['nearWall']['estimatedYPlus99']) < 1e-14
    assert preparation['normalization']['directExperimentalTotalTemperatureGate'] == 'pending-reference'
    assert preparation['normalization']['experimentalTotalTemperatureK'] is None
    assert preparation['solverRuns'] == 0 and not preparation['physicsQualified']
    assert all(sha(Path(path)) == expected for path, expected in preparation['sourceSHA256'].items())
    result = dict(passed=True, thermodynamicCases=len(cases), rejectedInputs=refused,
                  maximumThermodynamicRelativeError=maximum_error, allCoordinatesAndCellsEqual=True,
                  independentVolumeRelativeError=volume_error, maximumProfileLocationErrorM=maximum_location_error,
                  sourceStable=True, physicsQualified=False, elapsedSeconds=time.monotonic()-started)
    report_path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    print(json.dumps(verify(parser.parse_args().folder)))
