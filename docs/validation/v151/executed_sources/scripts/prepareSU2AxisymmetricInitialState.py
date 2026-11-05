"""为原公开轴对称域生成边界相容守恒初场，并调用既有受限运行器。"""
import argparse
import json
import os
from pathlib import Path

import numpy as np

from prepareSU2AxisymmetricBenchmark import ROOT, sha, save_json
from prepareSU2CheckpointContinuation import config
from runSU2AxisymmetricStudy import run
from su2AxisymmetricAnalysis import mesh_data
from su2SSTInterface import state, wall_geometry


def segment_projection(points, wall_points, segments):
    """分批求点到原闭线段的精确投影，保留脚点、线段编号和插值系数。"""
    start = wall_points[segments[:, 0]]
    edge = wall_points[segments[:, 1]] - start
    norm2 = np.sum(edge**2, axis=1)
    assert min(norm2) > 0
    distances = np.empty(len(points))
    indices = np.empty(len(points), dtype=int)
    weights = np.empty(len(points))
    feet = np.empty_like(points)
    for begin in range(0, len(points), 512):
        stop = min(begin + 512, len(points))
        delta = points[begin:stop, None] - start
        weight = np.clip(np.sum(delta*edge, axis=2)/norm2, 0, 1)
        displacement = delta - weight[:, :, None]*edge
        distance2 = np.sum(displacement**2, axis=2)
        chosen = np.argmin(distance2, axis=1)
        row = np.arange(stop-begin)
        indices[begin:stop] = chosen
        weights[begin:stop] = weight[row, chosen]
        distances[begin:stop] = np.sqrt(distance2[row, chosen])
        feet[begin:stop] = start[chosen]+weight[row, chosen, None]*edge[chosen]
    return distances, indices, weights, feet


def initial_state(points, markers, wall_ids, wall_lengths, conditions):
    """构造C1速度/总焓、固定压力及正SST初值，壁节点优先于共享入口节点。"""
    c = conditions
    distance, index, weight, feet = segment_projection(points, points, markers['wall'])
    inlet = np.unique(markers['inlet'])
    top = np.unique(markers['top'])
    length = np.maximum(np.max(points[inlet, 0])-feet[:, 0], 0)
    thickness = np.zeros(len(points))
    positive = length > 0
    thickness[positive] = .37*length[positive]/(c['unitReynoldsPerM']*length[positive])**.2
    eta = np.ones(len(points))
    np.divide(distance, thickness, out=eta, where=positive)
    eta = np.clip(eta, 0, 1)
    blend = eta**2*(3-2*eta)
    open_ids = np.setdiff1d(np.union1d(inlet, top), wall_ids)
    blend[open_ids] = 1
    blend[wall_ids] = 0
    cp = c['gamma']*c['gasConstant']/(c['gamma']-1)
    speed = c['speedMS']
    temperature = c['wallTemperatureK']+(c['temperatureK']-c['wallTemperatureK'])*blend
    temperature += speed**2/(2*cp)*blend*(1-blend)
    temperature[open_ids] = c['temperatureK']
    velocity = np.column_stack((-speed*blend, np.zeros(len(points))))
    density = c['pressurePa']/(c['gasConstant']*temperature)
    kinf = 1.5*(speed*c['turbulenceIntensity'])**2
    muinf = 1.716e-5*(c['temperatureK']/273.15)**1.5*383.55/(c['temperatureK']+110.4)
    winf = c['densityKgM3']*kinf/(muinf*c['eddyToLaminarRatio'])
    k = kinf*blend
    lengths = np.full(len(points), np.nan)
    lengths[wall_ids] = wall_lengths
    seg = markers['wall'][index]
    ell = (1-weight)*lengths[seg[:, 0]]+weight*lengths[seg[:, 1]]
    tw = c['wallTemperatureK']
    nuwall = 1.716e-5*(tw/273.15)**1.5*383.55/(tw+110.4)/(c['pressurePa']/(c['gasConstant']*tw))
    near = 60*nuwall/(.075*(distance+ell)**2)
    omega = np.exp((1-blend)*np.log(near)+blend*np.log(winf))
    omega[open_ids] = winf
    omega[wall_ids] = 60*nuwall/(.075*wall_lengths**2)
    energy = density*(c['gasConstant']*temperature/(c['gamma']-1)
                      + .5*np.sum(velocity**2, axis=1)+k)
    values = np.column_stack((density, density[:, None]*velocity, energy, k, omega))
    assert np.isfinite(values).all() and min(density) > 0 and min(omega) > 0
    assert min(temperature) > 0 and min(k) >= 0
    assert np.max(distance[wall_ids]) < 1e-12
    assert np.all(temperature[wall_ids] == tw) and np.all(velocity[wall_ids] == 0)
    assert np.all(temperature[open_ids] == c['temperatureK'])
    return values, dict(distance=distance, thickness=thickness, blend=blend,
                        temperature=temperature, feet=feet, openIDs=open_ids,
                        kInfinity=kinf, omegaInfinity=winf)


def prepare(source, folder, plan):
    """只变初场及明确的输出/预算，原网格、模型、边界和离散参数逐项比对。"""
    source = Path(source).resolve()
    folder = Path(folder).resolve()
    plan = Path(plan).resolve()
    folder.relative_to(ROOT/'outputs')
    assert not folder.exists() and plan.is_file()
    old_manifest = json.loads((source/'study_manifest.json').read_text(encoding='utf8'))
    old_receipt = json.loads((source/'engine_process.json').read_text(encoding='utf8'))
    assert old_receipt['status'] == 'terminal' and old_receipt['exitCode'] == 1
    assert all(old_receipt['inputsStable'].values())
    for path, expected in old_receipt['sourceSHA256'].items():
        assert sha(path) == expected, path
    c = old_manifest['benchmark']['conditions']
    points, cells, markers = mesh_data(source/'mesh.su2')
    assert len(points) == 16261 and len(cells) == 16000
    geometry = wall_geometry(source/'mesh.su2')
    values, detail = initial_state(points, markers, geometry['wallIDs'],
                                   geometry['nearestWallDistance'], c)
    current = config(source/'case.cfg')
    changes = dict(ITER='1200', OUTPUT_WRT_FREQ='100, 100', SCREEN_WRT_FREQ_INNER='100')
    current.update(changes)
    folder.mkdir(parents=True)
    os.link(str(source/'mesh.su2'), str(folder/'mesh.su2'))
    (folder/'case.cfg').write_bytes(('\n'.join(k+'= '+v for k, v in current.items())+'\n').encode('ascii'))
    columns = ['PointID', 'x', 'y', 'Density', 'Momentum_x', 'Momentum_y',
               'Energy', 'Turb_Kin_Energy', 'Omega']
    data = np.column_stack((np.arange(len(points)), points, values))
    np.savetxt(str(folder/'initial_flow.csv'), data, delimiter=',',
               fmt=['%d']+['%.17e']*8, header=','.join('"'+s+'"' for s in columns), comments='')
    a, _, nd, rho, velocity, thermal = state(folder/'initial_flow.csv', True)
    error = float(np.max(abs(thermal*(c['gamma']-1)/c['gasConstant']/detail['temperature']-1)))
    assert nd == 2 and error < 2e-12 and np.array_equal(a[:, 1:3], points)
    assert {k: v for k, v in current.items() if k not in changes} == {
        k: v for k, v in config(source/'case.cfg').items() if k not in changes}
    shared = np.intersect1d(markers['wall'], markers['inlet']).tolist()
    init = dict(kind='axisymmetric-C1-conservative-independent-initial-guess',
                initialGuessOnly=True, accepted=True, nodes=len(points), quadrilaterals=len(cells),
                method='exact line distance; C1 velocity/total enthalpy; positive SST',
                thicknessRule='0.37*l/(unitRe*l)^0.2; axial l from inlet to wall foot',
                maximumThicknessM=float(max(detail['thickness'])),
                temperatureRangeK=[float(min(detail['temperature'])), float(max(detail['temperature']))],
                givenWallTemperatureK=c['wallTemperatureK'], wallNoSlip=True,
                inletTopNonwallExact=True, sharedInletWallIDs=shared,
                sharedBoundaryPolicy='wall initial state; native inlet exterior flux retained',
                EOSRelativeError=error, kInfinity=detail['kInfinity'], omegaInfinity=detail['omegaInfinity'])
    save_json(folder/'initialization_manifest.json', init)
    dependencies = [Path(__file__).resolve(), ROOT/'scripts/su2AxisymmetricAnalysis.py',
                    ROOT/'scripts/su2BenchmarkAcceptance.py', ROOT/'scripts/auditSU2AxisymmetricStudy.py',
                    ROOT/'tests/testSU2AxisymmetricInitialization.py',
                    source/'study_manifest.json', source/'engine_process.json',
                    source/'case.cfg', source/'mesh.su2', folder/'initialization_manifest.json',
                    folder/'initial_flow.csv']
    manifest = dict(kind='NASA-ACC-SBLI-M7-axisymmetric-convergence-v1',
                    benchmark=old_manifest['benchmark'], sourceDirectory=str(source),
                    sourceLocalIteration=None, independentInitialization=init,
                    sourceOriginalStartupPassed=True, controlsChanged=changes,
                    samePhysicsAndNumerics=True, maximumSeconds=720, maximumRSSBytes=4*1024**3,
                    maximumOutputBytes=192*1024**2, sourceSHA256={str(p): sha(p) for p in dependencies},
                    plan=str(plan), planSHA256=sha(plan), physicsQualified=False, newQualifiedArea=0)
    save_json(folder/'study_manifest.json', manifest)
    print(json.dumps(init), flush=True)
    return manifest


def run_prepared(folder):
    """显式记录初场实现与验收源码的运行期稳定性，原生监督器拒绝重复任务。"""
    folder = Path(folder).resolve()
    manifest = json.loads((folder/'study_manifest.json').read_text(encoding='utf8'))
    for path, expected in manifest['sourceSHA256'].items():
        assert sha(path) == expected, path
    code = run(folder, manifest)
    stable = {path: sha(path) == expected for path, expected in manifest['sourceSHA256'].items()}
    save_json(folder/'initialization_run_binding.json', dict(exitCode=code, sourceStableAfter=stable,
              sourceSHA256=manifest['sourceSHA256']))
    return code if code else (0 if all(stable.values()) else 2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source')
    parser.add_argument('--folder', required=True)
    parser.add_argument('--plan')
    parser.add_argument('--run-prepared', action='store_true')
    args = parser.parse_args()
    if args.run_prepared:
        raise SystemExit(run_prepared(args.folder))
    prepare(args.source, args.folder, args.plan)
