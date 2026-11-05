"""新批处理状态、斜三角壁距初场和含k总焓的独立代数检查。"""
from pathlib import Path
import argparse
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from idealAirProperties import IdealAir
from prepareSU2IdealAirWallDistance import assemble
from prepareSU2WallDistanceState import surface_distance
from validateIdealAir import run, save, sha, write_rows, read_rows, relative


def main(output):
    output = Path(output).resolve()
    output.relative_to(ROOT/'outputs')
    output.mkdir(parents=True, exist_ok=False)
    core = IdealAir()
    temperatures = np.tile([60., 80., 226.5, 300., 1500., 2000.], 3)
    pressure = np.repeat([10., 1000., 100000.], 6)
    pt = core.state('PT', pressure, temperatures)
    ph = core.state('Ph', pressure, pt['h'])
    energy = core.state('rhoe', pt['rho'], pt['e'])
    errors = dict(enthalpyInverse=relative(ph['T'], temperatures),
                  energyInverse=relative(energy['T'], temperatures))
    assert max(errors.values()) < 2e-11
    rejected = []
    for label, args in [('lowT', ('PT', 1000, 59.)), ('highT', ('PT', 1000, 2001.)),
        ('pressure', ('PT', 0., 300.)), ('nan', ('PT', 1000, np.nan)),
        ('lowH', ('Ph', 1000, np.nextafter(pt['h'][0], -np.inf))),
        ('highH', ('Ph', 1000, np.nextafter(pt['h'][5], np.inf))),
        ('shape', ('PT', np.ones((2, 1)), np.ones(2))), ('empty', ('PT', [], [])),
        ('pair', ('bad', 1., 2.)), ('partialBatch', ('PT', 1000, [300., 2001., 300.]))]:
        try:
            core.state(*args)
        except ValueError:
            rejected.append(label)
    assert len(rejected) == 10
    assert core.state('PT', 1000., 300.)['T'][0] == 300.
    vertices = np.array([[0., 0, 0], [2, 0, .4], [0, 1, .3]])
    normal = np.cross(vertices[1]-vertices[0], vertices[2]-vertices[0])
    normal /= np.linalg.norm(normal)
    points = np.r_[vertices, [vertices.mean(axis=0)+z*normal for z in [.01, .13, .9]]]
    faces = np.array([[0, 1, 2]])
    distance, nearest, weights, _ = surface_distance(points, points, faces, .2)
    distance[:3] = 0.
    detail = dict(distance=distance, nearestFace=nearest, weights=weights,
                  wallIDs=np.arange(3), wallLengths=np.array([.001, .002, .003]), farIDs=np.array([5]))
    args = (.2, np.array([-1200., 100., 400.]), 1500., 230., 300., 25., 180., .0008)
    fields, thermo, h_error = assemble(core, detail, faces, *args)
    assert np.max(abs(distance[3:]-[.01, .13, .2])) < 2e-15
    assert np.all(fields[:3, 1:4] == 0) and np.all(fields[:3, 5] == 0)
    assert np.all(thermo['T'][:3] == 300.) and thermo['T'][5] == 230.
    assert np.max(abs(fields[5, 1:4]/fields[5, 0]-args[1])) < 2e-13
    assert np.max(abs(fields[:3, 6]/(60*.0008/(.075*detail['wallLengths']**2))-1)) < 2e-15
    perm = np.array([4, 0, 5, 2, 1, 3])
    inverse = np.argsort(perm)
    shifted = dict(detail, distance=distance[perm], nearestFace=nearest[perm], weights=weights[perm],
                   wallIDs=inverse[detail['wallIDs']], farIDs=inverse[detail['farIDs']])
    permuted, _, _ = assemble(core, shifted, inverse[faces], *args)
    assert relative(permuted[inverse], fields) < 2e-14
    try:
        assemble(core, dict(detail, farIDs=np.array([3])), faces, *args)
    except AssertionError:
        rejected.append('farfieldInsideBlend')
    assert len(rejected) == 11
    # CoolProp独立理想项同时验证新Ph入口的状态和非均匀混合场。
    query_rho = np.r_[pt['rho'], thermo['rho']]
    query_t = np.r_[pt['T'], thermo['T']]
    helper = ROOT/'outputs/v153_validation/native03/idealAirProbe.exe'
    assert sha(helper) == json.loads((helper.parent/'runtime.json').read_text())['sha256']
    write_rows(output/'oracle_input.csv', [('rhoT', float(r), float(t)) for r, t in zip(query_rho, query_t)])
    run(output, 'oracle', [str(helper), '--oracle', str(output/'oracle_input.csv'), str(output/'oracle.csv')], 30)
    oracle = read_rows(output/'oracle.csv')
    assert np.all(oracle['accepted'] == 1)
    for name in ('p', 'e', 'h', 'cp', 'sound'):
        errors['independent_'+name] = relative(np.r_[pt[name], thermo[name]], oracle[name])
    assert max(errors.values()) < 2e-10
    report = dict(passed=True, stateCount=len(query_t), rejected=rejected, errors=errors,
        mixedTotalEnthalpyError=h_error, permutationInvariant=True, exactWallAndFarfield=True,
        sourceSHA256=dict(core.locks, **{str(Path(__file__).resolve()): sha(Path(__file__).resolve()),
            str(ROOT/'scripts/prepareSU2IdealAirWallDistance.py'): sha(ROOT/'scripts/prepareSU2IdealAirWallDistance.py')}),
        physicsQualified=False, newQualifiedArea=0)
    save(output/'report.json', report)
    print(json.dumps({k: v for k, v in report.items() if k != 'sourceSHA256'}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    main(parser.parse_args().output)
