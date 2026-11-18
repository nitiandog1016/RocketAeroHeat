"""从绑定的存档逐原面比较工程交集与既有球帽候选；不改变任何输入。"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
FILES = {
    'outputs/v176_validation/main02/engineering_data.json': '9518d1d1f44abb83f60664a8c1d072d25797ce16d2182823900cfdf119110b1b',
    'outputs/v104_validation/baseline02/current_case.mat': 'b4889e42a28131c253a56f06472a58c9d2cb290d7c9937935fbb4b5d7349488f',
    'outputs/v255_validation/workbench01/coupled.json': '335cd47819b4c9a8db657107364f2e1f94f5cd133a6d852ecbc2db42ae34c1f4',
    'outputs/v253_validation/inlet02/inlet_profiles.npz': 'be5b68cb1ecb03dccc7622233c682a5531684c7a2509c8eb86ca239e998ca5c2',
}


def audit():
    data_path, mat_path, candidate_path, geometry_path = (ROOT / name for name in FILES)
    for name, expected in FILES.items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    data = json.loads(data_path.read_bytes())
    candidate = json.loads(candidate_path.read_bytes())
    ledger = data['evidenceLedger']
    assert ledger['sourceFrame']['sha256'] == FILES[mat_path.relative_to(ROOT).as_posix()]
    condition = data['resultCondition']
    assert [condition[k] for k in ('Mach', 'altitudeM', 'alphaDeg', 'betaDeg', 'wallUniformValueK')] == [5, 30000, 20, 0, 300]
    area = np.asarray(data['faceArea'])
    assert len(area) == 19752 and np.all(np.isfinite(area) & (area > 0))
    with np.load(geometry_path, allow_pickle=False) as source:
        for key, other in [('vertices','originalVertices'), ('faces','originalFaces'), ('faceArea','originalFaceArea')]:
            assert np.array_equal(data[key], source[other]), key
    total = math.fsum(area)
    assert abs(total - 1.9646292565715553) < 1e-11

    def stats(mask):
        n = int(np.count_nonzero(mask))
        value = math.fsum(area[mask])
        return {'faces': n, 'areaM2': value, 'percentOfOriginalArea': 100 * value / total}

    fields = {f['key']: f for f in data['fields']}
    bound = {f['key']: f for f in ledger['fields']}
    masks = {}
    states = {}
    for key in ('pressure', 'heatFlux', 'edgeTemperature', 'recoveryTemperature', 'wallInput', 'wallSolved'):
        status = np.asarray(fields[key]['status'])
        values = np.asarray(fields[key]['values'], dtype=float)
        assert status.shape == area.shape and np.isin(status, [0, 1, 2]).all()
        assert np.array_equal(status != 0, np.isfinite(values)), key
        masks[key] = status == 1
        states[key] = {name: stats(status == i) for i, name in enumerate(('unknown','computed','given'))}
    p, q = masks['pressure'], masks['heatFlux']
    both = p & q
    classes = {'bothComputed': both, 'onlyPressureComputed': p & ~q,
               'onlyHeatComputed': ~p & q, 'neitherComputed': ~p & ~q}
    groups = {k: stats(m) for k, m in classes.items()}
    assert sum(v['faces'] for v in groups.values()) == len(area)
    assert abs(sum(v['areaM2'] for v in groups.values()) - total) < 1e-11
    assert abs(ledger['joint']['pressureHeatFluxEngineeringArea'] - groups['bothComputed']['areaM2']) < 1e-11
    passed = both & (np.asarray(bound['pressure']['qualification']) == 2) & (np.asarray(bound['heatFlux']['qualification']) == 2)
    assert abs(stats(passed)['areaM2'] - ledger['joint']['pressureHeatFluxQualifiedArea']) < 1e-11
    region = np.asarray(data['region'])
    reasons = np.asarray(data['missingReason'])
    ids = np.asarray(candidate['faceIDs'])
    assert np.issubdtype(ids.dtype,np.integer) and len(ids) == len(set(ids.tolist())) == 1306
    assert min(ids) >= 1 and max(ids) <= len(area)
    cap = np.zeros(len(area), dtype=bool)
    cap[ids - 1] = True  # v255 导出器的 faceIDs 从 1 起算，面序已逐项对照。
    assert abs(stats(cap)['areaM2'] - candidate['identity']['candidateAreaM2']) < 1e-11
    assert candidate['identity']['physicalQualification'] is False
    by_reason = {str(name): stats((reasons == name) & ~q) for name in set(reasons[~q])}
    return dict(inputs=FILES, originalAreaM2=total, faceCount=len(area), condition=condition,
                fields=states, jointGroups=groups, jointQualified=stats(passed),
                jointGapByRegion={str(r): stats((region == r) & ~both) for r in sorted(set(region))},
                heatGapByReason=dict(sorted(by_reason.items(), key=lambda pair:-pair[1]['areaM2'])),
                candidateOverlap={k: stats(cap & m) for k, m in classes.items()},
                candidatePhysicalQualification=False, newSolvedAreaM2=0)


if __name__ == '__main__':
    print(json.dumps(audit(), ensure_ascii=False, indent=2, allow_nan=False))
