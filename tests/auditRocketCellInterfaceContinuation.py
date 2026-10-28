"""核对真实当前候选的全部单元/阻挡分区与集合，并以公开解抽查可检查界面。"""
from pathlib import Path
import json
import os
import sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from riemannReservoirReference import ROOT, sha, exact, PUBLIC_BINDINGS
from auditRocketCellInterfaceLookahead import primitive, vector


def main(folder):
    folder = Path(folder).resolve(); step = folder/'continuous01/step01'
    direction = step/'active_direction.mat'; prefix = step/'linear_checkpoints/lookahead_free_prefix.mat'
    context = loadmat(direction, simplify_cells=True, variable_names=['rows', 'sourceBinding', 'preparationFile', 'preparationSHA256', 'geometryBinding', 'directionReport'])
    p = json.loads(Path(context['preparationFile']).read_text('utf8')); assert sha(context['preparationFile']) == context['preparationSHA256']
    assert sha(context['sourceBinding']['sourceFile']) == context['sourceBinding']['sourceFileSHA256'] == p['sourceBinding']['sourceFileSHA256']
    assert sha(context['geometryBinding']['file']) == context['geometryBinding']['sha256']
    s = loadmat(p['artifacts']['cache']['file'], simplify_cells=True, variable_names=['Q', 'gas', 'scale'])
    cut = loadmat(context['geometryBinding']['file'], simplify_cells=True, variable_names=['domain'])['domain']
    saved = loadmat(prefix, simplify_cells=True); selection = context['directionReport']['lookaheadSelection']
    assert selection['selectionPolicy'] == 'original-and-predicted-cell-and-interface-faces'
    Q = s['Q']; scale = s['scale']; gas = s['gas']; low, high = gas['temperatureRange']
    candidate = Q + selection['targetAlpha']*saved['free']*scale
    rho = candidate[:, 0]; velocity = candidate[:, 1:4]/rho[:, None]
    T = (candidate[:, 4]-.5*rho*np.sum(velocity*velocity, axis=1))/(rho*gas['R']/(gas['gamma']-1))
    assert np.all(rho > 0) and np.all(T > 0)
    bad = np.flatnonzero((T < low) | (T > high))+1
    np.testing.assert_array_equal(bad, vector(selection['predictedViolatingCellIDs'], int))
    f = cut['openFaces']; minus = vector(f['cellMinus'], int); plus = vector(f['cellPlus'], int)
    ids = np.flatnonzero((f['area'] > 0) & (minus > 0) & (plus > 0))
    valid = (T >= low) & (T <= high); eligible = valid[minus[ids]-1] & valid[plus[ids]-1]
    assessed = ids[eligible]+1; blocked = ids[~eligible]+1; preview = selection['interfacePreview']
    np.testing.assert_array_equal(assessed, vector(preview['assessedFaceIDs'], int))
    np.testing.assert_array_equal(blocked, vector(preview['blockedFaceIDs'], int))
    np.testing.assert_array_equal(bad, vector(preview['blockedByOutOfRangeCellIDs'], int))
    assert bool(preview['allInteriorFacesAssessed']) == (not len(blocked)) and preview['allInputFacesAccountedFor']
    original = vector(saved['originalRows']['faceIDs'], int); required = vector(selection['requiredOriginalStateFaceIDs'], int)
    interfaces = vector(selection['predictedViolatingInterfaceFaceIDs'], int); rows = context['rows']
    selected = vector(rows['faceIDs'], int); np.testing.assert_array_equal(selected, np.union1d(np.union1d(original, required), interfaces))
    positions = np.searchsorted(selected, original)
    assert (saved['originalRows']['H']-rows['H'][positions]).nnz == 0
    np.testing.assert_array_equal(vector(saved['originalRows']['temperatureBoundSigns']), vector(rows['temperatureBoundSigns'])[positions])
    regions = vector(rows['sampleRegion'], int)[np.searchsorted(selected, required)]
    assert set(regions) <= {1, 4}
    np.testing.assert_array_equal(np.where(regions == 1, minus[required-1], plus[required-1]), bad)
    np.testing.assert_array_equal(vector(rows['ownerCV'], int), minus[selected-1])
    np.testing.assert_array_equal(vector(rows['neighbourCV'], int), plus[selected-1])
    np.testing.assert_array_equal(np.asarray(rows['QL']).reshape(-1, 5), Q[minus[selected-1]-1])
    np.testing.assert_array_equal(np.asarray(rows['QR']).reshape(-1, 5), Q[plus[selected-1]-1])
    # 真实大域不声称把849336个面全部再跑独立标量求根；覆盖所有报告越温面及确定性样本。
    margin = np.minimum(T-low, high-T)
    risk = np.minimum(margin[minus[assessed-1]-1], margin[plus[assessed-1]-1])
    samples = np.unique(np.r_[interfaces, assessed[np.argsort(risk, kind='stable')[:128]],
                              assessed[np.linspace(0, len(assessed)-1, 64, dtype=int)], np.intersect1d(assessed, selected)])
    checked = []
    for face in samples:
        sample, _, region = exact(primitive(candidate[minus[face-1]-1], gas), primitive(candidate[plus[face-1]-1], gas), f['normal'][face-1], gas)
        violating = sample[4] < low or sample[4] > high
        assert violating == (face in interfaces)
        checked.append(dict(faceID=int(face), publicTemperatureK=float(sample[4]), region=region, violating=bool(violating)))
    inputs = [direction, prefix, Path(context['preparationFile']), Path(p['artifacts']['cache']['file']), Path(context['sourceBinding']['sourceFile']), Path(context['geometryBinding']['file'])]
    bindings = {**PUBLIC_BINDINGS, **{path.relative_to(ROOT).as_posix(): sha(path) for path in inputs}}
    out = dict(status='passed-real-all-CV-partition-selected-union-original-H-and-public-interface-samples',
               cellCount=len(Q), interiorFaceCount=len(ids), predictedViolatingCells=len(bad), assessedInterfaceCount=len(assessed),
               blockedInterfaceCount=len(blocked), originalRows=len(original), selectedRows=len(selected), predictedViolatingInterfaces=len(interfaces),
               allCVTemperaturesAndInterfaceInputPartitionIndependentlyRecomputed=True,
               allAssessableRiemannStatesIndependentlyRecomputed=False, publicSampleCount=len(checked), publicSamples=checked,
               originalRowsExactlyRetained=True, inputSHA256=bindings, scriptSHA256=sha(__file__),
               newActualNSEvaluations=0, newFlowAcceptances=0, physicalSolution=False, newProductionSolvedAreaM2=0)
    target = folder/'independent_cell_interface_continuation.json'; assert not target.exists()
    target.write_bytes((json.dumps(out, indent=2)+'\n').encode('utf8'))
    print(json.dumps({k: v for k, v in out.items() if k not in ('publicSamples', 'inputSHA256')}))


if __name__ == '__main__':
    main(sys.argv[1])
