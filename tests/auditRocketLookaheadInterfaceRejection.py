"""公开压力曲线复算原内核定位的越温界面，不重放完整NS。"""
from pathlib import Path
import json, os, re, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from riemannReservoirReference import sha, exact, PUBLIC_BINDINGS
from auditRocketTemperatureProjection import row_reference

folder = Path(sys.argv[1]).resolve(); file = folder/'interface_rejection.mat'; digest = sha(file)
s = loadmat(file, simplify_cells=True); report = s['report']; inputs = s['rejectedInputs']; gas = inputs['gas']
cv = gas['R']/(gas['gamma']-1); checks = []
def primitive(q):
    rho = q[0]; u = q[1:4]/rho
    return np.r_[rho, u, (q[4]-.5*rho*(u@u))/(rho*cv)]
messages = inputs['messages']; messages = [messages] if isinstance(messages, str) else list(messages)
ids = np.atleast_1d(inputs['faceIDs'])
for face, ql, qr, normal, message in zip(ids, np.atleast_2d(inputs['QL']), np.atleast_2d(inputs['QR']), np.atleast_2d(inputs['normal']), messages):
    left, right = primitive(ql), primitive(qr)
    assert 80 <= left[4] <= 2000 and 80 <= right[4] <= 2000
    sample, _, region = exact(left, right, normal, gas)
    recorded = float(re.search(r'T=([^ ]+) K', message).group(1))
    assert sample[4] < 80 or sample[4] > 2000
    assert abs(sample[4]-recorded) <= 1e-8*max(1, abs(recorded))
    checks.append(dict(faceID=int(face), candidateSampleTemperatureK=float(sample[4]), candidateSampleRegion=int(region),
                       originalKernelRefusedSampleTemperatureK=recorded, publicDifferenceK=float(abs(sample[4]-recorded))))
assert not report['temperatureDomainRelaxed'] and report['newActualNSEvaluations'] == report['newFlowAcceptances'] == 0
public = row_reference(s['baseRows'], 3e-7)
assert sha(file) == digest
out = dict(status='passed-public-reference-for-original-kernel-rejected-interface-and-current-base-H',
           rejectedInterfaces=checks, currentBaseTemperatureDerivatives=public,
           originalKernelCoveredAllInteriorFaces=report['completeInteriorFaceCount'],
           allInteriorInterfacesIndependentlyRecomputed=False, newActualNSEvaluations=0, newFlowAcceptances=0,
           physicalSolution=False, newProductionSolvedAreaM2=0,
           inputSHA256={str(file): digest, **PUBLIC_BINDINGS}, scriptSHA256=sha(__file__))
target = folder/'independent_interface_rejection.json'; assert not target.exists()
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], rejectedInterfaces=checks)))
