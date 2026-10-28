# -*- coding: utf-8 -*-
"""独立复核已验方向入口、原NS算术、当前重新检查及八类损坏来源。"""
from pathlib import Path
import hashlib
import json
import os
import sys

for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
case, kind = sys.argv[1:3]
assert case.replace('_', '').isalnum() and kind in ('unit', 'real')
target = VAL / (case + '_independent_validated_probe.json')
assert not target.exists()
bindings = {}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def pathof(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def read(path, names=None):
    path = pathof(path)
    bindings[path.relative_to(ROOT).as_posix()] = sha(path)
    return loadmat(path, simplify_cells=True, variable_names=names)


def equal(a, b):
    if isinstance(a, dict) or isinstance(b, dict):
        return isinstance(a, dict) and isinstance(b, dict) and a.keys() == b.keys() and all(equal(a[k], b[k]) for k in a)
    return np.array_equal(np.asarray(a), np.asarray(b))


def rows(value):
    return [value] if isinstance(value, dict) else list(np.atleast_1d(value))


filename = 'validated_probe_unit_evidence.mat' if kind == 'unit' else 'wall_velocity_krylov_evidence.mat'
s = read(VAL / case / filename)
r = s['result']
logs = [a for a in rows(r['linearLog']) if a.get('initialGuessUsed', False)]
assert len(logs) == r['steadyInitialGuessModelApplications'] == 1
a = logs[0]
g = a['initialGuessEvidence']
source = g['sourceReport']
assert source['sourceVariant'] == 'validated-probe'
p = read(source['sourceFile'])
pr = p['report']
assert sha(pathof(source['sourceFile'])) == source['sourceFileSHA256']
assert equal(p['Q'], r['initialQ']) and equal(p['scale'], r['residualScale'])
assert equal(p['baseRate'], r['initialRate'])
assert equal(p['gas'], s['gas']) and equal(p['boundary'], s['boundary'])
assert source['sourceStateAndGeometryVerified'] and source['originalNSProbeArithmeticVerified']
assert not source['oldDerivativeReused'] and not source['originalNSAcceptanceClaimed']
assert pr['originalNSDirectionChecked'] and pr['originalNSDirectionTargetMet']
assert pr['originalNSDirectionTarget'] == .1 and not pr['physicalStateChanged'] and pr['newFlowAcceptances'] == 0
basefile = pathof(p['sourceBinding']['sourceFile'])
assert sha(basefile) == p['sourceBinding']['sourceFileSHA256']
bindings[basefile.relative_to(ROOT).as_posix()] = sha(basefile)
coarsefile = pathof(p['coarseFile'])
assert sha(coarsefile) == pr['coarseFileSHA256'] == source['sourceCoarseFileSHA256']
coarse = read(coarsefile, ['coarse'])['coarse']
assert equal(coarse['baseQ'], r['initialQ']) and equal(coarse['scale'], r['residualScale'])
assert equal(coarse['gas'], s['gas']) and equal(coarse['boundary'], s['boundary'])
assert equal(coarse['sourceIdentity'], p['sourceIdentity'])
if kind == 'unit':
    geometry = s['cut']
else:
    geometry = read(ROOT / 'outputs/v80_validation/domain02/shell_domain_evidence.mat', ['domain'])['domain']
assert equal(coarse['space']['geometry'], geometry)
for name, digest in p['sourceIdentity'].items():
    candidates = list((ROOT / 'internalFlow').rglob(name + '.m'))
    assert len(candidates) == 1 and sha(candidates[0]) == digest
    bindings[candidates[0].relative_to(ROOT).as_posix()] = digest
del coarse, geometry

rhs = r['initialRate'] / r['residualScale']
np.testing.assert_array_equal(p['rhs'], rhs)
e = p['nsEvidence']
np.testing.assert_array_equal(e['directionU'], p['totalDirection'])
assert e['fullThreeDimensionalResidualCalled'] and e['scheme'] in ('forward', 'central')
h = float(e['actualScalarStep'])
assert h > 0 and np.isfinite(h)
if e['scheme'] == 'forward':
    np.testing.assert_array_equal(e['rateMinusOrBase'], r['initialRate'])
factor = 2 if e['scheme'] == 'central' else 1
jv = (e['ratePlus'] - e['rateMinusOrBase']) / r['residualScale'] / (factor * h)
np.testing.assert_allclose(jv, e['value'], rtol=3e-13, atol=1e-12)
probe_defect = float(np.linalg.norm(rhs + jv) / np.linalg.norm(rhs))
probe_slope = float(np.sum(rhs * jv) / np.sum(rhs * rhs))
assert probe_defect <= .1 and probe_slope < 0 and probe_slope <= -1 + probe_defect + 1e-12
assert abs(probe_defect - pr['originalNSTotalDirectionDefect']) < 1e-12
for sign in ([1, -1] if factor == 2 else [1]):
    q = r['initialQ'] + sign * h * p['totalDirection'] * r['residualScale']
    rho = q[:, 0]
    energy = q[:, 4] - .5 * np.sum(q[:, 1:4] ** 2, axis=1) / rho
    temperature = energy / (rho * s['gas']['R'] / (s['gas']['gamma'] - 1))
    assert np.all(rho > 0) and np.all(energy > 0)
    assert np.all(temperature >= 80) and np.all(temperature <= 2000)

np.testing.assert_array_equal(g['directionU'], p['totalDirection'] * source['directionScale'])
np.testing.assert_array_equal(g['correctionRhs'], rhs + g['currentModelValue'])
model_defect = float(np.linalg.norm(g['correctionRhs']) / np.linalg.norm(rhs))
assert abs(model_defect - a['initialGuessModelRelativeDefect']) < 1e-12
assert a['gmresSkipped'] and model_defect <= a['target'] == .1
assert np.all(g['correctionDirectionU'] == 0)
fresh = a['steadyDirectionEvidence']
np.testing.assert_array_equal(fresh['directionU'], g['directionU'])
fresh_factor = 2 if fresh['scheme'] == 'central' else 1
fresh_jv = (fresh['ratePlus'] - fresh['rateMinusOrBase']) / r['residualScale'] / (fresh_factor * fresh['actualScalarStep'])
np.testing.assert_allclose(fresh_jv, fresh['value'], rtol=3e-13, atol=1e-12)
fresh_defect = float(np.linalg.norm(rhs + fresh_jv) / np.linalg.norm(rhs))
assert fresh_defect <= .1 and abs(fresh_defect - a['fullSystemRelativeDefect']) < 1e-12
assert fresh['fullThreeDimensionalResidualCalled'] and fresh['actualOriginalResidualEvaluations'] >= 1

bad_checks = []
if kind == 'unit':
    assert np.all(s['report']['rejected']) and s['report']['freshMainNSCalls'] == r['residualEvaluations'] == 3
    assert r['acceptedNonlinearIterations'] == 1
    for index in range(1, 9):
        b = read(VAL / case / ('bad_probe_%02d.mat' % index))
        tests = [not equal(b['Q'], r['initialQ']), not equal(b['scale'], r['residualScale']),
                 not b['report']['originalNSDirectionTargetMet'], not equal(b['totalDirection'], b['nsEvidence']['directionU']),
                 not equal(b['nsEvidence']['value'], e['value']), not equal(b['sourceIdentity'], p['sourceIdentity']),
                 sha(pathof(b['sourceBinding']['sourceFile'])) != b['sourceBinding']['sourceFileSHA256'],
                 not equal(b['boundary'], s['boundary'])]
        assert tests[index - 1]
        bad_checks.append(index)
assert not r['physicalSolution'] and r['newProductionSolvedAreaM2'] == 0
assert all(sha(ROOT / path) == digest for path, digest in bindings.items())
out = dict(status='passed-independent-validated-probe-identity-arithmetic-current-fresh-NS-and-guard-evidence',
           case=case, kind=kind, processID=os.getpid(), savedProbeDefect=probe_defect,
           currentModelDefect=model_defect, freshOriginalNSDefect=fresh_defect,
           sourceGuardEvidenceChecked=bad_checks, acceptedFlowUpdates=int(r['acceptedNonlinearIterations']),
           completeIndependentNSSolver=False, productionQualified=False, newProductionSolvedAreaM2=0,
           inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2) + '\n', encoding='utf8')
print(json.dumps(out))
