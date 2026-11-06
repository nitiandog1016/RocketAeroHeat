"""独立按原节点/原三角面复算新包导出，确认方法与未知壁温身份。"""
import csv
from pathlib import Path
import sys

import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from auditSU2SSTState import vtk_tokens
from su2ThermalRuntime import read, sha
from validateIdealAir import save

base = ROOT/'outputs/v166_validation'
source = ROOT/'outputs/v157_validation/fixture_workbench/mesh_input.mat'
mesh = loadmat(str(source), squeeze_me=True)
V, F = mesh['V'], mesh['F'].astype(int)-1
xyz, fields = vtk_tokens(base/'small01/surface_flow.vtk')
distance, indices = cKDTree(V).query(xyz)
assert max(distance) < 1e-12 and len(set(indices)) == len(V)
with (base/'main01/fields.csv').open(encoding='utf8') as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == len(F) == 8
errors = {}
for name, native in [('pressureRaw_Pa', 'Pressure'), ('heatFluxRaw_W_m2', 'Heat_Flux'),
                      ('fluidNodeTemperatureMean_K', 'Temperature')]:
    values = np.empty(len(V)); values[indices] = fields[native][:, 0]
    expected = values[F].mean(axis=1)
    actual = np.array([float(row[name]) for row in rows])
    errors[name] = float(max(abs(actual-expected)/np.maximum(1, abs(expected))))
    assert errors[name] < 1e-12
assert all(float(row['givenWallTemperature_K']) == 300 and np.isnan(float(row['solvedWallTemperature_K'])) for row in rows)
with (base/'main01/section.csv').open(encoding='utf-8-sig') as stream:
    section = list(csv.reader(stream))
assert len(section) == 5
assert read(base/'main01/section.csv.json')['field'] == 'pressureRaw'
for row in section[1:]:
    face = int(row[2])-1
    assert abs(float(row[10])-float(rows[face]['pressureRaw_Pa'])) < 1e-12*max(1, abs(float(row[10])))
    assert '17位有效数字ASCII' in row[12] and '同包冷重启' in row[12] and '非无中断等价' in row[12]
for name in ('fields.csv.json', 'section.csv.json'):
    metadata = read(base/'main01'/name)
    assert metadata['algorithmRuntime']['flavor'] == 'rah-ideal-air-admissible-v1'
    assert not metadata['productionQualified']
    assert metadata['convergenceAudit']['checkpointIteration'] == 2
    assert metadata['convergenceAudit']['completedUpdatesBeforeRun'] == 4
    assert metadata['convergenceAudit']['cumulativeCheckpointIteration'] == 6
save(base/'export_check.json', dict(passed=True, originalFaces=8, sectionSegments=4,
    fieldRelativeErrors=errors, methodAndUnknownTwVerified=True,
    inputs={str(p): sha(p) for p in (source, base/'small01/surface_flow.vtk',
        base/'main01/fields.csv', base/'main01/section.csv', Path(__file__).resolve())}, newQualifiedArea=0))
print(errors)
