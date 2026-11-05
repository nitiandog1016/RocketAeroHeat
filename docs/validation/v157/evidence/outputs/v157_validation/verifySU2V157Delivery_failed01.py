"""独立复核v157实际主流程的原面/周向CSV、阶段能量、容量及取消终态。"""
from pathlib import Path
import csv
import json
import sys

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from auditSU2SSTState import vtk_tokens
from runEmbeddedSU2 import atomic_write_json
from su2ThermalRuntime import read, sha
from su2IdealAirCheckpoints import output_plan, inventory


def relative(actual, expected):
    return float(np.max(abs(actual-expected)/np.maximum(1,abs(expected))))


def main():
    output = ROOT/'outputs/v157_validation'
    ui = read(output/'main03/report.json')
    assert ui['passed'] and read(output/'review01/report.json')['nativeMouseObserved']
    assert read(output/'process_tree/report.json')['passed']
    folder = Path(ui['engineFolder'])
    receipt = read(folder/'engine_process.json')
    pipeline = read(folder/'workflow/pipeline_process.json')
    assert receipt['exitCode'] == pipeline['exitCode'] == 0
    assert all(receipt['inputsStable'].values()) and all(pipeline['inputsStable'].values())
    for name, digest in receipt['outputSHA256'].items():
        assert sha(folder/name) == digest
    saved = inventory(folder, receipt)
    assert [r['localIteration'] for r in saved['complete']] == [0,2,3] and saved['normalScheduleComplete']
    geometry = loadmat(str(output/'fixture_workbench/hybrid_volume.mat'))
    wall, faces = geometry['sourceVertices'], geometry['sourceFaces'].astype(int)-1
    xyz, values = vtk_tokens(folder/'surface_flow.vtk')
    order = np.array([int(np.argmin(np.linalg.norm(xyz-p,axis=1))) for p in wall])
    assert np.array_equal(xyz[order], wall) and len(set(order)) == len(wall)
    fields = np.genfromtxt(str(output/'main03/fields.csv'), delimiter=',', names=True)
    errors = {}
    for target, key in [('pressureRaw_Pa','Pressure'),('heatFluxRaw_W_m2','Heat_Flux'),
                        ('fluidNodeTemperatureMean_K','Temperature')]:
        expected = np.mean(values[key][order][faces,:,0], axis=1)
        errors[target] = relative(fields[target], expected)
    vectors = np.cross(wall[faces[:,1]]-wall[faces[:,0]],wall[faces[:,2]]-wall[faces[:,0]])/2
    area = np.linalg.norm(vectors,axis=1)
    errors['area'] = relative(fields['area_m2'], area)
    scalar = np.zeros(len(wall)); vector = np.zeros_like(wall)
    for ids, av, a in zip(faces, vectors, area):
        scalar[ids] += a/3; vector[ids] += av/3
    factor = np.linalg.norm(vector,axis=1)/scalar
    expected_load = np.mean((values['Heat_Flux'][order,0]*factor)[faces],axis=1)
    errors['nativeLoad'] = relative(fields['nativeLoadDensity_W_m2'],expected_load)
    assert max(errors.values()) < 1e-12
    assert np.isnan(fields['solvedWallTemperature_K']).all() and np.all(fields['givenWallTemperature_K']==300)
    assert np.all(fields['productionQualified']==0) and np.array_equal(fields['originalFaceID'],np.arange(1,9))
    with (output/'main03/section.csv').open(encoding='utf-8-sig') as stream:
        rows = list(csv.reader(stream))[1:]
    assert len(rows) == 4
    for row in rows:
        face = int(row[2])-1
        assert abs(float(row[0])-.25) < 1e-14
        assert abs(float(row[10])-fields['heatFluxRaw_W_m2'][face]) < 1e-8
        assert 'OTH' in row[12] and 'Roe' in row[12] and 'SST-2003m' in row[12]
        assert abs(np.linalg.norm(np.array(row[7:9],float)-np.array(row[5:7],float))-float(row[9])) < 1e-14
    snapshot = loadmat(str(folder/'input_snapshot.mat'), squeeze_me=True, struct_as_record=False)
    plan = output_plan(24,6,8,4,2,32*1024**2)
    assert snapshot['budget'].estimatedBytes == plan['estimatedMaximumBytes']
    assert receipt['thermodynamicModel']['checkpointContract'] == plan
    index = read(folder/'ideal_air_review.json'); thermal = read(folder/index['auditPath'])
    native = np.genfromtxt(str(folder/'restart_flow_000003.csv'), delimiter=',', names=True)
    oracle = np.genfromtxt(str(folder/index['thermoPath']), delimiter=',', names=True)
    speed2 = sum((native['Momentum_'+a]/native['Density'])**2 for a in ('x','y','z'))
    energy = native['Energy']/native['Density']-.5*speed2-native['Flow_Primitive_K']
    phase = relative(energy,oracle['e']); assert phase < 2e-10
    cancel = ui['cancellation']; cr = read(Path(cancel['folder'])/'engine_process.json')
    assert cr['cancelled'] and all(cr['inputsStable'].values())
    ci = inventory(Path(cancel['folder']),cr)
    assert [r['localIteration'] for r in ci['complete']] == [0,20]
    assert read(Path(cancel['folder'])/'workflow/pipeline_process.json')['exitCode'] == 125
    matlab = {}
    for name in ('main01','main02','main03','review01'):
        record = read(output/(name+'_process.json'))
        assert record['status'] == 'terminal' and all(record['sourceStableAfter'].values())
        matlab[name] = dict(exitCode=record['exitCode'],sourceFiles=len(record['sourceSHA256']),
                           startedUTC=record['startedUTC'],endedUTC=record['endedUTC'])
    assert not read(Path(read(output/'main02/engine_folder.json')['folder'])/'workflow/pipeline_process.json')['engineStarted']
    report = dict(passed=True, mainTask=str(folder), nativeSeconds=receipt['elapsedSeconds'],
        savedLocalIterations=[0,2,3], actualSpeed=ui['actualFreestreamSpeed'], originalFaces=8,sectionSegments=4,
        independentFieldErrors=errors, independentPhaseEnergyError=phase, thermalErrors=thermal['errors'],
        cancellation=dict(nativeExitCode=cr['exitCode'],pipelineExitCode=125, seconds=cr['elapsedSeconds'],
                          savedLocalIterations=[0,20],abnormalInstantStateAvailable=ci['abnormalInstantStateAvailable']),
        matlab=matlab, actualMouseObserved=True, newPhysicalQualification=0, originalResidualGatePassed=False)
    atomic_write_json(output/'delivery_report.json',report)
    print(json.dumps(report,ensure_ascii=False))


if __name__ == '__main__':
    main()
