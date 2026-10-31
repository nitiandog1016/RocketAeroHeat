"""独立核对v107同一MAT、三份界面CSV、面积/热负荷及实际进程终态。"""
from pathlib import Path
from scipy.io import loadmat
import numpy as np
import csv
import json
import hashlib
import argparse

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    args = parser.parse_args()
    source = Path(args.source).resolve()
    output = ROOT / 'outputs/v107_validation'
    result_file = source / 'mapped_research.mat'
    r = loadmat(str(result_file), squeeze_me=True, struct_as_record=False)['result']
    assert not r.productionQualified and not r.eligibleForSurfaceCoverage
    n = len(r.faceArea)
    assert n == 2448
    triangles = r.V[r.F.astype(int) - 1]
    area = np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0], triangles[:,2]-triangles[:,0]),axis=1)/2
    assert np.allclose(area, r.faceArea, rtol=1e-13, atol=1e-18)
    comparisons = dict(area_m2='faceArea', pressureRaw_Pa='pressureRaw', heatFluxRaw_W_m2='heatFluxRaw',
                       nativeLoadDensity_W_m2='conservativeHeatLoadDensityRaw',
                       fluidNodeTemperatureMean_K='fluidNodeTemperatureMeanRaw')
    errors = {}
    inputs = [result_file, source/'history.csv']
    for name in ('workbench01','workbench02','replay01'):
        file = output/name/'research_snapshot.csv'
        inputs += [file, Path(str(file)+'.json')]
        with file.open(encoding='utf-8-sig') as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows)==n
        for column, field in comparisons.items():
            values = np.array([float(row[column]) for row in rows])
            expected = getattr(r, field)
            error = float(np.max(np.abs(values-expected))/max(1,np.max(np.abs(expected))))
            assert error<1e-12, (file,column,error)
            errors[name+'/'+column]=error
        assert [int(row['originalFaceID']) for row in rows] == list(range(1,n+1))
        assert all(float(row['givenWallTemperature_K'])==r.wallTemperatureInput for row in rows)
        assert all(row['solvedWallTemperature_K'] in ('NaN','nan','') for row in rows)
        assert all(row['productionQualified'] in ('0','false','False') for row in rows)
        meta=json.loads(Path(str(file)+'.json').read_text(encoding='utf8'))
        assert not meta['productionQualified'] and not meta['eligibleForSurfaceCoverage']
        assert meta['receipt']['request']['jobToken'] == r.receipt.request.jobToken
        assert meta['heatFluxAndHeatLoadAreDistinct']
    raw_power=float(area @ r.heatFluxRaw)
    native_power=float(area @ r.conservativeHeatLoadDensityRaw)
    with (source/'history.csv').open() as stream:
        rows=list(csv.reader(stream,skipinitialspace=True))
    history={key.strip():float(value) for key,value in zip(rows[0],rows[-1])}
    assert abs(native_power-history['HF'])/abs(history['HF'])<2e-8
    assert history['Inner_Iter']==1 and history['rms[RhoE]']>0  # 明确未收敛
    assert abs(raw_power-native_power)>1 and abs(raw_power-r.targetHeatPowerRaw)<1e-8
    matlab=[]
    for file in sorted(output.glob('*_process.json')):
        state=json.loads(file.read_text())
        assert state['status']=='terminal' and all(state['sourceStableAfter'].values())
        matlab.append(dict(name=file.name,exitCode=state['exitCode']))
        inputs.append(file)
    assert len(matlab)==4 and sum(item['exitCode']==0 for item in matlab)==3
    report=dict(originalFaces=n,areaM2=float(area.sum()),csvRelativeErrors=errors,matlabJobs=matlab,
                actualIterationIndex=history['Inner_Iter'],energyLogResidual=history['rms[RhoE]'],
                rawPhysicalAreaPowerW=raw_power,nativeLoadPowerW=native_power,solverMonitorW=history['HF'],
                areaConventionDifferenceFraction=(raw_power-native_power)/native_power,
                physicsQualified=False,newQualifiedRocketArea=0,
                inputSHA256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    dest=ROOT/'docs/validation/v107/independent_audit.json'
    dest.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({key:report[key] for key in ('originalFaces','areaM2','rawPhysicalAreaPowerW','nativeLoadPowerW','newQualifiedRocketArea')}))


if __name__=='__main__':
    main()
