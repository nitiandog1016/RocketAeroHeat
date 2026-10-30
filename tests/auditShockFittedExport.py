"""逐值核对实际MAT到三维工作台的映射，未知区域不得被数学测试值填充。"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.io import loadmat
from auditThermalPressureSurfaceFlow import sha

ROOT=Path(__file__).resolve().parents[1]


def audit(evidence,folder,target):
    assert not target.exists()
    s=loadmat(evidence,simplify_cells=True);d=json.loads((folder/'ale-data.json').read_text(encoding='utf8'))
    assert d['source']['evidenceSHA256']==sha(evidence)
    assert not d['productionQualified'] and d['newProductionSolvedAreaM2']==0
    assert (folder/'ale-data.js').read_bytes()==b'window.ALE_STUDY='+(folder/'ale-data.json').read_bytes()+b';\n'
    np.testing.assert_array_equal(d['originalGeometry']['vertices'],s['geometry']['V'])
    np.testing.assert_array_equal(d['originalGeometry']['faces'],s['geometry']['F']-1)
    p=s['pistons'][0];mesh=p['mesh'];bound=np.asarray(mesh['boundaryFaces'],int)-1;owner=np.asarray(mesh['owner'],int)[bound]-1
    benchmark=d['cases'][0];body=d['cases'][1];fields=['pressurePa','temperatureK','speedMS','wallPressurePa','heatFluxWm2','wallTemperatureK']
    np.testing.assert_array_equal(benchmark['faces'],mesh['faces'][bound]-1)
    np.testing.assert_array_equal(benchmark['boundaryKind'],p['boundary']['kind'])
    count=0
    for i,frame in enumerate(benchmark['frames']):
        np.testing.assert_array_equal(frame['position'],p['nodes0' if i==0 else 'nodes1'])
        assert frame['timeS']==(0 if i==0 else p['times'][-1])
        actual=p['jump']['state'] if i==0 else p['lastStep']['state']
        expected={'pressurePa':np.full(len(bound),actual['pressure']) if i==0 else actual['pressure'][owner],
                  'temperatureK':np.full(len(bound),actual['temperature']) if i==0 else actual['temperature'][owner],
                  'speedMS':np.full(len(bound),np.linalg.norm(actual['velocity'])) if i==0 else np.linalg.norm(actual['velocity'][owner],axis=1),
                  'wallPressurePa':np.full(len(bound),np.nan) if i==0 else p['lastStep']['wallPressurePa'],
                  'heatFluxWm2':np.full(len(bound),np.nan),'wallTemperatureK':np.full(len(bound),np.nan)}
        for key in fields:
            numeric=np.array([np.nan if v is None else v for v in frame[key]])
            np.testing.assert_allclose(numeric,expected[key],rtol=3e-15,atol=0,equal_nan=True);count+=len(numeric)
    shell=s['bodyCheck']['shell'];mesh=shell['mesh'];bound=np.asarray(mesh['boundaryFaces'],int)-1
    np.testing.assert_array_equal(body['faces'],mesh['faces'][bound]-1)
    np.testing.assert_array_equal(body['sourceFaceIDs'],shell['cellOriginalFaceIDs'][np.asarray(mesh['owner'],int)[bound]-1])
    for i,frame in enumerate(body['frames']):
        np.testing.assert_array_equal(frame['position'],shell['nodes'] if i==0 else s['bodyCheck']['nodes1'])
        for key in fields:assert all(v is None for v in frame[key]);count+=len(frame[key])
    assets={}
    for name in ['ale-study.js','ale-study.css']:
        assert (folder/name).read_bytes()==(ROOT/'web'/name).read_bytes();assets[name]=sha(folder/name)
    for name in ['geometry.js','geometry.css','ale-study-geometry.js']:
        assert (folder/'geometry'/name).read_bytes()==(ROOT/'web'/name).read_bytes();assets['geometry/'+name]=sha(folder/'geometry'/name)
    result=dict(status='passed-all-exported-fields-original-VF-face-bindings-and-unknowns',evidenceSHA256=sha(evidence),
        dataSHA256=sha(folder/'ale-data.json'),checkedFieldValues=count,originalFaces=len(d['originalGeometry']['faces']),
        currentUIAssetSHA256=assets,productionQualified=False,newProductionSolvedAreaM2=0)
    target.write_bytes((json.dumps(result,ensure_ascii=False,indent=2)+'\n').encode('utf8'));print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);p.add_argument('folder',type=Path);p.add_argument('target',type=Path)
    a=p.parse_args();audit(a.evidence,a.folder,a.target)
