"""独立核对实际NASA7结果到三维工作台的逐值映射、热状态身份及原始资源。"""
from pathlib import Path
import argparse
import json
import re

import numpy as np
from scipy.io import loadmat

from auditThermalPressureSurfaceFlow import sha, rows

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('folder');parser.add_argument('evidence');parser.add_argument('process')
    parser.add_argument('--reference-view');parser.add_argument('--report');args=parser.parse_args()
    folder=Path(args.folder).resolve();folder.relative_to(ROOT)
    evidence=Path(args.evidence).resolve();evidence.relative_to(ROOT)
    target=Path(args.report).resolve() if args.report else folder/'independent_thermal_surface_flow_study.json'
    target.relative_to(ROOT);assert not target.exists();target.parent.mkdir(parents=True,exist_ok=True)
    data=loadmat(evidence,simplify_cells=True);cases=rows(data['cases']);assert data['report']['cases']==6
    web=folder/'workbench';raw=(web/'pressure-surface-flow-data.json').read_bytes();payload=json.loads(raw)
    assert (web/'pressure-surface-flow-data.js').read_bytes()==b'window.PRESSURE_SURFACE_FLOW_STUDY='+raw+b';\n'
    assert payload['sourceEvidenceSHA256']==sha(evidence) and not payload['productionQualified']
    geom_text=(web/'geometry/rocket-geometry-data.js').read_text(encoding='utf8');prefix='window.ROCKET_GEOMETRY_CASES='
    assert geom_text.startswith(prefix) and geom_text.endswith(';');geometry=json.loads(geom_text[len(prefix):-1])
    mappings=[[2],[3],[4,5]];stations=0;segments=0;checks=[]
    for i,(case,model,source_ids) in enumerate(zip(payload['cases'],geometry,mappings)):
        assert case['pressureIsInput'] and not case['thermodynamicStateIsInput']
        assert case['thermodynamicStateIsConditionallyRestored'] and case['thermodynamicConvention']=='continuous-project-cp-v1'
        assert not case['heatFluxSolved'] and not case['wallTemperatureSolved'] and not case['productionQualified']
        assert case['newProductionSolvedAreaM2']==0 and '反算' in case['methodScope']
        assert model['surfaceFlowLines']['thermodynamicStateIsConditionallyRestored']
        assert '条件反算' in model['coordinateNote'] and 'Tr/Tw/热流未求' in model['coordinateNote']
        V=np.asarray(model['vertices']);F=np.asarray(model['faces'],dtype=int)
        assert V.shape[1]==F.shape[1]==3 and np.all(np.isfinite(V)) and F.min()>=0 and F.max()<len(V)
        edges=np.vstack([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]);unique,ids,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
        oriented=np.bincount(ids,weights=np.sign(edges[:,1]-edges[:,0]));volume=np.sum(V[F[:,0]]*np.cross(V[F[:,1]],V[F[:,2]]))/6
        if i==2:
            assert np.all(count==2) and np.all(oriented==0) and volume>0 and len(V)-len(unique)+len(F)==2
            assert case['geometry']['closed'] and abs(V[:,0].max())<1e-14
            assert '程序内部' in case['pressureSourceRole'] and 'Taylor' in case['pressureSourceRole']
        else:
            assert not case['geometry']['closed'] and np.all(V[:,2]==0)
            assert '验证输入' in case['pressureSourceRole']
        assert len(case['tubes'])==len(case['curves'])==len(source_ids)
        assert model['surfaceFlowLines']['curves']==case['curves']
        for j,source_id in enumerate(source_ids):
            source=cases[source_id]['result'];tube=case['tubes'][j];curve=case['curves'][j]
            assert np.all(source['validStation']) and tube['labelCenter']==source['labelCenter']
            for name in ['arcLengthM','metric','metricRatio','metricDifferenceIndicator']:
                np.testing.assert_array_equal(tube[name],source[name])
            for name in ['position','tangent']:np.testing.assert_array_equal(tube[name],source['center'][name])
            for name,key in [('pressureInputPa','pressure'),('temperatureRecoveredK','temperature'),('densityRecoveredKgM3','density'),('speedRecoveredMS','speed')]:
                np.testing.assert_array_equal(tube[name],source['center'][key])
            assert not any(k in tube for k in ['temperatureInputK','densityInputKgM3','speedInputMS','heatFlux','wallSolvedTemperature'])
            np.testing.assert_array_equal(curve['points'],source['center']['position']);np.testing.assert_array_equal(curve['normals'],source['center']['normal'])
            assert curve['sourceTubeIndexOneBased']==j+1
            if i==2:
                p=np.asarray(tube['position']);r=np.linalg.norm(p[:,1:],axis=1)
                np.testing.assert_allclose(r,-p[:,0]*np.tan(np.pi/12),rtol=0,atol=2e-12)
                assert np.all(p[:,0]>-.1/np.tan(np.pi/12))
            stations+=len(tube['arcLengthM']);segments+=len(tube['arcLengthM'])-1
        assert case['summary']['validStations']==sum(len(t['arcLengthM']) for t in case['tubes'])
        overlay=model['surfaceFlowLines'];az=float(overlay['initialCameraAzimuthRad']);el=float(overlay['initialCameraElevationRad'])
        assert model['initialCameraAzimuthRad']==az and model['initialCameraElevationRad']==el
        direction=np.array([np.cos(el)*np.cos(az),np.cos(el)*np.sin(az),np.sin(el)])
        visible=[int(np.count_nonzero(np.asarray(curve['normals'])@direction>0)) for curve in case['curves']]
        assert all(n>0 for n in visible),'默认视角不能把全部可选流管藏在背面。'
        checks.append(dict(model=case['id'],sourceCases=source_ids,closed=bool(case['geometry']['closed']),
                           tubes=len(source_ids),initialFrontFacingStationsPerTube=visible))
    assert len(geometry)==len(payload['cases'])==3 and stations==124 and segments==120
    for page in [web/'index.html',web/'geometry/index.html']:
        for name,stamp in re.findall(r'(?:src|href)="([^"?]+)\?v=([0-9a-f]+)"',page.read_text(encoding='utf8')):
            assert sha(page.parent/name).startswith(stamp)
    name,stamp=payload['geometryURL'].split('?v=');assert sha(web/name).startswith(stamp)
    assert 'href="pressure-surface-flow-data.json"' in (web/'index.html').read_text(encoding='utf8')
    process_path=Path(args.process).resolve();process=json.loads(process_path.read_text(encoding='utf-8-sig'))
    assert process['status']=='terminal' and process['exitCode']==0 and all(process['sourceStableAfter'].values())
    assert evidence.relative_to(ROOT).as_posix() in process['explicitInputs']
    bindings=[evidence,process_path,web/'pressure-surface-flow-data.json',web/'geometry/rocket-geometry-data.js',folder/'thermal_study_export.json']
    for rel in ['regionalFlow/surface/exportThermalPressureSurfaceFlowStudy.m','exportRocketTestGeometryView.m',
                'web/pressure-surface-study.js','web/pressure-surface-study.html','web/geometry.js']:
        assert sha(ROOT/rel)==process['sourceSHA256'][rel];bindings.append(ROOT/rel)
    unchanged=False
    if args.reference_view:
        previous=Path(args.reference_view).resolve();previous.relative_to(ROOT)
        older=json.loads((previous/'workbench/pressure-surface-flow-data.json').read_bytes())
        assert older['cases']==payload['cases'] and older['sourceEvidenceSHA256']==payload['sourceEvidenceSHA256']
        bindings.append(previous/'workbench/pressure-surface-flow-data.json');unchanged=True
    report=dict(status='passed-independent-NASA7-state-roles-actual-tube-mapping-topology-and-assets',models=3,
                conditionalTubes=4,actualStations=stations,visibleCurveSegments=segments,checks=checks,
                previousNumericalAndCasePayloadUnchanged=unchanged,restoredTeNotMislabelledAsTrOrTw=True,
                inputSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in bindings},scriptSHA256=sha(__file__),
                browserOrGPUVerified=False,independentExperimentalValidation=False,newProductionSolvedAreaM2=0)
    target.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps({k:report[k] for k in ['status','models','conditionalTubes','actualStations','previousNumericalAndCasePayloadUnchanged']}))


if __name__=='__main__':main()
