"""核对闭合经典箭体的条件流线、原始站点、可视化负载与输入/求解身份。"""
from pathlib import Path
import hashlib
import argparse
import json
import re
import sys

import numpy as np
from scipy.io import loadmat

from auditPressureDrivenSurfaceFlow import audit_case, rows, sha

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('folder');parser.add_argument('process',nargs='?')
    parser.add_argument('--report',help='另存新复核报告，不覆盖历史数值或审核记录。');args=parser.parse_args()
    folder=Path(args.folder).resolve();folder.relative_to(ROOT)
    target=Path(args.report).resolve() if args.report else folder/'independent_surface_flow_study.json'
    target.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
    assert not target.exists();source=folder/'pressure_surface_study_evidence.mat';d=loadmat(source,simplify_cells=True)
    web=folder/'workbench';raw=(web/'pressure-surface-flow-data.json').read_bytes();payload=json.loads(raw)
    assert 'href="pressure-surface-flow-data.json"' in (web/'index.html').read_text(encoding='utf8')
    assert (web/'pressure-surface-flow-data.js').read_bytes()==b'window.PRESSURE_SURFACE_FLOW_STUDY='+raw+b';\n'
    assert payload['sourceEvidenceSHA256']==sha(source) and not payload['productionQualified']
    geometry_raw=(web/'geometry/rocket-geometry-data.js').read_text(encoding='utf8')
    prefix='window.ROCKET_GEOMETRY_CASES=';assert geometry_raw.startswith(prefix) and geometry_raw.endswith(';')
    geometry=json.loads(geometry_raw[len(prefix):-1]);checks=[];total_stations=0;total_segments=0
    for saved,case,model in zip(rows(d['saved']),payload['cases'],geometry):
        assert case['id']==saved['kind'] and case['pressureIsInput'] and case['thermodynamicStateIsInput']
        assert not case['heatFluxSolved'] and not case['wallTemperatureSolved'] and not case['productionQualified']
        assert case['newProductionSolvedAreaM2']==0
        V=np.asarray(saved['mesh']['V']);F=np.asarray(saved['mesh']['F'],dtype=int)-1
        np.testing.assert_array_equal(model['vertices'],V);np.testing.assert_array_equal(model['faces'],F)
        edges=np.vstack([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]);unique,ids,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
        orientation=np.bincount(ids,weights=np.sign(edges[:,1]-edges[:,0]));volume=np.sum(V[F[:,0]]*np.cross(V[F[:,1]],V[F[:,2]]))/6
        if saved['kind']!='vortex':
            assert np.all(count==2) and np.all(orientation==0) and volume>0
            assert len(V)-len(unique)+len(F)==2
            assert case['geometry']['closed'] and case['geometry']['boundaryEdges']==case['geometry']['nonManifoldEdges']==0
        else:
            assert not case['geometry']['closed'] and np.any(count==1)
        overlay=model['surfaceFlowLines'];assert overlay['kind']=='pressure-surface-flow-research' and not overlay['productionQualified']
        results=rows(saved['results']);assert len(results)==len(case['tubes'])==len(case['curves'])==len(overlay['curves'])
        for j,(result,tube,curve) in enumerate(zip(results,case['tubes'],case['curves'])):
            checks.append(audit_case(result,saved['specification']));assert np.all(result['validStation'])
            center=result['center'];np.testing.assert_array_equal(tube['arcLengthM'],result['arcLengthM'])
            for name in ['metric','metricRatio','metricDifferenceIndicator']:
                np.testing.assert_array_equal(tube[name],result[name])
            for name in ['position','tangent']:
                np.testing.assert_array_equal(tube[name],center[name])
            for key,source_key in [('pressureInputPa','pressure'),('temperatureInputK','temperature'),('densityInputKgM3','density'),('speedInputMS','speed')]:
                np.testing.assert_array_equal(tube[key],center[source_key])
            assert curve['sourceTubeIndexOneBased']==j+1
            np.testing.assert_array_equal(curve['points'],center['position']);np.testing.assert_array_equal(curve['normals'],center['normal'])
            assert overlay['curves'][j]==curve
            total_stations+=len(tube['arcLengthM']);total_segments+=len(curve['points'])-1
        assert case['summary']['validStations']==sum(len(t['arcLengthM']) for t in case['tubes'])
        assert not any(name in case for name in ['heatFlux','wallSolvedTemperature'])
    assert len(checks)==31 and len(geometry)==4
    # 按真正文件字节核对所有缓存URL，不把旧嵌套资源当成新版本。
    for page in [web/'index.html',web/'geometry/index.html']:
        for name,short in re.findall(r'(?:src|href)="([^"?]+)\?v=([0-9a-f]+)"',page.read_text(encoding='utf8')):
            assert sha(page.parent/name).startswith(short),(page,name)
    name,short=payload['geometryURL'].split('?v=');assert sha(web/name).startswith(short)
    process_path=Path(args.process).resolve() if args.process else folder.parent/'surfaceflowview94_process.json'
    process=json.loads(process_path.read_text())
    assert process['status']=='terminal' and process['exitCode']==0 and all(process['sourceStableAfter'].values())
    export=json.loads((folder/'study_export.json').read_text());assert export['legacyGeometryThreeArgumentContractPassed']
    files=[source,process_path,folder/'study_export.json',web/'pressure-surface-flow-data.json',web/'geometry/rocket-geometry-data.js']
    source_files=['tests/studyPressureDrivenSurfaceFlow.m','tests/pressureSurfaceBenchmarkProblem.m','exportRocketTestGeometryView.m',
                  'regionalFlow/surface/marchPressureDrivenSurfaceFlow.m','regionalFlow/surface/pressureSurfaceFlowRHS.m','regionalFlow/surface/pressureSurfaceFlowOptions.m']
    for name in source_files:
        assert sha(ROOT/name)==process['sourceSHA256'][name];files.append(ROOT/name)
    report=dict(status='passed-independent-classical-body-topology-all-surface-rays-input-roles-and-web-binding',
                models=4,closedClassicalBodies=3,conditionalTubes=len(checks),conditionalRays=5*len(checks),
                actualStations=total_stations,visibleCurveSegments=total_segments,checks=checks,
                pressureAndThermalFieldsAreInputs=True,heatFluxAndSolvedWallTemperatureUnknown=True,
                inputSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in files},scriptSHA256=sha(__file__),
                referenceHelperSHA256=sha(ROOT/'tests/auditPressureDrivenSurfaceFlow.py'),
                independentExperimentalValidation=False,productionQualified=False,newProductionSolvedAreaM2=0)
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:report[k] for k in ['status','models','closedClassicalBodies','conditionalTubes','conditionalRays','actualStations','visibleCurveSegments']}))


if __name__=='__main__':main()
