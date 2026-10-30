"""逐值核对耦合研究视图、原面映射及真实浏览器下载；不授予壁面求解资格。"""

from pathlib import Path
import argparse
import csv
import json
import shutil
import numpy as np
from scipy.io import loadmat
from auditAxisymmetricShockLayerPatch import plain, graph
from auditThermalPressureSurfaceFlow import sha


def hit(V, F, point, direction):
    """模块一：独立射线—原三角面求交，保留并列交点中最小的原编号。"""
    v0=V[F[:,0]];e1=V[F[:,1]]-v0;e2=V[F[:,2]]-v0
    p=np.cross(np.broadcast_to(direction,e2.shape),e2);det=np.sum(e1*p,axis=1)
    usable=abs(det)>64*np.finfo(float).eps*np.linalg.norm(np.cross(e1,e2),axis=1)
    r=point-v0;q=np.cross(r,e1);u=np.full(len(F),np.nan);v=u.copy();t=u.copy()
    u[usable]=np.sum(r[usable]*p[usable],axis=1)/det[usable]
    v[usable]=(q[usable]@direction)/det[usable];t[usable]=np.sum(e2[usable]*q[usable],axis=1)/det[usable]
    good=usable&(u>=-1e-9)&(v>=-1e-9)&(u+v<=1+1e-9)&(t>128*np.finfo(float).eps*np.linalg.norm(np.ptp(V,axis=0)))
    if not np.any(good):return None
    nearest=np.min(t[good]);same=good&(abs(t-nearest)<=1e-9*np.linalg.norm(np.ptp(V,axis=0)))
    index=int(np.flatnonzero(same)[0]);return index+1,t[index],np.array([1-u[index]-v[index],u[index],v[index]])


def audit(evidence,folder,target,csv_file=None,browser_file=None,previous_file=None):
    """模块二：所有保存站点和三角面逐项核对，未知值不得变成数字零。"""
    assert not target.exists();s=plain(loadmat(evidence,simplify_cells=True))
    data=json.loads((folder/'shock-layer-data.json').read_text(encoding='utf-8'))
    assert data['kind']=='shock-layer-study-v1' and not data['productionQualified']
    assert not data['validation']['fullScopeValidation'] and not data['validation']['allIndependentIntegrationsPassed']
    assert data['sources']['coupledEvidence']['sha256']==sha(evidence)
    V=np.asarray(s['geometry']['V'],dtype=float);F=np.asarray(s['geometry']['F'],dtype=int)-1
    np.testing.assert_array_equal(data['geometry']['vertices'],V);np.testing.assert_array_equal(data['geometry']['faces'],F)
    focus=data['geometry']['focus'];sub=np.asarray(focus['vertices'])[np.asarray(focus['faces'],dtype=int)]
    np.testing.assert_array_equal(sub,V[F[np.asarray(focus['originalFaceIDs'],dtype=int)-1]])
    stations=0;ray_errors=[]
    for index,c in enumerate(data['cases']):
        mode=index//2;node=int(c['id'].split('_n')[-1]);where=int(np.flatnonzero(np.asarray(s['nodes'])==node)[0]);p=s['cases'][mode,where]
        keep=np.unique(np.r_[np.arange(0,len(p['eta']),20),len(p['eta'])-1]).astype(int)
        frame=p['shape']['frame'];z=np.asarray(p['radiusSquaredM2']);f,ig,_,kb,_=graph(p['shape'],z)
        for j,col in enumerate(c['columns']):
            np.testing.assert_array_equal(col['eta'],p['eta'][keep]);np.testing.assert_array_equal(col['position'],p['position'][keep,j,:])
            for name in ['normalDistanceM','pressurePa','temperatureK','densityKgM3','entropyJPerKgK','jacobianRatio']:
                np.testing.assert_array_equal(col[name],p[name][keep,j])
            np.testing.assert_allclose(col['speedMS'],np.hypot(p['streamwiseSpeedMS'][keep,j],p['normalSpeedMS'][keep,j]),rtol=3e-16,atol=2e-13)
            assert col['methods']==p['methods'] and not col['productionQualified']
            for field in ['heatFluxWPerM2','wallTemperatureInputK','wallTemperatureSolvedK','edgeTemperatureK','recoveryTemperatureK']:
                assert col[field]==[]
            origin=p['shape']['origin']+f[j]*frame[:,0]+np.sqrt(z[j])*frame[:,1]
            normal=ig[j]*frame[:,0]-np.sqrt(z[j])*kb[j]*frame[:,1];actual=hit(V,F,origin,normal)
            if actual is None:
                assert col['originalFaceID']==[] and col['bodyResidualM'] is None
            else:
                face,distance,bary=actual;assert col['originalFaceID']==face
                np.testing.assert_allclose(col['barycentric'],bary,rtol=0,atol=3e-11)
                error=abs(col['bodyResidualM']-(p['normalDistanceM'][-1,j]-distance));assert error<3e-12;ray_errors.append(float(error))
            stations+=len(keep)
    csv_rows=0
    if csv_file:
        rows=list(csv.DictReader(csv_file.open(encoding='utf-8-sig',newline='')));col=data['cases'][3]['columns'][7]
        assert len(rows)==len(col['eta'])
        mapping={'eta':'eta','normal_distance_m':'normalDistanceM','layer_pressure_Pa':'pressurePa','layer_temperature_K':'temperatureK',
                 'layer_density_kg_m3':'densityKgM3','layer_speed_m_s':'speedMS','entropy_J_kg_K':'entropyJPerKgK','normal_jacobian_ratio':'jacobianRatio'}
        for i,row in enumerate(rows):
            for key,name in mapping.items():assert float(row[key])==col[name][i]
            for axis,name in enumerate(['x_m','y_m','z_m']):assert float(row[name])==col['position'][i][axis]
            for key in ['heat_flux_W_m2','wall_temperature_input_K','wall_temperature_solved_K','edge_temperature_K','recovery_temperature_K']:assert row[key]==''
            assert row['pressure_method']==col['methods']['pressure'] and row['heat_flux_status']=='unknown';csv_rows+=1
        shutil.copyfile(str(csv_file),str(target.parent/'browser/actual_download_column8.csv'))
    if browser_file:
        observed=json.loads(browser_file.read_text(encoding='utf-8'));assert observed['caseValue']=='3' and observed['columnValue']=='7' and observed['stationValue']=='8'
        assert '原三角面号\n18849' in observed['methods'] and '热流 / W/m²\n未求解' in observed['probe']
    if previous_file:
        previous=json.loads(previous_file.read_text(encoding='utf-8'));assert previous['cases']==data['cases']
    report=dict(status='passed-exact-coupled-view-data-original-faces-and-actual-CSV-audit',evidenceSHA256=sha(evidence),
                dataSHA256=sha(folder/'shock-layer-data.json'),stations=stations,originalFaces=len(F),
                originalRayMappings=len(ray_errors),maximumRayDistanceDifferenceM=max(ray_errors),actualCSVrows=csv_rows,
                view01AndView02PhysicalNumbersUnchanged=bool(previous_file),productionQualified=False,newProductionSolvedAreaM2=0)
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('evidence',type=Path);parser.add_argument('folder',type=Path);parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--csv',type=Path);parser.add_argument('--browser',type=Path)
    parser.add_argument('--previous',type=Path)
    a=parser.parse_args();audit(a.evidence,a.folder,a.report,a.csv,a.browser,a.previous)
