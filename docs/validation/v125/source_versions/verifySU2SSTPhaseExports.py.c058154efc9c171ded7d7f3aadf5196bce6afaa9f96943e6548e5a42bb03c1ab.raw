"""从原生壁点和原三角面独立复算v125 CSV/交线；不调用MATLAB映射实现。"""
from pathlib import Path
import csv,json,sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from auditSU2SSTState import vtk_tokens
from prepareSU2CheckpointContinuation import sha,read


def main():
    base=ROOT/'outputs/v125_validation';view=base/'view02';case=base/'phase01'
    s=loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'),squeeze_me=True,struct_as_record=False)['s']
    V=s.V;F=s.F.astype(int)-1
    av=np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]])/2;area=np.linalg.norm(av,axis=1)
    xyz,native=vtk_tokens(case/'surface_flow.vtk');d,ids=cKDTree(V).query(xyz)
    assert max(d)<1e-6 and len(set(ids))==len(V)
    with (view/'fields.csv').open() as f:rows=list(csv.DictReader(f))
    col=lambda key:np.array([float(r[key]) for r in rows])
    assert len(rows)==19752 and np.array_equal(col('originalFaceID'),np.arange(1,len(F)+1))
    errors={}
    for output,key in [('pressureRaw_Pa','Pressure'),('heatFluxRaw_W_m2','Heat_Flux'),('fluidNodeTemperatureMean_K','Temperature')]:
        n=np.empty(len(V));n[ids]=native[key][:,0];expected=n[F].mean(axis=1)
        error=float(max(abs(col(output)-expected)/np.maximum(abs(expected),1)))
        assert error<1e-12;errors[output]=error
    assert max(abs(col('area_m2')/area-1))<1e-12
    assert np.all(col('givenWallTemperature_K')==300) and np.all(np.isnan(col('solvedWallTemperature_K')))
    assert not np.any(col('productionQualified')) and not np.any(col('heatFluxIsBoundaryInput'))
    with (view/'section.csv').open(encoding='utf-8-sig') as f:section=list(csv.reader(f))
    assert len(section)==93
    endpointError=0.;qerr=0.
    for row in section[1:]:
        x=float(row[0]);face=int(row[2])-1;assert x==-1.5
        a=np.array([x,float(row[5]),float(row[6])]);b=np.array([x,float(row[7]),float(row[8])]);tri=V[F[face]]
        assert abs(np.linalg.norm(a-b)-float(row[9]))<1e-12
        for point in [a,b]:
            matrix=(tri[1:]-tri[0]).T;uv=np.linalg.lstsq(matrix,point-tri[0],rcond=None)[0]
            dist=np.linalg.norm(matrix@uv-(point-tri[0]));endpointError=max(endpointError,float(dist))
            assert min(uv)>-1e-10 and sum(uv)<1+1e-10 and dist<1e-10
        qerr=max(qerr,abs(float(row[10])-col('heatFluxRaw_W_m2')[face])/max(abs(float(row[10])),1))
        assert 'SST-2003m' in row[12] and '阶段诊断' in row[12] and '未授予' in row[11]
    assert qerr<1e-12
    for name in ['fields.csv.json','section.csv.json']:
        m=read(view/name);a=m['convergenceAudit']
        assert a['checkpointIdentityConfirmed'] and a['checkpointIteration']==0 and a['sourceModelLocalIteration']==99
        assert not a['checkpointResidualPassed'] and a['turbulenceDiagnostic']['sourcePostAuditPassed'] is False
        assert '阶段诊断' in m['algorithmRuntime']['method'] and not m['productionQualified']
    phases=read(base/'phase_diagnosis.json');audit=read(case/'sst_state_audit.json');ui=read(view/'report.json')
    assert phases['phaseProvedForThisUpdate'] and not phases['nativeHighPrecisionTemperatureWritten']
    assert audit['diagnosticOnly'] and audit['last25HeatPowerRelativeChange'] is None
    assert not audit['screeningGates']['last25HeatPowerBelow2Percent'] and ui['oldFailedFieldStillRejected']
    report=dict(passed=True,originalFaces=len(F),originalAreaM2=float(sum(area)),sectionSegments=92,
        rawFieldRelativeErrors=errors,sectionHeatRelativeError=qerr,maximumEndpointTriangleDistanceM=endpointError,
        exactSavedOutputSHA256={str(p.relative_to(ROOT)):sha(p) for p in [case/'engine_process.json',view/'fields.csv',view/'section.csv',view/'fields.csv.json',view/'section.csv.json']},
        sourceSHA256={str(Path(__file__).resolve()):sha(Path(__file__).resolve())},
        physicsQualified=False,newQualifiedArea=0,solverRuns=0)
    target=base/'export_check.json';assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if not k.endswith('SHA256')},indent=2))

if __name__=='__main__':main()
