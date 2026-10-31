"""独立用原生壁点和原三角面核对两个SST数值方案的工作台导出。"""
from pathlib import Path
import csv,json,sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from auditSU2SSTState import vtk_tokens
from prepareSU2CheckpointContinuation import read,sha


def main():
    base=ROOT/'outputs/v127_validation';view=base/'workbench01';target=base/'export_check.json';assert not target.exists()
    s=loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'),squeeze_me=True,struct_as_record=False)['s'];V=s.V;F=s.F.astype(int)-1
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    reports={};locks={}
    for name in ['control','cap05']:
        folder=base/name;output=view/name;xyz,native=vtk_tokens(folder/'surface_flow.vtk');d,ids=cKDTree(V).query(xyz)
        assert max(d)<1e-6 and len(set(ids))==len(V)
        with (output/'fields.csv').open() as f:rows=list(csv.DictReader(f))
        col=lambda key:np.array([float(row[key]) for row in rows])
        assert len(rows)==19752 and np.array_equal(col('originalFaceID'),np.arange(1,len(F)+1))
        errors={}
        for key,field in [('pressureRaw_Pa','Pressure'),('heatFluxRaw_W_m2','Heat_Flux'),('fluidNodeTemperatureMean_K','Temperature')]:
            n=np.empty(len(V));n[ids]=native[field][:,0];expected=n[F].mean(axis=1)
            error=float(max(abs(col(key)-expected)/np.maximum(abs(expected),1)));assert error<1e-12;errors[key]=error
        assert max(abs(col('area_m2')/area-1))<1e-12
        assert np.all(col('givenWallTemperature_K')==300) and np.all(np.isnan(col('solvedWallTemperature_K')))
        assert not np.any(col('productionQualified')) and not np.any(col('heatFluxIsBoundaryInput'))
        with (output/'section.csv').open(encoding='utf-8-sig') as f:section=list(csv.reader(f))
        assert len(section)==93;distance=0.;heatError=0.
        for row in section[1:]:
            x=float(row[0]);face=int(row[2])-1;assert x==-1.5
            a=np.array([x,float(row[5]),float(row[6])]);b=np.array([x,float(row[7]),float(row[8])]);tri=V[F[face]]
            assert abs(np.linalg.norm(a-b)-float(row[9]))<1e-12
            for point in [a,b]:
                matrix=(tri[1:]-tri[0]).T;uv=np.linalg.lstsq(matrix,point-tri[0],rcond=None)[0]
                error=float(np.linalg.norm(matrix@uv-(point-tri[0])));distance=max(distance,error)
                assert min(uv)>-1e-10 and sum(uv)<1+1e-10 and error<1e-10
            heatError=max(heatError,abs(float(row[10])-col('heatFluxRaw_W_m2')[face])/max(abs(float(row[10])),1))
            assert '同源数值方案' in row[12] and '未授予' in row[11]
            assert ('上限0.5候选' if name=='cap05' else '上限5控制') in row[12]
        assert heatError<1e-12
        for name2 in ['fields.csv.json','section.csv.json']:
            m=read(output/name2);a=m['convergenceAudit']
            assert a['checkpointIdentityConfirmed'] and a['checkpointIteration']==149 and a['sourceModelLocalIteration']==49
            assert a['numericalStudy']['variant']==name and not m['productionQualified']
            assert m['algorithmRuntime']['physicsModel']['numericalStudy']['variant']==name
        reports[name]=dict(originalFaces=len(F),originalAreaM2=float(sum(area)),sectionSegments=92,
            fieldRelativeErrors=errors,sectionHeatRelativeError=heatError,maximumEndpointTriangleDistanceM=distance)
        for p in [folder/'engine_process.json',output/'fields.csv',output/'section.csv',output/'fields.csv.json',output/'section.csv.json']:locks[str(p)]=sha(p)
    m=read(base/'sst_numerical_workbench01_process.json');assert m['status']=='terminal' and m['exitCode']==0 and all(m['sourceStableAfter'].values())
    ui=read(view/'report.json');assert ui['wrongWallBoundaryRejected'] and ui['oldSSTRestartCompatible'] and ui['solverRuns']==0
    report=dict(passed=True,variants=reports,physicsQualified=False,newQualifiedArea=0,solverRuns=0,
        sourceSHA256=dict(locks,**{str(Path(__file__).resolve()):sha(Path(__file__).resolve())}))
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(reports,indent=2))


if __name__=='__main__':main()
