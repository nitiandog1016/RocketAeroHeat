"""独立核对启动研究的原三角面、算法标签和截面导出；不启动求解。"""
import csv,json
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2SSTFluxStudy import ROOT,BASE,read,sha,write
from auditSU2SSTState import vtk_tokens

def main():
    target=BASE/'export_review.json';assert not target.exists()
    report=read(BASE/'workbench01/report.json');items={};locks={}
    mesh=loadmat(str(ROOT/'outputs/v145_validation/mesh/hybrid_volume.mat'),variable_names=['sourceVertices','sourceFaces'])
    V=mesh['sourceVertices'];F=mesh['sourceFaces'].astype(np.int64)-1
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    for name,ui in report['variants'].items():
        folder=BASE/name;export=BASE/'workbench01'/name
        r=read(folder/'engine_process.json');assert sha(folder/'surface_flow.vtk')==r['outputSHA256']['surface_flow.vtk']
        xyz,fields=vtk_tokens(folder/'surface_flow.vtk');d,ids=cKDTree(V).query(xyz)
        assert max(d)<1e-6 and len(set(ids))==len(V)
        with (export/'fields.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
        col=lambda n:np.asarray([float(row[n]) for row in rows]);assert len(rows)==len(F)==19752
        assert np.array_equal(col('originalFaceID'),np.arange(1,len(F)+1)) and max(abs(col('area_m2')/area-1))<1e-12
        errors={};faces={}
        for field,key in [('pressureRaw_Pa','Pressure'),('heatFluxRaw_W_m2','Heat_Flux'),('fluidNodeTemperatureMean_K','Temperature')]:
            nodes=np.empty(len(V));nodes[ids]=fields[key][:,0];faces[key]=nodes[F].mean(axis=1)
            errors[field]=float(max(abs(col(field)-faces[key])/np.maximum(abs(faces[key]),1)));assert errors[field]<1e-12
        assert np.all(col('givenWallTemperature_K')==300) and np.all(np.isnan(col('solvedWallTemperature_K')))
        assert not np.any(col('productionQualified')) and not np.any(col('heatFluxIsBoundaryInput'))
        with (export/'section.csv').open(encoding='utf-8-sig') as f:section=list(csv.reader(f))
        assert len(section)==93;maximum=0.
        phrase='AUSM+-up2' if name=='up2' else 'AUSM'
        for row in section[1:]:
            i=int(row[2])-1;assert float(row[0])==-1.5
            maximum=max(maximum,abs(float(row[10])-faces['Heat_Flux'][i])/max(abs(faces['Heat_Flux'][i]),1))
            assert 'SST-2003m' in row[12] and phrase in row[12]
            for point in [np.array([-1.5,float(row[5]),float(row[6])]),np.array([-1.5,float(row[7]),float(row[8])])]:
                tri=V[F[i]];uv=np.linalg.lstsq((tri[1:]-tri[0]).T,point-tri[0],rcond=None)[0]
                assert np.linalg.norm((tri[1:]-tri[0]).T@uv-(point-tri[0]))<1e-10 and min(uv)>-1e-10 and sum(uv)<1+1e-10
        assert maximum<1e-12 and ui['webgl']['webglError']==0
        items[name]=dict(fieldRelativeErrors=errors,sectionHeatRelativeError=maximum,algorithmIdentityVerified=True,originalFaces=len(F))
        for p in [folder/'engine_process.json',folder/'surface_flow.vtk',export/'fields.csv',export/'section.csv']:locks[str(p)]=sha(p)
    result=dict(passed=True,variants=items,solverRuns=0,physicsQualified=False,newQualifiedArea=0,sourceSHA256=locks)
    write(target,result);print(json.dumps(result['variants']))
if __name__=='__main__':main()
