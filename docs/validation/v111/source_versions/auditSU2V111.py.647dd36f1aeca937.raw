"""独立核对控制体构造、同初态重构/CFL配对及全部原面p/q，失败筛选仍保留。"""
from pathlib import Path
from types import SimpleNamespace
import hashlib
import json
import re
import sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from analyzeSU2DualVolumes import dual_volumes,metrics
from auditSU2V110 import vtk_tokens


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def geometry_tests():
    points=[np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float),
            np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,0,1],[0,1,1]],float)]
    for x,kind,total in zip(points,[12,13],[1,.5]):
        cells=np.empty(1,object);cells[0]=np.arange(len(x))+1
        v=SimpleNamespace(vertices=x,cells=cells,cellTypes=np.array([kind]),cellVolume=np.array([total]),fluidVolume=total)
        vol,_=dual_volumes(v);assert np.max(abs(vol-total/len(x)))<1e-14
    return True


def main():
    out=ROOT/'outputs/v111_validation';assert not (out/'comparison.json').exists();geometry_tests()
    s=loadmat(str(ROOT/'outputs/v108_validation/normal01/normal_volume.mat'),squeeze_me=True,struct_as_record=False);v=s['volume'];proj=s['projection']
    vol=np.load(str(out/'dual_volume.npy'));xyz,f=vtk_tokens(ROOT/'outputs/v110_validation/probe02/volume_flow.vtk')
    d,ids=cKDTree(v.vertices).query(xyz);assert max(d)<1e-6 and len(set(ids))==len(ids)
    order=np.argsort(ids);base={k:(a[:,0] if a.shape[1]==1 else a)[order] for k,a in f.items()}
    xyz,f=vtk_tokens(out/'first_probe01/volume_flow.vtk');d,ids=cKDTree(v.vertices).query(xyz);assert max(d)<1e-6 and len(set(ids))==len(ids)
    order=np.argsort(ids);first={k:(a[:,0] if a.shape[1]==1 else a)[order] for k,a in f.items()}
    nw=len(np.unique(np.concatenate(v.wall)));a=metrics(v.vertices,base,vol,nw);b=metrics(v.vertices,first,vol,nw)
    ratio={k:np.sqrt(np.array(b[k]['sum'])/a[k]['sum']).tolist() for k in a}
    def warning(path):
        text=path.read_text();m=re.findall(r'(\d+) reconstructed states for upwinding are non-physical',text)
        return int(m[-1]) if m else 0
    original=ROOT/'outputs/v109_validation/rocket_seed01';initial=np.loadtxt(str(original/'restart_flow.csv'),delimiter=',',skiprows=1)
    # 远场法向用原几何面面积向量合成，分类是未收敛初态的特征诊断。
    normals=np.zeros_like(v.vertices)
    for face in v.farfield:
        ids=np.asarray(face,int)-1;q=v.vertices[ids];area=np.sum(np.cross(q[1:-1]-q[0],q[2:]-q[0]),axis=0)/2
        normals[ids]+=area/len(ids)
    far=np.flatnonzero(np.linalg.norm(normals,axis=1)>0);normal=normals[far]/np.linalg.norm(normals[far],axis=1)[:,None]
    rho=initial[far,4];velocity=initial[far,5:8]/rho[:,None];e=initial[far,8]/rho-(velocity**2).sum(axis=1)/2
    temperature=e*.4/287.05287;mn=(velocity*normal).sum(axis=1)/np.sqrt(1.4*287.05287*temperature)
    boundary=dict(nodes=len(far),minNormalMach=float(mn.min()),maxNormalMach=float(mn.max()),
                  supersonicOut=int(np.sum(mn>=1)),subsonicOut=int(np.sum((mn>=0)&(mn<1))),
                  subsonicIn=int(np.sum((mn<0)&(mn>-1))),supersonicIn=int(np.sum(mn<=-1)))
    runs={};mapped={};paths=[Path(__file__).resolve(),ROOT/'scripts/analyzeSU2DualVolumes.py',out/'dual_volume.npy',original/'restart_flow.csv']
    for name in ('first_probe01','control01','cfl01'):
        folder=out/name;r=json.loads((folder/'engine_process.json').read_text());assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        for file,h in r['outputSHA256'].items():assert sha(folder/file)==h
        assert sha(folder/'initial_flow.csv')==sha(original/'restart_flow.csv')
        rest=np.loadtxt(str(folder/'restart_flow.csv'),delimiter=',',skiprows=1);rho=rest[:,4];vel=rest[:,5:8]/rho[:,None];e=rest[:,8]/rho-(vel**2).sum(axis=1)/2
        assert len(rest)==len(v.vertices) and np.all(np.isfinite(rest)) and min(rho)>0 and min(e)>0
        h=np.atleast_1d(np.genfromtxt(str(folder/'history.csv'),delimiter=',',names=True));assert np.array_equal(h['Inner_Iter'],np.arange(len(h)))
        local=int(h['Inner_Iter'][-1]);assert local==(0 if name=='first_probe01' else 599)
        assert sha(folder/('restart_flow_%06d.csv'%local))==sha(folder/'restart_flow.csv')
        assert sha(folder/('surface_flow_%06d.vtk'%local))==sha(folder/'surface_flow.vtk')
        x,w=vtk_tokens(folder/'surface_flow.vtk');dist,node=cKDTree(v.vertices).query(x);assert len(set(node))==len(node) and max(dist)<1e-6
        values=np.zeros((len(v.vertices),3));values[node]=np.column_stack([w[k][:,0] for k in ('Pressure','Heat_Flux','Temperature')]);face=proj.weights@values
        native=float(proj.engineNodeArea@values[:,1]);assert abs(native/h['HF'][-1]-1)<2e-6
        mapped[name]=face
        runs[name]=dict(elapsedSeconds=r['elapsedSeconds'],iterations=len(h),nativeCheckpointLocalIteration=local,
                        residuals={k:float(h[k][-1]) for k in ['rmsRho','rmsRhoU','rmsRhoV','rmsRhoW','rmsRhoE']},
                        nativeHeatLoad=native,physicalAreaHeatPower=float(proj.targetFaceArea@face[:,1]),
                        nonphysicalReconstructionWarningCount=warning(folder/'engine_stdout.log'),
                        minimumDensity=float(rho.min()),minimumInternalEnergy=float(e.min()),
                        finalCFL=dict(minimum=float(h['Min_CFL'][-1]),maximum=float(h['Max_CFL'][-1]),average=float(h['Avg_CFL'][-1])),
                        originalMappedFaces=len(face),physicsQualified=False)
        paths += [folder/n for n in ('engine_process.json','case.cfg','history.csv','surface_flow.vtk','restart_flow.csv')]
    c0=(out/'control01/case.cfg').read_text();c1=(out/'cfl01/case.cfg').read_text()
    assert c0.replace('CFL_ADAPT_PARAM= (0.5, 1.1, 0.05, 5.0)','CFL_ADAPT_PARAM= (0.5, 1.1, 0.05, 50.0)')==c1
    improvement={k:runs['control01']['residuals'][k]-runs['cfl01']['residuals'][k] for k in ('rmsRho','rmsRhoE')}
    face0=mapped['control01'];face1=mapped['cfl01'];area=proj.targetFaceArea
    differences={name:dict(areaWeightedRelativeL1=float(area@abs(face1[:,i]-face0[:,i])/(area@abs(face0[:,i]))),
                           maximumAbsolute=float(max(abs(face1[:,i]-face0[:,i])))) for i,name in enumerate(('pressure','heatFlux','fluidNodeTemperature'))}
    report=dict(analyticCubeAndPrismPassed=True,geometry=json.loads((out/'normalization.json').read_text())['geometry'],
                firstOrderToOriginalNormRatios=ratio,originalReconstructionWarningCount=warning(ROOT/'outputs/v110_validation/probe02/engine_stdout.log'),
                firstOrderWarmupRejected=True,farfieldInitialCharacteristicCounts=boundary,runs=runs,
                cflCap50ImprovementLog10=improvement,cflCap50ScreenPassed=all(x>=.3 for x in improvement.values()),
                cflCap50OriginalFaceDifferences=differences,defaultCFLChanged=False,physicsQualified=False,newQualifiedArea=0,
                inputSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in paths})
    (out/'comparison.json').write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps({k:v for k,v in report.items() if k!='inputSHA256'},indent=2))


if __name__=='__main__':main()
