"""独立核对两网格原面/原生梯度/保存身份及运行前筛选，不赋予物理解资格。"""
from pathlib import Path
import json,hashlib,csv,sys,re
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v120_validation'
sys.path.insert(0,str(ROOT/'tests'));from auditSU2V110 import vtk_tokens
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def main():
    target=OUT/'independent_audit.json';assert not target.exists()
    records={};surfaces={};paths=[Path(__file__).resolve(),ROOT/'tests/auditSU2V110.py',OUT/'initial_audit.json']
    for name in ['control','candidate']:
        folder=OUT/'study01'/name;ui=OUT/'review01';receipt=read(folder/'engine_process.json')
        assert receipt['status']=='terminal' and receipt['exitCode']==0 and not receipt['timedOut'] and not receipt['cancelled']
        assert all(receipt['inputsStable'].values()) and receipt['request']['threads']==2
        for p,h in receipt['inputSHA256'].items():assert sha(Path(p))==h,p
        for p,h in receipt['outputSHA256'].items():assert sha(folder/p)==h,p
        r=loadmat(str(ui/(name+'_result.mat')),squeeze_me=True,struct_as_record=False)['r']
        V=r.V;F=r.F.astype(int)-1;n=len(V)
        areaVector=np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]])/2
        area=np.linalg.norm(areaVector,axis=1);normal=np.zeros_like(V)
        for j in range(3):np.add.at(normal,F[:,j],areaVector/3)
        assert len(F)==19752 and abs(sum(area)-1.9646292565715553)<1e-12
        def wall(path):
            xyz,fields=vtk_tokens(path);d,ids=cKDTree(V).query(xyz)
            assert max(d)<1e-6 and len(set(ids))==n
            values=np.zeros((n,3));values[ids]=np.column_stack([fields[k][:,0] for k in ('Pressure','Heat_Flux','Temperature')])
            return values,values[F].mean(axis=1),fields,ids
        node,face,f,ids=wall(folder/'surface_flow.vtk');_,mid,_,_=wall(folder/'surface_flow_000100.vtk');surfaces[name]=face
        hist=np.atleast_1d(np.genfromtxt(str(folder/'history.csv'),names=True,delimiter=','))
        assert np.array_equal(hist['Inner_Iter'],np.arange(200))
        for basename,suffix in [('restart_flow','csv'),('surface_flow','vtk')]:
            assert sha(folder/(basename+'.'+suffix))==sha(folder/(basename+'_000199.'+suffix))
        a=r.convergenceAudit;assert a.checkpointIdentityConfirmed and a.checkpointIteration==199 and not r.productionQualified
        if name=='control':assert a.cumulativeCheckpointIteration==799
        else:assert not hasattr(a,'cumulativeCheckpointIteration') and a.sourceCumulativeIteration==599
        meshPath=ROOT/'outputs/v118_validation/actual01/mesh/hybrid_volume.mat' if name=='control' else OUT/'mesh32/hybrid_volume.mat'
        mesh=loadmat(str(meshPath));points=mesh['vertices'];assert np.array_equal(points[:n],V)
        rest=np.loadtxt(str(folder/'restart_flow.csv'),delimiter=',',skiprows=1)
        assert np.array_equal(rest[:,0],np.arange(len(points))) and max(abs(rest[:,1:4]-points).ravel())<=1e-14
        rho=rest[:,4];vel=rest[:,5:8]/rho[:,None];e=rest[:,8]/rho-.5*np.sum(vel*vel,axis=1)
        assert np.isfinite(rest).all() and min(rho)>0 and min(e)>0
        T=e*.4/287.05287;p=.4*rho*e
        eos=max(float(max(abs(p[:n]/node[:,0]-1))),float(max(abs(T[:n]/node[:,2]-1))));assert eos<2e-6
        # 按真实拓扑选最大法向余弦邻点，并重算官方热梯度的法向投影距离。
        neighbors=[{i+n} for i in range(n)]
        for i,j,k in F:neighbors[i].update((j,k));neighbors[j].update((i,k));neighbors[k].update((i,j))
        near=[]
        for i,js in enumerate(neighbors):
            js=np.array(sorted(js));delta=points[js]-V[i]
            near.append(js[np.argmax(delta@normal[i]/np.linalg.norm(delta,axis=1))])
        near=np.asarray(near);assert np.array_equal(near,np.arange(n,2*n))
        unit=normal/np.linalg.norm(normal,axis=1)[:,None];delta=points[near]-V
        length=np.linalg.norm(delta,axis=1);distance=abs(np.sum(delta*unit,axis=1));assert min(distance/length)>.05
        mu=1.458e-6*T[:n]**1.5/(T[:n]+110.4)
        q=mu*(1.4*287.05287/.4)/.72*(T[near]-300)/distance
        gradientError=float(max(abs(q-node[:,1])/np.maximum(1,abs(q))));assert gradientError<2e-6
        flight=read(folder/'case_manifest.json')['flight'];rhoinf=flight['pressure']/(287.05287*flight['temperature'])
        speed=flight['Mach']*np.sqrt(1.4*287.05287*flight['temperature']);qinf=.5*rhoinf*speed**2
        shear=np.linalg.norm(f['Skin_Friction_Coefficient'],axis=1)*qinf
        estimate=length[ids]*np.sqrt(f['Density'][:,0]*shear)/f['Laminar_Viscosity'][:,0]
        yp=np.zeros(n);yp[ids]=estimate;assert np.all(f['Y_Plus']==0)
        native=float(np.linalg.norm(normal,axis=1)@node[:,1]);assert abs(native-hist['HF'][-1])/max(abs(native),1)<2e-6
        table=np.genfromtxt(str(ui/(name+'_fields.csv')),names=True,delimiter=',')
        assert np.array_equal(table['originalFaceID'],np.arange(1,len(F)+1))
        for j,(col,key) in enumerate(zip(['pressureRaw_Pa','heatFluxRaw_W_m2','fluidNodeTemperatureMean_K'],['pressureRaw','heatFluxRaw','fluidNodeTemperatureMeanRaw'])):
            assert np.allclose(face[:,j],table[col],rtol=1e-12,atol=1e-8) and np.allclose(face[:,j],getattr(r,key),rtol=1e-12,atol=1e-8)
        assert (table['givenWallTemperature_K']==300).all() and np.isnan(table['solvedWallTemperature_K']).all() and not table['productionQualified'].any()
        with (ui/(name+'_section.csv')).open(encoding='utf-8-sig') as stream:rows=list(csv.reader(stream))[1:]
        for row in rows:
            assert abs(float(row[0])+1.5)<1e-14 and abs(float(row[10])-face[int(row[2])-1,1])<1e-8*max(1,abs(float(row[10])))
        log=(folder/'engine_stdout.log').read_text();assert 'All volume elements are correctly oriented.' in log and 'All surface elements are correctly oriented.' in log
        warnings=re.findall(r'Warning: (\d+) reconstructed states',log)
        records[name]=dict(engineSeconds=receipt['elapsedSeconds'],nodes=len(points),originalFaces=len(F),areaM2=float(sum(area)),
            localCheckpointIteration=199,cumulativeCheckpointIteration=799 if name=='control' else None,sourceCumulativeIteration=599,
            finalResiduals=[float(hist[k][-1]) for k in ['rmsRho','rmsRhoE']],originalResidualThresholdPassed=bool(hist['rmsRho'][-1]<-8 and hist['rmsRhoE'][-1]<-8),
            withinRunHeatRelativeL1Change=float(area@abs(face[:,1]-mid[:,1])/(area@abs(mid[:,1]))),
            withinRunPressureRelativeL1Change=float(area@abs(face[:,0]-mid[:,0])/(area@abs(mid[:,0]))),
            last50PowerRelativeChange=float((hist['HF'][-1]-hist['HF'][-50])/max(abs(hist['HF'][-50]),1)),
            minimumDensity=float(min(rho)),minimumInternalEnergy=float(min(e)),temperatureRangeK=[float(min(T)),float(max(T))],
            nativeHeatPowerW=native,originalFaceHeatPowerW=float(area@face[:,1]),wallEOSRelativeError=eos,
            nativeGradientMaximumRelativeError=gradientError,firstNormalDistanceRangeM=[float(min(length)),float(max(length))],
            estimatedWallUnitsPercentiles=dict(zip(['p50','p95','p99','max'],np.percentile(estimate,[50,95,99,100]).tolist())),
            faceMeanWallUnitAbove1AreaFraction=float(sum(area[yp[F].mean(axis=1)>1])/sum(area)),
            nativeYPlusAllZero=True,finalReconstructionWarnings=int(warnings[-1]) if warnings else 0,sectionRows=len(rows),physicsQualified=False,newQualifiedArea=0)
        paths += [folder/'engine_process.json',meshPath,ui/(name+'_result.mat')]
    a=surfaces['control'];b=surfaces['candidate']
    qdiff=float(area@abs(b[:,1]-a[:,1])/(area@abs(a[:,1])));pdiff=float(area@abs(b[:,0]-a[:,0])/(area@abs(a[:,0])))
    gates=dict(geometryPositiveStateProvenanceOriginalArea=True,candidateWallUnitP99Below1=records['candidate']['estimatedWallUnitsPercentiles']['p99']<1,
        candidateEnergyLogNoMoreThanControlPlus03=records['candidate']['finalResiduals'][1]<=records['control']['finalResiduals'][1]+.3,
        candidateWithinRunHeatL1Below10Percent=records['candidate']['withinRunHeatRelativeL1Change']<.1,pairedFinalHeatL1Below10Percent=qdiff<.1)
    report=dict(evidenceChecksPassed=True,records=records,pairedFinalHeatRelativeL1=qdiff,pairedFinalPressureRelativeL1=pdiff,
        predeclaredScreeningGates=gates,allScreeningGatesPassed=all(gates.values()),productionDefaultPromoted=False,
        physicsQualified=False,newQualifiedArea=0,wallUnitEstimateMeaning='unconverged laminar wall shear and actual first normal distance; not validated turbulent y+',
        meshComparisonMeaning='two unconverged fields; not grid independence or discretization error',
        sourceSHA256={str(p):sha(p) for p in paths})
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='sourceSHA256'},indent=2))
if __name__=='__main__':main()
