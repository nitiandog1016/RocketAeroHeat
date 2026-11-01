"""按事前门槛比较固定通量；不同离散残差只作推进筛选，不授予热流资格。"""
from pathlib import Path
import json,numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2SSTFluxStudy import ROOT,BASE,PLAN,read,sha,write,config,validate_source,source_path
from compareSU2SSTStartupStudy import snapshot,RESIDUALS
from auditSU2SSTState import vtk_tokens

def gates(source,control,candidate):
    gain={k:control['last10MeanLogResidual'][k]-candidate['last10MeanLogResidual'][k] for k in RESIDUALS}
    checks=dict(interfaces=control['interfacePassed'] and candidate['interfacePassed'],
        densityImproved03=gain['rms[Rho]']>=.3,energyImproved03=gain['rms[RhoE]']>=.3,
        enthalpyExcessNoWorseThanSource=candidate['enthalpyExcess']<=source['enthalpyExcess'],
        enthalpyExcessReducedVersusControl=candidate['enthalpyExcess']<=.9*control['enthalpyExcess'],
        heatL1NotWorse=candidate['heatRelativeL1Change']<=control['heatRelativeL1Change'],
        heatPowerChangeNotWorse=abs(candidate['last25HeatPowerRelativeChange'])<=abs(control['last25HeatPowerRelativeChange']),
        noNewOmegaFloor=candidate['omegaFloor']<=source['omegaFloor'])
    return checks,gain

def main():
    target=BASE/'flux_comparison.json';assert not target.exists()
    geometry=ROOT/'outputs/v145_validation/mesh/hybrid_volume.mat'
    mesh=loadmat(str(geometry),variable_names=['sourceVertices','sourceFaces']);V=mesh['sourceVertices'];F=mesh['sourceFaces'].astype(np.int64)-1
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    assert len(F)==19752 and abs(sum(area)-1.9646292565715553)<1e-12
    variants={};faces={};locks={};source=source_path('control')
    for label,folder in [('source99',source),('control',BASE/'control'),('up2',BASE/'up2')]:
        if label!='source99':
            validate_source(folder,config(folder/'case.cfg'),read(folder/'sst_initialization.json'))
            assert read(folder/'engine_process.json')['physicsModel']['numericalStudy']['variant']==label
        v=snapshot(folder,len(V));assert v['localFinalIteration']==(99 if label=='source99' else 49)
        v['last10MeanLogResidual']={k:float(np.mean(v['sevenResidualHistory'][k][-10:])) for k in RESIDUALS}
        v['flowFlux']=config(folder/'case.cfg')['CONV_NUM_METHOD_FLOW'];variants[label]=v
        xyz,fields=vtk_tokens(folder/'surface_flow.vtk');d,ids=cKDTree(V).query(xyz)
        assert max(d)<1e-6 and len(ids)==len(set(ids))==len(V)
        ff={}
        for key in ['Pressure','Heat_Flux','Temperature']:
            nodes=np.empty(len(V));nodes[ids]=fields[key][:,0];ff[key]=nodes[F].mean(axis=1)
        faces[label]=ff;v['originalFaceHeatPowerW']=float(area@ff['Heat_Flux'])
        bins=np.linspace(min(V[:,0]),max(V[:,0]),9);xc=V[F].mean(axis=1)[:,0];regions=[]
        for i in range(8):
            mask=(xc>=bins[i])&(xc<(bins[i+1] if i<7 else np.inf))
            regions.append(dict(xIntervalM=bins[i:i+2].tolist(),areaM2=float(sum(area[mask])),heatPowerW=float(area[mask]@ff['Heat_Flux'][mask])))
        v['axialAreaBins']=regions
        for name in ['engine_process.json','sst_state_audit.json','case.cfg','history.csv','restart_flow.csv','surface_flow.vtk']:
            locks[str(folder/name)]=sha(folder/name)
    assert sha(BASE/'control/initial_flow.csv')==sha(BASE/'up2/initial_flow.csv')==sha(source/'restart_flow.csv')
    checks,gain=gates(variants['source99'],variants['control'],variants['up2'])
    physical={label:dict(originalResiduals=v['originalResidualPassed'],heatL1Below10Percent=v['heatRelativeL1Change']<.1,
        heatPowerBelow2Percent=abs(v['last25HeatPowerRelativeChange'])<.02,yPlusP99Below1=v['yPlusP99']<=1) for label,v in variants.items()}
    result=dict(variants=variants,screeningGates=checks,last10MeanLogResidualGains=gain,
        candidateEligibleForFurtherStudy=all(checks.values()),unchangedPhysicalGates=physical,
        physicsQualified=False,newQualifiedArea=0,defaultChanged=False,newSpatialRHSRecorded=False,
        originalFaces=len(F),originalAreaM2=float(sum(area)),sourceLocalIteration=99,newLocalIteration=49,
        pairedFaceRelativeL1Difference={key:float(area@abs(faces['up2'][key]-faces['control'][key])/(area@abs(faces['control'][key]))) for key in faces['control']},
        sourceSHA256=dict(locks,**{str(p):sha(p) for p in [PLAN,Path(__file__).resolve(),geometry,ROOT/'scripts/compareSU2SSTStartupStudy.py']}))
    write(target,result);print(json.dumps(dict(gates=checks,gain=gain,eligible=result['candidateEligibleForFurtherStudy'],physical=physical)))

if __name__=='__main__':main()
