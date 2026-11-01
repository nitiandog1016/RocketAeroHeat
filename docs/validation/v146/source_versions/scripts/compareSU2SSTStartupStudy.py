"""审核同源一阶启动对照；筛选只决定是否恢复二阶，绝不授予物理资格。"""
import csv, json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2SSTStartupStudy import ROOT, BASE, PLAN, original_source, read, sha, config, write, validate_source
from su2SSTInterface import state
from auditSU2SSTState import vtk_tokens

RESIDUALS=['rms[Rho]','rms[RhoE]','rms[k]','rms[w]']

def restore_gates(source, control, candidate):
    """事前固定的启动接受条件，原物理门槛在报告中另列且不降低。"""
    gain={k:control['last20MeanLogResidual'][k]-candidate['last20MeanLogResidual'][k] for k in RESIDUALS}
    gates=dict(interfaces=control['interfacePassed'] and candidate['interfacePassed'],
        energyImproved03=gain['rms[RhoE]']>=.3,
        otherResidualsNotWorse03=min(gain[k] for k in ['rms[Rho]','rms[k]','rms[w]'])>=-.3,
        kFloorHalvedVersusSource=candidate['kFloorOffWall']<=.5*source['kFloorOffWall'],
        kFloorHalvedVersusControl=candidate['kFloorOffWall']<=.5*control['kFloorOffWall'],
        noOmegaFloor=candidate['omegaFloor']==0,
        enthalpyExcessReducedVersusSource=candidate['enthalpyExcess']<=.75*source['enthalpyExcess'],
        enthalpyExcessReducedVersusControl=candidate['enthalpyExcess']<=.75*control['enthalpyExcess'],
        heatL1NotWorse=candidate['heatRelativeL1Change']<=control['heatRelativeL1Change'],
        heatPowerChangeNotWorse=abs(candidate['last25HeatPowerRelativeChange'])<=abs(control['last25HeatPowerRelativeChange']))
    return gates,gain

def snapshot(folder, wall_count):
    """读取全部原生节点；流动primitive总焓不混入更新后的k。"""
    r=read(folder/'engine_process.json');audit=read(folder/'sst_state_audit.json')
    assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
    assert audit['interfacePassed'] and sha(folder/'restart_flow.csv')==r['outputSHA256']['restart_flow.csv']
    a,keys,nd,_,_,_=state(folder/'restart_flow.csv',True,True);assert nd==3
    k=a[:,keys.index('Turb_Kin_Energy')];w=a[:,keys.index('Omega')];T=a[:,keys.index('Temperature')]
    velocity=np.column_stack([a[:,keys.index('Velocity_'+q)] for q in 'xyz'])
    c=config(folder/'case.cfg');gamma=float(c['GAMMA_VALUE']);R=float(c['GAS_CONSTANT']);Ti=float(c['FREESTREAM_TEMPERATURE'])
    cp=gamma*R/(gamma-1);h=cp*T+.5*np.sum(velocity**2,axis=1)
    hinf=cp*Ti+.5*float(c['MACH_NUMBER'])**2*gamma*R*Ti
    kfloor=np.flatnonzero(k[wall_count:]<=1e-10*(1+8*np.finfo(float).eps))+wall_count
    ext=int(np.argmax(h));tj=int(np.argmax(T))
    with (folder/'history.csv').open() as f:rows=list(csv.reader(f,skipinitialspace=True))
    hist=dict(zip([x.strip() for x in rows[0]],np.asarray(rows[1:],float).T))
    assert np.array_equal(hist['Inner_Iter'],np.arange(len(rows)-1))
    seven={key:values.tolist() for key,values in hist.items() if key.startswith('rms[')};assert len(seven)==7
    return dict(localFinalIteration=int(hist['Inner_Iter'][-1]),elapsedSeconds=r['elapsedSeconds'],
        interfacePassed=audit['interfacePassed'],kFloorOffWall=len(kfloor),
        kFloorXRangeM=[float(min(a[kfloor,1])),float(max(a[kfloor,1]))] if len(kfloor) else None,
        omegaFloor=int(sum(w<=1e-4*(1+8*np.finfo(float).eps))),
        maxTemperatureK=float(T[tj]),maxTemperatureXYZM=a[tj,1:4].tolist(),
        maximumPrimitiveEnthalpyRatio=float(h[ext]/hinf),enthalpyExcess=max(float(h[ext]/hinf)-1,0),
        enthalpyExtremeNodeZeroBased=ext,enthalpyExtremeXYZM=a[ext,1:4].tolist(),
        last20MeanLogResidual={key:float(np.mean(hist[key][-20:])) for key in RESIDUALS},
        finalLogResidual={key:float(hist[key][-1]) for key in RESIDUALS},sevenResidualHistory=seven,
        additionalHistory={key:v.tolist() for key,v in hist.items() if key in ['HF','Max CFL','Avg CFL','Min DT','Max DT','LinSolRes','LinSolResTurb']},
        heatRelativeL1Change=audit.get('withinRunHeatRelativeL1Change'),last25HeatPowerRelativeChange=audit['last25HeatPowerRelativeChange'],
        yPlusP99=audit['yPlusPercentiles']['p99'],originalResidualPassed=audit['originalResidualThresholdPassed'],
        nativeHeatPowerW=audit['nativeHeatPowerW'])

def main():
    target=BASE/'startup_comparison.json';assert not target.exists()
    mesh=loadmat(str(ROOT/'outputs/v145_validation/mesh/hybrid_volume.mat'),variable_names=['sourceVertices','sourceFaces'])
    V=mesh['sourceVertices'];F=mesh['sourceFaces'].astype(np.int64)-1
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    assert len(F)==19752 and abs(sum(area)-1.9646292565715553)<1e-12
    variants={};faces={};locks={};source=original_source()
    for label,folder in [('source149',source),('control',BASE/'control'),('first_order',BASE/'first_order')]:
        if label!='source149':validate_source(folder,config(folder/'case.cfg'),read(folder/'sst_initialization.json'))
        variants[label]=snapshot(folder,len(V))
        assert variants[label]['localFinalIteration']==(149 if label=='source149' else 99)
        xyz,fields=vtk_tokens(folder/'surface_flow.vtk');d,ids=cKDTree(V).query(xyz)
        assert max(d)<1e-6 and len(ids)==len(V) and len(set(ids))==len(V)
        ff={}
        for key in ['Pressure','Heat_Flux','Temperature']:
            nodes=np.empty(len(V));nodes[ids]=fields[key][:,0];ff[key]=nodes[F].mean(axis=1)
        faces[label]=ff
        variants[label]['originalFaceHeatPowerW']=float(area@ff['Heat_Flux'])
        bins=np.linspace(min(V[:,0]),max(V[:,0]),9);xc=V[F].mean(axis=1)[:,0];regions=[]
        for i in range(8):
            mask=(xc>=bins[i])&(xc<(bins[i+1] if i<7 else np.inf))
            regions.append(dict(xIntervalM=bins[i:i+2].tolist(),areaM2=float(sum(area[mask])),
                heatPowerW=float(area[mask]@ff['Heat_Flux'][mask])))
        variants[label]['axialAreaBins']=regions
        for name in ['engine_process.json','sst_state_audit.json','case.cfg','history.csv','restart_flow.csv','surface_flow.vtk']:
            locks[str(folder/name)]=sha(folder/name)
    assert sha(BASE/'control/initial_flow.csv')==sha(BASE/'first_order/initial_flow.csv')==sha(source/'restart_flow.csv')
    gates,gain=restore_gates(variants['source149'],variants['control'],variants['first_order'])
    physical={label:dict(originalResiduals=v['originalResidualPassed'],heatL1Below10Percent=v['heatRelativeL1Change'] is not None and v['heatRelativeL1Change']<.1,
        heatPowerBelow2Percent=abs(v['last25HeatPowerRelativeChange'])<.02,yPlusP99Below1=v['yPlusP99']<=1)
        for label,v in variants.items()}
    result=dict(variants=variants,screeningGates=gates,meanLogResidualGains=gain,restoreAllowed=all(gates.values()),
        unchangedPhysicalGates=physical,physicsQualified=False,newQualifiedArea=0,defaultChanged=False,
        originalFaces=len(F),originalAreaM2=float(sum(area)),sourceLocalIteration=149,newLocalIteration=99,
        pairedFaceRelativeL1Difference={key:float(area@abs(faces['first_order'][key]-faces['control'][key])/(area@abs(faces['control'][key]))) for key in faces['control']},
        sourceSHA256=dict(locks,**{str(p):sha(p) for p in [PLAN,Path(__file__).resolve()]}))
    write(target,result)
    print(json.dumps(dict(gates=gates,gain=gain,restoreAllowed=result['restoreAllowed'],physical=physical),indent=2))

if __name__=='__main__':main()
