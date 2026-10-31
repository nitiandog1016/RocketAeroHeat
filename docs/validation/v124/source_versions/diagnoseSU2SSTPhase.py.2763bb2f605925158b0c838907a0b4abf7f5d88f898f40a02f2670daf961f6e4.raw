"""只读已保存SST场：量化分步状态证据缺项，保留原热梯度失败，绝不重算/放行。"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.spatial import cKDTree
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from su2SSTInterface import state,wall_geometry
from auditSU2V110 import vtk_tokens
from prepareSU2CheckpointContinuation import sha,read


def main():
    out=ROOT/'outputs/v124_validation';target=out/'split_phase_diagnosis.json';assert not target.exists()
    p=Path(read(out/'workbench03/engine_folder.json')['folder']);r=read(p/'engine_process.json')
    assert r['exitCode']==0 and all(r['inputsStable'].values())
    for name in ['mesh.su2','initial_flow.csv','case.cfg']:assert sha(p/name)==r['inputSHA256'][str(p/name)]
    for name in ['restart_flow.csv','surface_flow.vtk','restart_flow_000050.csv','surface_flow_000050.vtk','history.csv']:
        assert sha(p/name)==r['outputSHA256'][name]
    g=wall_geometry(p/'mesh.su2');a,cols,nd,rho,u,e=state(p/'restart_flow.csv',True)
    xyz,f=vtk_tokens(p/'surface_flow.vtk');d,ids=cKDTree(a[:,1:4]).query(xyz);assert max(d)<1e-6
    order={v:i for i,v in enumerate(g['wallIDs'])};ix=np.array([order[j] for j in ids]);near=g['normalNeighbor'][ix]
    T=e*.4/287.05287;mu=1.458e-6*T[ids]**1.5/(T[ids]+110.4);cp=1.4*287.05287/.4
    q=mu*cp/.72*(T[near]-300)/g['thermalDistance'][ix];delta=q-f['Heat_Flux'][:,0]
    err=abs(delta)/np.maximum(abs(q),1);worst=np.argsort(err)[-10:]
    freestream=read(p/'case_manifest.json')['flight'];qinf=.5*freestream['pressure']/(287.05287*freestream['temperature'])*(5*np.sqrt(1.4*287.05287*freestream['temperature']))**2
    shear=np.linalg.norm(f['Skin_Friction_Coefficient'],axis=1)*qinf
    yp=g['nearestWallDistance'][ix]*np.sqrt(rho[ids]*shear)/mu
    yError=float(max(abs(yp-f['Y_Plus'][:,0])/np.maximum(abs(yp),1e-20)));assert yError<2e-6
    expectedOmega=60*mu/(rho[ids]*.075*g['nearestWallDistance'][ix]**2)
    s=loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'),squeeze_me=True,struct_as_record=False)['s'];V=s.V;F=s.F.astype(int)-1
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    _,vi=cKDTree(V).query(xyz);assert len(set(vi))==len(V)
    def faces(values):
        n=np.zeros(len(V));n[vi]=values;return n[F].mean(axis=1)
    face=faces(f['Heat_Flux'][:,0]);yface=faces(yp)
    midxyz,m=vtk_tokens(p/'surface_flow_000050.vtk');assert np.array_equal(midxyz,xyz);mid=faces(m['Heat_Flux'][:,0])
    hist=np.genfromtxt(str(p/'history.csv'),names=True,delimiter=',');native=float(g['nodeArea'][ix]@f['Heat_Flux'][:,0])
    assert abs(native-hist['HF'][-1])/max(abs(native),1)<2e-6
    closure=float(abs(native-g['nodeArea'][ix]@q)/max(abs(native),1))
    sources=[Path(__file__).resolve(),ROOT/'scripts/su2SSTInterface.py',p/'engine_process.json',p/'workflow/pipeline_process.json',
             ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src/iteration/CFluidIteration.cpp',
             ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src/integration/CMultiGridIntegration.cpp']
    report=dict(nativeExitCode=r['exitCode'],nativeElapsedSeconds=r['elapsedSeconds'],pipelineExitCode=read(p/'workflow/pipeline_process.json')['exitCode'],
        actualSourceLocalIteration=199,upstreamBeforeMeshTransferIteration=599,newModelLocalIteration=99,
        originalFaces=len(F),areaM2=float(sum(area)),nodes=len(a),nativeHeatPowerW=native,originalFaceHeatPowerW=float(area@face),
        finalResiduals=[float(hist[k][-1]) for k in ['rmsRho','rmsRhoE','rmsk','rmsw']],
        temperatureRangeK=[float(min(T)),float(max(T))],minimumDensity=float(min(rho)),minimumThermalInternalEnergy=float(min(e)),
        originalGradientThreshold=2e-6,maximumGradientRelativeError=float(max(err)),gradientGatePassed=False,
        maximumAbsoluteHeatDifferenceWm2=float(max(abs(delta))),sameStateReconstructedPowerRelativeDifference=closure,
        maximumEquivalentNeighborTemperatureDifferenceK=float(max(abs(delta*g['thermalDistance'][ix]/(mu*cp/.72)))),
        worst=[dict(nodeID=int(ids[j]),normalNeighborID=int(near[j]),qSaved=float(f['Heat_Flux'][j,0]),qFromSavedEnergyAndK=float(q[j]),
                    relativeError=float(err[j]),neighborK=float(a[near[j],-2])) for j in worst],
        yPlusRelativeError=yError,yPlusPercentiles=dict(zip(['p50','p95','p99','max'],np.percentile(yp,[50,95,99,100]).tolist())),
        faceMeanYPlusAbove1AreaFraction=float(sum(area[yface>1])/sum(area)),
        omegaWallRelativeDifference=float(max(abs(a[ids,-1]/expectedOmega-1))),
        withinRunHeatRelativeL1Change=float(area@abs(face-mid)/(area@abs(mid))),
        last25HeatPowerRelativeChange=float((hist['HF'][-1]-hist['HF'][-26])/abs(hist['HF'][-26])),
        savedLocalIterations=[50,99],previousTurbulence98Available=False,initialOutput0Available=False,
        sourceMechanism='flow update -> recompute primitive/heat -> turbulence update; same printed iteration does not guarantee same k phase',
        limitation='phase mismatch supported by fixed source ordering; missing predecessor98/precise native neighbor temperature prevents exact final heat reconstruction',
        fullResultAccepted=False,fullDomainInterfacePassed=False,physicsQualified=False,newQualifiedArea=0,
        sourceSHA256={str(p):sha(p) for p in sources})
    assert max(err)>=2e-6
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','worst']},indent=2))

if __name__=='__main__':main()
