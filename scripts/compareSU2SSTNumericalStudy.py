"""从具名原生场和完整历史审核v127事前门槛，筛选失败不自动重跑。"""
from pathlib import Path
import json,csv,re
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2CheckpointContinuation import ROOT,read,sha
from auditSU2SSTState import vtk_tokens


def main():
    base=ROOT/'outputs/v127_validation';target=base/'comparison.json';assert not target.exists()
    surface=loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'),squeeze_me=True,struct_as_record=False)['s']
    V=surface.V;F=surface.F.astype(int)-1;area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    x=V[F].mean(axis=1)[:,0];bins=np.linspace(min(V[:,0]),max(V[:,0]),9)
    variants={};heat={};locks={}
    residuals=['rms[Rho]','rms[RhoE]','rms[k]','rms[w]']
    for name in ['control','cap05']:
        folder=base/name;r=read(folder/'engine_process.json');audit=read(folder/'sst_state_audit.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        assert r['physicsModel']['numericalStudy']['variant']==name and audit['interfacePassed']
        with (folder/'history.csv').open() as f:
            reader=csv.DictReader(f,skipinitialspace=True);h={k.strip():[] for k in reader.fieldnames}
            for row in reader:
                for k,v in row.items():h[k.strip()].append(float(v))
        h={k:np.array(v) for k,v in h.items()};assert np.array_equal(h['Inner_Iter'],np.arange(150))
        required=residuals+['Min CFL','Max CFL','Avg CFL','Min DT','Max DT','LinSolRes','LinSolResTurb','Linear_Solver_Iterations','LinSolIterTurb','HF']
        assert all(np.isfinite(h[k]).all() for k in required)
        xy,q=vtk_tokens(folder/'surface_flow.vtk');dist,ids=cKDTree(V).query(xy);assert max(dist)<1e-6
        nodal=np.empty(len(V));nodal[ids]=q['Heat_Flux'][:,0];heat[name]=nodal[F].mean(axis=1)
        fallbacks=list(map(int,re.findall(r'Warning: (\d+) reconstructed states',(folder/'engine_stdout.log').read_text())))
        assert fallbacks
        regions=[]
        for j in range(8):
            mask=(x>=bins[j])&(x<(bins[j+1] if j<7 else np.inf))
            regions.append(dict(xIntervalM=bins[j:j+2].tolist(),areaM2=float(sum(area[mask])),
                meanHeatWm2=float(area[mask]@heat[name][mask]/sum(area[mask]))))
        variants[name]=dict(elapsedSeconds=r['elapsedSeconds'],initialSHA256=sha(folder/'initial_flow.csv'),
            last25MeanLogResidual={k:float(np.mean(h[k][-25:])) for k in residuals},
            firstLogResidual={k:float(h[k][0]) for k in residuals},lastLogResidual={k:float(h[k][-1]) for k in residuals},
            maximumActualCFL=float(max(h['Max CFL'])),finalCFL={k:float(h[k][-1]) for k in ['Min CFL','Max CFL','Avg CFL']},
            maximumCFLDropIterations=(np.where(h['Max CFL'][1:]<.5*h['Max CFL'][:-1])[0]+1).tolist(),
            sampledHistory={str(i):{k:float(h[k][i]) for k in residuals+['Max CFL','Avg CFL','LinSolRes','LinSolResTurb','HF']} for i in [0,16,25,49,75,99,100,124,149]},
            timestepRangeS=[float(min(h['Min DT'])),float(max(h['Max DT']))],
            worstLinearLogResidual={k:float(max(h[k])) for k in ['LinSolRes','LinSolResTurb']},
            last25MeanLinearLogResidual={k:float(np.mean(h[k][-25:])) for k in ['LinSolRes','LinSolResTurb']},
            linearIterationMaximum={k:int(max(h[k])) for k in ['Linear_Solver_Iterations','LinSolIterTurb']},
            fallbackWarningCounts=fallbacks,finalFallbackCount=fallbacks[-1],
            heatRelativeL1Change=audit['withinRunHeatRelativeL1Change'],last25HeatPowerRelativeChange=audit['last25HeatPowerRelativeChange'],
            nativeHeatPowerW=audit['nativeHeatPowerW'],originalFaceHeatPowerW=float(area@heat[name]),
            negativeHeatAreaFraction=float(sum(area[heat[name]<0])/sum(area)),axialRegions=regions,
            yPlusP99=audit['yPlusPercentiles']['p99'],originalResidualPassed=audit['originalResidualThresholdPassed'])
        for p in [folder/'engine_process.json',folder/'history.csv',folder/'sst_state_audit.json',folder/'surface_flow.vtk']:
            locks[str(p)]=sha(p)
    c=variants['control'];n=variants['cap05'];assert c['initialSHA256']==n['initialSHA256']
    gains={k:c['last25MeanLogResidual'][k]-n['last25MeanLogResidual'][k] for k in residuals}
    gates=dict(actualCFLBound=n['maximumActualCFL']<=.5+1e-12,
        densityMeanImproved03=gains['rms[Rho]']>=.3,energyMeanImproved03=gains['rms[RhoE]']>=.3,
        turbulenceNotWorse03=min(gains[k] for k in residuals[2:])>=-.3,
        heatL1Below10Percent=n['heatRelativeL1Change']<=.1,heatPowerBelow2Percent=abs(n['last25HeatPowerRelativeChange'])<=.02,
        yPlusP99Below1=n['yPlusP99']<=1,fallbackBelow110Percent=n['finalFallbackCount']<=1.1*c['finalFallbackCount'])
    report=dict(variants=variants,meanResidualGainsLog=gains,screeningGates=gates,allScreeningGatesPassed=all(gates.values()),
        pairedFaceHeatRelativeL1Difference=float(area@abs(heat['cap05']-heat['control'])/(area@abs(heat['control']))),
        originalFaces=len(F),originalAreaM2=float(sum(area)),physicsQualified=False,newQualifiedArea=0,
        defaultChanged=False,automaticContinuationAllowed=False,sourceSHA256=dict(locks,**{str(Path(__file__).resolve()):sha(Path(__file__).resolve())}))
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['variants','sourceSHA256']},indent=2))


if __name__=='__main__':main()
