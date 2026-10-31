"""只读v126具名场与历史：分开全局稳定筛选、残差反弹和表面变化来源。"""
from pathlib import Path
import sys,json,csv,re
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from auditSU2SSTState import vtk_tokens


def main():
    base=ROOT/'outputs/v126_validation';folder=Path(read(base/'main01/engine_folder.json')['folder'])
    receipt=read(folder/'engine_process.json');assert receipt['status']=='terminal' and receipt['exitCode']==0
    with (folder/'history.csv').open() as f:
        reader=csv.reader(f,skipinitialspace=True);keys=[s.strip() for s in next(reader)]
        h=np.array([[float(v) for v in row] for row in reader])
    assert np.array_equal(h[:,keys.index('Inner_Iter')],np.arange(50))
    columns=['rms[Rho]','rms[RhoE]','rms[k]','rms[w]','HF']
    records={str(i):{k:float(h[i,keys.index(k)]) for k in columns} for i in [0,24,25,27,49]}
    s=loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'),squeeze_me=True,struct_as_record=False)['s']
    V=s.V;F=s.F.astype(int)-1;centers=V[F].mean(axis=1)
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    q=[]
    for step in [25,49]:
        p=folder/('surface_flow_%06d.vtk'%step);assert sha(p)==receipt['outputSHA256'][p.name]
        xyz,fields=vtk_tokens(p);d,ids=cKDTree(V).query(xyz);assert max(d)<1e-6 and len(set(ids))==len(V)
        nodal=np.empty(len(V));nodal[ids]=fields['Heat_Flux'][:,0];q.append(nodal[F].mean(axis=1))
    change=area*abs(q[1]-q[0]);bins=np.linspace(min(V[:,0]),max(V[:,0]),9);regions=[]
    for j in range(8):
        mask=(centers[:,0]>=bins[j])&(centers[:,0]<(bins[j+1] if j<7 else np.inf))
        regions.append(dict(xIntervalM=bins[j:j+2].tolist(),faces=int(sum(mask)),areaM2=float(sum(area[mask])),
            fractionOfAbsoluteHeatChange=float(sum(change[mask])/sum(change)),
            meanHeat25Wm2=float(area[mask]@q[0][mask]/sum(area[mask])),
            meanHeat49Wm2=float(area[mask]@q[1][mask]/sum(area[mask]))))
    fallbacks=re.findall(r'Warning: (\d+) reconstructed states', (folder/'engine_stdout.log').read_text())
    c=config(folder/'case.cfg');a=read(base/'actual_surface_audit.json')
    report=dict(history=records,energyMinimumLocalIteration=int(np.argmin(h[:,keys.index('rms[RhoE]')])),
        energyResidualReboundLog=float(h[-1,keys.index('rms[RhoE]')]-min(h[:,keys.index('rms[RhoE]')])),
        densityResidualReboundLog=float(h[-1,keys.index('rms[Rho]')]-min(h[:,keys.index('rms[Rho]')])),
        fallbackWarningCounts=[int(x) for x in fallbacks],fallbackSpatialLocationsAvailable=False,
        actualCFLHistoryAvailable=False,configuredCFL={k:v for k,v in c.items() if 'CFL' in k},
        axialHeatChange=regions,negativeHeatAreaFraction=float(sum(area[q[1]<0])/sum(area)),
        withinRunHeatRelativeL1Change=a['withinRunHeatRelativeL1Change'],
        last25HeatPowerRelativeChange=a['last25HeatPowerRelativeChange'],
        interpretation='Residual rebound and reconstruction fallback remain; axial bins are diagnostics, not separation labels. No CFL causal proof or new solver run.',
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),folder/'history.csv',folder/'engine_stdout.log',folder/'case.cfg',folder/'engine_process.json']},
        physicsQualified=False,newQualifiedArea=0,solverRuns=0)
    target=base/'multistep_diagnosis.json';assert not target.exists()
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','history','axialHeatChange']},indent=2))


if __name__=='__main__':main()
