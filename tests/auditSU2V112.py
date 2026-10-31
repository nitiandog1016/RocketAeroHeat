"""独立核对多重网格筛选、原面值及真实超时UI的CSV，不重复启动求解。"""
from pathlib import Path
import csv
import hashlib
import json
import re
import sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from auditSU2V110 import vtk_tokens


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    out=ROOT/'outputs/v112_validation';assert not (out/'independent_audit.json').exists()
    s=loadmat(str(ROOT/'outputs/v108_validation/normal01/normal_volume.mat'),squeeze_me=True,struct_as_record=False)
    v=s['volume'];proj=s['projection'];mapped={};runs={};paths=[Path(__file__).resolve(),ROOT/'tests/auditSU2V110.py',ROOT/'outputs/v108_validation/normal01/normal_volume.mat']
    for name,relative in [('control','outputs/v111_validation/control01'),('multigrid','outputs/v112_validation/mg01')]:
        p=ROOT/relative;receipt=json.loads((p/'engine_process.json').read_text(encoding='utf-8'))
        assert receipt['status']=='terminal' and receipt['exitCode']==0 and all(receipt['inputsStable'].values())
        for file,h in receipt['outputSHA256'].items():assert sha(p/file)==h
        h=np.atleast_1d(np.genfromtxt(str(p/'history.csv'),delimiter=',',names=True));assert np.array_equal(h['Inner_Iter'],np.arange(600))
        rest=np.loadtxt(str(p/'restart_flow.csv'),delimiter=',',skiprows=1)
        assert np.array_equal(rest[:,0],np.arange(len(v.vertices))) and np.all(np.isfinite(rest))
        rho=rest[:,4];e=rest[:,8]/rho-(rest[:,5:8]/rho[:,None]) .__pow__(2).sum(axis=1)/2
        assert min(rho)>0 and min(e)>0
        x,f=vtk_tokens(p/'surface_flow.vtk');distance,ids=cKDTree(v.vertices).query(x)
        assert max(distance)<1e-6 and len(set(ids))==len(ids)
        nodal=np.zeros((len(v.vertices),3));nodal[ids]=np.column_stack([f[k][:,0] for k in ('Pressure','Heat_Flux','Temperature')])
        values=proj.weights@nodal;mapped[name]=values
        native=float(proj.engineNodeArea@nodal[:,1]);assert abs(native/h['HF'][-1]-1)<2e-6
        linsol=h['Linear_Solver_Iterations'];counts={str(int(i)):int(sum(linsol==i)) for i in np.unique(linsol)}
        runs[name]=dict(seconds=receipt['elapsedSeconds'],rhoLog=float(h['rmsRho'][-1]),energyLog=float(h['rmsRhoE'][-1]),
                        nativePower=native,physicalAreaPower=float(proj.targetFaceArea@values[:,1]),
                        linearReturnedIndexCounts=counts,linearResidualLogMaximum=float(h['LinSolRes'].max()),
                        minimumDensity=float(min(rho)),minimumInternalEnergy=float(min(e)),
                        last250HeatPowerRelativeChange=float(h['HF'][-1]/h['HF'][-251]-1),originalFaces=len(values))
        paths += [p/name for name in ('engine_process.json','case.cfg','history.csv','surface_flow.vtk','restart_flow.csv','initial_flow.csv')]
    p=out/'mg01';assert sha(p/'initial_flow.csv')==sha(ROOT/'outputs/v111_validation/control01/initial_flow.csv')
    cfg0=(ROOT/'outputs/v111_validation/control01/case.cfg').read_text();cfg1=(p/'case.cfg').read_text()
    base=cfg1.split('MGCYCLE=')[0].replace('MGLEVEL= 1','MGLEVEL= 0');assert base==cfg0
    log=(p/'engine_stdout.log').read_text();assert 'Initialize Jacobian structure (Navier-Stokes). MG level: 1.' in log
    coarse=re.search(r'\|\s*1\|\s*(\d+)\|\s*1/7.59\|',log);assert coarse
    control=mapped['control'];mg=mapped['multigrid'];area=proj.targetFaceArea
    differences={name:float(area@abs(mg[:,i]-control[:,i])/(area@abs(control[:,i]))) for i,name in enumerate(['pressure','heatFlux','fluidNodeTemperature'])}
    improvement={key:runs['control'][key]-runs['multigrid'][key] for key in ('rhoLog','energyLog')}
    saved=loadmat(str(out/'review01/mg_mapped.mat'),squeeze_me=True,struct_as_record=False)['result']
    for i,key in enumerate(['pressureRaw','heatFluxRaw','fluidNodeTemperatureMeanRaw']):assert np.allclose(mg[:,i],getattr(saved,key),rtol=2e-12,atol=1e-9)
    ui=out/'ui02';r=loadmat(str(ui/'reopened.mat'),squeeze_me=True,struct_as_record=False)['r']
    table=np.genfromtxt(str(ui/'fields.csv'),delimiter=',',names=True)
    assert np.array_equal(table['originalFaceID'],np.arange(1,2449)) and np.all(table['productionQualified']==0)
    assert np.all(np.isnan(table['solvedWallTemperature_K'])) and np.all(table['givenWallTemperature_K']==300)
    for col,key in [('pressureRaw_Pa','pressureRaw'),('heatFluxRaw_W_m2','heatFluxRaw'),('fluidNodeTemperatureMean_K','fluidNodeTemperatureMeanRaw')]:
        assert np.allclose(table[col],getattr(r,key),rtol=1e-12,atol=1e-8)
    with (ui/'section.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.reader(f))[1:]
    for row in rows:assert abs(float(row[10])-r.heatFluxRaw[int(row[2])-1])<1e-8
    assert len(rows)==48
    meta=json.loads((ui/'fields.csv.json').read_text(encoding='utf-8'));a=meta['convergenceAudit']
    assert a['checkpointIteration']==100 and a['terminalIteration']==103 and not a['productionQualified']
    paths += [out/'review01/mg_mapped.mat',ui/'reopened.mat',ui/'fields.csv',ui/'section.csv',ui/'fields.csv.json',ui/'section.csv.json']
    report=dict(runs=runs,coarseControlVolumes=int(coarse.group(1)),relativeRuntime=runs['multigrid']['seconds']/runs['control']['seconds'],
                logImprovement=improvement,multigridScreenPassed=all(x>=.3 for x in improvement.values()),
                originalFaceAreaWeightedRelativeL1=differences,defaultMultigridStillZero=True,
                matlabMappingIndependentlyMatched=True,uiFullFaceRows=2448,uiSectionRows=48,
                uiCheckpointIteration=100,uiTerminalIteration=103,physicsQualified=False,newQualifiedArea=0,
                inputSHA256={p.relative_to(ROOT).as_posix():sha(p) for p in paths})
    (out/'independent_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('logImprovement','relativeRuntime','originalFaceAreaWeightedRelativeL1','multigridScreenPassed','newQualifiedArea')}))


if __name__=='__main__':main()
