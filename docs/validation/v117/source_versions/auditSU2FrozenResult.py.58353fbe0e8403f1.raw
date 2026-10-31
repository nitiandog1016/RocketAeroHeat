"""独立重算冻结节点恢复、EOS、原面/CSV和局部具名迭代；恢复不是精度认证。"""
from pathlib import Path
import json,hashlib,csv,re,sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v117_validation'
sys.path.insert(0,str(ROOT/'tests'));from auditSU2V110 import vtk_tokens
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def main():
    target=OUT/'independent_audit.json';assert not target.exists()
    gate=read(OUT/'diagnostic_gate.json')
    cases={};paths=[]
    names=['diagnostic','apply'] if gate['passed'] else ['diagnostic']
    if not gate['passed']:assert not (OUT/'apply').exists(),'启动门槛失败后不得运行候选'
    for name in names:
        folder=OUT/name;rec=read(folder/'engine_process.json')
        assert rec['status']=='terminal' and rec['exitCode']==0 and not rec['timedOut']
        assert all(rec['inputsStable'].values()) and rec['frozenEnergyMode']==name
        for p,h in rec['inputSHA256'].items():assert sha(Path(p))==h,p
        for p,h in rec['outputSHA256'].items():assert sha(folder/p)==h,p
        steps=2 if name=='diagnostic' else 600
        for base,ext in [('restart_flow','csv'),('surface_flow','vtk')]:
            assert sha(folder/(base+'.'+ext))==sha(folder/(base+'_'+str(steps-1).zfill(6)+'.'+ext))
        d=np.loadtxt(str(folder/'restart_flow.csv'),delimiter=',',skiprows=1)
        rho=d[:,4];e=d[:,8]/rho-.5*np.sum((d[:,5:8]/rho[:,None])**2,axis=1)
        T=e*.4/287.05287;p=e*.4*rho
        h=np.atleast_1d(np.genfromtxt(str(folder/'history.csv'),delimiter=',',names=True))
        assert np.array_equal(h['Inner_Iter'],np.arange(steps)) and (h['HF']==0).all()
        counts=np.atleast_1d(np.genfromtxt(str(folder/'rah_frozen_counts.csv'),delimiter=',',names=True))
        assert set(counts['iteration'])==set(range(steps))
        applied=int(sum(counts['applied']));assert (applied==0)==(name=='diagnostic')
        init=np.loadtxt(str(folder/'initial.csv'),delimiter=',',skiprows=1)
        nodes=[]
        for j in [29214,31040,32866,34692]:
            nodes.append(dict(originalNodeID=j,temperatureK=float(T[j]),conservativeChanged=not np.array_equal(init[j,4:9],d[j,4:9])))
        log=(folder/'engine_stdout.log').read_text(encoding='utf-8',errors='replace')
        warnings=re.findall(r'Warning: (\d+) reconstructed states',log)
        cases[name]=dict(seconds=rec['elapsedSeconds'],iterations=steps,checkpointIteration=steps-1,
            finite=bool(np.isfinite(d).all()),minimumRho=float(min(rho)),temperatureRangeK=[float(min(T)),float(max(T))],
            residualRho=float(h['rmsRho'][-1]),residualEnergy=float(h['rmsRhoE'][-1]),
            residualGatePassed=bool(h['rmsRho'][-1]<=-8 and h['rmsRhoE'][-1]<=-8),
            lastReportedReconstructedStates=int(warnings[-1]) if warnings else 0,appliedUpdates=applied,
            targetNodes=nodes,physicsQualified=False,newQualifiedArea=0)
        paths += [folder/n for n in ['engine_process.json','history.csv','rah_frozen_counts.csv','rah_frozen_directions.csv','restart_flow.csv','surface_flow.vtk']]
        if name!=names[-1]:continue
        ui=OUT/'review01';r=loadmat(str(ui/'result.mat'),squeeze_me=True,struct_as_record=False)['r'];V=r.V;F=r.F.astype(int)-1
        x,fields=vtk_tokens(folder/'surface_flow.vtk');distance,ids=cKDTree(V).query(x)
        assert max(distance)<1e-6 and len(np.unique(ids))==len(V)
        v=np.zeros((len(V),3));v[ids]=np.column_stack([fields[k][:,0] for k in ['Pressure','Heat_Flux','Temperature']]);expected=v[F].mean(axis=1)
        table=np.genfromtxt(str(ui/'fields.csv'),delimiter=',',names=True)
        for j,k in enumerate(['pressureRaw_Pa','heatFluxRaw_W_m2','fluidNodeTemperatureMean_K']):assert np.allclose(table[k],expected[:,j],rtol=1e-12,atol=1e-8)
        assert len(table)==3648 and np.array_equal(table['originalFaceID'],np.arange(1,3649))
        assert table['heatFluxIsBoundaryInput'].all() and not table['productionQualified'].any()
        assert np.isnan(table['givenWallTemperature_K']).all() and np.isnan(table['solvedWallTemperature_K']).all()
        distance,node=cKDTree(d[:,1:4]).query(V);assert max(distance)<1e-14
        eos=max(float(max(abs(p[node]/v[:,0]-1))),float(max(abs(T[node]/v[:,2]-1))));assert eos<2e-6
        with (ui/'section.csv').open(encoding='utf-8-sig') as f:rows=list(csv.reader(f))[1:]
        assert len(rows)==32
        for row in rows:
            assert abs(float(row[0])-.158496)<1e-14 and abs(float(row[10])-expected[int(row[2])-1,0])<1e-8
        for file in ['fields.csv.json','section.csv.json']:
            meta=read(ui/file);assert meta['algorithmRuntime']['mode']==name and not meta['algorithmRuntime']['physicsQualified']
        original=loadmat(str(ROOT/'outputs/v116_validation/review02/control_result.mat'),squeeze_me=True,struct_as_record=False)['r']
        area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
        cases[name].update(originalFaces=len(F),sectionRows=len(rows),EOSrelativeError=eos,
            pressureDifferenceToControl=float(area@abs(expected[:,0]-original.pressureRaw)/(area@abs(original.pressureRaw))),
            temperatureDifferenceToControl=float(area@abs(expected[:,2]-original.fluidNodeTemperatureMeanRaw)/(area@abs(original.fluidNodeTemperatureMeanRaw))))
        paths += [ui/n for n in ['result.mat','report.json','fields.csv','section.csv','fields.csv.json','section.csv.json']]
    control=read(ROOT/'outputs/v116_validation/initial_comparison.json')['control']
    result=dict(passed=True,candidateAuthorized=gate['passed'],candidateRun=('apply' in cases),cases=cases,priorControl=control,sourceFiles={str(p):sha(p) for p in paths+[Path(__file__),ROOT/'tests/auditSU2V110.py']},
        physicsQualified=False,newQualifiedArea=0,meaning='mechanism and field identity audit; not validation or discretization error')
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(cases,indent=2))
if __name__=='__main__':main()
