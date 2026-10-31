"""独立审核真实19752原面、保存身份、EOS与三个任意截面；不将全值当合格覆盖。"""
from pathlib import Path
import json,hashlib,csv,sys,re
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v118_validation';ENGINE=OUT/'actual01/budget_run';UI=OUT/'review01'
sys.path.insert(0,str(ROOT/'tests'));from auditSU2V110 import vtk_tokens
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def main():
    target=OUT/'independent_audit.json';assert not target.exists()
    receipt=read(ENGINE/'engine_process.json');assert receipt['status']=='terminal'
    assert (receipt['exitCode']==0 or receipt['timedOut']) and all(receipt['inputsStable'].values())
    for p,h in receipt['inputSHA256'].items():assert sha(Path(p))==h,p
    for p,h in receipt['outputSHA256'].items():assert sha(ENGINE/p)==h,p
    r=loadmat(str(UI/'result.mat'),squeeze_me=True,struct_as_record=False)['r'];V=r.V;F=r.F.astype(int)-1
    x,fields=vtk_tokens(ENGINE/'surface_flow.vtk');dist,ids=cKDTree(V).query(x)
    assert max(dist)<1e-6 and len(np.unique(ids))==len(V)
    values=np.zeros((len(V),3));values[ids]=np.column_stack([fields[k][:,0] for k in ['Pressure','Heat_Flux','Temperature']])
    expected=values[F].mean(axis=1)
    table=np.genfromtxt(str(UI/'fields.csv'),delimiter=',',names=True);assert len(table)==19752
    assert np.array_equal(table['originalFaceID'],np.arange(1,19753))
    for i,col in enumerate(['pressureRaw_Pa','heatFluxRaw_W_m2','fluidNodeTemperatureMean_K']):
        assert np.allclose(table[col],expected[:,i],rtol=1e-12,atol=1e-8)
    assert not table['productionQualified'].any() and not table['heatFluxIsBoundaryInput'].any()
    assert (table['givenWallTemperature_K']==300).all() and np.isnan(table['solvedWallTemperature_K']).all()
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    assert abs(sum(area)-1.9646292565715553)<1e-12
    assert abs(area@expected[:,1]-r.targetHeatPowerRaw)<1e-6*max(1,abs(r.targetHeatPowerRaw))
    av=np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]])/2;nav=np.zeros_like(V)
    for j in range(3):np.add.at(nav,F[:,j],av/3)
    native=np.linalg.norm(nav,axis=1)@values[:,1]
    h=np.atleast_1d(np.genfromtxt(str(ENGINE/'history.csv'),names=True,delimiter=','));last=int(h['Inner_Iter'][-1])
    step=int(r.convergenceAudit.checkpointIteration)
    assert np.array_equal(h['Inner_Iter'],np.arange(last+1)) and step<=last
    assert sha(ENGINE/'restart_flow.csv')==sha(ENGINE/('restart_flow_'+str(step).zfill(6)+'.csv'))
    assert sha(ENGINE/'surface_flow.vtk')==sha(ENGINE/('surface_flow_'+str(step).zfill(6)+'.vtk'))
    assert abs(native-h['HF'][step])/max(1,abs(h['HF'][step]))<2e-6
    state=np.loadtxt(str(ENGINE/'restart_flow.csv'),skiprows=1,delimiter=',');rho=state[:,4]
    e=state[:,8]/rho-.5*np.sum((state[:,5:8]/rho[:,None])**2,axis=1);T=e*.4/287.05287;p=e*rho*.4
    assert np.isfinite(state).all() and min(rho)>0 and min(e)>0
    dist,nodes=cKDTree(state[:,1:4]).query(V);assert max(dist)<1e-14
    eos=max(float(max(abs(p[nodes]/values[:,0]-1))),float(max(abs(T[nodes]/values[:,2]-1))));assert eos<2e-6
    sections=[]
    for k,(xx,j) in enumerate([(-.2,0),(-1.5,1),(-2.65,2)],1):
        file=UI/('section'+str(k)+'.csv')
        with file.open(encoding='utf-8-sig') as f:rows=list(csv.reader(f))[1:]
        assert len(rows)>2
        for row in rows:
            assert abs(float(row[0])-xx)<1e-14
            assert abs(float(row[10])-expected[int(row[2])-1,j])<1e-8*max(1,abs(expected[int(row[2])-1,j]))
        meta=read(Path(str(file)+'.json'));assert meta['algorithmRuntime']['flavor']=='official'
        assert meta['convergenceAudit']['checkpointIteration']==step
        sections.append(dict(X=xx,rows=len(rows),variableIndex=j))
    text=(ENGINE/'engine_stdout.log').read_text(encoding='utf-8',errors='replace')
    warnings=re.findall(r'Warning: (\d+) reconstructed states',text)
    report=dict(passed=True,physicsQualified=False,newQualifiedArea=0,originalFaces=len(F),areaM2=float(sum(area)),
        engineSeconds=receipt['elapsedSeconds'],timedOut=receipt['timedOut'],exitCode=receipt['exitCode'],
        terminalIteration=last,checkpointIteration=step,terminalResiduals=[float(h['rmsRho'][-1]),float(h['rmsRhoE'][-1])],
        checkpointResiduals=[float(h['rmsRho'][step]),float(h['rmsRhoE'][step])],
        minimumDensity=float(min(rho)),temperatureRangeK=[float(min(T)),float(max(T))],maximumWallEOSRelativeError=eos,
        fullAreaHeatPowerRawW=float(area@expected[:,1]),nativeHeatPowerRawW=float(native),sections=sections,
        pressureRangePa=[float(min(expected[:,0])),float(max(expected[:,0]))],heatFluxRangeWm2=[float(min(expected[:,1])),float(max(expected[:,1]))],
        lastReportedNonphysicalReconstructions=int(warnings[-1]) if warnings else 0,
        savedStateNotNecessarilyTerminal=(step!=last),
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'tests/auditSU2V110.py',UI/'result.mat',ENGINE/'engine_process.json',ENGINE/'history.csv']})
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
