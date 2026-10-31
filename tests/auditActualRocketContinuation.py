"""独立核对实际续算原面/累计身份、正状态及运行前已固定的继续门槛。"""
from pathlib import Path
import json,hashlib,csv,sys,re
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v119_validation';ENGINE=OUT/'continuation01'
sys.path.insert(0,str(ROOT/'tests'));from auditSU2V110 import vtk_tokens

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))

def main():
    target=OUT/'independent_audit.json';assert not target.exists()
    receipt=read(ENGINE/'engine_process.json');assert receipt['status']=='terminal' and all(receipt['inputsStable'].values())
    for p,h in receipt['inputSHA256'].items():assert sha(Path(p))==h,p
    for p,h in receipt['outputSHA256'].items():assert sha(ENGINE/p)==h,p
    source=ROOT/'outputs/v118_validation/actual01/budget_run';lineage=read(ENGINE/'continuation.json')
    for p,h in lineage['inputSourceHashes'].items():assert sha(Path(p))==h,p
    assert sha(ENGINE/'initial_flow.csv')==sha(source/'restart_flow_000299.csv')
    assert lineage['completedUpdatesBeforeRun']==300
    r=loadmat(str(OUT/'review01/result.mat'),squeeze_me=True,struct_as_record=False)['r'];V=r.V;F=r.F.astype(int)-1
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    assert len(F)==19752 and abs(sum(area)-1.9646292565715553)<1e-12
    def wall(path):
        xyz,f=vtk_tokens(path);d,ids=cKDTree(V).query(xyz);assert max(d)<1e-6 and len(set(ids))==len(V)
        v=np.zeros((len(V),3));v[ids]=np.column_stack([f[k][:,0] for k in ['Pressure','Heat_Flux','Temperature']])
        return v,v[F].mean(axis=1),f
    values,faces,fields=wall(ENGINE/'surface_flow.vtk')
    table=np.genfromtxt(str(OUT/'review01/fields.csv'),names=True,delimiter=',')
    assert len(table)==19752 and np.array_equal(table['originalFaceID'],np.arange(1,19753))
    for j,col in enumerate(['pressureRaw_Pa','heatFluxRaw_W_m2','fluidNodeTemperatureMean_K']):
        assert np.allclose(table[col],faces[:,j],rtol=1e-12,atol=1e-8)
    assert (table['givenWallTemperature_K']==300).all() and np.isnan(table['solvedWallTemperature_K']).all()
    assert not table['productionQualified'].any()
    h=np.atleast_1d(np.genfromtxt(str(ENGINE/'history.csv'),delimiter=',',names=True));step=int(r.convergenceAudit.checkpointIteration)
    assert np.array_equal(h['Inner_Iter'],np.arange(len(h)))
    assert sha(ENGINE/'restart_flow.csv')==sha(ENGINE/('restart_flow_%06d.csv'%step))
    assert sha(ENGINE/'surface_flow.vtk')==sha(ENGINE/('surface_flow_%06d.vtk'%step))
    assert r.convergenceAudit.cumulativeCheckpointIteration==step+300
    av=np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]])/2;nav=np.zeros_like(V)
    for j in range(3):np.add.at(nav,F[:,j],av/3)
    native=np.linalg.norm(nav,axis=1)@values[:,1]
    assert abs(native-h['HF'][step])/max(abs(h['HF'][step]),1)<2e-6
    state=np.loadtxt(str(ENGINE/'restart_flow.csv'),delimiter=',',skiprows=1);rho=state[:,4]
    velocity=state[:,5:8]/rho[:,None];e=state[:,8]/rho-.5*np.sum(velocity**2,axis=1)
    positive=bool(np.isfinite(state).all() and min(rho)>0 and min(e)>0);assert positive
    T=e*.4/287.05287;p=e*rho*.4;d,ids=cKDTree(state[:,1:4]).query(V);assert max(d)<1e-14
    eos=max(float(max(abs(p[ids]/values[:,0]-1))),float(max(abs(T[ids]/values[:,2]-1))));assert eos<2e-6
    for name in ['fields.csv.json','section.csv.json']:
        m=read(OUT/'review01'/name);assert m['convergenceAudit']['cumulativeCheckpointIteration']==step+300
        assert m['convergenceAudit']['completedUpdatesBeforeRun']==300
    with (OUT/'review01/section.csv').open(encoding='utf-8-sig') as f:rows=list(csv.reader(f))[1:]
    for row in rows:assert abs(float(row[0])+1.5)<1e-14 and abs(float(row[10])-faces[int(row[2])-1,1])<1e-8*max(1,abs(float(row[10])))
    baseline=np.genfromtxt(str(source/'history.csv'),names=True,delimiter=',')
    improvement=[float(baseline[k][-1]-h[k][-1]) for k in ['rmsRho','rmsRhoE']]
    _,mid,_=wall(ENGINE/'surface_flow_000150.vtk')
    qchange=float(area@abs(faces[:,1]-mid[:,1])/(area@abs(mid[:,1])))
    pchange=float(area@abs(faces[:,0]-mid[:,0])/(area@abs(mid[:,0])))
    tail=float((h['HF'][-1]-h['HF'][-50])/max(abs(h['HF'][-50]),1))
    gates=dict(normalTerminal=receipt['exitCode']==0 and not receipt['timedOut'] and not receipt['cancelled'],
               densityImprovesAtLeast03=improvement[0]>=.3,energyImprovesAtLeast03=improvement[1]>=.3,
               heatDistributionChangeBelow1643875Percent=qchange<.1643875,
               last50HeatPowerChangeBelow10Percent=abs(tail)<.10,positiveState=positive,originalFullArea=True)
    warning=re.findall(r'Warning: (\d+) reconstructed states',(ENGINE/'engine_stdout.log').read_text())
    yp=fields['Y_Plus'][:,0]
    report=dict(passed=True,physicsQualified=False,newQualifiedArea=0,originalFaces=len(F),areaM2=float(sum(area)),
                localCheckpointIteration=step,cumulativeCheckpointIteration=step+300,
                engineSeconds=receipt['elapsedSeconds'],successfulChainSeconds=lineage['priorSuccessfulEngineSeconds']+receipt['elapsedSeconds'],
                discardedV118AttemptSeconds=434.156,finalResiduals=[float(h[k][-1]) for k in ['rmsRho','rmsRhoE']],
                residualLogImprovement=improvement,withinRunPressureRelativeL1Change=pchange,withinRunHeatRelativeL1Change=qchange,
                last50NativePowerRelativeChange=tail,predeclaredNumericalGates=gates,allNumericalGatesPassed=all(gates.values()),
                originalResidualThresholdPassed=bool(h['rmsRho'][-1]<-8 and h['rmsRhoE'][-1]<-8),
                minimumDensity=float(min(rho)),temperatureRangeK=[float(min(T)),float(max(T))],wallEOSRelativeError=eos,
                reportedYPlusRange=[float(min(yp)),float(max(yp))],reportedYPlusAbove1Nodes=int(sum(yp>1)),
                yPlusSemantics='native laminar starting field; not turbulent-mesh qualification',
                finalReportedReconstructionWarnings=int(warning[-1]) if warning else 0,
                nativeHeatPowerW=float(native),originalFaceHeatPowerW=float(area@faces[:,1]),sectionRows=len(rows),
                sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'tests/auditSU2V110.py',ENGINE/'engine_process.json',OUT/'review01/result.mat']})
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
