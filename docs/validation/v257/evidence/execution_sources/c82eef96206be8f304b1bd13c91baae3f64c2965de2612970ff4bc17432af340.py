"""模块职责：独立核对KINSOL实际返回态、原方程残差和热壁量；不再次运行求解器。"""
import csv,json,time
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from auditFullResidualV256 import matrix
from phaseAwareViscousAir import PhaseAwareViscousAir
BASE=ROOT/'outputs/v257_validation';CASE=BASE/'rocket01';CLOSED=False


def vector(path):
    with path.open('rb') as f:
        assert f.read(8)==b'RAH256V1';n=int(np.frombuffer(f.read(8),dtype='<u8')[0])
    assert path.stat().st_size==16+n*8
    return np.memmap(path,mode='r',offset=16,dtype='<f8',shape=(n,))


def audit():
    start=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    receipt=read(BASE/'process_control/report.json');native=read(CASE/'native_completion.json')
    assert receipt['passed'] and receipt['inputsStable'] and receipt['status']=='terminal' and not receipt['liveGroupAfterCleanup']
    assert native['normalReturn'] and native['KINSOLExecuted'] and native['nativeNonlinearAttempts']<=3
    assert native['initialStateAndResidualBitwiseV256'] and not native['physicalQualification']
    assert sha(CASE/'initial_state.bin')==sha(ROOT/'outputs/v256_validation/rocket02/initial_state.bin')
    assert sha(CASE/'initial_full_residual.bin')==sha(ROOT/'outputs/v256_validation/rocket02/base_full_residual.bin')
    rows=list(csv.DictReader((CASE/'frames.csv').open(encoding='utf8')))
    faces=np.loadtxt(ROOT/'outputs/v250_validation/rocket01/original_faces.txt');air=PhaseAwareViscousAir()
    assert faces.shape==(19752,5) and np.array_equal(faces[:,0],np.arange(1,19753))
    reports=[];previous=None;changes=0;surface={}
    for row in rows:
        label=row['frame'];s=matrix(CASE/(label+'_state.bin'),b'RAH256S1',19);r=matrix(CASE/(label+'_full_residual.bin'),b'RAH256R1',7)
        x=vector(CASE/(label+'_unknowns.bin'));assert len(s)==335726 and len(x)==1974844
        assert np.isfinite(s).all() and np.isfinite(r).all() and np.isfinite(x).all()
        assert np.array_equal(np.sort(s[:,0]),np.arange(len(s))) and np.array_equal(s[:,0],r[:,0])
        wall=s[:,17]==1;assert wall.sum()==9878 and np.all(s[wall,6:9]==0) and np.all(s[wall,10]==0)
        keep=np.ones((len(s),6),dtype=bool);keep[np.ix_(wall,[1,2,3,5])]=False
        assert np.array_equal(np.asarray(s[:,5:11])[keep],x)
        if previous is not None:changes+=not np.array_equal(previous,x)
        previous=x
        gas=air.evaluate(s[:,12],s[:,11]);transport=air.interpolated_transport(s[:,12],s[:,11]);rho=s[:,5]
        kinetic=np.sum(s[:,6:9]**2,axis=1)/(2*rho*rho);energy=rho*(gas['internalEnergyJkg']+kinetic)
        errors=dict(density=float(np.max(abs(gas['densityKgm3']/rho-1))),physicalEnergy=float(np.max(abs(energy-s[:,9])/np.maximum(1,abs(s[:,9])))),
            mu=float(np.max(abs(transport['dynamicViscosityPas']/s[:,14]-1))),conductivity=float(np.max(abs(transport['conductivityWmK']/s[:,13]-1))))
        nu=s[:,10];chi=nu/(s[:,14]/rho);mut=np.zeros(len(s));positive=nu>=0
        mut[positive]=rho[positive]*nu[positive]*chi[positive]**3/(chi[positive]**3+7.1**3)
        errors['SAEddyViscosity']=float(np.max(abs(mut-s[:,15])/np.maximum(1e-30,np.maximum(abs(mut),abs(s[:,15])))))
        errors['sameStateFlowSAEddyViscosity']=float(np.max(abs(s[:,16]-s[:,15])/np.maximum(1e-30,np.maximum(abs(s[:,16]),abs(s[:,15])))))
        assert max(errors.values())<=2e-10,errors
        inverse=np.empty(len(s),dtype=int);inverse[s[:,0].astype(int)]=np.arange(len(s));node=inverse[faces[:,1:4].astype(int)]
        assert np.all(wall[node]);xyz=s[node,1:4];area=np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)/2
        areaError=float(np.max(abs(area-faces[:,4])));assert areaError<1e-14 and abs(area.sum()-1.9646292565715553)<2e-13
        w=np.genfromtxt(CASE/(label+'_wall.csv'),delimiter=',',names=True);wi=inverse[w['globalID'].astype(int)];wj=inverse[w['neighborID'].astype(int)]
        assert len(w)==9878 and np.array_equal(np.sort(wi),np.where(wall)[0]) and np.all(w['distanceM']>0)
        for key,col in [('pressurePa',12),('wallFluidK',11),('k',13),('mut',16)]:assert np.array_equal(w[key],s[wi,col]),key
        assert np.array_equal(w['neighborK'],s[wj,11]);assert np.max(abs(w['Cp']/gas['cpJkgK'][wi]-1))<1e-12
        q=(s[wi,13]+gas['cpJkgK'][wi]*s[wi,16]/w['Prt'])*(s[wj,11]-300)/w['distanceM']
        qError=float(np.max(abs(q-w['qIntoWallWm2'])/np.maximum(1,abs(q))));assert qError<2e-12
        qnode=np.full(len(s),np.nan);qnode[wi]=q;qs=np.mean(qnode[node],axis=1);ps=np.mean(s[node,12],axis=1);tf=np.mean(s[node,11],axis=1)
        rms=np.sqrt(np.mean(r[:,1:]**2,axis=0));claimed=np.array([float(row[k]) for k in ['mass','momX','momY','momZ','physicalEnergy','SA']])
        assert np.max(abs(rms/claimed-1))<1e-11
        reports.append(dict(frame=label,nativeAttempt=int(row['nativeAttempt']),observedCurrentChanges=int(row['observedCurrentChanges']),
            physicalChecks=errors,physicalSixEquationRMS=rms.tolist(),physicalSixEquationSum=np.sum(r[:,1:],axis=0).tolist(),
            allSixResidualsBelow1eMinus8=bool(np.all(rms<=1e-8)),nativeReconstructionWarnings=int(row['nonphysicalRecon']),
            originalFaceAreaAbsoluteError=areaError,wallQCheckRelativeError=qError,originalFaceSignedPowerW=float(qs@faces[:,4]),
            researchFaceQRangeWm2=[float(qs.min()),float(qs.max())],negativeHeatFaces=int(np.count_nonzero(qs<0)),
            wallFluidTemperatureRangeK=[float(s[wall,11].min()),float(s[wall,11].max())]))
        if label in ['initial','terminal']:
            for key,value in [('pressurePa',ps),('qConvGivenTwWm2',qs),('wallFluidTemperatureK',tf),('physicalEnergyResidual',np.mean(r[node,5],axis=1))]:surface[label+'_'+key]=value
    assert changes==native['observedCurrentStateChanges']
    assert np.array_equal(vector(CASE/'returned_unknowns.bin'),vector(CASE/'terminal_unknowns.bin'))
    gates=list(csv.DictReader((CASE/'linear_gate.csv').open(encoding='utf8')))
    for row in gates:
        expected=int(row['flag'])==0 and 0<=float(row['reportedResidual'])<=float(row['target'])
        assert bool(int(row['passed']))==expected
    refused=[row for row in gates if row['passed']=='0'];assert bool(refused)==native['linearGateRefused']
    scales=np.genfromtxt(CASE/'scaling.csv',delimiter=',',names=True)
    terminal=matrix(CASE/'terminal_full_residual.bin',b'RAH256R1',7)
    scaledNorm=float(np.linalg.norm((np.asarray(terminal[:,1:])/scales['residualCharacteristic'])[keep]))
    # 返回状态和KINSOL内部范数各自保存，不用其退出码替代真正更新或资格。
    surface.update(faceIDs=faces[:,0].astype(int),nodeIDs=faces[:,1:4].astype(int),areaM2=faces[:,4]);np.savez_compressed(out/'surface_values.npz',**surface)
    report=dict(passed=True,frames=reports,actualObservedStateChanges=changes,nativeSummary=native,
        initialBitwiseV256=True,returnedStateMatchesSavedTerminal=True,terminalScaledL2Norm=scaledNorm,
        rawResidualRatios=(np.array(reports[-1]['physicalSixEquationRMS'])/reports[0]['physicalSixEquationRMS']).tolist(),
        linearGates=gates,originalFaces=19752,originalAreaM2=1.9646292565715553,
        physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,
        givenWallTemperatureK=300,solvedSolidWallTemperatureK=None,actualProfilesFilled=0,
        sameStageWallThermalFormulaVerified=True,completeBoundaryConservationVerified=False,independentLinearDirectionAccuracyProved=False,
        heatStabilityVerified=False,meshAccuracyVerified=False,experimentalValidationPassed=False,
        inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in CASE.glob('*') if p.is_file()},elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ['passed','actualObservedStateChanges','rawResidualRatios','elapsedSeconds']}))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v257 fixed coupled audit closed; follow NEXT_DECISION')
    audit()
