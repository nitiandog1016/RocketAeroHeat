"""模块职责：独立读取完整状态/残差、同阶段物性与SA刷新、原壁及真实拒绝回执。"""
import json,re,time
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseAwareViscousAir import PhaseAwareViscousAir
BASE=ROOT/'outputs/v256_validation';CASE=BASE/'rocket02';CLOSED=False


def matrix(path,magic,cols):
    with open(path,'rb') as f:
        assert f.read(8)==magic;n,k=np.frombuffer(f.read(16),dtype='<u8')
    assert int(k)==cols and path.stat().st_size==24+int(n)*cols*8
    return np.memmap(path,mode='r',offset=24,dtype='<f8',shape=(int(n),cols))


def audit():
    started=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    parent=read(BASE/'process_control2/report.json');n=read(CASE/'native_completion.json')
    assert parent['passed'] and parent['inputsStable'] and parent['status']=='terminal' and not parent['liveGroupAfterCleanup']
    assert n['normalReturn'] and n['actualTimeUpdates']==n['nonlinearIterations']==0
    state=matrix(CASE/'initial_state.bin',b'RAH256S1',19);end=matrix(CASE/'terminal_state.bin',b'RAH256S1',19)
    r=matrix(CASE/'base_full_residual.bin',b'RAH256R1',7);p=matrix(CASE/'probe_full_residual.bin',b'RAH256R1',7);rr=matrix(CASE/'restored_full_residual.bin',b'RAH256R1',7)
    assert len(state)==335726 and np.array_equal(np.sort(state[:,0]),np.arange(len(state)))
    assert sha(CASE/'initial_state.bin')==sha(CASE/'terminal_state.bin')
    assert sha(CASE/'base_full_residual.bin')==sha(CASE/'restored_full_residual.bin')
    assert np.array_equal(state[:,0],r[:,0]) and np.array_equal(r[:,0],p[:,0])
    assert np.isfinite(state).all() and np.isfinite(r).all() and np.isfinite(p).all()
    wall=state[:,17]==1;assert wall.sum()==9878 and np.all(state[wall,6:9]==0) and np.all(state[wall,10]==0)
    assert 6*len(state)-4*wall.sum()==n['unknowns']==1974844 and not n['wallEnergyDOFEliminated']
    changed=int(np.count_nonzero(r[:,1:]!=p[:,1:]));assert changed==n['probeChangedResidualEntries']==157
    relative=float(np.max(abs(p[:,1:]-r[:,1:])/np.maximum(1,abs(r[:,1:]))));assert abs(relative-n['probeMaximumScaledDifference'])<1e-15
    # 所有节点的同阶段量热/输运核对；不是只读返回温度就认定状态合法。
    air=PhaseAwareViscousAir();gas=air.evaluate(state[:,12],state[:,11]);transport=air.interpolated_transport(state[:,12],state[:,11])
    rho=state[:,5];kinetic=np.sum(state[:,6:9]**2,axis=1)/(2*rho*rho)
    expectedEnergy=rho*(gas['internalEnergyJkg']+kinetic)
    checks=dict(density=float(np.max(abs(gas['densityKgm3']/rho-1))),
        physicalEnergy=float(np.max(abs(expectedEnergy-state[:,9])/np.maximum(1,abs(state[:,9])))),
        mu=float(np.max(abs(transport['dynamicViscosityPas']/state[:,14]-1))),
        conductivity=float(np.max(abs(transport['conductivityWmK']/state[:,13]-1))))
    nu=state[:,10];chi=nu/(state[:,14]/rho);mut=np.zeros(len(state));positive=nu>=0
    mut[positive]=rho[positive]*nu[positive]*chi[positive]**3/(chi[positive]**3+7.1**3)
    checks['SAEddyViscosity']=float(np.max(abs(mut-state[:,15])/np.maximum(1e-30,np.maximum(abs(mut),abs(state[:,15])))))
    checks['sameStateFlowSAEddyViscosity']=float(np.max(abs(state[:,16]-state[:,15])/np.maximum(1e-30,np.maximum(abs(state[:,16]),abs(state[:,15])))))
    assert max(checks.values())<=2e-10,checks
    # 完整原面ID和真实面积，按体网格全局ID恢复坐标与墙节点；不缩小分母。
    faces=np.loadtxt(ROOT/'outputs/v250_validation/rocket01/original_faces.txt');assert faces.shape==(19752,5) and np.array_equal(faces[:,0],np.arange(1,19753))
    inverse=np.empty(len(state),dtype=int);inverse[state[:,0].astype(int)]=np.arange(len(state));node=inverse[faces[:,1:4].astype(int)]
    assert np.all(wall[node]);xyz=state[node,1:4];a=np.linalg.norm(xyz[:,1]-xyz[:,0],axis=1);b=np.linalg.norm(xyz[:,2]-xyz[:,1],axis=1);c=np.linalg.norm(xyz[:,0]-xyz[:,2],axis=1);semi=(a+b+c)/2
    areas=np.sqrt(semi*(semi-a)*(semi-b)*(semi-c));areaError=float(np.max(abs(areas-faces[:,4])));assert areaError<1e-14 and abs(areas.sum()-1.9646292565715553)<2e-13
    calls=(CASE/'calls.csv').read_text(encoding='utf8').splitlines();assert len(calls)==6 and calls[4]=='4,invalid_density,1,NONPOSITIVE_DENSITY,0'
    warning=re.findall(r'Warning: (\d+) reconstructed states for upwinding are non-physical', (BASE/'process_control2/stdout.txt').read_text(encoding='utf8',errors='replace'))
    rms=np.sqrt(np.mean(r[:,1:]**2,axis=0));summation=np.sum(r[:,1:],axis=0)
    np.savez_compressed(out/'surface_values.npz',faceIDs=faces[:,0].astype(int),nodeIDs=faces[:,1:4].astype(int),areaM2=faces[:,4],
        pressurePa=np.mean(state[node,12],axis=1),wallFluidTemperatureK=np.mean(state[node,11],axis=1),
        massResidual=np.mean(r[node,1],axis=1),physicalEnergyResidual=np.mean(r[node,5],axis=1),SAResidual=np.mean(r[node,6],axis=1))
    report=dict(passed=True,nodes=335726,originalFaces=19752,originalAreaM2=1.9646292565715553,wallNodes=9878,unknowns=n['unknowns'],
        completeStateBitwiseRestored=True,completeResidualBitwiseRestored=True,physicalStateChecks=checks,originalFaceAreaAbsoluteError=areaError,
        probeChangedResidualEntries=changed,probeMaximumScaledDifference=relative,recoverableInvalidTrials=1,invalidTrialChangedState=False,
        physicalSixEquationRMS=rms.tolist(),physicalSixEquationSum=summation.tolist(),allSixResidualsBelow1eMinus8=bool(np.all(rms<=1e-8)),
        wallFluidTemperatureRangeK=[float(state[wall,11].min()),float(state[wall,11].max())],givenWallTemperatureK=300,solvedSolidWallTemperatureK=None,
        nativeReconstructionWarningCounts=[int(x) for x in warning],reconstructionWarningIsNotCellStateRollback=True,
        originalFluxFormulaIndependentlyVerified=False,completeBoundaryConservationVerified=False,nonlinearIterations=0,timeUpdates=0,
        physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,KINSOLExecuted=False,
        inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in [CASE/'initial_state.bin',CASE/'terminal_state.bin',CASE/'base_full_residual.bin',CASE/'probe_full_residual.bin',CASE/'restored_full_residual.bin',CASE/'native_completion.json']},
        elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ['passed','physicalStateChecks','physicalSixEquationRMS','nativeReconstructionWarningCounts','elapsedSeconds']}))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v256 fixed residual audit closed; follow NEXT_DECISION')
    audit()
