"""模块职责：独立分析固定方向的真实差分、近似矩阵与ILU效果；不运行或替换求解内核。"""
import json,time
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from auditFullResidualV256 import matrix
from auditKinsolCoupledV257 import vector
from phaseAwareViscousAir import PhaseAwareViscousAir
BASE=ROOT/'outputs/v258_validation';CASE=BASE/'rocket01';PREV=ROOT/'outputs/v257_validation/rocket01';CLOSED=True


def relative(a,b):
    numerator=float(np.linalg.norm(a-b));denominator=float(np.linalg.norm(b))
    return dict(absoluteL2=numerator,referenceL2=denominator,relativeL2=numerator/denominator if denominator else None)


def physical(values):
    a=np.array(values,copy=True);a[:,4]-=300000*a[:,0];return a


def audit():
    start=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    receipt=read(BASE/'process_control/report.json');native=read(CASE/'native_completion.json');step=read(CASE/'step.json')
    assert receipt['passed'] and receipt['inputsStable'] and receipt['status']=='terminal' and not receipt['liveGroupAfterCleanup']
    assert native['normalReturn'] and native['fullPhysicalEvaluations']==4 and native['nonlinearUpdates']==native['timeUpdates']==0
    assert native['initialBitwiseV257'] and native['terminalIsRestoredInitial'] and not native['KINSOLExecuted']
    for suffix in ['state.bin','full_residual.bin']:
        assert sha(CASE/('initial_'+suffix))==sha(PREV/('initial_'+suffix))==sha(CASE/('terminal_restored_'+suffix))
    states={k:matrix(CASE/(k+'_state.bin'),b'RAH256S1',19) for k in ['initial','sigma','half']}
    residuals={k:matrix(CASE/(k+'_full_residual.bin'),b'RAH256R1',7) for k in ['initial','sigma','half']}
    s=states['initial'];wall=s[:,17]==1;assert len(s)==335726 and wall.sum()==9878
    mask=np.ones((len(s),6),bool);mask[np.ix_(wall,[1,2,3,5])]=False
    x0=vector(PREV/'initial_unknowns.bin');trial=vector(PREV/'linear_refused_trial.bin');d=vector(CASE/'saved_direction.bin')
    assert np.array_equal(trial-x0,d) and np.array_equal(np.asarray(s[:,5:11])[mask],x0)
    scales=np.genfromtxt(PREV/'scaling.csv',delimiter=',',names=True);char=scales['variableCharacteristic'];fchar=scales['residualCharacteristic']
    component=np.broadcast_to(np.arange(6),mask.shape)[mask];us=1/char[component];uds=x0*us;dds=d*us
    dot=float(uds@dds);vv=float(dds@dds);l1=float(np.abs(dds).sum());sigma=(1 if dot>=0 else -1)*np.sqrt(np.finfo(float).eps)*max(abs(dot),l1)/vv
    assert abs(sigma/step['sigma']-1)<1e-10 and abs(vv/step['weightedDirectionSquaredNorm']-1)<1e-10
    sigma=step['sigma'];air=PhaseAwareViscousAir();stateReports={}
    for label in ['sigma','half']:
        t=states[label];r=residuals[label];assert np.array_equal(t[:,:5],s[:,:5]) and np.array_equal(t[:,17:],s[:,17:])
        assert np.array_equal(t[:,0],r[:,0]);assert np.isfinite(t).all() and np.isfinite(r).all()
        assert np.all(t[wall,6:9]==0) and np.all(t[wall,10]==0)
        gas=air.evaluate(t[:,12],t[:,11]);transport=air.interpolated_transport(t[:,12],t[:,11]);rho=t[:,5]
        kinetic=np.sum(t[:,6:9]**2,axis=1)/(2*rho**2);energy=rho*(gas['internalEnergyJkg']+kinetic)
        checks=dict(rho=float(np.max(abs(gas['densityKgm3']/rho-1))),energy=float(np.max(abs(energy-t[:,9])/np.maximum(1,abs(t[:,9])))),
            mu=float(np.max(abs(transport['dynamicViscosityPas']/t[:,14]-1))),k=float(np.max(abs(transport['conductivityWmK']/t[:,13]-1))))
        nu=t[:,10];chi=nu/(t[:,14]/rho);mut=np.zeros(len(t));positive=nu>=0;mut[positive]=rho[positive]*nu[positive]*chi[positive]**3/(chi[positive]**3+7.1**3)
        checks['SA']=float(np.max(abs(mut-t[:,15])/np.maximum(1e-30,np.maximum(abs(mut),abs(t[:,15])))))
        checks['sameStateSAFlow']=float(np.max(abs(t[:,15]-t[:,16])/np.maximum(1e-30,np.maximum(abs(t[:,15]),abs(t[:,16])))))
        assert max(checks.values())<2e-10
        intended=x0+(sigma if label=='sigma' else sigma/2)*d;actual=np.asarray(t[:,5:11])[mask]
        mappingError=float(np.max(abs(actual-intended)/np.maximum(1,abs(intended))));assert mappingError<=1e-10
        stateReports[label]=dict(physicalChecks=checks,maximumCoordinateRoundtripError=mappingError,
            relativeActualPerturbationError=relative(actual-x0,intended-x0),temperatureRangeK=[float(t[:,11].min()),float(t[:,11].max())])
    block=matrix(CASE/'direction_matrix.bin',b'RAH258B1',13);ilu=matrix(CASE/'ilu_matrix.bin',b'RAH258B1',13)
    assert np.array_equal(block[:,0],s[:,0]) and np.array_equal(ilu[:,0],s[:,0]) and np.isfinite(block).all() and np.isfinite(ilu).all()
    nativeD=block[:,1:7];Ad=physical(block[:,7:13]);z=physical(ilu[:,1:7]);Az=physical(ilu[:,7:13]);roundtrip=physical(nativeD)[mask]
    mappingError=float(np.max(abs(roundtrip-d)/np.maximum(1,abs(d))));assert mappingError<1e-10
    strong={label:float(np.max(abs(values[np.ix_(wall,[1,2,3,5])]))) for label,values in [('direction',nativeD),('A_direction',block[:,7:13]),('ILU_solution',ilu[:,1:7]),('A_ILU_solution',ilu[:,7:13])]}
    # 强约束量应为零；极小非零也保留原值，不能以投影回填消除问题。
    strongExact=all(v==0 for v in strong.values())
    f0=np.asarray(residuals['initial'][:,1:]);J=(residuals['sigma'][:,1:]-f0)/sigma;half=(residuals['half'][:,1:]-f0)/(sigma/2)
    fd=[relative(J[:,j],half[:,j]) for j in range(6)];approx=[relative(Ad[:,j],half[:,j]) for j in range(6)]
    fdScaled=relative((J/fchar)[mask],(half/fchar)[mask]);approxScaled=relative((Ad/fchar)[mask],(half/fchar)[mask])
    # 只有固定方向与两步长这一局部结论，不转授全Jacobian或原线性求解资格。
    consistent=all((r['relativeL2']<=1e-3 if r['relativeL2'] is not None else r['absoluteL2']==0) for r in fd+[fdScaled])
    originalNativeF=f0.copy();originalNativeF[:,4]+=300000*originalNativeF[:,0]
    iluNative=[relative(ilu[:,7+j],-originalNativeF[:,j]) for j in range(6)]
    iluPhysical=[relative(Az[:,j],-f0[:,j]) for j in range(6)]
    iluScaled=relative((Az/fchar)[mask],(-f0/fchar)[mask])
    regions={}
    for name,selected in [('wall',wall),('interior',~wall)]:
        regions[name]=dict(finiteDifference=[relative(J[selected,j],half[selected,j]) for j in range(6)],
            approximateMatrix=[relative(Ad[selected,j],half[selected,j]) for j in range(6)])
    faces=np.loadtxt(ROOT/'outputs/v250_validation/rocket01/original_faces.txt');inverse=np.empty(len(s),int);inverse[s[:,0].astype(int)]=np.arange(len(s));nodes=inverse[faces[:,1:4].astype(int)]
    assert faces.shape==(19752,5) and np.all(wall[nodes]) and abs(faces[:,4].sum()-1.9646292565715553)<2e-13
    np.savez_compressed(out/'surface_values.npz',faceIDs=faces[:,0].astype(int),areaM2=faces[:,4],
        energyTrueDirection=np.mean(half[nodes,4],axis=1),energyApproxDirection=np.mean(Ad[nodes,4],axis=1),
        energyDirectionDifference=np.mean((Ad-half)[nodes,4],axis=1))
    report=dict(passed=True,diagnosticCompleted=True,baseCommit='29d7f95a2a01ba61cdf3db46799e88fab96f94ca',nativeSummary=native,
        sigma=sigma,independentSigmaRelativeError=abs((1 if dot>=0 else -1)*np.sqrt(np.finfo(float).eps)*max(abs(dot),l1)/vv/sigma-1),
        probeStates=stateReports,maximumDirectionCoordinateRoundtripError=mappingError,strongRowsMaximumAbs=strong,strongRowsExactlyZero=strongExact,
        finiteDifferencePerPhysicalEquation=fd,finiteDifferenceScaled=fdScaled,directionConsistencyThreshold=1e-3,fixedDirectionHalfStepConsistent=consistent,
        approximateMatrixPerPhysicalEquation=approx,approximateMatrixScaled=approxScaled,significantApproximationDifferenceThreshold=.1,
        significantApproximationDifference=bool(approxScaled['relativeL2']>.1),regionDiagnostics=regions,
        ILUApproxMatrixResidualPerNativeEquation=iluNative,ILUApproxMatrixResidualPerPhysicalEquation=iluPhysical,ILUApproxMatrixResidualScaled=iluScaled,
        sourceSpecificInterfaceBugProved=False,entireJacobianAccuracyProved=False,uniqueNumericalFailureCauseProved=False,
        initialAndRestoredStateResidualBitwiseV257=True,nonlinearUpdates=0,newQualifiedAreaM2=0,physicalQualification=False,profileReadyForCFD=False,
        originalFaces=19752,originalAreaM2=1.9646292565715553,
        inputSHA256={p.name:sha(p) for p in CASE.glob('*') if p.is_file()},elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);plot(out,report)
    print(json.dumps({k:report[k] for k in ['passed','fixedDirectionHalfStepConsistent','finiteDifferenceScaled','approximateMatrixScaled','ILUApproxMatrixResidualScaled','strongRowsMaximumAbs','elapsedSeconds']}))


def plot(out,r):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels=['Mass','Mom X','Mom Y','Mom Z','Energy','SA'];x=np.arange(6)
    fig,axes=plt.subplots(2,1,figsize=(10,6.7),layout='constrained')
    for ax,key,title,threshold in [(axes[0],'finiteDifferencePerPhysicalEquation','Same saved direction: sigma versus half sigma',1e-3),
        (axes[1],'approximateMatrixPerPhysicalEquation','Native block matrix versus full coupled F directional response',.1)]:
        values=[v['relativeL2'] for v in r[key]];ax.semilogy(x,values,'o-',color='#32597a');ax.axhline(threshold,color='#a44d43',linestyle='--',label=f'Diagnostic threshold {threshold:g}')
        ax.set(xticks=x,xticklabels=labels,ylabel='Relative L2 difference',title=title);ax.grid(alpha=.2);ax.legend()
    fig.suptitle('v258 | fixed rejected trial direction | diagnostic only | no accepted flow update')
    fig.savefig(out/'direction_consistency.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v258 fixed direction audit closed; follow NEXT_DECISION')
    audit()
