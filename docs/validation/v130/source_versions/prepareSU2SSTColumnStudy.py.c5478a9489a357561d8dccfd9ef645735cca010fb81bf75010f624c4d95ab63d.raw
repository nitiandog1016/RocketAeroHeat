"""完整箭体近壁k/omega扩散—耗散初猜；保持热内能，原SST和数值参数不变。"""
from pathlib import Path
import argparse,copy,json,os
import numpy as np
from scipy.linalg import solve_banded
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from prepareSU2SSTCheckpoint import OUTPUT_CONTROL,validate_source_artifacts
from prepareSU2SSTCase import set_value
from su2SSTInterface import state

KIND='same-mesh-SST-wall-column-initialization-study'
STUDY='v130-near-wall-diffusion-destruction-initial-guess'
SOURCE=ROOT/'outputs/v127_validation/control'
PLAN=ROOT/'docs/validation/v130/PLAN.md'


def column_solution(y,nu,k0,kend,w0,wend,kfloor=0.):
    """无对流/生产的近壁极限；正边界的M矩阵离散保持内部正值。"""
    y=np.asarray(y,float);assert len(y)>=3 and y[0]==0 and np.all(np.diff(y)>0)
    assert nu>0 and k0>=kfloor>=0 and kend>=kfloor and kend>0 and w0>0 and wend>0
    length=y[-1];x=y/length;scale=nu/(.075*length**2)
    lo=2/(np.diff(x)[:-1]*(x[2:]-x[:-2]));hi=2/(np.diff(x)[1:]*(x[2:]-x[:-2]));diag=lo+hi
    w=np.linspace(w0/scale,wend/scale,len(y));rhs=np.zeros(len(y)-2)
    rhs[0]+=lo[0]*w[0];rhs[-1]+=hi[-1]*w[-1]
    def residual(z):return lo*(z[1:-1]-z[:-2])+hi*(z[1:-1]-z[2:])+z[1:-1]**2
    def local_error(z):
        scale=np.maximum(z[1:-1]**2,lo*abs(z[1:-1]-z[:-2])+hi*abs(z[2:]-z[1:-1]))
        return float(max(abs(residual(z))/scale))
    for it in range(80):
        f=residual(w);den=max(float(max(diag*w[1:-1])),float(max(abs(rhs))),1.)
        err=float(max(abs(f))/den)
        # 壁值跨多个数量级，整列尺度不能替代每个节点的方程检查。
        localError=local_error(w)
        if err<2e-12 and localError<2e-12:break
        band=np.zeros((3,len(diag)));band[1]=diag+2*w[1:-1];band[0,1:]=-hi[:-1];band[2,:-1]=-lo[1:]
        delta=solve_banded((1,1),band,-f);step=1.
        while step>=2**-30:
            trial=w.copy();trial[1:-1]+=step*delta
            if min(trial)>0 and local_error(trial)<localError:break
            step*=.5
        if step<2**-30:
            # 归一化分母随大幅更新变化时，改用本次固定权重的Newton下降范数。
            weights=np.maximum(w[1:-1]**2,lo*abs(w[1:-1]-w[:-2])+hi*abs(w[2:]-w[1:-1]))
            step=1.
            while step>=2**-30:
                trial=w.copy();trial[1:-1]+=step*delta
                if min(trial)>0 and max(abs(residual(trial))/weights)<localError:break
                step*=.5
        assert step>=2**-30,'正值Newton未找到下降步'
        w=trial
    else:raise AssertionError('近壁初猜未在固定80次非线性迭代内收敛')
    assert err<1e-10
    # k方程是固定omega下的线性反应—扩散，不采用经验y平方补值。
    band=np.zeros((3,len(diag)));band[1]=diag+(.09/.075)*w[1:-1];band[0,1:]=-hi[:-1];band[2,:-1]=-lo[1:]
    b=np.zeros(len(diag));b[0]+=lo[0]*k0;b[-1]+=hi[-1]*kend
    kval=solve_banded((1,1),band,b);active=kval<kfloor
    # 保留SU2原下限，不在后处理直接裁剪而破坏相邻自由方程。
    for _ in range(80):
        constrained=band.copy();rhs=b.copy()
        for j in np.flatnonzero(active):
            constrained[1,j]=1.;rhs[j]=kfloor
            if j<len(diag)-1:constrained[0,j+1]=0.
            if j>0:constrained[2,j-1]=0.
        kval=solve_banded((1,1),constrained,rhs);kval[active]=kfloor
        if np.any((kval<kfloor)&~active):
            active[np.argmin(np.where(active,np.inf,kval-kfloor))]=True;continue
        k=np.r_[k0,kval,kend]
        gradient=diag*kval-lo*k[:-2]-hi*k[2:]+(.09/.075)*w[1:-1]*kval
        if np.any(active&(gradient<-1e-12*max(float(max(abs(b))),1.))):
            active[np.argmin(np.where(active,gradient,np.inf))]=False;continue
        break
    else:raise AssertionError('k下限互补问题未收敛')
    k=np.r_[k0,kval,kend]
    fk=diag*k[1:-1]-lo*k[:-2]-hi*k[2:]+(.09/.075)*w[1:-1]*k[1:-1]
    defect=np.where(active,np.minimum(fk,0),fk)
    ke=float(max(abs(defect))/max(float(max(abs(b))),float(max(diag*k[1:-1])),1.))
    assert min(k[1:])>0 and ke<1e-10
    result=w*scale;result[0]=w0;result[-1]=wend
    return k,result,err,ke,it


def reconstruct(a,keys,c):
    """按已有32层列和壁剪固定选择近壁范围，不按下限值选择修复节点。"""
    nw=9878;layers=32;assert len(a)==349206
    xyz=a[:(layers+1)*nw,1:4].reshape(layers+1,nw,3)
    heights=np.linalg.norm(xyz-xyz[0],axis=2)
    mr=read(ROOT/'outputs/v120_validation/mesh32/mesh_report.json')
    assert mr['normalLayers']==layers and max(abs(heights-np.array(mr['normalHeightsM'])[:,None]).ravel())<2e-12
    rho=a[:nw,keys.index('Density')];mu=a[:nw,keys.index('Laminar_Viscosity')];nu=mu/rho
    R=float(c['GAS_CONSTANT']);T=float(c['FREESTREAM_TEMPERATURE']);gamma=float(c['GAMMA_VALUE'])
    rhoi=float(c['FREESTREAM_PRESSURE'])/(R*T);U=float(c['MACH_NUMBER'])*np.sqrt(gamma*R*T)
    cf=np.linalg.norm(a[:nw,[keys.index('Skin_Friction_Coefficient_'+v) for v in 'xyz']],axis=1)
    yp=heights*np.sqrt(.5*rhoi*U**2*cf/rho)[None,:]/nu[None,:]
    k=a[:,keys.index('Turb_Kin_Energy')];w=a[:,keys.index('Omega')]
    ids=[];ks=[];ws=[];anchors=[];we=[];ke=[];its=[]
    for j in range(nw):
        crossing=np.flatnonzero(yp[1:13,j]>=1)
        end=int(crossing[0]+1) if len(crossing) else 12
        anchors.append(end)
        if end<2:continue
        ix=np.arange(end+1)*nw+j
        kk,ww,ew,ek,it=column_solution(heights[:end+1,j],nu[j],k[ix[0]],k[ix[-1]],w[ix[0]],w[ix[-1]],1e-10)
        ids.extend(ix[1:-1]);ks.extend(kk[1:-1]);ws.extend(ww[1:-1]);we.append(ew);ke.append(ek);its.append(it)
    ids=np.array(ids,int);order=np.argsort(ids);ids=ids[order];ks=np.array(ks)[order];ws=np.array(ws)[order]
    assert len(np.unique(ids))==len(ids) and min(ids)>=nw and max(ids)<12*nw
    assert min(ks)>=1e-10 and min(ws)>1e-4
    report=dict(nodesChanged=int(len(ids)),wallColumns=nw,columnsChanged=int(sum(np.array(anchors)>=2)),
        anchorLayerCounts={str(i):anchors.count(i) for i in sorted(set(anchors))},
        maximumOmegaEquationRelativeError=max(we),maximumKEquationRelativeError=max(ke),maximumNewtonIterations=max(its),
        newKRange=[float(min(ks)),float(max(ks))],newOmegaRange=[float(min(ws)),float(max(ws))],
        kNativeLowerBound=1e-10,kLowerBoundActiveNodes=int(sum(ks==1e-10)),
        kEquationContract='free-node equation and active-node nonnegative residual complementarity',
        initialGuessOnly=True,physicsQualified=False,sourcePrimitiveColumnsRetainedForProvenanceOnly=True)
    return ids,ks,ws,np.array(anchors),report


def validate_initial(folder,source,init):
    """逐字段绑定初场变换，保持rho、动量、未改节点和热内能。"""
    a,keys,_,rho,_,energy=state(source/'restart_flow.csv',True,True)
    b,newkeys,_,_,_,newenergy=state(folder/'initial_flow.csv',True,True);assert keys==newkeys
    report=read(folder/'column_initialization.json');npz=folder/'column_values.npz'
    assert sha(folder/'column_initialization.json')==init['columnReportSHA256'] and sha(npz)==report['valuesSHA256']
    values=np.load(str(npz));ids=values['nodeIDs'];kcol=keys.index('Turb_Kin_Energy');wcol=keys.index('Omega');ecol=keys.index('Energy')
    assert len(ids)==report['nodesChanged'] and np.all(np.diff(ids)>0)
    assert np.array_equal(b[ids,kcol],values['newK']) and np.array_equal(b[ids,wcol],values['newOmega'])
    for j in range(len(keys)):
        if j not in [kcol,wcol,ecol]:assert np.array_equal(a[:,j],b[:,j]),keys[j]
    off=np.ones(len(a),bool);off[ids]=False;assert np.array_equal(a[off],b[off])
    expected=a[ids,ecol]+rho[ids]*(b[ids,kcol]-a[ids,kcol]);assert np.array_equal(b[ids,ecol],expected)
    assert max(abs(newenergy/energy-1))<2e-12
    # 按源几何/剪切重算选择和BVP，杜绝仅篡改清单使任意正数过关。
    ri,rk,rw,ra,rr=reconstruct(a,keys,config(source/'case.cfg'))
    assert np.array_equal(ri,ids) and np.array_equal(rk,values['newK']) and np.array_equal(rw,values['newOmega'])
    assert np.array_equal(ra,values['anchorLayers'])
    assert all(report[k]==rr[k] for k in rr)


def validate_source(folder,c,init):
    source=Path(init['sourceDirectory']).resolve()
    assert source==SOURCE and init['kind']==KIND and init['studyID']==STUDY
    assert init['variant'] in ['control','column'] and init['sourceLocalIteration']==149 and init['localIterationStartsAt']==0
    assert init['sameMesh'] and init['samePhysicsAndNumerics'] and not init['physicsQualified']
    assert init['changedField']==('k-omega-and-compatible-rhoE-in-near-wall-columns' if init['variant']=='column' else 'none')
    a=config(source/'case.cfg')
    assert {k:v for k,v in c.items() if k not in OUTPUT_CONTROL}=={k:v for k,v in a.items() if k not in OUTPUT_CONTROL}
    assert c['ITER']=='100' and c['OUTPUT_WRT_FREQ']=='50' and c['VOLUME_OUTPUT']=='(COORDINATES, SOLUTION, PRIMITIVE)'
    locked=validate_source_artifacts(folder,init,validate_initial if init['variant']=='column' else None)
    assert sha(PLAN)==init['planSHA256'];locked[str(PLAN)]=sha(PLAN);locked[str(Path(__file__).resolve())]=sha(Path(__file__).resolve())
    if init['variant']=='column':
        boundplan=ROOT/'docs/validation/v130/NATIVE_FLOOR.md';assert sha(boundplan)==init['nativeFloorPlanSHA256'];locked[str(boundplan)]=sha(boundplan)
        localplan=ROOT/'docs/validation/v130/LOCAL_RESIDUAL.md';assert sha(localplan)==init['localResidualPlanSHA256'];locked[str(localplan)]=sha(localplan)
        for name in ['column_initialization.json','column_values.npz']:locked[str(folder/name)]=sha(folder/name)
    return locked


def prepare(variant):
    folder=ROOT/'outputs/v130_validation'/('column_final' if variant=='column' else variant);assert not folder.exists();folder.mkdir(parents=True)
    os.link(str(SOURCE/'mesh.su2'),str(folder/'mesh.su2'));text=(SOURCE/'case.cfg').read_text()
    for key,value in dict(ITER=100,OUTPUT_WRT_FREQ=50,SCREEN_WRT_FREQ_INNER=10).items():text=set_value(text,key,value)
    (folder/'case.cfg').write_text(text,encoding='utf-8')
    if variant=='control':os.link(str(SOURCE/'restart_flow.csv'),str(folder/'initial_flow.csv'))
    else:
        a,keys,_,rho,_,e=state(SOURCE/'restart_flow.csv',True,True)
        ids,kk,ww,anchors,report=reconstruct(a,keys,config(SOURCE/'case.cfg'))
        kcol=keys.index('Turb_Kin_Energy');wcol=keys.index('Omega');ecol=keys.index('Energy')
        ne=a[ids,ecol]+rho[ids]*(kk-a[ids,kcol]);mapping={int(i):(k,w,E) for i,k,w,E in zip(ids,kk,ww,ne)}
        with (SOURCE/'restart_flow.csv').open('rb') as old,(folder/'initial_flow.csv').open('wb') as new:
            new.write(next(old))
            for i,line in enumerate(old):
                if i in mapping:
                    k,w,E=mapping[i];parts=line.split(b',')
                    for j,v in [(kcol,k),(wcol,w),(ecol,E)]:parts[j]=(' %.17e'%v).encode('ascii')
                    line=b','.join(parts)
                new.write(line)
        np.savez_compressed(str(folder/'column_values.npz'),nodeIDs=ids,newK=kk,newOmega=ww,anchorLayers=anchors)
        _,_,_,_,_,newenergy=state(folder/'initial_flow.csv',True,True)
        report.update(valuesSHA256=sha(folder/'column_values.npz'),sourceSHA256=sha(SOURCE/'restart_flow.csv'),
            thermalInternalEnergyMaxRelativeDifference=float(max(abs(newenergy/e-1))),
            rhoEnergyAdjustmentL1=float(np.sum(abs(ne-a[ids,ecol]))),rhoEnergyAdjustmentIsInitialStateNotPhysicalHeating=True)
        assert report['thermalInternalEnergyMaxRelativeDifference']<2e-12
        (folder/'column_initialization.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    locks=['engine_process.json','sst_state_audit.json','case.cfg','mesh.su2','sst_contract.json','case_manifest.json',
        'restart_flow.csv','surface_flow.vtk','history.csv','sst_initialization.json']
    init=dict(kind=KIND,studyID=STUDY,variant=variant,changedField='k-omega-and-compatible-rhoE-in-near-wall-columns' if variant=='column' else 'none',
        sourceDirectory=str(SOURCE),sourceLocalIteration=149,sourceWasDiagnostic=False,localIterationStartsAt=0,
        sourceStateHasNativePrimitive=True,sameMesh=True,samePhysicsAndNumerics=True,
        restartCFLReset=True,uninterruptedTrajectoryEquivalent=False,physicsQualified=False,
        initialStateSHA256=sha(folder/'initial_flow.csv'),sourceSHA256={str(SOURCE/k):sha(SOURCE/k) for k in locks},planSHA256=sha(PLAN))
    mr=ROOT/'outputs/v120_validation/mesh32/mesh_report.json';init['sourceSHA256'][str(mr)]=sha(mr)
    if variant=='column':
        init['columnReportSHA256']=sha(folder/'column_initialization.json')
        init['nativeFloorPlanSHA256']=sha(ROOT/'docs/validation/v130/NATIVE_FLOOR.md')
        init['localResidualPlanSHA256']=sha(ROOT/'docs/validation/v130/LOCAL_RESIDUAL.md')
    validate_source(folder,config(folder/'case.cfg'),init)
    (folder/'sst_initialization.json').write_text(json.dumps(init,indent=2)+'\n',encoding='utf-8')
    contract=copy.deepcopy(read(SOURCE/'sst_contract.json'));contract.pop('numericalStudy',None)
    contract.update(initializationSHA256=sha(folder/'sst_initialization.json'),configurationSHA256=sha(folder/'case.cfg'),
        initializationStudy={k:init[k] for k in ['studyID','variant','changedField']})
    (folder/'sst_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    m=copy.deepcopy(read(SOURCE/'case_manifest.json'));m['flight'].update(iterations=100,outputEvery=50,namedCheckpoints=True)
    m['sourceRunControlsFromConfig']={k:config(SOURCE/'case.cfg')[k] for k in ['ITER','OUTPUT_WRT_FREQ']}
    m.update(initialization=init,sst=contract,sstContractSHA256=sha(folder/'sst_contract.json'))
    (folder/'case_manifest.json').write_text(json.dumps(m,indent=2)+'\n',encoding='utf-8')
    req=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=360,dryRun=False,
        sourceFiles=['scripts/prepareSU2SSTColumnStudy.py','docs/validation/v130/PLAN.md'])
    (folder/'engine_request.json').write_text(json.dumps(req,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(variant=variant,folder=str(folder),newLocalLastIteration=99)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('variant',choices=['control','column']);prepare(p.parse_args().variant)
