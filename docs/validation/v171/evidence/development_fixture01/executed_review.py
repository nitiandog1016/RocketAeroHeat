"""模块用途：只读评估实际层流保存帧发展；分开研究统计、严格热审核和物理资格。"""
import argparse
import csv
from pathlib import Path
import sys
import time
import numpy as np

from analyzeSU2DualVolumes import FACES
from idealAirProperties import IdealAir
from naca3344Model import R,SHOULDER
from nacaSpatialWallFlux import geometry,evaluate,read_temperature_table,temperature_at,segment_integrals
from su2NACALaminarInterface import mesh,load_state
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import save

RESIDUALS=('rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]','rms[RhoE]')


def wls_operators(points,cells,selected):
    """模块一：真实单元棱边上的加权最小二乘，不声称重放原生历史梯度。"""
    selected=np.unique(np.asarray(selected,int));mask=np.zeros(len(points),bool);mask[selected]=True
    neighbors={int(i):set() for i in selected}
    for kind in (12,13):
        ids=np.array([row[1:] for row in cells if row[0]==kind],int)
        assert ids.ndim==2 and len(ids)
        ids=ids[np.any(mask[ids],axis=1)]
        edges=sorted({tuple(sorted((a,b))) for f in FACES[kind] for a,b in zip(f,f[1:]+f[:1])})
        for a,b in edges:
            for i,j in ids[:,[a,b]]:
                if mask[i]:neighbors[int(i)].add(int(j))
                if mask[j]:neighbors[int(j)].add(int(i))
    operators={};maximum_condition=0.;moment_error=0.
    for i,adjacent in neighbors.items():
        ids=np.array(sorted(adjacent),int);delta=points[ids]-points[i]
        assert len(ids)>=3 and np.all(np.sum(delta*delta,axis=1)>0)
        weighted=delta/np.sum(delta*delta,axis=1)[:,None]
        matrix=weighted.T@delta;condition=float(np.linalg.cond(matrix))
        assert np.isfinite(condition) and condition<1e14,'局部WLS退化，不能构造物理通量诊断'
        op=np.linalg.solve(matrix,weighted.T);err=float(np.max(abs(op@delta-np.eye(3))))
        assert err<2e-8,'WLS仿射再现未通过'
        operators[i]=(ids,op);maximum_condition=max(maximum_condition,condition);moment_error=max(moment_error,err)
    return operators,dict(nodes=len(selected),maximumCondition=maximum_condition,affineMomentMaximumError=moment_error,
        method='cell-edge inverse-squared-distance WLS diagnostic; not native stored gradient')


def gradients(operators,values):
    grad=np.zeros((len(values),values.shape[1],3))
    for i,(ids,op) in operators.items():grad[i]=(op@(values[ids]-values[i])).T
    return grad


def polygon_flux(points,rows,flux):
    """模块二：有向多边形按面心三角片积分节点线性物理通量。"""
    result=[]
    for size in (3,4):
        ids=np.array([r[1:] for r in rows if len(r)==size+1],int)
        if not len(ids):continue
        p=points[ids];center=p.mean(axis=1);f=flux[ids];fc=f.mean(axis=1);total=np.zeros(len(ids))
        for j in range(size):
            k=(j+1)%size
            area_vector=np.cross(p[:,k]-p[:,j],center-p[:,j])/2
            total+=np.einsum('ij,ij->i',area_vector,(f[:,j]+f[:,k]+fc)/3)
        result.extend(total.tolist())
    return np.asarray(result)


def boundary_geometry(points,markers):
    """完整外边界方向/闭合的只读量；不是流动收敛检查。"""
    area_sum=0.;area_vector=np.zeros(3);volume=0.
    for rows in markers.values():
        for size in (3,4):
            ids=np.array([r[1:] for r in rows if len(r)==size+1],int)
            if not len(ids):continue
            p=points[ids];center=p.mean(axis=1)
            for j in range(size):
                av=np.cross(p[:,(j+1)%size]-p[:,j],center-p[:,j])/2
                area_vector+=av.sum(axis=0);area_sum+=float(np.linalg.norm(av,axis=1).sum())
                volume+=float(np.einsum('ij,ij->i',av,center).sum()/3)
    closure=float(np.linalg.norm(area_vector)/area_sum)
    assert volume>0 and closure<2e-10
    return dict(outwardBoundaryVolumeM3=volume,areaVectorRelativeClosure=closure,totalBoundaryAreaM2=area_sum)


def physical_balance(points,markers,operators,values,thermo,wall_power):
    """模块三：层流Newton应力/Pr导热的独立物理边界通量，原生离散通量仍未知。"""
    rho=values['Density'];u=np.column_stack([values['Momentum_'+a]/rho for a in 'xyz'])
    grad=gradients(operators,np.column_stack((u,values['Temperature'])))
    mu=values['Laminar_Viscosity'];g=grad[:,:3]
    stress=mu[:,None,None]*(g+g.transpose(0,2,1))
    stress-=(2*mu*np.trace(g,axis1=1,axis2=2)/3)[:,None,None]*np.eye(3)
    conductivity=mu*thermo['cp']/.72
    mass=rho[:,None]*u
    energy=(values['Energy']+values['Pressure'])[:,None]*u-np.einsum('nij,nj->ni',stress,u)-conductivity[:,None]*grad[:,3]
    by_marker={};all_mass=[];all_energy=[]
    for tag,rows in markers.items():
        if tag=='wall':continue
        mf=polygon_flux(points,rows,mass);ef=polygon_flux(points,rows,energy)
        by_marker[tag]=dict(massNetKgS=float(mf.sum()),massGrossKgS=float(abs(mf).sum()),energyNetW=float(ef.sum()),
            energyGrossW=float(abs(ef).sum()))
        all_mass.extend(mf);all_energy.extend(ef)
    mf=np.asarray(all_mass);ef=np.asarray(all_energy)
    mass_in=float(np.maximum(-mf,0).sum());energy_in=float(np.maximum(-ef,0).sum())
    netm=float(mf.sum());nete=float(ef.sum()+wall_power)
    return dict(byMarker=by_marker,massNetKgS=netm,massIncomingKgS=mass_in,
        massRelativeImbalance=abs(netm)/mass_in if mass_in>0 else None,
        exteriorEnergyNetW=float(ef.sum()),wallPowerOutOfFluidW=float(wall_power),energyNetW=nete,energyIncomingW=energy_in,
        energyRelativeImbalance=abs(nete)/energy_in if energy_in>0 else None,
        method='independent cell-edge WLS; laminar Newtonian stress, Pr=.72 conduction; face-center linear flux quadrature',
        allExteriorMarkersIncluded=True,symmetryLeakageReported=True,energyReference='OTH rhoE+p; no SST k or eddy viscosity',
        wallTerm='native weak-wall qIntoWall times native projected nodal area; impermeable wall assumed',
        nativeDiscreteBoundaryFluxAvailable=False,nativeConservationProven=False,physicalQualified=False)


def l1_change(area,first,last):
    """模块四：有符号量取绝对差，分母为首态绝对积分；零分母保持未定义。"""
    area=np.asarray(area,float);first=np.asarray(first,float);last=np.asarray(last,float)
    assert area.shape==first.shape==last.shape and np.all(area>0) and np.isfinite(area+first+last).all()
    denominator=float(area@abs(first));numerator=float(area@abs(last-first))
    return dict(numerator=numerator,referenceAbsoluteIntegral=denominator,
        relativeL1=numerator/denominator if denominator>0 else None,zeroReference=denominator==0)


def history_statistics(path):
    """保持真实有限历史位置，末20窗不足时不伪造收敛统计。"""
    with Path(path).open() as f:raw=list(csv.reader(f,skipinitialspace=True))
    assert raw;keys=[k.strip() for k in raw[0]];position=keys.index('Inner_Iter')
    rows=[];rejected=0
    for row in raw[1:]:
        try:
            a=np.array(row,float)
            if len(a)!=len(keys) or not np.isfinite(a).all() or a[position]!=int(a[position]):raise ValueError()
            rows.append(a)
        except (ValueError,OverflowError):rejected+=1
    data=np.asarray(rows);assert data.ndim==2 and len(data)>0
    steps=data[:,position].astype(int);assert np.array_equal(steps,np.arange(len(steps)))
    windows={};power=data[:,keys.index('HF')]
    for name,start,stop in [('first20',0,20),('previous20',len(data)-40,len(data)-20),('last20',len(data)-20,len(data))]:
        if start<0 or stop>len(data):windows[name]=None;continue
        sample=data[start:stop];p=power[start:stop];den=float(np.mean(abs(p)))
        windows[name]=dict(firstLocalIteration=int(steps[start]),lastLocalIteration=int(steps[stop-1]),
            meanLogResidual={n:float(sample[:,keys.index(n)].mean()) for n in RESIDUALS},
            heatPowerMeanW=float(p.mean()),heatPowerRangeW=[float(p.min()),float(p.max())],
            heatPowerRelativeRange=float(np.ptp(p)/den) if den>0 else None,zeroPowerScale=den==0)
    delta={}
    for ref in ('first20','previous20'):
        delta[ref]=({n:windows['last20']['meanLogResidual'][n]-windows[ref]['meanLogResidual'][n] for n in RESIDUALS}
                    if windows['last20'] is not None and windows[ref] is not None else None)
    final={n:float(data[-1,keys.index(n)]) for n in RESIDUALS}
    return dict(finiteHistoryRows=len(data),lastFiniteLocalIteration=int(steps[-1]),rejectedNonfiniteOrIncompleteRows=rejected,
        windows=windows,lastWindowMinusReferenceLog=delta,finalLogResidual=final,
        originalDensityEnergyMinus8Passed=bool(max(final[n] for n in ('rms[Rho]','rms[RhoE]'))<=-8),
        historyPowerSign='positive into wall',historyUsedAsPhysicalQualification=False)


def main(folder,audit_folder,out):
    """模块五：只分析终态已保存完整帧；异常瞬间、缺失审计和实验资格不补造。"""
    started=time.monotonic();folder=Path(folder).resolve();out=Path(out).resolve();audit_folder=Path(audit_folder).resolve()
    receipt=read(folder/'engine_process.json');manifest=read(folder/'case_manifest.json')
    assert receipt['status']=='terminal' and all(receipt['inputsStable'].values()),'真实任务未终态或输入已变化'
    out.mkdir(parents=True,exist_ok=False);(out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    normal=bool(receipt['exitCode']==0 and receipt['checkpoints']['normalScheduleComplete'] and not any(
        receipt.get(k,False) for k in ('timedOut','cancelled','outputBudgetExceeded','supervisorError','memoryBudgetExceeded')))
    entries=receipt['checkpoints']['complete'];assert entries,'没有可分析的完整保存帧'
    input_paths=[folder/name for name in ('engine_process.json','case_manifest.json','mesh.su2','geometry.json','mesh_arrays.npz',
                                        'wall_temperature.csv','case.cfg','history.csv')]
    input_paths += [Path(__file__).resolve(),ROOT/'scripts/nacaSpatialWallFlux.py',ROOT/'scripts/nacaWallFlux.py',
                   ROOT/'scripts/su2NACALaminarInterface.py',ROOT/'scripts/analyzeSU2DualVolumes.py']
    locks={str(p):sha(p) for p in input_paths}
    audit_path=audit_folder/'report.json';audit=None
    if audit_path.exists():
        audit=read(audit_path);locks[str(audit_path)]=sha(audit_path)
        assert all(audit.get('inputsStable',{}).values()),'已有热审核来源变化'
    audit_frames={f['localIteration']:f for f in audit.get('frames',[])} if audit else {}
    p,cells,markers=mesh(folder/'mesh.su2');meta=read(folder/'geometry.json')
    with np.load(str(folder/'mesh_arrays.npz')) as a:
        assert np.array_equal(a['points'],p);labels=a['nodeIndices'];ix=a['columnNodeIDs']
    node_s=np.asarray(meta['s'])[labels[:,1]];g=geometry(p,cells,markers,node_s);ids=g['ids']
    exterior=np.unique([i for tag,rows in markers.items() if tag!='wall' for row in rows for i in row[1:]])
    operators,wls_check=wls_operators(p,cells,exterior);boundary=boundary_geometry(p,markers)
    table=read_temperature_table(folder/'wall_temperature.csv');tw=temperature_at(p[ids,0],table)
    core=IdealAir();locks.update(core.locks);c=manifest['publicBenchmark']['conditions']
    frames=[];fields=[];measurements=[];sector=meta['azimuth'][-1]-meta['azimuth'][0]
    for pair in entries:
        path=folder/pair['volume'];assert sha(path)==pair['volumeSHA256'];locks[str(path)]=sha(path)
        _,v=load_state(path);assert len(v['Density'])==len(p)
        rho=v['Density'];u=np.column_stack([v['Momentum_'+axis]/rho for axis in 'xyz'])
        thermo=core.state('PT',v['Pressure'],v['Temperature']);flux=evaluate(g,v['Temperature'],thermo['cp']*v['Laminar_Viscosity']/.72,tw)
        q=flux['qIntoWall'];pressure=v['Pressure'][ids];step=pair['localIteration'];audited=audit_frames.get(step)
        cf=np.column_stack([v['Skin_Friction_Coefficient_'+axis][ids] for axis in 'xyz'])
        shear=np.linalg.norm(cf,axis=1)*.5*c['density_kg_m3']*c['speed_m_s']**2
        yp=g['euclideanDistance']*np.sqrt(rho[ids]*shear)/v['Laminar_Viscosity'][ids]
        delta=v['Temperature'][ids]-tw;area=g['area']
        regions=segment_integrals(g,flux,[[0,SHOULDER],[SHOULDER,float(max(meta['s']))]],sector)
        for name,region in zip(('hemisphere','cylinder'),regions):
            region['name']=name
            lower,upper=region['boundsM'];selected=(node_s[ids]>=lower-1e-14)&(node_s[ids]<=upper+1e-14)
            region.update(wallNodesIncludingSharedEndpoints=int(selected.sum()),pressureNodeRangePa=[float(pressure[selected].min()),float(pressure[selected].max())],
                signedHeatNodeRangeWm2=[float(q[selected].min()),float(q[selected].max())],
                givenTwFluidMaximumDifferenceK=float(abs(delta[selected]).max()),independentFirstNeighborYPlus99=float(np.quantile(yp[selected],.99)))
        measurement=None
        if audited and audited.get('measurementJSON'):
            mp=audit_folder/audited['measurementJSON'];assert sha(mp)==audited['measurementSHA256'];locks[str(mp)]=sha(mp)
            m=read(mp);measurement=dict(localIteration=step,status=m.get('status'),acceptedEdgeDiagnostics=m.get('edgeColumnAcceptedDiagnostics',0),
                edgeColumnStatusCounts=m.get('edgeColumnStatusCounts',{}),all20CenterIdentitiesClosed=m.get('all20CenterIdentitiesClosed',False),
                originalTransportTableEquivalent=m.get('originalTransportTableEquivalent',False),experimentalQualification=False,segments=[])
            for item in m.get('segments',[]):
                measurement['segments'].append({key:item.get(key) for key in ('row','boundsM','meanQByNominalAreaWm2',
                    'powerIntoWallSectorW','discreteVsNominalAreaRelativeDifference','strictLocalAbscissaIdentityClosed','measuredParameter')})
            measurements.append(measurement)
        frame=dict(localIteration=step,successfulUpdates=pair['successfulUpdates'],strictThermalAuditPassed=audited.get('interfacePassed') if audited else None,
            auditStatus='bound-to-existing-report' if audited else 'strict-thermal-audit-not-available',
            wallFluidTemperatureRangeK=[float(v['Temperature'][ids].min()),float(v['Temperature'][ids].max())],
            givenTwRangeK=[float(tw.min()),float(tw.max())],wallTemperatureMaximumDifferenceK=float(abs(delta).max()),
            wallTemperatureProjectedAreaRMSDifferenceK=float(np.sqrt(area@(delta*delta)/area.sum())),
            pressureNodeRangePa=[float(pressure.min()),float(pressure.max())],signedHeatNodeRangeWm2=[float(q.min()),float(q.max())],
            positiveIntoWallProjectedAreaFraction=float(area[q>0].sum()/area.sum()),negativeIntoWallProjectedAreaFraction=float(area[q<0].sum()/area.sum()),
            nativeProjectedAreaSectorM2=float(area.sum()),nativeNodalPowerIntoWallSectorW=flux['powerIntoWall'],
            domainTemperatureRangeK=[float(v['Temperature'].min()),float(v['Temperature'].max())],
            maximumWallSpeedMS=float(np.linalg.norm(u[ids],axis=1).max()),
            independentFirstNeighborYPlus99=float(np.quantile(yp,.99)),independentFirstNeighborYPlusMaximum=float(yp.max()),
            actualYPlusThresholdPassed=bool(np.quantile(yp,.99)<=1),nativeYPlusUsedForQualification=False,
            physicalBoundaryDiagnostic=physical_balance(p,markers,operators,v,thermo,flux['powerIntoWall']),
            regions=regions,measuredSegmentStatisticsAvailable=measurement is not None,physicalQualified=False,newQualifiedAreaM2=0)
        frames.append(frame);fields.append(dict(p=pressure.copy(),q=q.copy()));save(out/('frame_%06d.json'%step),frame)
    changes=[]
    for i in range(1,len(frames)):
        changes.append(dict(fromLocalIteration=frames[i-1]['localIteration'],toLocalIteration=frames[i]['localIteration'],
            successfulUpdateDifference=frames[i]['localIteration']-frames[i-1]['localIteration'],
            fields={name:l1_change(g['area'],fields[i-1][name],fields[i][name]) for name in ('p','q')}))
    first_to_last={name:l1_change(g['area'],fields[0][name],fields[-1][name]) for name in ('p','q')} if len(fields)>1 else None
    history=history_statistics(folder/'history.csv')
    assert history['lastFiniteLocalIteration']==receipt['checkpoints']['lastFiniteHistoryLocalIteration']
    stable={name:sha(Path(name))==h for name,h in locks.items()};assert all(stable.values())
    result=dict(reviewCompleted=True,normalTerminal=normal,processExitCode=receipt['exitCode'],
        processStatus=dict((key,receipt.get(key,False)) for key in ('timedOut','cancelled','outputBudgetExceeded','supervisorError','memoryBudgetExceeded')),
        completeSavedLocalIterations=[f['localIteration'] for f in frames],lastCompleteLocalIteration=frames[-1]['localIteration'],
        lastFiniteHistoryLocalIteration=history['lastFiniteLocalIteration'],lastHistoryHasFullState=frames[-1]['localIteration']==history['lastFiniteLocalIteration'],
        abnormalInstantStateAvailable=False,history=history,frames=frames,consecutiveSurfaceChanges=changes,
        firstToLastSurfaceChange=first_to_last,surfaceChangeArea='native projected nodal wall area; not smooth analytic area',
        boundaryGeometry=boundary,wlsGeometry=wls_check,measurementFrames=measurements,
        strictThermalAuditAvailable=audit is not None,strictAllSavedThermalChecksPassed=bool(audit and all(f['strictThermalAuditPassed'] for f in frames)),
        heatStabilityQualified=False,gridAndDomainQualified=False,experimentalQualification=False,physicalQualified=False,newQualifiedAreaM2=0,
        inputSHA256=locks,inputsStable=stable,elapsedSeconds=time.monotonic()-started,newCFDRuns=0,automaticAdditionalTasks=0)
    save(out/'report.json',result)
    print({key:result[key] for key in ('normalTerminal','completeSavedLocalIterations','lastFiniteHistoryLocalIteration','elapsedSeconds')})
    return result


def synthetic(out):
    """解析六面体仿射真值：验证梯度与有向物理通量，另测零分母和符号。"""
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    p=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
    # 加第二个远离的棱柱给混合类型分支独立真值；各自完整封闭，不声称同一连通求解域。
    wedge=np.array([[3,0,0],[4,0,0],[3,1,0],[3,0,1],[4,0,1],[3,1,1]],float)
    points=np.vstack((p,wedge));cells=[[12]+list(range(8)),[13]+list(range(8,14))]
    faces=[[9 if len(f)==4 else 5]+[i for i in f] for f in FACES[12]]
    faces.extend([[9 if len(f)==4 else 5]+[i+8 for i in f] for f in FACES[13]])
    markers={'outer':faces}
    op,wls=wls_operators(points,cells,np.arange(len(points)))
    velocity=np.column_stack((10+2*points[:,0],1-3*points[:,1],2+points[:,2]))
    temp=300+points@np.array([1.,2.,3.]);values=np.column_stack((velocity,temp))
    grad=gradients(op,values);wanted=np.vstack((np.diag([2.,-3.,1.]),[1,2,3]))
    error=float(np.max(abs(grad-wanted)));assert error<1e-12
    v=dict(Density=np.full(len(points),2.),Pressure=np.full(len(points),20.),Energy=np.full(len(points),100.),
           Temperature=temp,Laminar_Viscosity=np.full(len(points),.01))
    v.update({'Momentum_'+axis:2*velocity[:,j] for j,axis in enumerate('xyz')})
    balance=physical_balance(points,markers,op,v,dict(cp=np.full(len(points),1000.)),0.)
    # div(u)=0，常应力应变功为0.28 W/m³，体积1+1/2，线性T的导热净通量0。
    mass_error=abs(balance['massNetKgS']);energy_error=abs(balance['energyNetW']-(-.28*1.5))
    assert mass_error<1e-12 and energy_error<2e-11
    geo=boundary_geometry(points,markers);assert abs(geo['outwardBoundaryVolumeM3']-1.5)<1e-14
    sign=l1_change([1,2],[-1,1],[1,-1]);assert sign['relativeL1']==2.
    zero=l1_change([1,2],[0,0],[1,2]);assert zero['relativeL1'] is None and zero['zeroReference']
    report=dict(passed=True,gradientAffineMaximumError=error,massAnalyticAbsoluteError=mass_error,
        energyAnalyticAbsoluteError=energy_error,expectedEnergyNetW=-.42,mixedClosedVolumeM3=geo['outwardBoundaryVolumeM3'],
        signChangingRelativeL1=sign['relativeL1'],zeroDenominatorKeptUndefined=True,wls=wls,
        sourceSHA256=sha(Path(__file__).resolve()),syntheticThermodynamicInputsAreFluxAlgebraFixture=True,newPhysicalQualification=False)
    save(out/'report.json',report);print(report)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder');parser.add_argument('--audit')
    parser.add_argument('--output',required=True);parser.add_argument('--synthetic',action='store_true');args=parser.parse_args()
    if args.synthetic:synthetic(args.output)
    else:
        assert args.folder and args.audit
        main(args.folder,args.audit,args.output)
