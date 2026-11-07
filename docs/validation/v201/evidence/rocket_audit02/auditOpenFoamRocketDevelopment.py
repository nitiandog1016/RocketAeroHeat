"""一次只读核对真实候选与完整早期场；保存失败而不以研究值冒充外缘。"""
from pathlib import Path
import json,re,shutil,time,csv
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from readOpenFoamRocketState import read_frame,BASE,read_json
from openFoamAirThermo import evaluate
from openFoamFields import internal_field,read_state
from openFoamConservation import boundary_values,conservative,flux_and_residual,native_mass_flux
from auditOpenFoamDevelopment import log_progress


def initial_update(case,folder,g,dt):
    """模块一：由准确均匀初态和原生壁投影独立核对首步五守恒量更新。"""
    n=len(g['cells']);ni=len(g['neighbor']);p=internal_field(case/'0/p',n);T=internal_field(case/'0/T',n);u=internal_field(case/'0/U',n,True)
    gas=evaluate(p,T);f=dict(p=p,T=T,U=u,rho=gas['densityKgm3'],e=gas['internalEnergyJkg'])
    ub=boundary_values(case/'0','U',u,g)
    # 所有初态内部面两侧相同，边界两侧取同一原生边界状态；Kurganov退化为物理通量。
    uf=np.r_[u[g['owner'][:ni]],ub];pf=p[g['owner']];rf=f['rho'][g['owner']];ef=f['e'][g['owner']]
    phiv=np.einsum('ij,ij->i',uf,g['sf']);mass=rf*phiv
    momentum=mass[:,None]*uf+pf[:,None]*g['sf']
    energy=(rf*(ef+np.einsum('ij,ij->i',uf,uf)/2)+pf)*phiv
    flux=np.column_stack([mass,momentum,energy]);summed=np.zeros((n,5))
    np.add.at(summed,g['owner'],flux);np.add.at(summed,g['neighbor'],-flux[:ni])
    expected=conservative(f)-dt*summed/g['volume'][:,None];actual=conservative(read_state(folder,n))
    scale=np.array([f['rho'][0],*([f['rho'][0]*np.linalg.norm(u[0])]*3),f['rho'][0]*(abs(f['e'][0])+np.linalg.norm(u[0])**2)])
    difference=np.max(abs(actual-expected)/scale,axis=0)
    phi_error=float(np.max(abs(native_mass_flux(folder,g)-mass))/(scale[0]*np.linalg.norm(u[0])*g['area'].max()))
    integrated=np.sum((actual-conservative(f))*g['volume'][:,None],axis=0)+dt*flux[ni:].sum(axis=0)
    balance=integrated/(scale*g['volume'].sum())
    return dict(conservativeUpdateMaxScaled=difference.tolist(),nativeLaggedMassFluxMaxScaled=phi_error,
        discreteGlobalBalanceScaled=balance.tolist(),passed=bool(difference.max()<2e-10 and phi_error<2e-10 and abs(balance).max()<2e-10)),scale


def main():
    start=time.monotonic();out=BASE/'rocket_audit02';out.mkdir(exist_ok=False)
    for name in ['auditOpenFoamRocketDevelopment.py','readOpenFoamRocketState.py','openFoamConservation.py']:
        shutil.copyfile(ROOT/'scripts'/name,out/name)
    run=read_json(BASE/'rocket_run01/report.json');process=run['process'];case=BASE/'rocket01'
    assert not run['passed'] and process['stopReason']=='scientific-guard' and process['guardReason']['token']=='out of temperature range'
    assert all(process['inputStable'].values()) and not process['liveGroupAfterCleanup'] and process['logExportPassed']
    text=(BASE/'rocket_run01/native/stdout.txt').read_text(encoding='utf-8');progress=log_progress(text,inspect_gaps=True)
    first_dt=float(re.search(r'^deltaT = ([\deE.+-]+)',text,re.M).group(1))
    current=None;step=0;warnings=[]
    for line in text.splitlines():
        match=re.match(r'^Time = ([\deE.+-]+)s',line)
        if match:current=float(match.group(1));step+=1
        match=re.search(r'out of temperature range 200 -> 2500;\s+T = ([\deE.+-]+)',line)
        if match:warnings.append(dict(announcedStep=step,timeSeconds=current,requestedTemperatureK=float(match.group(1))))
    assert warnings
    names=sorted(run['saved'],key=float);assert len(names)==1
    meta,fields,g,m,ids,owners=read_frame(names[0],float(names[0]),diagnostic=True)
    update,scale=initial_update(case,case/names[0],g,first_dt)
    f,flux,residual=flux_and_residual(case/names[0],g)
    budget=read_json(BASE/'budget01/report.json');length=budget['lengthReferenceM'];speed=budget['speedReferenceMps'];volume=g['volume']
    norms=np.sum(abs(residual)*volume[:,None],axis=0)/volume.sum()/scale*length/speed
    integrated=(residual*volume[:,None]).sum(axis=0);boundary_flux=flux[len(g['neighbor']):].sum(axis=0)
    closure=(integrated-boundary_flux)/(scale*speed*length**2)
    h0inf=float(evaluate(1197.0279007383137,226.50908361133003)['enthalpyJkg'])+speed**2/2
    gas=evaluate(f['p'],f['T']);h0=gas['enthalpyJkg']+np.einsum('ij,ij->i',f['U'],f['U'])/2
    all_dirs=[]
    for d in case.iterdir():
        if d.is_dir():
            try:t=float(d.name)
            except ValueError:continue
            all_dirs.append(dict(name=d.name,timeSeconds=t,files=sorted(p.name for p in d.iterdir() if p.is_file())))
    # 模块二：准确异常时序/可读状态范围，不从警告逆算请求推断最终全场最低温度。
    result=dict(auditPassed=True,planAcceptance=False,candidateClosed=True,normalCompletion=False,
        progress=progress,warningCount=len(warnings),firstWarning=warnings[0],
        minimumRequestedTemperatureK=min(w['requestedTemperatureK'] for w in warnings),
        requestedTemperatureMeaning='native inverse-iteration request, not complete-field minimum or unique failure cause',
        savedDirectories=all_dirs,savedState=meta,firstNativeDtSeconds=first_dt,
        initialCourantRelativeError=abs(run['initialMaxCo']/budget['initialMaxCoAtDt1e7']-1),
        earlyUpdate=update,residualVolumeL1Scaled=norms.tolist(),residualThreshold=1e-8,
        residualScale=dict(referenceLengthM=length,referenceSpeedMps=speed,conservedScale=scale.tolist(),definition='volume L1 divergence / conservedScale * L/U'),
        boundaryFlux=boundary_flux.tolist(),boundaryVsVolumeClosureScaled=closure.tolist(),
        totalEnthalpyRelativeMax=float(abs(h0/h0inf-1).max()),
        wallPressureRangePa=[float(fields['pressurePa'].min()),float(fields['pressurePa'].max())],
        wallFluidTemperatureRangeK=[float(fields['fluidTemperatureK'].min()),float(fields['fluidTemperatureK'].max())],
        areaStability=None,residualWindowImprovement=None,numericallyConverged=False,attachedEdgeQualified=False,
        physicalQualification=False,newQualifiedHeatAreaM2=0,processSeconds=process['elapsedSeconds'],entrySeconds=run['elapsedSeconds'],
        sourceReportSHA256=sha(BASE/'rocket_run01/report.json'),elapsedSeconds=time.monotonic()-start)
    assert update['passed'] and abs(closure).max()<2e-10
    np.savez(out/'original_wall_state.npz',faceIDs=m['originalFaceIDs'],areaM2=g['area'][ids],wallCentersM=m['wallCentersM'],
        samplePositionsM=m['samplePositionsM'],nativeOwnerCellsZeroBased=owners,**fields)
    save(out/'state_metadata.json',meta);save(out/'report.json',result);save(out/'warning_requests.json',warnings)
    with (out/'original_wall.csv').open('w',newline='',encoding='utf-8') as stream:
        w=csv.writer(stream);w.writerow(['faceID','x_m','y_m','z_m','area_m2','p_Pa','Tfluid_K','Ux_mps','Uy_mps','Uz_mps','native_es_Jkg','engineering_h0_Jkg','q_Wm2','solid_Tw_K'])
        for i in range(len(ids)):w.writerow([i+1,*m['wallCentersM'][i],g['area'][ids[i]],fields['pressurePa'][i],fields['fluidTemperatureK'][i],*fields['velocityMps'][i],fields['nativeSensibleInternalEnergyJkg'][i],fields['totalEnthalpyJkg'][i],'',''])
    print(json.dumps({k:result[k] for k in ['auditPassed','progress','firstWarning','minimumRequestedTemperatureK','earlyUpdate','residualVolumeL1Scaled','elapsedSeconds']}))


if __name__=='__main__':main()
