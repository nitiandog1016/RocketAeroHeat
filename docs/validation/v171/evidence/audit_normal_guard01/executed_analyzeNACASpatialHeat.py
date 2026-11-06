"""模块用途：空间Tw层流保存态的严格热审核、原电极段功率和原节点外缘诊断。"""
import argparse
import csv
from pathlib import Path
import time
import numpy as np

from idealAirProperties import IdealAir
from naca3344Model import D,R,SHOULDER
from nacaBoundaryLayerEdge import extract_edge
from nacaMeasurementSegments import segments,area,reduce_segment
from nacaSpatialWallFlux import geometry,evaluate,read_temperature_table,temperature_at,segment_integrals,verify_applied
from su2NACALaminarInterface import mesh,load_state,inlet_state_check
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import save,paper_cp


def column_indices(meta,points,arrays_path=None):
    """模块一：只使用原节点列，不以最近邻插值掩盖测量中心缺失。"""
    shape=(len(meta['offsets']),len(meta['s']),len(meta['azimuth']))
    if arrays_path is not None and Path(arrays_path).exists():
        with np.load(str(arrays_path)) as a:
            assert np.array_equal(a['points'],points)
            labels=a['nodeIndices'].copy();ix=a['columnNodeIDs'].copy()
    else:
        labels=np.array(meta['nodeIndices'],int);ix=np.full(shape,-1,int)
        for node,(j,i,k) in enumerate(labels):ix[j,i,k]=node
        ix[:,0,:]=ix[:,0,:1]
    assert ix.shape==shape and labels.shape==(len(points),3) and ix.min()>=0 and ix.max()<len(points)
    assert len(np.unique(ix))==len(points)
    for node,(j,i,k) in enumerate(labels):assert ix[j,i,k]==node
    return labels,ix


def column_edge(meta,points,values,thermo,ix,station_s,azimuth_index,source,state_kind='computed-solution'):
    """模块二：把同一真实法向列送到受限外缘提取器，边缘定义仍可拒绝。"""
    s=np.asarray(meta['s']);matches=np.flatnonzero(abs(s-station_s)<1e-13)
    if len(matches)!=1:
        return dict(accepted=False,status='measurement-center-not-an-original-column',
            requestedStationM=float(station_s),closestOriginalDistanceM=float(np.min(abs(s-station_s))),
            physicalQualification=False,experimentalQualification=False,freestreamFallbackUsed=False)
    i=int(matches[0]);k=int(azimuth_index);ids=ix[:,i,k]
    phi=meta['azimuth'][k];theta=min(s[i]/R,np.pi/2)
    normal=np.array([np.cos(theta),np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi)])
    tangent=np.array([-np.sin(theta),np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi)])
    rho=values['Density'][ids]
    velocity=np.column_stack([values['Momentum_'+axis][ids]/rho for axis in 'xyz'])
    profile=dict(distanceM=np.asarray(meta['offsets']),coordinatesM=points[ids],velocityMS=velocity,
        pressurePa=values['Pressure'][ids],temperatureK=values['Temperature'][ids],densityKgM3=rho,
        enthalpyJkg=thermo['h'][ids],specificHeatJkgK=thermo['cp'][ids],
        viscosityPas=values['Laminar_Viscosity'][ids],
        conductivityWmK=values['Laminar_Viscosity'][ids]*thermo['cp'][ids]/.72)
    result=extract_edge(profile,normal,tangent,float(s[i]),
        dict(stateKind=state_kind,flowStateSource=source,normalColumnSource='original-mesh-column-node-IDs',
             thermodynamicModel='Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1'))
    result.update(originalNodeIDs=ids.tolist(),sIndex=i,azimuthIndex=k,azimuthRadians=float(phi))
    return result


def measurement_analysis(meta,p,values,thermo,ix,g,flux,given_temperature,source):
    """模块三：实际对偶功率按20原段分别积分，中心/原物性缺口不会消失。"""
    rows,summary=segments()
    bounds=np.array([[float(r[name]) for name in ('s_left_m','s_right_m')] for r in rows])
    if bounds.max()>np.max(g['patchS'])+1e-14:
        return dict(status='wall-does-not-cover-original-measurement-segments',experimentalQualification=False)
    actual=segment_integrals(g,flux,bounds,meta['azimuth'][-1]-meta['azimuth'][0])
    out=[];edge_counts={};accepted=0
    for row,item in zip(rows,actual):
        a,b=item['boundsM'];nominal=float(area(a,b));power=item['symmetryScaledPowerIntoWallFullW']
        item.update(row=int(row['heat_row_1based']),nominalFullAreaM2=nominal,
            meanQByNominalAreaWm2=power/nominal,
            discreteVsNominalAreaRelativeDifference=item['discreteGeometricAreaSectorM2']*item['fullRevolutionScale']/nominal-1,
            strictLocalAbscissaIdentityClosed=row['strict_local_abscissa_identity_closed']=='True',
            measuredParameter=float(row['original_Nu_s_over_sqrt_Re_s']),experimentalQualification=False)
        centers=dict(plot=float(row['original_plot_s_over_D'])*D,
            oppositeTemperatureStation=float(row['temperature_station_s_over_D'])*D,
            voltageIntervalMidpoint=float(row['voltage_interval_arc_midpoint_s_over_D'])*D)
        reduced=dict(row=item['row'],meanQIntoWallWm2=item['meanQByNominalAreaWm2'],centerCandidatesS_m=centers,
            strictLocalAbscissaIdentityClosed=item['strictLocalAbscissaIdentityClosed'],measuredParameter=item['measuredParameter'])
        details={}
        for name,center in centers.items():
            # 所有原周向列逐一留痕；不挑选更接近实验的一条。
            columns=[]
            for az in range(len(meta['azimuth'])):
                result=column_edge(meta,p,values,thermo,ix,center,az,source)
                edge_counts[result['status']]=edge_counts.get(result['status'],0)+1
                if result['accepted']:
                    accepted+=1;i=result['sIndex'];node=int(ix[0,i,az]);tw=float(given_temperature(node))
                    result['givenTwAtSameOriginalCenterK']=tw
                    if item['meanQByNominalAreaWm2']<0:
                        try:
                            reduction=reduce_segment(reduced,name,tw,result['edge'],result['method'])
                            reduction.update(discreteWallPowerMapped=True,geometryPowerMethod=item['method'],
                                symmetryAssumption='quarter-domain mirrored to full revolution; quarter variation separately retained')
                            result['thermalReduction']=reduction
                        except AssertionError as exc:
                            result['thermalReduction']=dict(status='reduction-domain-rejected',reason=str(exc),experimentalQualification=False)
                    else:result['thermalReduction']=dict(status='heated-experiment-sign-not-satisfied',experimentalQualification=False)
                columns.append(result)
            details[name]=columns
        item['centerColumns']=details;out.append(item)
    return dict(status='discrete-segment-power-mapped',segments=out,
        nativeNodalPowerDecomposition=True,originalTransportTableEquivalent=False,
        all20CenterIdentitiesClosed=summary['strict_all20_local_abscissa_identity_closed'],
        edgeColumnAcceptedDiagnostics=accepted,edgeColumnStatusCounts=edge_counts,
        overlappingSegmentPowersAreNotTotal=True,experimentalQualification=False,newQualifiedArea_m2=0)


def audit(folder,out,saved_frames_only=False):
    """模块四：严格审核保存帧；异常任务须显式选择研究模式，绝不授予正常完成。"""
    started=time.monotonic();folder=Path(folder).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    receipt=read(folder/'engine_process.json');manifest=read(folder/'case_manifest.json')
    from nacaSavedFrameStudy import select_saved_frames
    checkpoints,identity=select_saved_frames(folder,receipt,saved_frames_only)
    files=[folder/name for name in ('engine_process.json','case_manifest.json','mesh.su2','geometry.json','case.cfg','history.csv')]
    files += [Path(__file__).resolve(),ROOT/'scripts/nacaSpatialWallFlux.py',ROOT/'scripts/nacaBoundaryLayerEdge.py',
              ROOT/'scripts/nacaMeasurementSegments.py',ROOT/'scripts/su2NACALaminarInterface.py',
              ROOT/'scripts/nacaSavedFrameStudy.py',ROOT/'scripts/su2NACASpatialInterface.py']
    locks={str(path):sha(path) for path in files}
    for path in files:
        if path.suffix=='.py':(out/('executed_'+path.name)).write_bytes(path.read_bytes())
    p,cells,markers=mesh(folder/'mesh.su2');meta=read(folder/'geometry.json')
    arrays_path=folder/'mesh_arrays.npz'
    if arrays_path.exists():locks[str(arrays_path)]=sha(arrays_path)
    labels,ix=column_indices(meta,p,arrays_path);node_s=np.array(meta['s'])[labels[:,1]]
    g=geometry(p,cells,markers,node_s)
    from su2IdealAirInterface import settings
    cfg=settings((folder/'case.cfg').read_text());assert float(cfg['PRANDTL_LAM'])==.72
    assert cfg['KIND_TURB_MODEL']=='NONE' and cfg['AXISYMMETRIC']=='NO'
    c=manifest['publicBenchmark']['conditions'];boundary=inlet_state_check((folder/'case.cfg').read_text(),c,manifest['expectedVelocity'])
    table_path=folder/'wall_temperature.csv';applied=None
    if table_path.exists():
        table=read_temperature_table(table_path);locks[str(table_path)]=sha(table_path)
        given=lambda node:temperature_at(p[node,0],table)
        tw=given(g['ids']);applied_path=folder/'rah_wall_temperature_applied.csv'
        applied=verify_applied(applied_path,p,g,tw);locks[str(applied_path)]=sha(applied_path)
    else:
        values=[token.strip() for token in cfg['MARKER_ISOTHERMAL'].strip('() ').split(',')]
        assert len(values)==2 and values[0]=='wall'
        scalar=float(values[1]);given=lambda node:np.full(np.asarray(node).shape,scalar)
        tw=given(g['ids'])
    with (folder/'history.csv').open() as f:
        history={int(float(row['Inner_Iter'])):row for row in
            [{k.strip():v for k,v in r.items()} for r in csv.DictReader(f,skipinitialspace=True)]}
    core=IdealAir();locks.update(core.locks);frames=[]
    for pair in checkpoints['complete']:
        path=folder/pair['volume'];assert sha(path)==pair['volumeSHA256'];locks[str(path)]=sha(path)
        data,v=load_state(path);assert len(data)==len(p) and np.max(abs(data[:,1:4]-p))<2e-13
        rho=v['Density'];u=np.column_stack([v['Momentum_'+axis]/rho for axis in 'xyz'])
        thermo=core.state('PT',v['Pressure'],v['Temperature']);thermal=v['Energy']/rho-.5*np.sum(u*u,axis=1)
        energy=float(np.max(abs(thermal-thermo['e'])/np.maximum(1,abs(thermo['e']))))
        eos=float(np.max(abs(rho/thermo['rho']-1)));cp=float(np.max(abs(paper_cp(v['Temperature'])/thermo['cp']-1)))
        viscosity=1.716e-5*(v['Temperature']/273.15)**1.5*(273.15+110.4)/(v['Temperature']+110.4)
        mu_error=float(np.max(abs(v['Laminar_Viscosity']/viscosity-1)))
        flux=evaluate(g,v['Temperature'],thermo['cp']*v['Laminar_Viscosity']/.72,tw)
        ids=g['ids'];q=flux['qIntoWall'];q_error=float(np.max(abs(q-v['Heat_Flux'][ids])/np.maximum(1,abs(q))))
        step=pair['localIteration'];h=history[step];native_power=float(h['HF']);power_error=abs(flux['powerIntoWall']-native_power)/max(1,abs(flux['powerIntoWall']))
        cf=np.column_stack([v['Skin_Friction_Coefficient_'+axis][ids] for axis in 'xyz'])
        shear=np.linalg.norm(cf,axis=1)*.5*c['density_kg_m3']*c['speed_m_s']**2
        # 层流原生Y_Plus为零；用真实首邻距离、同阶段壁rho/mu与实际Cf另算尺度。
        yp=g['euclideanDistance']*np.sqrt(rho[ids]*shear)/v['Laminar_Viscosity'][ids]
        filename='wall_spatial_%06d.csv'%step
        columns=['PointID','X_m','Y_m','Z_m','normalNeighbor','givenTw_K','wallFluidT_K','normalDistance_m',
                 'euclideanFirstDistance_m','normalDistanceFloorTriggered','nativeProjectedArea_m2','qIntoWall_Wm2',
                 'nodalPowerIntoWall_W','independentFirstNeighborYPlus','method','physicalQualified']
        with (out/filename).open('w',encoding='utf8',newline='') as f:
            writer=csv.writer(f);writer.writerow(columns)
            writer.writerows([[int(node),*p[node],int(near),float(twi),float(v['Temperature'][node]),float(dn),float(de),bool(floor),
                float(ar),float(qi),float(power),float(yplus),flux['method'],False] for node,near,twi,dn,de,floor,ar,qi,power,yplus in zip(
                    ids,g['normalNeighbor'],tw,g['distance'],g['euclideanDistance'],g['distanceFloorTriggered'],g['area'],q,flux['nodalPowerIntoWall'],yp)])
        measured=measurement_analysis(meta,p,v,thermo,ix,g,flux,given,str(path)+'#'+sha(path))
        measured_file='measurement_%06d.json'%step;save(out/measured_file,measured)
        residuals={key:float(h[key]) for key in ('rms[Rho]','rms[RhoE]')}
        passed=max(energy,eos,mu_error)<2e-11 and cp<2e-8 and q_error<2e-10 and power_error<2e-9
        frames.append(dict(localIteration=step,interfacePassed=bool(passed),energyMaximumRelativeError=energy,
            eosMaximumRelativeError=eos,paperCpMaximumRelativeError=cp,viscosityMaximumRelativeError=mu_error,
            nativeHeatMaximumScaledError=q_error,powerScaledDifference=power_error,
            nativeNodalPowerIntoWallW=flux['powerIntoWall'],nativeHistoryPowerW=native_power,
            geometricAreaSectorM2=g['discreteGeometricArea'],nativeProjectedAreaSectorM2=g['nativeProjectedArea'],
            givenTwMaximumFluidDifferenceK=float(np.max(abs(tw-v['Temperature'][ids]))),
            independentFirstNeighborYPlus99=float(np.percentile(yp,99)),independentFirstNeighborYPlusMaximum=float(yp.max()),
            independentFirstNeighborYPlusDefinition='Euclidean original normal-neighbor distance * sqrt(rho_wall*abs(tau_native))/mu_wall',
            normalDistanceFloorTriggered=int(g['distanceFloorTriggered'].sum()),residuals=residuals,
            originalResidualGatePassed=bool(max(residuals.values())<=-8),heatStabilityQualified=False,
            experimentalQualification=False,wallCSV=filename,wallCSV_SHA256=sha(out/filename),
            measurementJSON=measured_file,measurementSHA256=sha(out/measured_file),
            acceptedEdgeDiagnostics=measured.get('edgeColumnAcceptedDiagnostics',0)))
    stable={name:sha(Path(name))==value for name,value in locks.items()};assert all(stable.values())
    saved_passed=bool(frames and all(f['interfacePassed'] for f in frames))
    identity.update(thermalInterfaceAudited=True,savedFrameInterfacePassed=saved_passed)
    result=dict(interfacePassed=bool(identity['normalRunCompleted'] and saved_passed),
        savedFrameInterfacePassed=saved_passed,auditScope=identity['auditScope'],
        normalRunCompleted=identity['normalRunCompleted'],savedFrameIdentity=identity,frames=frames,
        appliedTemperatureCheck=applied,inletStateCheck=boundary,inputSHA256=locks,inputsStable=stable,
        stage='same-flow-update-laminar-no-SST',nativeLaminarYPlusIsQualification=False,
        givenTwIsSolvedSolidTemperature=False,gridAndDomainQualified=False,experimentalQualification=False,
        newQualifiedArea_m2=0,elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',result);print(dict(interfacePassed=result['interfacePassed'],
        savedFrameInterfacePassed=saved_passed,normalRunCompleted=result['normalRunCompleted'],frames=len(frames)))
    assert saved_passed, '严格空间壁热接口未过；保存失败，不重算CFD'
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',required=True);parser.add_argument('--output',required=True)
    parser.add_argument('--saved-frames-only',action='store_true',help='显式审核已终态异常任务的完整帧；不补造异常瞬间或正常完成')
    args=parser.parse_args();audit(args.folder,args.output,args.saved_frames_only)
