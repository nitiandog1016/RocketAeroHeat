"""只读M6两份已存场，区分压力偏置、壁热变化和未闭合入口；不运行求解器。"""
from pathlib import Path
import csv,json
import numpy as np
from prepareSU2CheckpointContinuation import ROOT,read,sha
from compareSU2DNSPlate import reference,viscosity,station_bracket
from su2SSTInterface import state


def main():
    folder=ROOT/'outputs/v129_validation/dns01/plate'
    out=ROOT/'outputs/v130_validation/dns_pressure_review.json';assert not out.exists()
    engine=read(folder/'engine_process.json');comparison=read(folder.parent/'comparison03/comparison.json')
    ref,_=reference();station=comparison['stationM'];cases={};locks={}
    saved={}
    with (folder/'history.csv').open() as f:rows=list(csv.reader(f,skipinitialspace=True))
    h=dict(zip([x.strip() for x in rows[0]],np.array(rows[1:],float).T))
    for step in [5000,9999]:
        path=folder/('restart_flow_%06d.csv'%step)
        assert sha(path)==engine['outputSHA256'][path.name];locks[str(path)]=sha(path)
        a,keys,nd,rho,uv,_=state(path,True,True);assert nd==2
        ss=np.unique(-a[:,1]);ys=np.unique(a[:,2]);profiles=[];wall=[]
        for s in ss:
            ids=np.flatnonzero(-a[:,1]==s);ids=ids[np.argsort(a[ids,2])];assert np.array_equal(a[ids,2],ys)
            U=-uv[ids,0]/ref['speedMs'];theta=np.trapz(rho[ids]/ref['densityKgM3']*U*(1-U),ys)
            ret=ref['densityKgM3']*ref['speedMs']*theta/viscosity(ref['temperatureK'])
            iw=ids[0];p=a[ids,keys.index('Pressure')]/ref['pressurePa']
            profiles.append(np.column_stack([p,a[ids,keys.index('Temperature')]/ref['temperatureK'],U,uv[ids,1]/ref['speedMs']]))
            wall.append([s,ret,a[iw,keys.index('Heat_Flux')],p[0]])
        profiles=np.array(profiles);wall=np.array(wall);saved[step]=(wall,profiles)
        i=int(np.searchsorted(ss,station))-1;f=(station-ss[i])/(ss[i+1]-ss[i]);fixed=(1-f)*profiles[i]+f*profiles[i+1]
        j,g=station_bracket(ss,wall[:,1],ref['ReTheta'],.35);matched=(1-g)*wall[j]+g*wall[j+1]
        samples=[]
        for y in [0,.0036,.01,.02,.04]:
            value=[float(np.interp(y,ys,fixed[:,k])) for k in range(4)]
            samples.append(dict(yM=y,pressureRatio=value[0],temperatureRatio=value[1],streamwiseSpeedRatio=value[2],normalSpeedRatio=value[3]))
        use=(ss>=.05*.35)&(ss<=.85*.35)
        cases[str(step)]=dict(fixedStationM=station,fixedStationWallHeatWm2=float((1-f)*wall[i,2]+f*wall[i+1,2]),
            fixedStationProfiles=samples,independentlyMatchedStationM=float(matched[0]),
            matchedHeatWm2=float(matched[2]),matchedWallPressureRatio=float(matched[3]),
            interiorWallPressureRatioRange=[float(min(wall[use,3])),float(max(wall[use,3]))],
            interiorTopPressureRatioRange=[float(min(profiles[use,-1,0])),float(max(profiles[use,-1,0]))],
            nativeLogResidual={k:float(h[k][step]) for k in ['rms[Rho]','rms[RhoE]','rms[k]','rms[w]']})
    a,b=saved[5000],saved[9999];use=a[0][:,0]>=0
    qa,qb=a[0][use,2],b[0][use,2];ss=a[0][use,0]
    qchange=float(np.trapz(abs(qb-qa),ss)/np.trapz(abs(qa),ss))
    assert abs(cases['9999']['matchedHeatWm2']/comparison['heatFluxIntoWallWm2']-1)<1e-12
    report=dict(savedFields=cases,streamwiseWallHeatRelativeL1Change=qchange,
        maximumPrimitivePressureRatioChange=float(max(abs(a[1][:,:,0]-b[1][:,:,0]).ravel())),
        exactDNSInletAvailable=False,causalBoundaryOrGridDiagnosisEstablished=False,solverRuns=0,
        originalQualificationChanged=False,physicsQualified=False,
        interpretation='Two saved fields show changes, not a separation of inlet, boundary, grid and unconverged errors.',
        sourceSHA256=dict(locks,**{str(Path(__file__).resolve()):sha(Path(__file__).resolve()),
                                 str(folder/'history.csv'):sha(folder/'history.csv')}))
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
