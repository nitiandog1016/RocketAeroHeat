"""完整箭体原生隐式右端与伪时间的空间分析；积分量和局部变化率分列。"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.spatial import cKDTree
from prepareSU2CheckpointContinuation import ROOT,read,sha,config
from su2SSTInterface import state,mesh,triangle_distances,SPATIAL_OUTPUT


def median_volumes(points,cells):
    """四面体V/4与棱柱中位对偶子体积；不以单元平均体积代替壁节点体积。"""
    volume=np.zeros(len(points));totals={};tetquality=np.full(len(points),np.inf)
    faces=[[0,2,1],[3,4,5],[0,1,4,3],[1,2,5,4],[2,0,3,5]]
    for typ in [10,13]:
        ids=np.array([r[1:] for r in cells if r[0]==typ],int);q=points[ids]
        if typ==10:
            v=abs(np.einsum('ij,ij->i',q[:,1]-q[:,0],np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0])))/6
            assert min(v)>0
            edge2=sum(np.sum((q[:,j]-q[:,i])**2,axis=1) for i in range(4) for j in range(i+1,4))
            quality=12*(3*v)**(2/3)/edge2
            for k in range(4):
                np.add.at(volume,ids[:,k],v/4);np.minimum.at(tetquality,ids[:,k],quality)
        else:
            center=q.mean(axis=1);v=np.zeros(len(ids))
            for face in faces:
                fc=q[:,face].mean(axis=1)
                for a,b in zip(face,face[1:]+face[:1]):
                    midpoint=(q[:,a]+q[:,b])/2
                    for k in [a,b]:
                        p=q[:,k];sub=abs(np.einsum('ij,ij->i',np.cross(midpoint-p,fc-p),center-p))/6
                        assert min(sub)>0;np.add.at(volume,ids[:,k],sub);v+=sub
        totals[str(typ)]=dict(cells=len(ids),volumeM3=float(sum(v)),minimumCellVolumeM3=float(min(v)))
    assert len(cells)==sum(v['cells'] for v in totals.values()) and min(volume)>0
    return volume,tetquality,totals


def analyze():
    folder=ROOT/'outputs/v128_validation/spatial01';out=ROOT/'outputs/v128_validation/spatial_analysis.json'
    assert not out.exists();receipt=read(folder/'engine_process.json');c=config(folder/'case.cfg')
    assert receipt['status']=='terminal' and receipt['exitCode']==0 and all(receipt['inputsStable'].values())
    assert receipt['physicsModel']['spatialOutputContract']==SPATIAL_OUTPUT
    for name in ['restart_flow.csv','history.csv']:assert sha(folder/name)==receipt['outputSHA256'][name]
    a,keys,nd,rho,uv,e=state(folder/'restart_flow.csv',True,True,True);assert nd==3
    points,cells,wall=mesh(folder/'mesh.su2')
    # 原生CSV用15位小数的科学记数输出；按节点ID逐点检查文本舍入界，不做最近邻重排。
    coordinateError=float(np.max(abs(points-a[:,1:4])/np.maximum(abs(points),np.finfo(float).tiny)))
    assert coordinateError<8*np.finfo(float).eps
    n=len(a);ids=np.arange(n);nw=len(set(j for f in wall for j in f[1:]));assert nw==9878
    assert set(j for f in wall for j in f[1:])==set(range(nw))
    mr=read(ROOT/'outputs/v120_validation/mesh32/mesh_report.json');assert mr['nodes']==n and mr['normalLayers']==32
    collar=(mr['normalLayers']+1)*nw;layer=np.where(ids<collar,ids//nw,-1)
    volume,quality,totals=median_volumes(points,cells)
    expected=read(folder/'case_manifest.json')['fluidVolume'];error=abs(volume.sum()/expected-1);assert error<2e-10
    with (folder/'history.csv').open() as f:rows=list(csv.reader(f,skipinitialspace=True))
    history=dict(zip([s.strip() for s in rows[0]],map(float,rows[-1])));assert history['Inner_Iter']==0
    names=['Density','Momentum_x','Momentum_y','Momentum_z','Energy','TKE','Omega']
    histnames=['rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]','rms[RhoE]','rms[k]','rms[w]']
    residual=np.column_stack([a[:,keys.index('Residual_'+s)] for s in names])
    logs=np.log10(np.sqrt(np.mean(residual**2,axis=0)))
    discrepancy=max(abs(logs-np.array([history[k] for k in histnames])));assert discrepancy<2e-8
    dt=a[:,keys.index('Delta_Time')];tdt=a[:,keys.index('Turb_Delta_Time')];cfl=a[:,keys.index('CFL')]
    assert min(dt)>0 and np.array_equal(dt,tdt)
    assert abs(min(dt)/history['Min DT']-1)<1e-8 and abs(max(dt)/history['Max DT']-1)<1e-8
    assert abs(min(cfl)/history['Min CFL']-1)<1e-8 and abs(max(cfl)/history['Max CFL']-1)<1e-8
    masks={'wall':layer==0,'layers1to4':(layer>=1)&(layer<=4),'layers5to16':(layer>=5)&(layer<=16),
        'layers17to31':(layer>=17)&(layer<=31),'interface32':layer==32,'tetraOnly':layer==-1}
    assert np.all(sum(m.astype(int) for m in masks.values())==1)
    # 仅几何分带；鼻前候选与尖端附近分列，不冒充激波/分离分类器。
    axial={'aheadOfNose':points[:,0]>max(points[:nw,0]),'noseQuarter':(points[:,0]>=-.31)&(points[:,0]<=max(points[:nw,0])),
        'midBody':(points[:,0]>=-1.5)&(points[:,0]<-.31),'aftBody':(points[:,0]>=min(points[:nw,0]))&(points[:,0]<-1.5),
        'downstream':points[:,0]<min(points[:nw,0])}
    assert np.all(sum(m.astype(int) for m in axial.values())==1)
    total=(residual**2).sum(axis=0);rates=residual**2/volume[:,None];rateTotal=rates.sum(axis=0)
    def groups(mm):return {name:dict(nodes=int(sum(m)),volumeM3=float(sum(volume[m])),
        integratedRHSSquareFraction=dict(zip(names,((residual[m]**2).sum(axis=0)/total).tolist())),
        volumeWeightedRateSquareFraction=dict(zip(names,(rates[m].sum(axis=0)/rateTotal).tolist())),
        minimumDeltaTimeSeconds=float(min(dt[m]))) for name,m in mm.items() if any(m)}
    tri=np.array([f[1:] for f in wall]);p0,p1,p2=(points[tri[:,j]] for j in range(3));ctr=(p0+p1+p2)/3
    radius=max(np.linalg.norm(points[tri]-ctr[:,None],axis=2).ravel());tree=cKDTree(ctr);vtree=cKDTree(points[:nw])
    def sample(i):
        near=vtree.query(points[i]);cand=tree.query_ball_point(points[i],radius+near[0]+1e-12)
        distance=min(triangle_distances(points[i],p0[cand],p1[cand],p2[cand]))
        return dict(node=int(i),xyzM=points[i].tolist(),normalLayer=int(layer[i]),wallDistanceM=float(distance),
            nearestWallVertex=int(near[1]),dualVolumeM3=float(volume[i]),deltaTimeSeconds=float(dt[i]),
            nextCFL=float(cfl[i]),temperatureK=float(a[i,keys.index('Temperature')]),
            minimumIncidentTetMeanRatio=None if not np.isfinite(quality[i]) else float(quality[i]),
            RHS=dict(zip(names,residual[i].tolist())))
    extremes={'minimumTimeStep':[sample(i) for i in np.argsort(dt)[:8]],
        'v127EnthalpyAndTemperatureNodes':[sample(i) for i in [286416,167659]]}
    for j,name in enumerate(names):extremes['largestRHS_'+name]=[sample(i) for i in np.argsort(abs(residual[:,j]))[-5:][::-1]]
    length=float(np.ptp(points[:nw,0]));R=float(c['GAS_CONSTANT']);T=float(c['FREESTREAM_TEMPERATURE']);gamma=float(c['GAMMA_VALUE'])
    rhoi=float(c['FREESTREAM_PRESSURE'])/(R*T);speed=float(c['MACH_NUMBER'])*np.sqrt(gamma*R*T);tref=length/speed
    ki=receipt['physicsModel']['kInfinity'];wi=receipt['physicsModel']['omegaInfinity'];Ei=R*T/(gamma-1)+speed**2/2+ki
    scales=rhoi*np.array([1,speed,speed,speed,Ei,ki,wi])
    normalization=dict(bodyLengthM=length,advectiveReferenceSeconds=float(tref),transportedVariableScales=scales.tolist(),
        dimensionlessVolumeWeightedRateRMS=(np.sqrt(rateTotal/volume.sum())*tref/scales).tolist(),
        definition='sqrt(sum(R_i^2/V_i)/sum(V_i))*(L/Uinf)/Qinf; diagnostic only; constrained rows and stage-lag retained')
    hotspots={name:dict(top1PercentRHSSquareFraction=float(np.sort(residual[:,j]**2)[-int(np.ceil(n*.01)):].sum()/total[j]),
        top10RHSSquareFraction=float(np.sort(residual[:,j]**2)[-10:].sum()/total[j])) for j,name in enumerate(names)}
    fields=['Delta_Time','CFL','Turb_Delta_Time','Turb_CFL'];stepRanges={s:[float(min(a[:,keys.index(s)])),float(max(a[:,keys.index(s)]))] for s in fields}
    result=dict(nodes=n,originalFaces=len(wall),originalAreaM2=mr['bodyArea'],sourceLocalIteration=149,newLocalIteration=0,
        sourceRuntime=receipt['runtimeFlavor'],nativeSeconds=receipt['elapsedSeconds'],recomputedLogRMS=dict(zip(names,logs.tolist())),
        maximumLogRMSDifference=float(discrepancy),meshVolumeRelativeDifference=float(error),cellTypes=totals,
        maximumCoordinateRelativeRounding=coordinateError,coordinateRoundingBound=8*np.finfo(float).eps,
        timestepRanges=stepRanges,layerGroups=groups(masks),axialGroups=groups(axial),hotspots=hotspots,
        perLayer=groups({'layer'+str(k):layer==k for k in range(33)}),extrema=extremes,normalization=normalization,
        diagnosticOnly=True,physicsQualified=False,newQualifiedArea=0,originalResidualGatePassed=False,
        semantics='Implicit integrated RHS before each subsystem update; dt used this update, CFL adapted for next update; not a new-state defect or physical error.',
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),folder/'engine_process.json',folder/'restart_flow.csv',
            folder/'history.csv',folder/'mesh.su2',ROOT/'scripts/su2SSTInterface.py',ROOT/'outputs/v120_validation/mesh32/mesh_report.json']})
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(maximumLogRMSDifference=float(discrepancy),groups=result['layerGroups'],hotspots=hotspots,
        minDT=extremes['minimumTimeStep'][0],energyPeak=extremes['largestRHS_Energy'][0]),indent=2))


if __name__=='__main__':analyze()
