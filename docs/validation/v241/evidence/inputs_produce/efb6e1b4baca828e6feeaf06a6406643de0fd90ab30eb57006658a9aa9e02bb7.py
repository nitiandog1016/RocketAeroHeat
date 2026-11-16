"""模块职责：直接读取内部前体完整场，生成四个原面种子的真实三维九线流管研究数据。"""
from pathlib import Path
import math
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseEulerRegionalVolume import read_tetra_mesh,topology
from tetraStreamtubes import TetraField,tube,LABELS

BASE=ROOT/'outputs/v241_validation';OLD=ROOT/'outputs/v240_validation/run01'


def produce():
    """模块一：绑定内部源身份、原面职责和原始坐标，不借用外部CFD。"""
    assert read(BASE/'test02/report.json')['passed']
    x,c,markers=read_tetra_mesh(OLD/'mesh.su2')
    a=np.loadtxt(str(OLD/'state_000200.csv'),delimiter=',',skiprows=1);a=a[np.argsort(a[:,0])]
    assert a.shape==(len(x),15) and np.array_equal(a[:,0],np.arange(len(x))) and np.array_equal(a[:,1:4],x)
    old=read(ROOT/'outputs/v240_validation/audit02/report.json')
    assert sha(OLD/'state_000200.csv')==next(r['sourceSHA256'] for r in old['savedFrameChecks'] if r['update']==200)
    assert not old['physicalQualification']
    field=TetraField(x,c,a[:,4:9]);topo=topology(x,c,markers)
    mapping=np.loadtxt(str(OLD/'original_faces.csv'),delimiter=',',skiprows=1)
    target=set(read(ROOT/'outputs/v228_validation/prepare03/region.json')['originalFaceIDs'])
    rows=np.flatnonzero(np.isin(mapping[:,0],list(target)));tri=mapping[:,1:4].astype(int)
    first_x=float(x[tri[rows],0].max());ring=rows[np.any(abs(x[tri[rows],0]-first_x)<1e-12,axis=1)]
    centers=x[tri[ring]].mean(axis=1);angle=np.arctan2(centers[:,2],centers[:,1])
    seeds=[int(ring[np.argmin(abs(np.angle(np.exp(1j*(angle-beta)))))] ) for beta in (0,np.pi/2,np.pi,3*np.pi/2)]
    assert len(set(seeds))==4
    owners={tuple(sorted(topo['faces'][i])):int(topo['owner'][i]) for i in topo['markers']['wall']}
    length=float(np.ptp(x[np.unique(markers['wall']),0]));ds=min(length/256,.002);steps=int(math.ceil(1.5*length/ds))
    out=BASE/'streamtubes01';out.mkdir(exist_ok=False)
    manifest=dict(schema='RAH.InternalVolumeStreamtubes.v1',status='unqualified-research',physicalQualification=False,
      edgeAccepted=False,heatProviderCalled=False,newHeatAreaM2=0,newQualifiedAreaM2=0,
      sourceBranch='phase-air-planar-forebody-euler-v240',sourceUpdate=200,originalFaces=19752,
      originalAreaM2=1.9646292565715553,engineeringTargetFaces=3128,sourceRegionFaces=4434,
      sourceIdentity='Euler wall-near volume samples; no boundary-layer-edge qualification',
      sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in [OLD/'mesh.su2',OLD/'state_000200.csv',OLD/'original_faces.csv',
        ROOT/'outputs/v228_validation/prepare03/region.json',ROOT/'data/thermo/phase_aware_dry_air_v1.json']},
      energyShiftJkg=300000,labelUnit='m',dsM=ds,requestedSteps=steps,firstTargetRingXM=first_x,tubes=[])
    for k,row in enumerate(seeds):
        owner=owners[tuple(sorted(tri[row]))];seed=x[c[owner]].mean(axis=0)
        _,_,u=field.sample(seed,owner);state=field.states(u);t=state['velocity']/state['speed']
        circum=np.array([0.,-seed[2],seed[1]]);beta=circum-t*np.dot(t,circum);beta/=np.linalg.norm(beta)
        eta=np.cross(t,beta)
        gradients=np.vstack([-field.inverse[owner].sum(axis=0),field.inverse[owner]])
        distance=float(np.min(.25/np.linalg.norm(gradients,axis=1)));delta=.2*distance
        for la,lb in LABELS:
            _,w=field.locate(seed+delta*(la*beta+lb*eta),owner);assert w.min()>0
        rays,metric=tube(field,seed,beta,eta,delta,ds,steps)
        arrays={}
        for j,r in enumerate(rays):
            for name in ('xyz','cell','weights','U','s'):arrays['ray%d_%s'%(j,name)]=r[name]
        for name,value in metric.items():arrays['metric_'+name]=np.asarray(value)
        path=out/('tube%02d.npz'%(k+1));np.savez_compressed(str(path),**arrays)
        record=dict(id=k+1,originalSeedFaceID=int(mapping[row,0]),ownerCell=owner,seedXYZ=seed.tolist(),
          seedLabelBeta=beta.tolist(),seedLabelEta=eta.tolist(),labelDeltaM=delta,ownerMinimumFaceDistanceM=distance,
          requestedStations=steps+1,rayStationCounts=[len(r['s']) for r in rays],rayStopReasons=[r['stopReason'] for r in rays],
          commonStations=metric['commonStations'],validMetricStations=int(sum(metric['valid'])),
          maxAbsHRelativeDrift=float(max(abs(metric['HRelativeDrift']))),
          maxAbsEntropyDriftJkgK=float(max(abs(metric['entropyDrift']))),
          massFluxRatioRange=[float(min(metric['massFluxRatio'])),float(max(metric['massFluxRatio']))],
          maxLabelIndicator=float(max(metric['labelIndicator'])),file=path.name,sha256=sha(path))
        manifest['tubes'].append(record);save(out/'manifest.json',manifest)
        print('completed tube',k+1,'stations',metric['commonStations'],flush=True)
        if k==0:
            fine,fm=tube(field,seed,beta,eta,delta,ds/2,steps*2)
            compare={};errors=[]
            for j,(coarse,refined) in enumerate(zip(rays,fine)):
                n=min(len(coarse['s']),(len(refined['s'])+1)//2)
                errors.append(float(np.max(np.linalg.norm(coarse['xyz'][:n]-refined['xyz'][::2][:n],axis=1))))
                for name in ('xyz','cell','weights','U','s'):compare['ray%d_%s'%(j,name)]=refined[name]
            for name,value in fm.items():compare['metric_'+name]=np.asarray(value)
            fp=out/'tube01_halfstep.npz';np.savez_compressed(str(fp),**compare)
            nc=min(metric['commonStations'],(fm['commonStations']+1)//2)
            manifest['halfStepCheck']=dict(tubeID=1,commonStations=nc,maxPositionDifferencesM=errors,
              maxHBetaDifference=float(np.max(abs(metric['hBeta'][:nc]-fm['hBeta'][::2][:nc]))),
              maxAreaDifference=float(np.max(abs(metric['signedArea'][:nc]-fm['signedArea'][::2][:nc]))),
              rayStopReasons=[r['stopReason'] for r in fine],rayStationCounts=[len(r['s']) for r in fine],
              dsM=ds/2,file=fp.name,sha256=sha(fp),spatialMeshOrCFDConvergence=False)
            save(out/'manifest.json',manifest)
    manifest['passed']=True;save(out/'manifest.json',manifest)


if __name__=='__main__':produce()
