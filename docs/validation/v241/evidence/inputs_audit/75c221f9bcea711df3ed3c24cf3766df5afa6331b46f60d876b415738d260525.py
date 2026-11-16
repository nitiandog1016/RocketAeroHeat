"""模块职责：独立重心解、二分热恢复和邻线度量审核；只读已保存轨迹，不重积分。"""
from pathlib import Path
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseEulerRegionalVolume import read_tetra_mesh

BASE=ROOT/'outputs/v241_validation';SOURCE=BASE/'streamtubes01'


def independent_states(u):
    """独立实现：直接Cp分段积分及全域二分内能根，不调用生产反解函数。"""
    d=read(ROOT/'data/thermo/phase_aware_dry_air_v1.json');t=np.asarray(d['temperatureK']);cp=np.asarray(d['cpJkgK']);R=d['gasConstantJkgK']
    a=np.diff(cp)/np.diff(t);b=cp[:-1]-a*t[:-1]
    ih=np.r_[0,np.cumsum(.5*(cp[1:]+cp[:-1])*np.diff(t))]
    iss=np.r_[0,np.cumsum(a*np.diff(t)+b*np.log(t[1:]/t[:-1]))]
    ref=int(np.flatnonzero(t==298.15)[0]);ih-=ih[ref];iss-=iss[ref]
    def thermal(temp):
        j=np.searchsorted(t[1:-1],temp);dt=temp-t[j]
        return ih[j]+cp[j]*dt+.5*a[j]*dt**2,iss[j]+a[j]*dt+b[j]*np.log(temp/t[j])
    rho=u[:,0];v=u[:,1:4]/rho[:,None];ke=np.sum(v*v,axis=1)/2;e=u[:,4]/rho-ke-300000
    assert np.all(rho>0) and np.isfinite(u).all()
    lo=np.full(len(u),t[0]);hi=np.full(len(u),t[-1]);elo=thermal(lo)[0]-R*lo;ehi=thermal(hi)[0]-R*hi
    assert np.all(e>=elo)&np.all(e<=ehi)
    for _ in range(52):
        mid=(lo+hi)/2;em=thermal(mid)[0]-R*mid;low=em<e;lo[low]=mid[low];hi[~low]=mid[~low]
    temp=(lo+hi)/2;h,s=thermal(temp);p=rho*R*temp
    upper=np.minimum(1e6,1e4*(temp/100)**3);cold=temp<=132.6312
    if np.any(cold):
        theta=1-temp[cold]/132.6312
        dew=3785020*np.exp(132.6312/temp[cold]*(-.1567266*theta**.5-5.539635*theta+.7567212*theta**2.5-3.514322*theta**4))
        upper[cold]=np.minimum(upper[cold],.1*dew)
    assert np.all((p>0)&(p<=upper))
    return dict(rho=rho,v=v,U=np.linalg.norm(v,axis=1),p=p,T=temp,H=h+ke,S=s-R*np.log(p/100000))


def segment_boundary_hits(points,x,markers):
    """独立几何：检查已接受中心线段是否穿过实际边界三角形内部，不能只看端点。"""
    faces=np.concatenate(list(markers.values()));p=x[faces];e1=p[:,1]-p[:,0];e2=p[:,2]-p[:,0]
    low=p.min(axis=1);high=p.max(axis=1);hits=[]
    for k,(start,end) in enumerate(zip(points[:-1],points[1:])):
        near=np.all(high>=np.minimum(start,end),axis=1)&np.all(low<=np.maximum(start,end),axis=1)
        idx=np.flatnonzero(near)
        if not len(idx):continue
        direction=end-start;cross=np.cross(direction,e2[idx]);det=np.einsum('ij,ij->i',e1[idx],cross)
        safe=abs(det)>1e-24;idx=idx[safe];cross=cross[safe];det=det[safe]
        if not len(idx):continue
        z=start-p[idx,0];u=np.einsum('ij,ij->i',z,cross)/det;cc=np.cross(z,e1[idx]);v=cc@direction/det
        tau=np.einsum('ij,ij->i',e2[idx],cc)/det
        inside=(u>1e-10)&(v>1e-10)&(u+v<1-1e-10)&(tau>1e-10)&(tau<1-1e-10)
        if np.any(inside):hits.append(dict(segment=k,boundaryRows=idx[inside].tolist(),fraction=tau[inside].tolist()))
    return hits


def audit():
    manifest=read(SOURCE/'manifest.json');assert manifest['passed'] and manifest['status']=='unqualified-research'
    assert not manifest['physicalQualification'] and not manifest['edgeAccepted'] and not manifest['heatProviderCalled']
    for name,h in manifest['sourceSHA256'].items():assert sha(ROOT/name)==h,name
    x,t,markers=read_tetra_mesh(ROOT/'outputs/v240_validation/run01/mesh.su2')
    raw=np.loadtxt(str(ROOT/'outputs/v240_validation/run01/state_000200.csv'),delimiter=',',skiprows=1);raw=raw[np.argsort(raw[:,0])];u=raw[:,4:9]
    mapping=np.loadtxt(str(ROOT/'outputs/v240_validation/run01/original_faces.csv'),delimiter=',',skiprows=1)
    target=set(read(ROOT/'outputs/v228_validation/prepare03/region.json')['originalFaceIDs'])
    records=[];maxw=maxu=maxm=0.;counts=0;thermos=[]
    sets=manifest['tubes']+[manifest['halfStepCheck']]
    for info in sets:
        path=SOURCE/info['file'];assert sha(path)==info['sha256'];z=np.load(str(path));common=int(z['metric_commonStations'])
        raystates=[]
        for j in range(9):
            p=z['ray%d_xyz'%j];ids=z['ray%d_cell'%j];weights=z['ray%d_weights'%j];uu=z['ray%d_U'%j];s=z['ray%d_s'%j]
            assert len(p)==info['rayStationCounts'][j] and np.isfinite(p).all()
            assert np.all(ids>=0) and np.all(ids<len(t))
            q=x[t[ids]];v=np.linalg.solve((q[:,1:]-q[:,0,None]).transpose(0,2,1),p-q[:,0])
            w=np.column_stack([1-v.sum(axis=1),v]);dw=float(abs(w-weights).max());maxw=max(maxw,dw);assert dw<5e-10
            assert w.min()>=-2e-11 and np.max(abs(np.sum(w[:,:,None]*q,axis=1)-p))<2e-12
            expected=np.einsum('ni,nij->nj',w,u[t[ids]])
            du=float(np.max(abs(expected-uu)/np.maximum(1,abs(expected))));maxu=max(maxu,du);assert du<2e-10
            state=independent_states(uu);raystates.append(state);counts+=len(p)
            ds=info.get('dsM',manifest['dsM']);assert np.array_equal(s,np.arange(len(p))*ds)
            assert np.all(np.linalg.norm(np.diff(p,axis=0),axis=1)<=ds*(1+1e-10))
        if 'id' in info:
            assert info['originalSeedFaceID'] in target
            face=mapping[mapping[:,0]==info['originalSeedFaceID'],1:4].astype(int)[0]
            assert set(face).issubset(set(t[info['ownerCell']]))
            assert np.array_equal(z['ray0_xyz'][0],x[t[info['ownerCell']]].mean(axis=0))
            delta=info['labelDeltaM'];center=raystates[0];direction=center['v'][:common]/center['U'][:common,None]
            jb=(z['ray3_xyz'][:common]-z['ray2_xyz'][:common])/delta
            je=(z['ray7_xyz'][:common]-z['ray6_xyz'][:common])/delta
            area=np.array([np.linalg.det(np.vstack([a,b,c])) for a,b,c in zip(direction,jb,je)])
            transverse=jb-direction*np.einsum('ij,ij->i',direction,jb)[:,None]
            flux=center['rho'][:common]*center['U'][:common]*area
            for name,expected in [('hBeta',np.linalg.norm(transverse,axis=1)),('signedArea',area),
              ('massFluxRatio',flux/flux[0]),('HRelativeDrift',(center['H'][:common]-center['H'][0])/abs(center['H'][0])),
              ('entropyDrift',center['S'][:common]-center['S'][0])]:
                diff=float(np.max(abs(expected-z['metric_'+name])/np.maximum(1,abs(expected))));maxm=max(maxm,diff);assert diff<2e-10,(name,diff)
            crossings=segment_boundary_hits(z['ray0_xyz'],x,markers);assert not crossings,crossings
            # 仅标识保存末点靠近的边界；原失败RK阶段未保存，不冒称实际首次穿壁点。
            last=z['ray0_xyz'][-1];distances={name:float(np.linalg.norm(x[np.unique(f)]-last,axis=1).min()) for name,f in markers.items()}
            records.append(dict(tubeID=info['id'],commonStations=common,centerStations=len(center['p']),
              lastCenterXYZ=last.tolist(),endpointNearestMarkerNodeDistancesM=distances,actualFailedRKStageUnknown=True,
              acceptedCenterSegmentsNoInteriorBoundaryCrossing=True))
            thermos.append((info['id'],center,common))
    out=BASE/'audit01';out.mkdir(exist_ok=False)
    header='station,sM,xM,yM,zM,cellID,rhoKgm3,Ux,Uy,Uz,speedMps,pPa,TK,HJkg,SJkgK,hBeta,hEta,signedArea,massFluxRatio,HRelativeDrift,entropyDrift,metricValid'
    for tubeid,state,n in thermos:
        z=np.load(str(SOURCE/('tube%02d.npz'%tubeid)))
        a=np.column_stack([np.arange(n),z['metric_s'],z['ray0_xyz'][:n],z['ray0_cell'][:n],state['rho'][:n],state['v'][:n],
          state['U'][:n],state['p'][:n],state['T'][:n],state['H'][:n],state['S'][:n]]+
          [z['metric_'+k] for k in ('hBeta','hEta','signedArea','massFluxRatio','HRelativeDrift','entropyDrift','valid')])
        np.savetxt(str(out/('tube%02d.csv'%tubeid)),a,delimiter=',',header=header,comments='',fmt='%.17g')
    save(out/'report.json',dict(passed=True,savedRayStationsAudited=counts,tubeCount=4,sourceStable=True,
      maxIndependentBarycentricDifference=maxw,maxIndependentFiveStateScaledDifference=maxu,maxIndependentMetricScaledDifference=maxm,
      independentThermalInversion='52-step bisection of independently integrated piecewise Cp',
      tubes=records,totalCommonStations=sum(r['commonStations'] for r in records),
      acceptedCenterPolylinesChecked=True,allNineContinuousCurvesProvedInside=False,
      originalFailedStageLocationsRecorded=False,qualificationFalsePreserved=True,physicalQualification=False,newQualifiedAreaM2=0,
      CFDUpdates=0,flowReevaluations=0,trajectoryReintegrations=0))


if __name__=='__main__':audit()
