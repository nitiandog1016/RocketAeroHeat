"""模块职责：原面共同环、已存实际剖面、横向标签与激波质量匹配合成为入口研究数据。"""
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq
from phaseAwareEntropyInlet import ShockMassModel,spherical_source_table,match_inlet_mass,integrate_metric_shock
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v253_validation';SOURCE=ROOT/'outputs/v252_validation/matlab01'


def common_ring(g,region):
    """模块一：共同边必须真实连接原鼻部与3128前体；不凭最近点选造入口。"""
    F=g['originalFaces'];nose=set(np.array(region['retainedNoseCandidateIDs'])-1);target=set(np.array(region['originalFaceIDs'])-1)
    def edges(ids):
        result={}
        for i in ids:
            for a,b in zip(F[i],np.roll(F[i],-1)):result.setdefault(tuple(sorted((int(a),int(b)))),[]).append(i)
        return result
    a=edges(nose);b=edges(target);shared=sorted(a.keys()&b.keys())
    assert shared and all(len(a[k])==len(b[k])==1 for k in shared)
    ids=np.unique(shared);V=g['originalVertices'];phi=np.arctan2(V[ids,2],V[ids,1]);ids=ids[np.argsort(phi)]
    pairs=np.column_stack([ids,np.roll(ids,-1)])
    assert {tuple(sorted(x)) for x in pairs}==set(shared)
    assert len(ids)==len(shared) and len(set(V[ids,0]))==1
    return pairs,np.array([a[tuple(sorted(x))][0] for x in pairs]),np.array([b[tuple(sorted(x))][0] for x in pairs])


def controls(model,table):
    """模块二：非恒定尺度解析质量、标签不变性和具体非法入口拒绝。"""
    s=np.linspace(0,.05,21);h=.01*(1+8*s);u=.9;j=model.shock(u)
    p=model.p+(j['pressurePa']-model.p)*(model.p0-model.p)/(model.normal['pressurePa']-model.p)
    initial={'sinShockAngle':u,'massCoordinate':.001,'sourceSHA256':'a'*64}
    actual=integrate_metric_shock(model,s,h,np.full(len(s),p),initial)
    expected=initial['massCoordinate']+model.rho*model.speed*u*.01*(s+4*s*s)
    error=float(max(abs(actual[1]-expected))/max(expected));angle_error=float(max(abs(actual[0]-u)))
    assert error<2e-5 and angle_error<2e-5
    altered=integrate_metric_shock(model,s,h*3,np.full(len(s),p),dict(initial,massCoordinate=.003))
    assert max(abs(altered[0]-actual[0]))<1e-10 and max(abs(altered[1]/3-actual[1]))<1e-10
    rejected=[]
    for name,fn in [
        ('missing-finite-source',lambda:integrate_metric_shock(model,s,h,np.full(len(s),p),{})),
        ('zero-finite-source',lambda:integrate_metric_shock(model,s,h,np.full(len(s),p),dict(initial,massCoordinate=0))),
        ('folded-label',lambda:match_inlet_mass(model,table,p,.1,.01,-1)),
        ('out-of-source-mass',lambda:match_inlet_mass(model,table,p,1e9,.01,1)),
        ('nonmonotone-stations',lambda:integrate_metric_shock(model,s[::-1],h,np.full(len(s),p),initial))]:
        try:fn()
        except ValueError:rejected.append(name)
    assert len(rejected)==5
    return dict(passed=True,analyticMassRelativeError=error,constantAngleAbsoluteError=angle_error,
        labelReparameterizationPassed=True,rejections=rejected,scope='new general metric and inlet identity only')


def produce():
    out=BASE/'inlet01';out.mkdir(exist_ok=False)
    e=read(ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json');c=e['condition'];center=np.array(e['sphereCenter'])
    model=ShockMassModel(c['pressurePa'],c['temperatureK'],c['speedMps']);R=e['geometryEstimate']['radius']
    assert abs(model.p0/e['stagnation']['pressurePa']-1)<1e-8
    table=spherical_source_table(model,R,np.deg2rad(c['maximumPolarAngleDeg']));save(out/'shock_mass_source.json',table)
    save(out/'controls.json',controls(model,table))
    g=np.load(str(ROOT/'outputs/v222_validation/geometry01/cache.npz'));region=read(ROOT/'outputs/v228_validation/prepare03/region.json')
    pairs,noseFaces,targetFaces=common_ring(g,region);V=g['originalVertices'];wall=(V[pairs[:,0]]+V[pairs[:,1]])/2
    phi=np.arctan2(wall[:,2],wall[:,1]);dphi=np.angle(np.exp(1j*(np.arctan2(V[pairs[:,1],2],V[pairs[:,1],1])-np.arctan2(V[pairs[:,0],2],V[pairs[:,0],1]))))
    assert np.all(dphi>0)
    J=(V[pairs[:,1]]-V[pairs[:,0]])/dphi[:,None];radial=wall-center;radius=np.linalg.norm(radial,axis=1);n=radial/radius[:,None]
    flow=np.array(c['flowDirection']);t=flow-n*(n@flow)[:,None];t/=np.linalg.norm(t,axis=1)[:,None];b=np.cross(n,t)
    metric=np.sum(J*b,axis=1);along=np.sum(J*t,axis=1);normalComponent=np.sum(J*n,axis=1)
    ea=np.array([0.,1.,0.]);eb=np.cross(-flow,ea);a=radial@ea;z=radial@eb
    label=(a*(J@eb)-z*(J@ea))/(a*a+z*z);theta=np.arccos(n@(-flow))
    assert np.all(metric>0) and np.all(label>0) and min(theta)>0 and max(theta)<np.pi/3
    pressure=model.p+(model.p0-model.p)*np.cos(theta)**2
    sourceEntropy=e['stagnation']['entropyJkgK']
    sourceTe=np.array([brentq(lambda T:float(model.air.evaluate(p,T)['entropyJkgK'])-sourceEntropy,200.,model.t0,xtol=1e-10) for p in pressure])
    sourceU=np.sqrt(2*(model.H-model.air.thermal(sourceTe)[1]))
    arrays=dict(originalVertices=V,originalFaces=g['originalFaces'],originalFaceArea=g['faceArea'],wallXYZ=wall,
        ringEdges=pairs,upstreamFaceIDs=noseFaces+1,downstreamFaceIDs=targetFaces+1,
        phiRad=phi,thetaRad=theta,metricM=metric,labelJacobian=label,alongLabelShiftM=along,
        normalLabelComponentM=normalComponent,unitNormal=n,unitTangent=t,ringLabelDerivative=J,
        sourcePressurePa=pressure,sourceEdgeTemperatureK=sourceTe,sourceEdgeSpeedMps=sourceU)
    rows=[];summaries={};steps=np.linspace(0,1,601)
    for branch in ('base','fine','history'):
        chain=read(SOURCE/(branch+'.json'));angles=np.array([s['thetaRad'] for s in chain['stations']]);profiles=chain['profiles']
        yy=[];temps=[];speeds=[];densities=[];states=[];records=[];massFlux=[];enthalpyFlux=[]
        for i,angle in enumerate(theta):
            j=int(np.searchsorted(angles,angle)-1);w=float((angle-angles[j])/(angles[j+1]-angles[j]));assert 0<=w<=1
            def blend(key):return (1-w)*np.array(profiles[j]['profile'][key])+w*np.array(profiles[j+1]['profile'][key])
            y=blend('physicalY');T=blend('temperature');u=blend('speed')
            assert len(y)==601 and np.all(np.diff(y)>0) and y[0]==0
            state=model.air.evaluate(pressure[i],T);rho=state['densityKgm3'];velocity=u[:,None]*t[i]
            mass=float(simpson(rho*u,x=y));energy=float(simpson(rho*u*(state['enthalpyJkg']+.5*u*u),x=y))
            matched=match_inlet_mass(model,table,float(pressure[i]),mass,float(metric[i]),float(label[i]))
            invariant=match_inlet_mass(model,table,float(pressure[i]),mass,float(metric[i]/2),float(label[i]/2))
            assert matched['temperatureK']==invariant['temperatureK']
            deltaT=matched['temperatureK']/sourceTe[i]-1;deltaU=matched['speedMps']/sourceU[i]-1
            flags={'sourceGeometryQualified':False,'sourcePhysicalQualification':False,'coupledEntropyProfileConverged':False,
                'crossflowMomentumSolved':False,'SAInletAvailable':False,'profileReadyForCFD':False}
            rec=dict(branch=branch,index=i,upstreamFaceID=int(noseFaces[i]+1),downstreamFaceID=int(targetFaces[i]+1),
                ringNodeIDs=(pairs[i]+1).tolist(),phiDeg=float(np.rad2deg(phi[i])),sourcePolarAngleDeg=float(np.rad2deg(angle)),
                lowerProfile=j,upperProfile=j+1,upperWeight=w,sourceChainSHA256=chain['sourceSHA256'],
                pressurePa=float(pressure[i]),sourceEdgeTemperatureK=float(sourceTe[i]),sourceEdgeSpeedMps=float(sourceU[i]),
                massFluxPerWidthKgMs=mass,enthalpyFluxPerWidthWm=energy,
                sourceFiniteProfileTopTemperatureK=float(T[-1]),sourceFiniteProfileTopSpeedMps=float(u[-1]),
                matchedEdge=matched,edgeTemperatureRelativeCorrection=float(deltaT),edgeSpeedRelativeCorrection=float(deltaU),
                correctionBelowOnePercent=bool(max(abs(deltaT),abs(deltaU))<=.01),wallTemperatureK=float(T[0]),
                pressureSource='existing rotated spherical-cap modified Newtonian pressure; unchanged upstream engineering source',
                metricSource='original shared-edge chord derivative projected onto n cross t; not ring radius',
                profileSource='v252 actual finite profiles; linear polar/normalized-normal-coordinate interpolation',
                correctionMethod='CR112180 65-69; actual profile mass integral and wind-azimuth label conversion',**flags)
            records.append(rec);rows.append(rec);yy.append(y);temps.append(T);speeds.append(u);densities.append(rho);massFlux.append(mass);enthalpyFlux.append(energy)
            conservative=np.column_stack([rho,rho[:,None]*velocity,rho*(state['internalEnergyJkg']+.5*u*u)])
            states.append(conservative)
        for key,value in [('yM',yy),('temperatureK',temps),('speedMps',speeds),('densityKgm3',densities),('physicalConservative',states)]:arrays[branch+'_'+key]=np.array(value)
        summaries[branch]=dict(inlets=len(records),normalSamples=601,
            minMaxEntropyCorrection=[min(r['matchedEdge']['entropyJkgK']-sourceEntropy for r in records),max(r['matchedEdge']['entropyJkgK']-sourceEntropy for r in records)],
            maximumEdgeTemperatureCorrection=max(abs(r['edgeTemperatureRelativeCorrection']) for r in records),
            maximumEdgeSpeedCorrection=max(abs(r['edgeSpeedRelativeCorrection']) for r in records),
            correctionsBelowOnePercent=sum(r['correctionBelowOnePercent'] for r in records),
            sourceFiniteProfileMassFluxPerWidthRange=[min(massFlux),max(massFlux)],profileReadyForCFD=False)
    np.savez_compressed(str(out/'inlet_profiles.npz'),**arrays)
    save(out/'inlets.json',rows)
    save(out/'report.json',dict(passed=True,originalFaces=19752,originalAreaM2=float(sum(g['faceArea'])),targetFaces=len(region['originalFaceIDs']),
        commonRingEdges=len(pairs),ringXM=float(wall[0,0]),sourceRadiusM=R,
        maximumRadialGeometryErrorM=float(max(abs(radius-R))),maximumNormalLabelComponentM=float(max(abs(normalComponent))),
        projectedMetricRangeM=[float(min(metric)),float(max(metric))],labelJacobianRange=[float(min(label)),float(max(label))],
        maximumMetricVsRadiusRelativeDifference=float(max(abs(metric/np.linalg.norm(wall[:,1:],axis=1)-1))),
        sourceEntropyJkgK=sourceEntropy,sourceTotalEnthalpyJkg=model.H,branches=summaries,
        sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in [SOURCE/(k+'.json') for k in ('base','fine','history')]},
        newQualifiedAreaM2=0,newHeatAreaM2=0,acceptedEdge=False,profileReadyForCFD=False,
        exactExistingViscousCutBoundaryProduced=False,threeDimensionalRANSInletProduced=False,
        outputRole='actual shared-ring profile and mass-matched entropy candidate; joint feedback not closed')))


if __name__=='__main__':produce()
