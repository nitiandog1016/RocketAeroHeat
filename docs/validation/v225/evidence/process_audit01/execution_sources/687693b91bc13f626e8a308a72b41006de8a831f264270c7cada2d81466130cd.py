"""模块职责：独立积分核对已保存波系；从实际端点恢复鼻部几何诊断，不重跑MATLAB。"""
from pathlib import Path
import json,math,time
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.io import loadmat
import psutil
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v225_validation';OUT=BASE/'audit01'

def main():
    start=time.monotonic();OUT.mkdir(exist_ok=False)
    data=read(ROOT/'data/thermo/phase_aware_dry_air_v1.json')
    x=np.array(data['temperatureK']);cp=np.array(data['cpJkgK']);R=data['gasConstantJkgK']
    slope=np.diff(cp)/np.diff(x)
    def caloric(t):
        assert 60<=t<=2000
        # 从参考温度直接分段积分，独立于MATLAB的累计表插值。
        lo=min(t,298.15);hi=max(t,298.15);sign=1 if t>=298.15 else -1
        points=np.r_[lo,x[(x>lo)&(x<hi)],hi];h=0.;s=0.
        for a,b in zip(points[:-1],points[1:]):
            j=min(len(slope)-1,np.searchsorted(x,.5*(a+b))-1);B=slope[j];A=cp[j]-B*x[j]
            h+=A*(b-a)+.5*B*(b*b-a*a);s+=A*math.log(b/a)+B*(b-a)
        C=float(np.interp(t,x,cp));return C,sign*h,sign*s,math.sqrt(C/(C-R)*R*t)
    def entropy_delta(t1,t2,p1,p2):return caloric(t2)[2]-caloric(t1)[2]-R*math.log(p2/p1)
    rows=read(BASE/'matlab02/waves.json');checks=[]
    for i,row in enumerate(rows):
        p1=row['upstreamPressurePa'];t1=row['upstreamTemperatureK'];u1=row['upstreamSpeedMps'];w=row['result']
        p2=w['pressurePa'];t2=w['temperatureK'];u2=w['speedMps'];C,h1,s1,a1=caloric(t1);_,h2,s2,a2=caloric(t2)
        energy=(h2+.5*u2*u2-h1-.5*u1*u1)/max(1,abs(h1)+.5*u1*u1)
        ent=entropy_delta(t1,t2,p1,p2);mass=None;momentum=None;angle=None;weak=None
        assert abs(energy)<1e-8
        if w['kind']=='expansion' and w['turnDeg']>0:
            h0=h1+.5*u1*u1
            def integrand(t):
                C,h,_,a=caloric(t);v2=2*(h0-h);return math.sqrt(v2/a**2-1)*C/v2
            # 按源Cp分段，端点显式保留；独立QUADPACK积分，不复用MATLAB节点或角度。
            points=np.r_[t2,x[(x>t2)&(x<t1)],t1]
            theta=sum(quad(integrand,a,b,epsabs=1e-12,epsrel=1e-11)[0] for a,b in zip(points[:-1],points[1:]))
            angle=theta-math.radians(w['turnDeg']);assert abs(angle)<1e-7 and abs(ent)<1e-7
        elif w['kind']=='normal' or w['turnDeg']>0:
            beta=math.radians(w['shockAngleDeg']);n1=u1*math.sin(beta);ut=u1*math.cos(beta)
            # 用独立质量方程恢复Un2，而非读取MATLAB的法向速度输出。
            n2=p1/(R*t1)*n1/(p2/(R*t2))
            mass=(p2/(R*t2)*n2-p1/(R*t1)*n1)/(p1/(R*t1)*n1)
            momentum=(p2+p2/(R*t2)*n2**2-p1-p1/(R*t1)*n1**2)/(p1+p1/(R*t1)*n1**2)
            assert abs(momentum)<1e-8 and abs(math.hypot(n2,ut)/u2-1)<1e-8 and ent>=-1e-8
            if w['kind']=='oblique':
                angle=beta-math.atan2(n2,ut)-math.radians(w['turnDeg']);assert abs(angle)<1e-7
            if 0<w['turnDeg']<1e-3:
                gamma=C/(C-R);M=u1/a1;linear=gamma*M*M/math.sqrt(M*M-1)*math.radians(w['turnDeg'])
                weak=abs((p2/p1-1)/linear-1);assert weak<2e-6
            else:
                # 独立压力比变量的Hugoniot求根，避免复制生产z坐标与二分实现。
                rho1=p1/(R*t1)
                def jump(pr):
                    v=n1-(p1*pr-p1)/(rho1*n1);temp=p1*pr*v/(rho1*n1*R)
                    if not 60<=temp<=2000 or v<=0:return None
                    return caloric(temp)[1]+.5*v*v-h1-.5*n1*n1
                upper=1+rho1*n1*n1/p1
                pts=np.linspace(1.00001,upper*(1-1e-8),80);bracket=None
                for aa,bb in zip(pts[:-1],pts[1:]):
                    fa=jump(aa);fb=jump(bb)
                    if fa is not None and fb is not None and fa>0 and fb<0:bracket=(aa,bb);break
                assert bracket is not None
                ratio=brentq(jump,*bracket,xtol=1e-12);assert abs(ratio/(p2/p1)-1)<1e-8
        else:assert p1==p2 and t1==t2 and u1==u2
        checks.append(dict(index=i,energyResidual=energy,entropyChangeJkgK=ent,massResidual=mass,
          momentumResidual=momentum,angleResidualRad=angle,weakLinearRelativeDifference=weak))
    path=read(BASE/'matlab02/path.json');blocked=read(BASE/'matlab02/blocked_path.json')
    assert all(path['valid']) and blocked['valid']==[True,False,False]
    assert blocked['pressurePa'][1:]==[None,None]
    for j in range(1,4):
        row=rows[7+j]
        assert row['upstreamPressurePa']==path['pressurePa'][j-1]
        assert row['upstreamTemperatureK']==path['temperatureK'][j-1]
        assert row['upstreamSpeedMps']==path['speedMps'][j-1]

    # sc.sectionX只列截面环，实际单鼻尖节点在环外，不能断言末环半径为零。
    mat=loadmat(str(ROOT/'outputs/v104_validation/baseline02/current_case.mat'),squeeze_me=True,struct_as_record=False)
    m=mat['current'];fs=mat['fs'];sc=m.surfaceCoord;V=m.V;F=m.F.astype(int)-1
    area=.5*np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)
    assert len(F)==19752 and abs(area.sum()-1.9646292565715553)<2e-13
    ids=np.flatnonzero(V[:,0]==V[:,0].max());assert len(ids)==1
    tip=V[ids[0]];last=np.asarray(sc.sectionNodeIDs[-1],dtype=int)-1;ring=V[last]
    center=np.asarray(sc.sectionCenterYZ[-1]);radius=np.linalg.norm(ring[:,1:]-center,axis=1).mean()
    dx=ring[:,0].mean()-tip[0];assert dx<0 and radius>0
    geometric=math.degrees(math.atan2(radius,-dx));wind=geometric+fs.alpha;lee=geometric-fs.alpha
    maximum=rows[3]['result']['maximumTurnDeg'];assert min(wind,lee)>maximum
    incident=np.any(F==ids[0],axis=1)
    rocket=dict(faces=len(F),areaM2=float(area.sum()),tipVertexID=int(ids[0]+1),tip=tip.tolist(),
      lastRingX=float(ring[:,0].mean()),lastRingRadiusM=float(radius),lastRingIsNotTip=True,
      tipAdjacentFaces=int(incident.sum()),tipAdjacentAreaM2=float(area[incident].sum()),
      meridianInclinationDeg=geometric,windwardFirstTurnDeg=wind,leewardFirstTurnDeg=lee,
      maximumWeakShockTurnDeg=maximum,geometryOnlyDiagnostic=True,oldEdgeValuesUsed=False,
      newCandidateHeatAreaM2=0,newQualifiedAreaM2=0,originalFaceStatesProduced=False,
      needs='detached blunt-nose/subsonic-to-sonic-line edge and boundary-layer inlet history; no downstream reset')
    save(OUT/'rocket_inlet.json',rocket)
    admission=read(BASE/'process_matlab02/admission.json');live=read(BASE/'process_matlab02/live.json')
    assert not psutil.pid_exists(live['pid'])
    assert not [p for p in psutil.process_iter(['name']) if (p.info['name'] or '').lower()=='matlab.exe']
    stable={n:sha(n)==h for n,h in admission['inputSHA256'].items()};assert all(stable.values())
    processObservation=dict(status='terminal-process-absent',observedPid=live['pid'],inputStable=stable,
      originalSupervisorReportPresent=(BASE/'process_matlab02/report.json').exists(),
      originalExitCode=None,originalPeakRSS=None,actualElapsedSeconds=None,
      MATLABFullTestPassed=False,knownFailure='test geometry assumption: final ring is the tip; assertion line 51',
      completedWaveArtifactsSHA256={n:sha(BASE/'matlab02'/n) for n in ('waves.json','path.json','blocked_path.json')},
      budgetConservativelyChargedSeconds=admission['limitSeconds'],automaticRestart=False,
      note='Saved output survives; missing supervisor completion is not reconstructed as a successful run.')
    save(OUT/'matlab02_terminal_observation.json',processObservation)
    save(OUT/'wave_checks.json',checks)
    result=dict(passed=True,scope='independent audit of saved wave/path outputs and corrected original tip geometry only',
      states=len(rows),maxEnergyResidual=max(abs(c['energyResidual']) for c in checks),
      maxAngleResidualRad=max(abs(c['angleResidualRad'] or 0) for c in checks),
      weakLinearRelativeDifference=checks[2]['weakLinearRelativeDifference'],
      exactUpstreamInheritance=True,downstreamBlockPreserved=True,rocket=rocket,
      physicalQualification=False,fullMatlabSuitePassed=False,elapsedSeconds=time.monotonic()-start)
    save(OUT/'report.json',result);print(json.dumps(result))

if __name__=='__main__':main()
