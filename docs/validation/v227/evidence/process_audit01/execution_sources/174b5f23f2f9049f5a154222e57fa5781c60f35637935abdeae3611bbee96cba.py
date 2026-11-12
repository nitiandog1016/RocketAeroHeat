"""模块职责：独立物理温度BVP、物理壁梯度/积分、原面映射与工作台审核。"""
from pathlib import Path
import sys,json,time,hashlib,subprocess
import numpy as np
from scipy.integrate import solve_bvp,quad
from scipy.optimize import brentq
from phaseAwareViscousAir import PhaseAwareViscousAir
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v227_validation'

def independent_profile(air,r):
    e=r['edge'];st=r['station'];Te=e['temperatureK'];p=e['pressurePa'];U=e['speedMps'];Tw=r['wallTemperature']
    cp60,h60,_=air.thermal(60.);offset=float(h60-60*cp60);H=e['totalEnthalpyJkg']-offset
    mu_e=float(np.interp(Te,air.T,air.mu));rho=p/(air.R*Te);beta=st['beta']
    K=np.sqrt(mu_e/rho/(2*st['velocityGradientPerSecond'])) if U==0 else np.sqrt(2*st['xi'])/(rho*U*st['metricM'])
    def props(y):
        T=y[3];cp,h,_=air.thermal(T);mu=np.interp(T,air.T,air.mu);k=np.interp(T,air.T,air.kappa)
        l=mu/mu_e*Te/T;Pr=cp*mu/k;g=(h+.5*U*U*y[1]**2-offset)/H
        return T,cp,mu,k,l,Pr,g
    def ode(x,y):
        T,cp,mu,k,l,Pr,g=props(y);F=y[1];m=y[2];n=y[4]
        fp=m/l;gp=Pr/l*(n-(1-1/Pr)*U*U/H*F*m)
        tp=Pr/(l*cp)*(H*n-U*U*F*m)
        return np.array([F,fp,-y[0]*fp+beta*(F*F-T/Te),tp,-y[0]*gp])
    def bc(a,b):return np.array([a[0],a[1],(a[3]-Tw)/Te,b[1]-1,(b[3]-Te)/Te])
    x=np.linspace(0,12,81);F=np.tanh(x);T=Tw+(Te-Tw)*F;fp=1-F*F
    y=np.array([np.log(np.cosh(x)),F,np.zeros(len(x)),T,np.zeros(len(x))])
    _,cp,mu,k,l,Pr,_=props(y);y[2]=l*fp;y[4]=l/Pr*cp/H*(Te-Tw)*fp+U*U/H*F*y[2]
    sol=solve_bvp(ode,bc,x,y,tol=1e-7,max_nodes=12000)
    assert sol.success,sol.message
    audit=np.linspace(0,12,1201);Y=sol.sol(audit);T,cp,mu,k,l,Pr,g=props(Y)
    air.evaluate(p,T);F=Y[1];fp=Y[2]/l;gp=Pr/l*(Y[4]-(1-1/Pr)*U*U/H*F*Y[2])
    fs=max(abs(fp));gs=max(.005*max(abs(gp)),1e-8)
    def defect(x):
        y=sol.sol(x);T,cp,mu,k,l,Pr,g=props(y)
        gp=Pr/l*(y[4]-(1-1/Pr)*U*U/H*y[1]*y[2])
        return np.max(np.array([abs(1-y[1])/.005,abs(1-g)/.005,abs(T/Te-1)/.005,
            abs(Te/T-1)/.005,abs(y[2]/l)/fs/.005,abs(gp)/gs]),axis=0)
    bad=np.flatnonzero(defect(audit)>1);assert len(bad) and audit[bad[-1]+1]<=11
    ee=brentq(lambda x:float(defect(x)-1),audit[bad[-1]],audit[bad[-1]+1],xtol=1e-12)
    delta=K*quad(lambda x:sol.sol(x)[3]/Te,0,ee,epsabs=1e-10,epsrel=1e-8)[0]
    mass=K*quad(lambda x:sol.sol(x)[1],0,ee,epsabs=1e-10,epsrel=1e-8)[0]
    theta=K*quad(lambda x:sol.sol(x)[1]*(1-sol.sol(x)[1]),0,ee,epsabs=1e-10,epsrel=1e-8)[0]
    wall=sol.sol(0);_,_,mu,k,_,_,_=props(wall);der=ode(0,wall)
    q=float(k*der[3]/(K*Tw/Te));tau=float(mu*U*der[1]/(K*Tw/Te))
    calc=dict(heatFlux=q,wallShearStress=tau,delta=delta,deltaStar=delta-mass,thetaPhysical=theta,massThickness=mass)
    errors={name:abs(value-r[name])/max(abs(value),1e-12) for name,value in calc.items()}
    for key,err in errors.items():assert err<(2e-3 if key in ('delta','deltaStar','massThickness') else 2e-4),(key,err)
    xp=np.array(r['profile']['eta']);actual=np.array(r['profile']['temperature']);py=sol.sol(xp)
    temp_error=float(max(abs(py[3]-actual))/Te);vel_error=float(max(abs(py[1]-np.array(r['profile']['states'])[1])))
    assert max(temp_error,vel_error)<2e-4
    return dict(thetaRad=st['thetaRad'],values=calc,relativeErrors=errors,temperatureScaledError=temp_error,
      velocityRatioError=vel_error,etaEdge=ee,nodes=len(sol.x),physicalTemperatureCoordinates=True)

def main(stage):
    started=time.monotonic();out=BASE/stage;out.mkdir(exist_ok=False)
    source=BASE/'matlab02';assert read(BASE/'process_matlab02/report.json')['passed']
    h=read(source/'heating.json');m=read(source/'rocket_workbench.json');e=read(ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json')
    if stage=='render01':plot(out,m,h);save(out/'report.json',dict(passed=True,numericalRerun=False));return
    air=PhaseAwareViscousAir();profiles=h['profiles'];angles=np.array(h['anglesRad'])
    indices=[0,int(np.argmin(abs(angles-np.deg2rad(30)))),len(angles)-1]
    refs=[independent_profile(air,profiles[j]) for j in indices]
    # xi独立自适应积分：不复用MATLAB的固定32点积分或存储密度。
    st0=e['stagnation'];R=e['geometryEstimate']['radius'];pinf=e['condition']['pressurePa'];p0=st0['pressurePa'];T0=st0['temperatureK']
    _,_,s0=air.thermal(T0);h0=st0['enthalpyJkg']
    def integrand(th):
        p=pinf+(p0-pinf)*np.cos(th)**2;target=s0+air.R*np.log(p/p0)
        T=brentq(lambda t:float(air.thermal(t)[2]-target),60,T0,xtol=1e-11)
        cp,hh,_=air.thermal(T);U=np.sqrt(2*(h0-hh));mu=np.interp(T,air.T,air.mu)
        return float(p/(air.R*T)*mu*U*(R*np.sin(th))**2*R)
    xi_errors=[]
    for j in indices[1:]:
        xi=quad(integrand,0,angles[j],epsabs=1e-16,epsrel=1e-8,limit=100)[0]
        err=abs(xi/h['stations'][j]['xi']-1);assert err<2e-6;xi_errors.append(err)
    valid=np.array(h['validFace'],bool);area=np.array(m['faceArea']);assert valid.sum()==1306 and len(valid)==19752
    assert abs(area.sum()-1.9646292565715553)<2e-13
    keys=['heatFlux','wallShearStress','delta','deltaStar','thetaPhysical','massThickness','skinFriction']
    for key in keys:
        actual=np.array(h[key],float);assert np.all(np.isnan(actual[~valid]))
        f=next(f for f in m['fields'] if f['key']==key);assert np.array_equal(np.array(f['values'],float)[valid],actual[valid])
        if key=='skinFriction':continue
        lower=np.array(h['lowerStation'],int)[valid]-1;upper=np.array(h['upperStation'],int)[valid]-1;w=np.array(h['upperWeight'])[valid]
        assert np.all(upper==lower+1) and np.all((w>=0)&(w<=1))
        values=np.array([r[key] for r in profiles]);mapped=(1-w)*values[lower]+w*values[upper]
        assert np.max(abs(mapped-actual[valid])/np.maximum(1,abs(mapped)))<2e-15
    power=float(np.sum(area[valid]*np.array(h['heatFlux'],float)[valid]));assert abs(power-h['summary']['signedPowerW'])<1e-8
    digest=hashlib.sha256(h['sourceText'].encode('utf8'))
    for key in keys+['validFace','wallTemperature','lowerStation','upperStation','upperWeight']:
        x=np.array(h[key],dtype='<f8');x.view('<u8')[np.isnan(x)]=0xfff8000000000000
        digest.update(np.array([len(x),1],dtype='<f8').tobytes());digest.update(x.tobytes(order='F'))
    assert digest.hexdigest()==h['sourceSHA256']
    run=subprocess.run(['node',str(ROOT/'tests/testPhaseNoseHeatingViewerV227.js'),str(source)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend.txt').write_bytes(run.stdout);assert run.returncode==0,run.stdout
    plot(out,m,h)
    report=dict(passed=True,independentProfiles=refs,xiRelativeErrors=xi_errors,heatCandidateFaces=int(valid.sum()),
      heatCandidateAreaM2=float(area[valid].sum()),signedPowerW=power,sourceSHA256=digest.hexdigest(),
      newQualifiedAreaM2=0,physicalQualification=False,frontendSHA256=sha(ROOT/'web/viewer.js'),
      browserInteractionTested=False,elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',report);print(json.dumps(report))

def plot(out,m,h):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.colors import Normalize
    V=np.array(m['vertices']);F=np.array(m['faces']);valid=np.array(h['validFace'],bool);q=np.array(h['heatFlux'],float)/1e6
    norm=Normalize(q[valid].min(),q[valid].max());cmap=plt.get_cmap('viridis');colors=np.tile([.78,.79,.80,1.],(len(F),1));colors[valid]=cmap(norm(q[valid]))
    fig=plt.figure(figsize=(14,5.3),facecolor='white')
    for j,zoom in enumerate((False,True)):
        ax=fig.add_subplot(1,3,j+1,projection='3d');mask=np.min(V[F,0],axis=1)>.025 if zoom else np.ones(len(F),bool)
        tri=V[F[mask]];ax.add_collection3d(Poly3DCollection(tri,facecolors=colors[mask],edgecolors='none',rasterized=True))
        pts=tri.reshape(-1,3);lo=pts.min(0);hi=pts.max(0)
        ax.set_xlim(lo[0],hi[0]);ax.set_ylim(lo[1],hi[1]);ax.set_zlim(lo[2],hi[2]);ax.set_box_aspect(hi-lo);ax.view_init(18,-65)
        if zoom:
            ax.set_xticks([.03,.05]);ax.set_yticks([-.03,0,.03]);ax.set_zticks([-.03,0,.03]);ax.tick_params(labelsize=8,pad=1)
            ax.set_xlabel('X [m]');ax.set_ylabel('Y [m]');ax.set_zlabel('Z [m]')
        else:ax.set_axis_off()
        ax.set_title('Original rocket; grey = unresolved' if not zoom else 'Nose heat candidate [MW/m2]',fontsize=10)
    ax=fig.add_subplot(1,3,3);ax.plot(np.rad2deg(h['anglesRad']),[r['heatFlux']/1e6 for r in h['profiles']],'-o',ms=3)
    ax.set_xlabel('Angle from shifted stagnation [deg]');ax.set_ylabel('Heat into wall [MW/m2]');ax.grid(alpha=.25)
    ax.set_title('Local-similarity BVP stations; Tw = 300 K',fontsize=10)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=fig.axes[1],shrink=.55,pad=.02,label='MW/m2')
    fig.suptitle('M5 / alpha20: engineering nose candidate, not full-surface validation',fontsize=12)
    fig.text(.5,.035,'1306 / 19752 faces; area 0.179089%; physical qualification remains 0',ha='center')
    fig.subplots_adjust(left=.02,right=.98,top=.87,bottom=.2,wspace=.4);fig.savefig(str(out/'rocket_nose_heat.png'),dpi=160);plt.close(fig)
if __name__=='__main__':main(sys.argv[1])
