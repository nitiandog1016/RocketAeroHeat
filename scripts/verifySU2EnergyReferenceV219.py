"""模块职责：用独立物理通量、能量账及热壁导数验证坐标接口，不调用CFD。"""
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib,json,time
import numpy as np
from phaseAwareDryAir import pressure_ceiling,SOURCE
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2EnergyReference import (PHYSICAL_REFERENCE as P,SOLVER_REFERENCE as S,SHIFT_JKG as C,
                               flow_packet,convert_flow,specific_energy,contract)

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v219_validation/energy01'

def scaled(actual,reference):
    return float(np.max(np.abs(np.asarray(actual)-reference)/np.maximum(1.,np.abs(reference))))

def main():
    start=time.monotonic();OUT.mkdir(exist_ok=False);air=PhaseAwareViscousAir()
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985'
    t=np.array([60,75,100,132.6312,180,226.50908361133003,298.15,300,600,1000,1500,2000.])
    p=.25*pressure_ceiling(t);v=np.c_[-np.linspace(30,1500,12),np.linspace(-25,50,12),np.linspace(5,500,12)]
    state=air.evaluate(p,t);rho=state['densityKgm3'];e=state['internalEnergyJkg'];h=state['enthalpyJkg']
    energy=e+np.sum(v*v,axis=1)/2
    U=np.c_[rho,rho[:,None]*v,rho*energy]
    packet=flow_packet(U,P,'state');shifted=convert_flow(packet,S)
    direct=np.c_[rho,rho[:,None]*v,rho*(e+C+np.sum(v*v,axis=1)/2)]
    errors={'conserved':scaled(shifted['values'],direct),'roundTrip':scaled(convert_flow(shifted,P)['values'],U)}
    assert np.array_equal(packet['values'],U)
    assert np.array_equal(convert_flow(shifted,S)['values'],shifted['values'])
    assert np.min(air.Cp-air.R)>0
    min_e=float(air.thermal(60.)[1]-air.R*60)
    assert min_e+C>0

    # 独立60位物理通量：从rho、速度、压力、热焓定义组装，不读接口结果。
    normals=np.array([[1.,0,0],[0,1.,0],[0,0,1.],[-.3,.4,.5]])
    f=np.zeros((12,4,5));fs=np.zeros_like(f)
    with localcontext() as ctx:
        ctx.prec=60
        D=lambda x:Decimal.from_float(float(x))
        for i in range(12):
            rv=D(rho[i]);pv=D(p[i]);vv=list(map(D,v[i]));hv=D(h[i])
            ke=sum(x*x for x in vv)/2
            for j,n in enumerate(normals):
                nn=list(map(D,n));vn=sum(x*y for x,y in zip(vv,nn));mass=rv*vn
                row=[mass]+[mass*x+pv*y for x,y in zip(vv,nn)]+[mass*(hv+ke)]
                f[i,j]=list(map(float,row))
                row[-1]=mass*(hv+ke+D(C));fs[i,j]=list(map(float,row))
    errors['physicalFlux']=scaled(convert_flow(flow_packet(f,P,'flux'),S)['values'],fs)
    # 同一闭合账中正负面有向组装；转换在求范数之前，禁止两个RMS相减。
    weights=np.array([1.,-.75,.125,-.5])
    residual=np.sum(f*weights[None,:,None],axis=1)
    solver_residual=np.sum(fs*weights[None,:,None],axis=1)
    errors['signedResidual']=scaled(convert_flow(flow_packet(solver_residual,S,'residual'),P)['values'],residual)
    counter=np.array([[1.,0,0,0,0],[-1.,0,0,0,1.]])
    counter_s=convert_flow(flow_packet(counter,P,'residual'),S)['values']
    false_norm=np.linalg.norm(counter_s[:,4])-C*np.linalg.norm(counter_s[:,0])
    assert abs(false_norm-np.linalg.norm(counter[:,4]))>.1

    # general EOS组合只检验能量参考代数，不实现或宣称测试了Roe数值通量。
    kappa=air.R/state['cvJkgK'];chi=air.R*t-kappa*e
    shifted_a2=(chi-kappa*C)+kappa*(h+C)
    errors['soundSpeedSquared']=scaled(shifted_a2,state['soundSpeedMps']**2)
    pure_viscous=np.array([[0.,2.,-3.,4.,-500.],[0.,-.1,.2,0.,45.]])
    assert np.array_equal(convert_flow(flow_packet(pure_viscous,P,'flux'),S)['values'],pure_viscous)

    inverse_error=0.;wall_error=0.
    for temp in [75.,180.,226.50908361133003,300.,600.,1500.]:
        pp=float(.1*pressure_ceiling(temp));st=air.evaluate(pp,temp);rr=float(st['densityKgm3'])
        vv=np.array([-700.,40.,200.]);E=float(st['internalEnergyJkg'])+C+np.dot(vv,vv)/2
        u=np.r_[rr,rr*vv,rr*E];cv=float(st['cvJkgK'])
        jac=np.r_[(np.dot(vv,vv)-E)/(rr*cv),-vv/(rr*cv),1/(rr*cv)]
        def decode(x):
            ee=x[4]/x[0]-np.dot(x[1:4],x[1:4])/(2*x[0]**2)-C
            # 量热反解本身不依赖压力；合法低压用于复用已有严格反解，然后检查真实p/T。
            tt=float(air.inverse(1.,ee,kind='e'))
            air.evaluate(x[0]*air.R*tt,tt)
            return tt
        inverse_error=max(inverse_error,abs(decode(u)-temp))
        for j in range(5):
            delta=abs(u[j])*1e-6;up=u.copy();um=u.copy();up[j]+=delta;um[j]-=delta
            numerical=(decode(up)-decode(um))/(2*delta)
            wall_error=max(wall_error,abs(numerical-jac[j])/max(1.,abs(jac[j])))
    assert inverse_error<=1e-7 and wall_error<=2e-6
    assert max(errors.values())<=2e-12,errors

    # 失配与非法输入拒绝；通量/残差质量为负仍合法。
    rejected=[]
    def reject(label,fn):
        try:fn()
        except ValueError:rejected.append(label)
        else:raise AssertionError('not rejected: '+label)
    valid=flow_packet([1,2,3,4,5],P,'state')
    for key,value in [('schema','old'),('gasModel','OTH'),('domain','unbounded'),('transport','Sutherland'),
                      ('reference','unknown'),('quantity','RMS'),('units','nondimensional')]:
        reject(key,lambda key=key,value=value:convert_flow(dict(valid,**{key:value}),S))
    reject('missingIdentity',lambda:convert_flow({'values':[1,2,3,4,5]},S))
    reject('unknownTarget',lambda:convert_flow(valid,'unknown'))
    for label,value in [('sixComponents',[1]*6),('scalar',1),('empty',[]),('NaN',[1,2,3,4,np.nan]),
                        ('Inf',[1,2,3,4,np.inf]),('zeroDensity',[0,2,3,4,5]),('negativeDensity',[-1,2,3,4,5])]:
        reject(label,lambda value=value:flow_packet(value,P,'state'))
    reject('overflow',lambda:convert_flow(flow_packet([1e308,0,0,0,0],P,'flux'),S))
    reject('specificNaN',lambda:specific_energy(np.nan,P,S))
    reject('specificEmpty',lambda:specific_energy([],P,S))
    reject('specificUnknown',lambda:specific_energy([1],'unknown',S))
    for value in [-1.,0.,1.]:
        for quantity in ('residual','flux','increment'):
            convert_flow(flow_packet([value,0,0,0,2],P,quantity),S)
    assert float(specific_energy(specific_energy(-200000.,P,S),S,P))==-200000.
    report=dict(passed=True,contract=contract(),temperatureCount=len(t),minimumCvJkgK=float(np.min(air.Cp-air.R)),
        minimumPhysicalEnergyJkg=min_e,minimumSolverEnergyJkg=min_e+C,scaledErrors=errors,
        temperatureInverseMaxAbsK=inverse_error,wallDerivativeScaledError=wall_error,
        rejectionCount=len(rejected),rejections=rejected,wrongRMSSubtractionCounterexample=float(false_norm),
        seconds=time.monotonic()-start,nativeQueries=0,CFDRuns=0,newQualifiedAreaM2=0.)
    (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:report[k] for k in ('passed','scaledErrors','wallDerivativeScaledError','rejectionCount','seconds')}))

if __name__=='__main__':main()
