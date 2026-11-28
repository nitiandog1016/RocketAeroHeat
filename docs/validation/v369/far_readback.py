"""在同一个低速半步的远界上，对照积分位置与特征状态对单侧输运读数的影响。"""
import argparse
from pathlib import Path
import numpy as np
from short_step import BASE, DOC, OLD, numbers, read, save, sha
from energy_budget import GAUSS

GAMMA = 1.4
GAS = 287.05287
CP = GAMMA*GAS/(GAMMA-1)


def transport(state, normal_area):
    """给定 SI 状态及旋转面积向量，输出质量与总焓的分项和带符号净值。"""
    rho,u,v,p = state.T
    mass = rho*(u*normal_area[:,0]+v*normal_area[:,1])
    enthalpy = GAMMA/(GAMMA-1)*p/rho + .5*(u*u+v*v)
    energy = mass*enthalpy
    return dict(inMassKgS=float(np.minimum(mass,0).sum()),outMassKgS=float(np.maximum(mass,0).sum()),
                netMassKgS=float(mass.sum()),netAdvectiveW=float(energy.sum()))


def calculate(xy, primitive, c):
    gamma1 = GAMMA-1
    rho0,p0,u0 = c["densityKgM3"],c["freePressurePa"],c["speedMps"]
    a0 = np.sqrt(GAMMA*p0/rho0)
    assert len(xy)==109 and xy[0,1]==xy[-1,1]==0
    edge = np.diff(xy,axis=0)
    vector = np.column_stack((edge[:,1],-edge[:,0]))
    # 先按双线性场的边界限制插值原始量。能量沿用原收支的静温插值。
    rho,u,v,p,T = primitive.T
    bulk = np.column_stack((rho,u,v,p))
    m_in=m_out=advection=0.
    for s in GAUSS:
        state=(1-s)*primitive[:-1]+s*primitive[1:]
        r=(1-s)*xy[:-1,1]+s*xy[1:,1]
        rho_s,us,vs,_,ts=state.T
        mass=np.pi*r*rho_s*(us*vector[:,0]+vs*vector[:,1])
        m_in+=np.minimum(mass,0).sum()
        m_out+=np.maximum(mass,0).sum()
        advection+=mass@(CP*ts+.5*(us**2+vs**2))
    gauss=dict(inMassKgS=float(m_in),outMassKgS=float(m_out),netMassKgS=float(m_in+m_out),
               netAdvectiveW=float(advection))
    dual=.5*np.vstack((vector[0],vector[:-1]+vector[1:],vector[-1]))
    area=dual*(2*np.pi*xy[:,1,None])
    unit=dual/np.linalg.norm(dual,axis=1)[:,None]
    assert np.all(np.sum((xy-[c["lengthM"]/2,0])*unit,axis=1)>0)
    free=np.tile([rho0,u0,0,p0],(len(xy),1))
    uniform=transport(free,area)
    assert abs(uniform["netMassKgS"])<1e-9 and abs(uniform["netAdvectiveW"])<1e-5
    # 与固定版 CEulerSolver 的亚声速静止网格分支对应，原生 Roe 边界通量另需读回。
    def characteristic(interior):
        den,ux,ur,pres=interior.T
        un=ux*unit[:,0]+ur*unit[:,1]
        unf=u0*unit[:,0]
        assert max(abs(unf))<a0 and max(abs(un)/np.sqrt(GAMMA*pres/den))<1
        plus=un+2*np.sqrt(GAMMA*pres/den)/gamma1
        minus=unf-2*a0/gamma1
        vn=.5*(plus+minus)
        sound=.25*gamma1*(plus-minus)
        assert min(sound)>0
        out=unf>0
        ux_char=np.where(out,ux,u0)+(vn-np.where(out,un,unf))*unit[:,0]
        ur_char=np.where(out,ur,0)+(vn-np.where(out,un,unf))*unit[:,1]
        entropy=np.where(out,den**GAMMA/pres,rho0**GAMMA/p0)
        density=(entropy*sound**2/GAMMA)**(1/gamma1)
        return np.column_stack((density,ux_char,ur_char,density*sound**2/GAMMA))
    np.testing.assert_allclose(characteristic(free),free,rtol=1e-12,atol=1e-10)
    return dict(fluidGauss=gauss,fluidNodal=transport(bulk,area),characteristicNodal=transport(characteristic(bulk),area),
                uniformNodal=uniform, nativeNumericalFluxMeasured=False)


def capture(out):
    folder=BASE/"pair01/half1"
    frame=folder/"restart_00001.csv"
    rec=read(folder/"step.json")
    assert rec["passed"] and sha(frame)==rec["frameSHA256"]=="7c783fbdee774c07729c86f6cc39daa778f84955a3111ef5de0fd10bea88f970"
    c=read(OLD/"condition.json")
    a=numbers(frame)
    end=slice(-109,None)
    scale=np.sqrt(c["freePressurePa"]/c["densityKgM3"])
    xy=np.column_stack((a["x"][end],a["y"][end]))
    fields=np.column_stack((a["Density"][end]*c["densityKgM3"],a["Velocity_x"][end]*scale,
                            a["Velocity_y"][end]*scale,a["Pressure"][end]*c["freePressurePa"],
                            a["Temperature"][end]*c["freeTemperatureK"]))
    r=calculate(xy,fields,c)
    b=rec["fieldBalance"]
    np.testing.assert_allclose(r["fluidGauss"]["netMassKgS"],b["farOutMassKgS"],atol=1e-9,rtol=0)
    np.testing.assert_allclose(r["fluidGauss"]["netAdvectiveW"],
                               b["farOutTotalEnergyW"]-b["farOutConductionW"]-b["farOutViscousWorkW"],atol=1e-5,rtol=0)
    r.update(frameSHA256=sha(frame),condition={k:c[k] for k in ("lengthM","densityKgM3","freePressurePa","speedMps")})
    out.mkdir(parents=True,exist_ok=False)
    with (out/"far_nodes.csv").open("x",encoding="ascii",newline="\n") as f:
        np.savetxt(f,np.column_stack((xy,fields)),fmt="%.17g",delimiter=",",comments="",header="x,r,rho,u,v,p,T")
    save(out/"far.json",r)
    print(r)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    args=p.parse_args()
    capture(args.out)
