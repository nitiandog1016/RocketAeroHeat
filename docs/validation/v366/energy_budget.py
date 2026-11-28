"""在完整旋转域用同一步前后总能量和步末边界场独立核对质量、热量收支。"""
import numpy as np

GAUSS=(.5-.5/np.sqrt(3),.5+.5/np.sqrt(3))


def patch(f,s,t):
    """四边形双线性插值及两个参数方向上的解析导数。"""
    a,b,c,d=f[:-1,:-1],f[:-1,1:],f[1:,:-1],f[1:,1:]
    return ((1-t)*((1-s)*a+s*b)+t*((1-s)*c+s*d),
            (1-t)*(b-a)+t*(d-c),(1-s)*(c-a)+s*(d-b))


def budget(before,after,dt,c):
    shape=tuple(c["shape"]);n=shape[1]
    p0,rho0,t0=c["freePressurePa"],c["densityKgM3"],c["freeTemperatureK"]
    v0=np.sqrt(p0/rho0);cp=1.4*287.05287/.4
    h0=cp*t0+.5*c["speedMps"]**2
    for key in ("x","y","PointID"):
        np.testing.assert_array_equal(before[key],after[key])
    xy=np.stack([after[k] for k in ("x","y")],axis=-1).reshape(*shape,2)
    rates=np.stack([(after["Density"]-before["Density"])*rho0/dt,
                    (after["Energy"]-before["Energy"])*p0/dt],axis=-1).reshape(*shape,2)
    stored=np.zeros(2);volume=0.
    for s in GAUSS:
        for t in GAUSS:
            x,xs,xt=patch(xy,s,t)
            jac=xs[:,:,1]*xt[:,:,0]-xs[:,:,0]*xt[:,:,1]
            assert min(jac.ravel())>0
            dv=.25*2*np.pi*x[:,:,1]*jac
            stored+=(dv[:,:,None]*patch(rates,s,t)[0]).sum(axis=(0,1))
            volume+=dv.sum()

    fields=np.stack((after["Density"]*rho0,after["Pressure"]*p0,
        after["Velocity_x"]*v0,after["Velocity_y"]*v0,after["Temperature"]*t0),axis=-1).reshape(*shape,5)
    mass=adv=conduct=work=wall=area=0.
    q=after["Heat_Flux"][:n]
    for s in GAUSS:
        # 在远场折线上取流体侧双线性梯度，包含黏性功和热传导。
        x,xs,xt=(v[0] for v in patch(xy[-2:],s,1.))
        val,vs,vt=(v[0] for v in patch(fields[-2:],s,1.))
        rho,p,u,v,T=val.T
        grad=np.linalg.solve(np.stack((xs,xt),axis=1),np.stack((vs,vt),axis=1))
        ux,vx,Tx=grad[:,0,2:].T;ur,vr,Tr=grad[:,1,2:].T
        vec=np.column_stack((xs[:,1],-xs[:,0]))*(np.pi*x[:,1,None])
        mu=1.716e-5*(T/273.15)**1.5*(273.15+110.4)/(T+110.4)
        div=ux+vr+v/x[:,1]
        tauxx=mu*(2*ux-2*div/3);taurr=mu*(2*vr-2*div/3);tauxr=mu*(ur+vx)
        m=rho*(vec[:,0]*u+vec[:,1]*v)
        mass+=m.sum()
        # 直接积分超额焓，避免把大幅出入流各自舍入后再作差。
        adv+=(m*(cp*T+.5*(u*u+v*v)-h0)).sum()
        conduct-=((mu*cp/.72)*(vec[:,0]*Tx+vec[:,1]*Tr)).sum()
        work-=(vec[:,0]*(tauxx*u+tauxr*v)+vec[:,1]*(tauxr*u+taurr*v)).sum()
        radius=(1-s)*xy[0,:-1,1]+s*xy[0,1:,1]
        da=np.pi*radius*np.linalg.norm(np.diff(xy[0],axis=0),axis=1)
        wall+=da@((1-s)*q[:-1]+s*q[1:]);area+=da.sum()
    # 用网格子午多边形的一阶面积矩独立检查旋转体积。
    def revolution(curve):
        end=np.roll(curve,-1,axis=0)
        cross=curve[:,0]*end[:,1]-end[:,0]*curve[:,1]
        return 2*np.pi*abs(np.sum((curve[:,1]+end[:,1])*cross)/6)
    np.testing.assert_allclose(volume,revolution(xy[-1])-revolution(xy[0]),rtol=1e-11)
    storage_mass,storage_energy=stored
    far_energy=adv+conduct+work+h0*mass
    mass_defect=storage_mass+mass
    energy_defect=storage_energy+far_energy+wall
    excess=storage_energy-h0*storage_mass+adv+conduct+work+wall
    np.testing.assert_allclose(energy_defect-h0*mass_defect,excess,atol=1e-6)
    return dict(volumeM3=float(volume),revolvedWallAreaM2=float(area),massRateKgS=float(storage_mass),
        farOutMassKgS=float(mass),massImbalanceKgS=float(mass_defect),storedTotalEnergyRateW=float(storage_energy),
        farOutTotalEnergyW=float(far_energy),farOutExcessAdvectionW=float(adv),farOutConductionW=float(conduct),
        farOutViscousWorkW=float(work),wallIntoSolidW=float(wall),totalEnergyImbalanceW=float(energy_defect),
        freeTotalEnthalpyJkg=float(h0),excessEnergyImbalanceW=float(excess),nativeFluxClosureProven=False)
