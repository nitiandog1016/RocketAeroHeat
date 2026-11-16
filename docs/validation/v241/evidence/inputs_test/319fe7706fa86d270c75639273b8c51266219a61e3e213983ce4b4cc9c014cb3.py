"""模块职责：用可解析三维场验证插值、方向、邻线投影及不跨域；不是实验热验证。"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from scipy.spatial import Delaunay
from scipy.integrate import solve_ivp
from tetraStreamtubes import TetraField,trace,tube,TraceStop
from runSU2RocketPhaseV222 import save


def run():
    # 这里只为解析测试建立小网格，生产器从不重新剖分真实源网格。
    x=np.array([[a,b,c] for a in (0.,.5,1.) for b in (-1.,0.,1.) for c in (-1.,0.,1.)])
    t=Delaunay(x).simplices.copy();det=np.linalg.det((x[t[:,1:]]-x[t[:,0,None]]).transpose(0,2,1))
    t=t[abs(det)>1e-14];det=det[abs(det)>1e-14];negative=det<0;t[negative,1],t[negative,2]=t[negative,2].copy(),t[negative,1].copy()
    # 常密度、线性转向速度；总能量线性，使五量P1插值具有独立精确参照。
    def exact(p):
        p=np.asarray(p);return np.column_stack([np.ones(len(p)),300+0*p[:,0],60*p[:,0],0*p[:,0],4e5+1000*p[:,2]])
    f=TetraField(x,t,exact(x));points=np.array([[.1,.2,.3],[.73,-.17,.43],[.5,0.,0.]])
    maximum=0
    for p in points:
        i,w,u=f.sample(p);maximum=max(maximum,float(abs(u-exact(p[None])[0]).max()));assert np.allclose(w@x[t[i]],p,atol=2e-14,rtol=0)
    assert maximum<1e-9
    r=trace(f,[.1,.2,.1],.007,90);q=r['xyz'];turn_error=float(abs(q[:,1]-(.2+.1*(q[:,0]**2-.1**2))).max());assert turn_error<1e-10
    oracle=solve_ivp(lambda s,p:np.array([300.,60*p[0],0.])/np.linalg.norm([300.,60*p[0],0.]),[0,r['s'][-1]],q[0],
                     method='DOP853',rtol=1e-12,atol=1e-13,t_eval=r['s'])
    oracle_error=float(abs(q-oracle.y.T).max());assert oracle_error<1e-10
    # 斜置入口标签带沿流向分量，|J|和真实横向投影必须不同。
    u=exact(x);u[:,2]=0;f2=TetraField(x,t,u)
    rays,m=tube(f2,[.1,0.,0.],[.4,1.,0.],[0.,.2,1.],.01,.01,50)
    assert np.max(abs(m['hBeta']-1))<1e-12 and np.max(abs(m['signedArea']-1))<1e-12
    assert np.all(np.linalg.norm(m['Jbeta'],axis=1)>m['hBeta'])
    assert np.max(abs(m['massFluxRatio']-1))<1e-12
    assert abs(m['HEta'][0]-1000)<1e-6 and np.all(m['valid'])
    stopped=trace(f2,[.9,0.,0.],.03,30);assert stopped['stopReason']=='TUBE_OUTSIDE_ORIGINAL_DOMAIN' and stopped['xyz'][-1,0]<=1
    failures=[]
    for name,call in [('outside',lambda:f.sample([-.1,0,0])),('nan',lambda:f.sample([np.nan,0,0])),
      ('negative-density',lambda:TetraField(x,t,-u)),('degenerate-label',lambda:tube(f2,[.1,0,0],[0,1,0],[0,1,0],.01,.01,2))]:
        try:call()
        except ValueError:failures.append(name)
        else:raise AssertionError(name)
    # 刚体旋转使轨迹/横向面积协变；避免只通过轴对齐实现。
    theta=.7;c,s=np.cos(theta),np.sin(theta);rot=np.array([[c,-s,0],[s,c,0],[0,0,1.]])
    uu=u.copy();uu[:,1:4]=u[:,1:4]@rot.T;fr=TetraField(x@rot.T,t,uu)
    rr,mm=tube(fr,np.array([.1,0,0])@rot.T,np.array([.4,1,0])@rot.T,np.array([0,.2,1])@rot.T,.01,.01,50)
    rotation_error=float(np.max(abs(rr[0]['xyz']-rays[0]['xyz']@rot.T)));assert rotation_error<1e-12
    assert np.max(abs(mm['signedArea']-m['signedArea']))<1e-11
    return dict(passed=True,affineFiveStateAbsoluteError=maximum,analyticTurningAbsoluteErrorM=turn_error,
       independentDOP853MaxDifferenceM=oracle_error,rigidRotationMaxDifferenceM=rotation_error,
       obliqueLabelProjectionPassed=True,entropyAndHNotForcedConstant=True,rejections=failures,
       outsidePrefixStations=len(stopped['s']),experimentalValidation=False,physicalQualification=False)


if __name__=='__main__':save(Path(sys.argv[1]),run())
