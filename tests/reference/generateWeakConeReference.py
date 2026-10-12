"""独立数值校核：SciPy DOP853 积分无量纲定比热 Taylor–Maccoll 方程。

用于检查 MATLAB 三变量变比热算法在 80~150 K 冻结比热域内的极限。
不调用 MATLAB，不读取 MATLAB 求解结果；这不是实验数据。
"""
from pathlib import Path
import json
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# 与项目材料模型相同的 cp(150 K) 定义；参考流动解完全独立计算。
t=150.0
low=np.array([[3.53100528,-1.23660987e-4,-5.02999433e-7,2.43530612e-9,-1.40881235e-12],
              [3.78245636,-2.99673416e-3,9.84730201e-6,-9.68129509e-9,3.24372837e-12],
              [2.5,0,0,0,0]])
cp=float(np.sum(np.array([.755,.232,.013])*8.31446261815324/
                np.array([.0280134,.0319988,.039948])*(low@t**np.arange(5))))
R=287.05; gamma=cp/(cp-R); M=4.0; T1=100.0
mu=np.arcsin(1/M); scale=np.sqrt((gamma-1)*M*M/(2+(gamma-1)*M*M))

def integrate(beta,cone,rtol):
    mn2=(M*np.sin(beta))**2
    density=(gamma+1)*mn2/((gamma-1)*mn2+2)
    initial=np.array([scale*np.cos(beta),-scale*np.sin(beta)/density])
    def rhs(theta,y):
        vr,vt=y; a=(gamma-1)/2*(1-vr*vr-vt*vt)
        return [vt,(vt*vt*vr-a*(2*vr+vt/np.tan(theta)))/(a-vt*vt)]
    def wall(theta,y):return y[1]
    wall.terminal=True
    sol=solve_ivp(rhs,(beta,cone),initial,method='DOP853',rtol=rtol,atol=rtol*.01,
                  first_step=min((beta-cone)/100,max((beta-mu)*.1,1e-12)),
                  max_step=(beta-cone)/40,events=wall)
    if not sol.success:raise RuntimeError(sol.message)
    residual=sol.t[-1]-cone if len(sol.t_events[0]) else sol.y[1,-1]
    return residual,sol.y[:,-1],initial,mn2

def solve(angle,rtol):
    cone=np.deg2rad(angle)
    grid=np.unique(np.r_[mu+np.logspace(-8,-1,60),np.linspace(mu+.1,np.pi/2-.001,60)])
    previous=None
    for beta in grid:
        if beta<=cone:continue
        residual=integrate(beta,cone,rtol)[0]
        if previous is not None and previous[1]*residual<=0:
            root=brentq(lambda b:integrate(b,cone,rtol)[0],previous[0],beta,xtol=2e-14)
            _,w,shock,mn2=integrate(root,cone,rtol)
            ratioT=(1-np.dot(w,w))/(1-np.dot(shock,shock))
            p2=(2*gamma*mn2-(gamma-1))/(gamma+1)
            pw=p2*ratioT**(gamma/(gamma-1))
            tw=T1*(1+(gamma-1)*M*M/2)*(1-np.dot(w,w))
            return [angle,float(np.rad2deg(root)),float(pw),float(tw)]
        previous=(beta,residual)
    raise RuntimeError('no weak bracket')

angles=[.5,1,2,3.75,4,5,10]
rows=np.array([solve(a,1e-11) for a in angles])
fine=np.array([solve(a,1e-12) for a in angles])
assert np.all(fine[:,3]<150)
assert np.max(np.abs(rows[:,1:]-fine[:,1:]))<1e-6
result={'evidenceType':'independent numerical constant-cp limit, not experiment',
        'solver':'SciPy DOP853 + Brent, normalized two-variable Taylor-Maccoll',
        'scipyVersion':scipy.__version__,'Mach':M,'temperature':T1,'R':R,'gamma':gamma,
        'specificHeat':cp,'columns':['coneAngleDeg','shockAngleDeg','pressureRatio','wallTemperature'],
        'rows':fine.tolist(),'refinementMaximumAbsoluteChange':float(np.max(abs(rows[:,1:]-fine[:,1:])))}
dest=Path(__file__).resolve().parents[2]/'docs/validation/weak_cone_reference.json'
dest.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
