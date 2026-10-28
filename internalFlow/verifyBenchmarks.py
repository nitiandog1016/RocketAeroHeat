# -*- coding: utf-8 -*-
"""独立解析已求MAT稳态：定比热平壁导热/Couette闭式关系及能量账。"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--folder', type=Path, required=True)
args = parser.parse_args()
data = loadmat(str(args.folder/'advance_test_evidence.mat'), squeeze_me=True, struct_as_record=False)
couette = data['couette']
thermal = data['thermal']
primitive = np.asarray(couette.finalSpatialInfo.primitive)
assert primitive.shape == (8,5)
eta = (np.arange(8)+.5)/8
gamma, gas_r, mu, pr, wall_t, speed, length = 1.4, 1., .5, .72, 300., 3., .06
cp = gamma*gas_r/(gamma-1)
conductivity = mu*cp/pr
amplitude = mu*speed**2/(2*conductivity)
exact_u = speed*eta
exact_t = wall_t+amplitude*eta*(1-eta)
offset = amplitude/(4*8**2)
u_error = np.max(np.abs(primitive[:,1]-exact_u))
t_grid_error = np.max(np.abs(primitive[:,4]-exact_t))
t_stop_error = np.max(np.abs(primitive[:,4]-(exact_t+offset)))
assert u_error < 2e-6 and t_stop_error < 3e-5
assert abs(t_grid_error-offset) < 3e-5
walls = np.atleast_1d(couette.wallFaces)
exact_heat = mu*speed**2/(2*length)
q_error = max(abs(float(w.heatFluxIntoWall)/exact_heat-1) for w in walls)
heat = sum(float(w.heatFluxIntoWall)*float(w.area) for w in walls)
power = -sum(np.dot(np.asarray(w.velocityInput),np.asarray(w.traction))*float(w.area) for w in walls)
assert q_error<1e-5 and abs(heat-power)/abs(power)<1e-5
assert abs(power/(mu*speed**2/length)-1)<1e-5
assert couette.steadyConverged and not couette.physicalSolution and not couette.eligibleForSurfaceCoverage
history = np.asarray(couette.history.residualNorm)
assert np.all(history[-8:]<=2e-6)
thermal_p = np.asarray(thermal.finalSpatialInfo.primitive)
thermal_eta=(np.arange(6)+.5)/6
thermal_error=np.max(np.abs(thermal_p[:,4]-(300+30*thermal_eta)))
thermal_exact=conductivity*30/length
thermal_q_error=max(abs(abs(float(w.heatFluxIntoWall))/thermal_exact-1) for w in np.atleast_1d(thermal.wallFaces))
assert thermal_error<3e-5 and thermal_q_error<1e-5
record=dict(oracle='independent Python reading actual MATLAB states; artificial constant-property gas',
            actualRocketSolved=False, experimentalValidation=False,
            couetteVelocityErrorMps=float(u_error),couetteTemperatureGridErrorK=float(t_grid_error),
            couettePredictedGridOffsetK=float(offset),couetteStoppingErrorK=float(t_stop_error),
            couetteWallHeatRelativeError=float(q_error),couetteHeatW=float(heat),couettePowerW=float(power),
            thermalTemperatureErrorK=float(thermal_error),thermalWallHeatRelativeError=float(thermal_q_error))
(args.folder/'independent_benchmarks.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig, axes=plt.subplots(2,2,figsize=(9,6),constrained_layout=True)
x=np.linspace(0,1,250)
axes[0,0].plot(x,speed*x,'k-',label='Analytical')
axes[0,0].plot(eta,primitive[:,1],'o',ms=4,color='#2b8397',label='Computed steady state')
axes[0,0].set(xlabel='y / L',ylabel='Velocity u / m s$^{-1}$',title='Couette velocity')
axes[0,0].legend(frameon=False)
axes[0,1].plot(x,amplitude*x*(1-x),'k-',label='Continuous analytical')
axes[0,1].plot(eta,primitive[:,4]-wall_t,'o',ms=4,color='#2b8397',label='Computed steady state')
axes[0,1].set(xlabel='y / L',ylabel='$T-T_w$ / K',title='Viscous heating; grid offset shown')
axes[0,1].legend(frameon=False)
axes[1,0].semilogy(np.arange(len(history)),history,color='#2b8397')
axes[1,0].axhline(2e-6,color='#555',ls='--',lw=.8)
axes[1,0].set(xlabel='Accepted RK2 step',ylabel='Scaled spatial residual / s$^{-1}$',title='Actual convergence history')
axes[1,1].bar(['Wall heat out','Moving-wall power in'],[heat,power],color=['#2b8397','#6c7b8b'],width=.5)
axes[1,1].set(ylabel='Power / W',title='Steady energy balance')
fig.suptitle('Numerical verification: artificial constant-property ideal gas; not rocket/experimental validation',fontsize=10)
fig.savefig(str(args.folder/'couette_verification.png'),dpi=160)
plt.close(fig)
print(json.dumps(record,indent=2))
