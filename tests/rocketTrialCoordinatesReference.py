"""独立恢复两种候选路径；只作核对，不裁状态，不调用MATLAB或新增NS。"""
import numpy as np


def trial_candidate(Q, value, scale, alpha, gas, coordinates='conservative'):
    scale = np.asarray(scale).reshape(1, 5)
    if coordinates == 'conservative':
        return Q+alpha*value*scale
    assert coordinates == 'primitive-tangent'
    rho = Q[:, 0]; velocity = Q[:, 1:4]/rho[:, None]
    cv = gas['R']/(gas['gamma']-1)
    temperature = (Q[:, 4]-.5*rho*np.sum(velocity**2, axis=1))/(rho*cv)
    direction = value*scale; drho = direction[:, 0]
    du = (direction[:, 1:4]-velocity*drho[:, None])/rho[:, None]
    dT = (direction[:, 4]-np.sum(velocity*direction[:, 1:4], axis=1)
          +(.5*np.sum(velocity**2, axis=1)-cv*temperature)*drho)/(rho*cv)
    r = rho+alpha*drho; u = velocity+alpha*du; T = temperature+alpha*dT
    return np.column_stack([r, r[:, None]*u, r*(cv*T+.5*np.sum(u**2, axis=1))])
