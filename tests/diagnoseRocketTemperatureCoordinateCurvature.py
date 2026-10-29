"""从已保存实际状态诊断冷单元及两条有限步路径；不外推新的物理解、不计算新NS。"""
from pathlib import Path
import os
import sys
import json
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from riemannReservoirReference import ROOT, sha
from rocketTrialCoordinatesReference import trial_candidate


def main(state_file, folder):
    state_file = Path(state_file).resolve(); folder = Path(folder).resolve()
    domain_file = ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat'
    target = folder/'current_coordinate_curvature.json'
    assert not target.exists()
    files = [state_file, domain_file]
    bindings = {p.relative_to(ROOT).as_posix(): sha(p) for p in files}
    s = loadmat(state_file, simplify_cells=True, variable_names=['result', 'gas'])
    r = s['result']; gas = s['gas']; cv = gas['R']/(gas['gamma']-1)
    cut = loadmat(domain_file, simplify_cells=True, variable_names=['domain'])['domain']
    def primitive(Q):
        rho = Q[:, 0]; u = Q[:, 1:4]/rho[:, None]
        return rho, u, (Q[:, 4]-.5*rho*np.sum(u*u, axis=1))/(rho*cv)
    def tangent(Q, dQ):
        rho, u, T = primitive(Q); drho = dQ[:, 0]
        du = (dQ[:, 1:4]-u*drho[:, None])/rho[:, None]
        dT = (dQ[:, 4]-np.sum(u*dQ[:, 1:4], axis=1)
              +(.5*np.sum(u*u, axis=1)-cv*T)*drho)/(rho*cv)
        return drho, du, dT
    rho, u, T = primitive(r['Q']); _, _, temperature_rate = tangent(r['Q'], r['rate'])
    nf = len(np.atleast_1d(cut['walls']['area']))
    cold = np.argsort(T)[:12]
    cold_cells = [dict(cellID=int(i+1), sourceFacet=int(i % nf+1), layer=int(i//nf+1),
                       centroidM=cut['cells']['centroid'][i].tolist(), temperatureK=float(T[i]),
                       rhoKgM3=float(rho[i]), Mach=float(np.linalg.norm(u[i])/np.sqrt(gas['gamma']*gas['R']*T[i])),
                       actualNSTemperatureRateKPerS=float(temperature_rate[i])) for i in cold]
    trace = r['trace'] if isinstance(r['trace'], dict) else list(r['trace'])[-1]
    Q = trace['initialQ']; value = trace['linearUpdateU']; alpha = float(trace['alpha']); scale = r['residualScale']
    drho, du, dT = tangent(Q, value*scale)
    initial_rho, initial_u, initial_T = primitive(Q)
    ratio = initial_rho/(initial_rho+alpha*drho)
    kinetic = alpha**2*ratio**2*np.sum(du*du, axis=1)/(2*cv)
    exact_temperature = initial_T+alpha*ratio*dT-kinetic
    conservative = trial_candidate(Q, value, scale, alpha, gas, 'conservative')
    _, _, conservative_T = primitive(conservative)
    identity_error = float(np.max(abs(exact_temperature-conservative_T))/max(1., np.max(abs(conservative_T))))
    assert identity_error < 3e-15
    coordinates = r['options'].get('trialStateCoordinates', 'conservative')
    actual = trial_candidate(Q, value, scale, alpha, gas, coordinates)
    update_error = float(np.max(abs(actual-r['Q']))/max(1., np.max(abs(r['Q']))))
    assert update_error < 3e-15
    i = int(cold[0])
    replay = dict(initialStateIsSavedAcceptedPredecessor=True, actualCoordinates=coordinates, alpha=alpha,
                  conservativeTemperatureIdentityRelativeDifference=identity_error, acceptedUpdateRelativeDifference=update_error,
                  coldestCellAtHistoricalStep=dict(cellID=i+1, initialTemperatureK=float(initial_T[i]),
                      linearTemperatureK=float(initial_T[i]+alpha*dT[i]),
                      conservativeTemperatureK=float(conservative_T[i]),
                      densityAmplificationRatio=float(ratio[i]), kineticCurvatureCoolingK=float(kinetic[i]),
                      alphaLinearTemperatureChangeK=float(alpha*dT[i])))
    assert all(sha(ROOT/p) == h for p, h in bindings.items())
    result = dict(status='passed-current-cold-cell-and-exact-conservative-curvature-diagnostic',
                  currentTemperatureRangeK=[float(T.min()), float(T.max())], coldCells=cold_cells,
                  lowerMarginCounts={str(m): int(np.count_nonzero(T-80 < m)) for m in (.05, .1, .5, 1, 5, 10)},
                  historicalAcceptedStepReplay=replay, inputSHA256=bindings, scriptSHA256=sha(__file__),
                  historicalDirectionUsedToAdvanceCurrentState=False, physicalEquilibriumBelow80KInferred=False,
                  newActualNSEvaluations=0, newFlowAcceptances=0, physicalSolution=False, newProductionSolvedAreaM2=0)
    folder.mkdir(parents=True, exist_ok=True)
    target.write_bytes((json.dumps(result, indent=2)+'\n').encode('utf8'))
    print(json.dumps(dict(status=result['status'], minimumTemperatureK=float(T.min()), coldest=cold_cells[0], replay=replay)))


if __name__ == '__main__':
    main(*sys.argv[1:])
