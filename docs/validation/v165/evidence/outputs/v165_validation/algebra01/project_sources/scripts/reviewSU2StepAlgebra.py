"""模块用途。复用真实拒绝状态，独立核验有限守恒线段及原生调用顺序。"""
import csv
import json
from pathlib import Path
import numpy as np
from idealAirProperties import IdealAir
from su2ThermalRuntime import ROOT, sha
from validateIdealAir import save


def thermal(u, k):
    return u[4]/u[0] - .5*np.dot(u[1:4]/u[0], u[1:4]/u[0]) - k


def main():
    out = ROOT/'outputs/v165_validation/algebra01'
    out.mkdir(parents=True, exist_ok=False)
    source = ROOT/'outputs/v164_validation/real01/admissibility.csv'
    core = IdealAir()
    rows = list(csv.DictReader(source.open()))
    results = []
    for row in rows:
        k = float(row['k'])
        old, trial = [np.array([float(row[p+'_'+v]) for v in ('rho','rhou','rhov','rhow','rhoE')])
                      for p in ('old','trial')]
        core.state('rhoe', old[0], thermal(old,k))
        rejected = False
        try:
            core.state('rhoe', trial[0], thermal(trial,k))
        except ValueError:
            rejected = True
        assert rejected
        theta = 1.0
        for attempt in range(1,33):
            theta *= .5
            accepted = old+theta*(trial-old)
            try:
                state = core.state('rhoe',accepted[0],thermal(accepted,k))
                break
            except ValueError:
                pass
        else:
            raise AssertionError('实际事件没有找到非零可接受步')
        assert np.linalg.norm(accepted-old)>0 and 60<=state['T'][0]<=2000
        # 对独立能量不等式复核，不以正内能取代固定物性上下限。
        bounds=core.state('PT', [101325,101325], [60,2000])['e']
        assert bounds[0] <= thermal(accepted,k) <= bounds[1]
        results.append(dict(event=int(row['event']), localRCM=int(row['local_point']),
            theta=theta, trials=attempt, acceptedTemperatureK=float(state['T'][0]),
            oldTemperatureK=float(core.state('rhoe',old[0],thermal(old,k))['T'][0]),
            acceptedState=accepted.tolist(), currentK=k,
            linearIncrementKnown=False, globalPointKnown=False))
    assert len(results)==14
    parent=ROOT/'outputs/v164_validation/build01'
    source_root=ROOT/'tools/runtime/su2-thermal-build/source'
    files=['SU2_CFD/src/solvers/CEulerSolver.cpp','SU2_CFD/src/solvers/CNSSolver.cpp',
        'SU2_CFD/include/solvers/CFVMFlowSolverBase.inl',
        'SU2_CFD/src/integration/CMultiGridIntegration.cpp',
        'SU2_CFD/src/integration/CSingleGridIntegration.cpp',
        'SU2_CFD/src/iteration/CFluidIteration.cpp', 'SU2_CFD/src/solvers/CTurbSolver.cpp']
    locked=json.loads((parent/'source_tree.json').read_text())
    sources={}
    for name in files:
        snap=parent/'executed_sources'/name
        src=snap if snap.exists() else source_root/name
        assert sha(src)==locked[str(source_root/name)], name
        target=out/'native_sources'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(src.read_bytes())
        sources[name]=dict(source=str(src),sha256=sha(src))
    for name in ('scripts/reviewSU2StepAlgebra.py','docs/validation/v165/PLAN.md'):
        target=out/'project_sources'/name
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((ROOT/name).read_bytes())
    save(out/'report.json',dict(passed=True,events=results,
        sourceSHA256={str(source):sha(source),**core.locks}, nativeSources=sources,
        minimumTheta=min(r['theta'] for r in results),maximumTheta=max(r['theta'] for r in results),
        originalGlobalConservationKnown=False,physicsQualified=False,newQualifiedArea=0))
    print(json.dumps(dict(passed=True,events=14,theta=[r['theta'] for r in results])))


if __name__=='__main__':
    main()
