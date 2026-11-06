"""模块用途。只读核对实际包的SST更新路径与v166实测能量对，不启动引擎。"""
import json
from pathlib import Path
import time
import csv
import numpy as np
from su2ThermalRuntime import ROOT, read, sha
from buildSU2IdealAirRuntime import SOURCE, BUILD
from idealAirProperties import IdealAir
from validateIdealAir import save

BASE=ROOT/'outputs/v167_validation'
PARENT=ROOT/'outputs/v164_validation/build01'


def main():
    started=time.monotonic(); out=BASE/'path01';out.mkdir(parents=True,exist_ok=False)
    tree=read(PARENT/'source_tree.json'); selected={}
    names=('include/solvers/CScalarSolver.inl','include/variables/CVariable.hpp',
           'src/solvers/CTurbSolver.cpp','src/solvers/CTurbSSTSolver.cpp',
           'src/integration/CSingleGridIntegration.cpp','src/iteration/CFluidIteration.cpp',
           'src/solvers/CNSSolver.cpp','src/variables/CNSVariable.cpp')
    locks={}
    for short in names:
        name='SU2_CFD/'+short; snapshot=PARENT/'executed_sources'/name
        actual=snapshot if snapshot.exists() else SOURCE/name
        assert sha(actual)==tree[str(SOURCE/name)]
        target=out/'exact_sources'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(actual.read_bytes())
        selected[short]=actual.read_text('utf8');locks[str(actual)]=sha(actual)
    scalar=selected['include/solvers/CScalarSolver.inl'];variable=selected['include/variables/CVariable.hpp']
    assert 'flowNodes->GetDensity(iPoint)' in scalar and 'density, density);' in scalar
    assert 'Solution_Old(iPoint,iVar)*Sol2Conservative_old + solution' in variable
    assert scalar.index('ComputeUnderRelaxationFactor(solver_container, config);') < scalar.index('nodes->AddClippedSolution',scalar.index('CompleteImplicitIteration'))
    single=selected['src/integration/CSingleGridIntegration.cpp']
    assert single.index('->Preprocessing(')<single.index('->Set_OldSolution();')<single.index('Time_Integration(')
    assert 'solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0)' in selected['src/solvers/CNSSolver.cpp']
    events_path=ROOT/'outputs/v166_validation/real02/admissibility.csv'
    rows=list(csv.DictReader(events_path.open()));a,b=rows[335],rows[337]
    cols=('rho','rhou','rhov','rhow','rhoE')
    u=np.array([float(a['old_'+n]) for n in cols]);v=np.array([float(b['old_'+n]) for n in cols])
    assert np.array_equal(u,v) and a['local_point']==b['local_point']=='292375'
    e0=u[4]/u[0]-.5*np.sum((u[1:4]/u[0])**2)
    core=IdealAir();emin=float(core.state('PT',101325,60)['e'][0]);oldk,newk=float(a['k']),float(b['k'])
    assert e0-oldk>=emin>e0-newk
    cap=(e0-emin-oldk)/(newk-oldk)
    assert 0<cap<1
    objects=read(ROOT/'outputs/v165_validation/build01/objects_after.json')
    actual={str(p.relative_to(BUILD)):sha(p) for p in BUILD.rglob('*.obj')}
    assert actual==objects and len(actual)==298
    parent_objects=read(PARENT/'objects_after.json')
    changed=[n for n in actual if actual[n]!=parent_objects[n]]
    assert len(changed)==2 and all('CNSVariable' in n or 'CNSSolver' in n for n in changed)
    locks[str(events_path)]=sha(events_path);locks.update(core.locks)
    save(out/'report.json',dict(passed=True,sourceSHA256=locks,sourceTreeSHA256=sha(PARENT/'source_tree.json'),
        update='clip((Solution_Old*rho+alpha*LinSysSol)/rho); same current rho in numerator and denominator',
        sharedKAndOmegaAlpha=True,oldSolutionSavedBeforeSSTSpaceIntegration=True,
        recordedPair=dict(firstEvent=336,lastEvent=338,thermalWithoutK=e0,minimumEnergy=emin,
            oldK=oldk,newK=newk,maximumAnalyticFractionOfObservedKChange=float(cap),
            actualSSTIncrementUnknown=True,globalPointUnknown=True),
        objectsRequiringParentRestoration=changed,nativeRuns=0,elapsedSeconds=time.monotonic()-started,
        physicsQualified=False,newQualifiedArea=0))
    (out/'executed_audit.py').write_bytes(Path(__file__).read_bytes())
    assert time.monotonic()-started<120
    print(json.dumps(dict(passed=True,seconds=time.monotonic()-started,observedPairMaximumFraction=cap)))


if __name__=='__main__':main()
