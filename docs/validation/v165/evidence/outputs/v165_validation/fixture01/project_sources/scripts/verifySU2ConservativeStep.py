"""模块用途。有限编译实际类夹具，并以独立物性/状态线段检查全部事件。"""
import csv
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from prepareSU2ConservativeStepPatch import ROOT,SOURCE,BUILD,FILE,overlay
from idealAirProperties import IdealAir
from reviewSU2StepAlgebra import thermal
from validateIdealAir import BIN,run,save,sha

OUT=ROOT/'outputs/v165_validation/fixture01'
CONFIG=ROOT/'outputs/v163_validation/small02/case.cfg'
EVENTS=ROOT/'outputs/v164_validation/real01/admissibility.csv'


def main():
    OUT.mkdir(parents=True,exist_ok=False)
    assert json.loads((ROOT/'outputs/v165_validation/algebra01/report.json').read_text())['passed']
    for name in ('scripts/verifySU2ConservativeStep.py','tests/testSU2ConservativeStep.cpp'):
        target=OUT/'project_sources'/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((ROOT/name).read_bytes())
    before={str(p.relative_to(BUILD)):sha(p) for p in BUILD.rglob('*.obj')}
    save(OUT/'objects_before.json',before)
    command=[str(BIN/'g++.exe'),'-std=c++17','-O3','-fopenmp','-ffast-math','-fno-finite-math-only',
        '-I'+str(SOURCE/'SU2_CFD/include'),'-I'+str(BUILD/'Common/src'),
        '-I'+str(SOURCE/'externals/eigen'),'-I'+str(SOURCE/'externals/mel'),
        str(ROOT/'tests/testSU2ConservativeStep.cpp'),str(SOURCE/FILE),'-Wl,--start-group',
        str(BUILD/'SU2_CFD/src/libSU2core.a'),str(BUILD/'Common/src/libSU2Common.a'),
        '-Wl,--end-group','-o',str(OUT/'fixture.exe')]
    save(OUT/'command.json',dict(compile=command,configSHA256=sha(CONFIG),eventsSHA256=sha(EVENTS)))
    with overlay(OUT):
        run(OUT,'compile',command,30)
        run(OUT,'test',[sys.executable,str(Path(__file__).resolve()),'--worker'],30)
    assert all(sha(BUILD/n)==v for n,v in before.items())
    result=json.loads((OUT/'results.json').read_text())
    rows=list(csv.DictReader((OUT/'conservative_steps.csv').open()))
    assert len(rows)==38 and [int(r['event']) for r in rows]==list(range(1,39))
    assert sum(r['old_valid']=='0' for r in rows)==1
    assert sum(r['old_valid']=='1' and float(r['theta'])==0 for r in rows)==1
    core=IdealAir();maximum=0;theta_values=[]
    for row in rows:
        if row['old_valid']=='0': continue
        old,trial,new=[np.array([float(row[p+'_'+v]) for v in ('rho','rhou','rhov','rhow','rhoE')])
                       for p in ('old','trial','accepted')]
        theta=float(row['theta']);k=float(row['k'])
        difference=float(np.max(abs(new-(old+theta*(trial-old)))/np.maximum(1,abs(new))))
        maximum=max(maximum,difference)
        state=core.state('rhoe',new[0],thermal(new,k))
        assert 60<=state['T'][0]<=2000 and difference<2e-13
        theta_values.append(theta)
    assert result['passed']
    save(OUT/'report.json',dict(passed=True,actualClass=result,events=len(rows),
        maximumSegmentError=maximum,minimumPositiveTheta=min(v for v in theta_values if v>0),
        sourceRestored=True,engineObjectsUnchanged=298,physicsQualified=False,newQualifiedArea=0))
    print(json.dumps(result))


if __name__=='__main__':
    if sys.argv[1:]==['--worker']:
        raise SystemExit(subprocess.call([str(OUT/'fixture.exe'),str(CONFIG),str(EVENTS),str(OUT/'results.json')],cwd=str(OUT)))
    main()
