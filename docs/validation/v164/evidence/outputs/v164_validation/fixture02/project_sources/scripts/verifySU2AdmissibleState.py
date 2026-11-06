"""模块用途。有限编译真实CNSVariable夹具，审核其状态恢复和线程安全记录。"""
import csv
import json
from pathlib import Path
import subprocess
import sys

from prepareSU2AdmissiblePatch import ROOT, SOURCE, BUILD, FILE, overlay
from validateIdealAir import BIN, run, save, sha

OUT = ROOT/'outputs/v164_validation/fixture02'


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    for name in ('scripts/verifySU2AdmissibleState.py','tests/testSU2AdmissibleState.cpp'):
        target=OUT/'project_sources'/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((ROOT/name).read_bytes())
    objects={str(p.relative_to(BUILD)):sha(p) for p in BUILD.rglob('*.obj')}
    save(OUT/'objects_before.json',objects)
    executable=OUT/'fixture.exe'
    config=ROOT/'outputs/v163_validation/small02/case.cfg'
    command=[str(BIN/'g++.exe'),'-std=c++17','-O3','-fopenmp','-ffast-math','-fno-finite-math-only',
        '-I'+str(SOURCE/'SU2_CFD/include'),'-I'+str(BUILD/'Common/src'),
        '-I'+str(SOURCE/'externals/eigen'),'-I'+str(SOURCE/'externals/mel'),
        str(ROOT/'tests/testSU2AdmissibleState.cpp'),str(SOURCE/FILE),'-Wl,--start-group',
        str(BUILD/'SU2_CFD/src/libSU2core.a'),str(BUILD/'Common/src/libSU2Common.a'),
        '-Wl,--end-group','-o',str(executable)]
    save(OUT/'command.json',dict(compile=command,configSHA256=sha(config)))
    with overlay(OUT):
        run(OUT,'compile',command,30)
        run(OUT,'test',[sys.executable,str(Path(__file__).resolve()),'--worker'],30)
    stable={name:sha(BUILD/name)==digest for name,digest in objects.items()}
    assert all(stable.values())
    save(OUT/'objects_unchanged.json',stable)
    result=json.loads((OUT/'results.json').read_text())
    assert result['passed']
    with (OUT/'admissibility.csv').open() as stream:
        rows=list(csv.DictReader(stream))
    assert len(rows)==21 and [int(row['event']) for row in rows]==list(range(1,22))
    assert sum(int(row['recovered']) for row in rows)==20
    assert sum(row['recovered']=='0' for row in rows)==1
    save(OUT/'report.json',dict(passed=True,realClassTests=result,traceRows=len(rows),
        localIndexOnly=True,sourceRestored=True,originalEngineObjectsUnchanged=298,
        executableSHA256=sha(executable),newQualifiedArea=0))
    print(json.dumps(result))


if __name__=='__main__':
    if sys.argv[1:]==['--worker']:
        raise SystemExit(subprocess.call([str(OUT/'fixture.exe'),str(ROOT/'outputs/v163_validation/small02/case.cfg'),
            str(OUT/'results.json')],cwd=str(OUT)))
    main()
