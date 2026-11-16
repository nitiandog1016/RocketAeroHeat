"""模块职责：监督v242具名准备、一次外围/体网格与独立审核，绝不启动CFD。"""
from pathlib import Path
import sys
import runSU2RocketPhaseV222 as guard
ROOT=guard.ROOT;BASE=ROOT/'outputs/v242_validation';DOC=ROOT/'docs/validation/v242'
CLOSED=True


def main(stage):
    if CLOSED:raise SystemExit('v242 fixed execution closed; use NEXT_DECISION')
    specs={'prepare':('scripts/prepareForebodyPLCV242.py',[],30,'offline'),
      'surface':('scripts/surfaceForebodyPLCV242.py',[],120,'mesh-surface'),
      'build':('scripts/buildForebodyPLCV242.py',[],60,'build'),
      'tetgen':('scripts/meshForebodyPLCV242.py',[],120,'mesh-volume'),
      'audit':('scripts/auditForebodyPLCV242.py',[],60,'audit')}
    script,args,seconds,kind=specs[stage];BASE.mkdir(exist_ok=True)
    paths=[Path(__file__).resolve(),ROOT/script,DOC/'PLAN.md',ROOT/'scripts/runSU2RocketPhaseV222.py',
      ROOT/'scripts/phaseEulerRegionalVolume.py',ROOT/'scripts/reviewSU2ExplicitQuality.py',
      ROOT/'outputs/v240_validation/run01/mesh.su2',ROOT/'outputs/v240_validation/run01/original_faces.csv',
      ROOT/'outputs/v240_validation/region01/provenance.npz',ROOT/'outputs/v240_validation/region01/report.json',
      ROOT/'outputs/v241_validation/buffer01/cavity.npz']
    if stage=='surface':paths += [BASE/'prepare01/wall.npz',BASE/'prepare01/report.json']+[ROOT/'tools/runtime/gmsh/4.15.2'/n for n in ('runtime.json','gmsh.py','gmsh-4.15.dll','LICENSE.txt')]
    if stage=='build':paths += [BASE/'surface01/report.json',ROOT/'docs/validation/v142/storage_and_sources.json',
      ROOT/'third_party/tetgen/single_shell_plc_v242.cpp',ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/g++.exe']+[
      ROOT/'tools/runtime/tetgen/1.5.1'/n for n in ('source/tetgen.h','tetgen.o','predicates.o')]
    if stage=='tetgen':paths += [BASE/'build01/runtime.json',BASE/'build01/single_shell_plc.exe',BASE/'surface01/plc.dat',BASE/'surface01/plc.npz',
      BASE/'surface01/report.json',ROOT/'scripts/tetgenFixedInterface.py',ROOT/'scripts/cgalTetrahedralContract.py',ROOT/'scripts/auditSU2HybridTopology.py',ROOT/'scripts/clipEulerForebody.py']
    if stage=='audit':paths += [BASE/'tetgen01'/n for n in ('mesh.npz','mesh.su2','report.json')]+[ROOT/'tools/runtime/gmsh/4.15.2'/n for n in ('runtime.json','gmsh.py','gmsh-4.15.dll','LICENSE.txt')]
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;size=guard.size_bytes
    guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
    snapshot=BASE/('inputs_'+stage);snapshot.mkdir(exist_ok=False)
    for p in paths:
        if p.suffix in ('.py','.md','.cpp','.hpp'):(snapshot/(guard.sha(p)+p.suffix)).write_bytes(p.read_bytes())
    return guard.supervise(stage,[sys.executable,str(ROOT/script)]+args,seconds,256,paths,BASE,kind)


if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
