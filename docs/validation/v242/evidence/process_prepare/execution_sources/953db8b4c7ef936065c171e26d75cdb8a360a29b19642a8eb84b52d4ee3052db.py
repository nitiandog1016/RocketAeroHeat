"""模块职责：监督v242具名准备、一次外围/体网格与独立审核，绝不启动CFD。"""
from pathlib import Path
import sys
import runSU2RocketPhaseV222 as guard
ROOT=guard.ROOT;BASE=ROOT/'outputs/v242_validation';DOC=ROOT/'docs/validation/v242'
CLOSED=False


def main(stage):
    if CLOSED:raise SystemExit('v242 fixed execution closed; use NEXT_DECISION')
    specs={'prepare':('scripts/prepareForebodyPLCV242.py',[],30,'offline')}
    script,args,seconds,kind=specs[stage];BASE.mkdir(exist_ok=True)
    paths=[Path(__file__).resolve(),ROOT/script,DOC/'PLAN.md',ROOT/'scripts/runSU2RocketPhaseV222.py',
      ROOT/'scripts/phaseEulerRegionalVolume.py',ROOT/'scripts/reviewSU2ExplicitQuality.py',
      ROOT/'outputs/v240_validation/run01/mesh.su2',ROOT/'outputs/v240_validation/run01/original_faces.csv',
      ROOT/'outputs/v240_validation/region01/provenance.npz',ROOT/'outputs/v240_validation/region01/report.json',
      ROOT/'outputs/v241_validation/buffer01/cavity.npz']
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;size=guard.size_bytes
    guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
    snapshot=BASE/('inputs_'+stage);snapshot.mkdir(exist_ok=False)
    for p in paths:
        if p.suffix in ('.py','.md','.cpp','.hpp'):(snapshot/(guard.sha(p)+p.suffix)).write_bytes(p.read_bytes())
    return guard.supervise(stage,[sys.executable,str(ROOT/script)]+args,seconds,256,paths,BASE,kind)


if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
