"""模块职责：登记并监督v241离线生产/独立审核；不启动CFD、构建或原生查询。"""
from pathlib import Path
import sys
import runSU2RocketPhaseV222 as guard
ROOT=guard.ROOT;BASE=ROOT/'outputs/v241_validation';DOC=ROOT/'docs/validation/v241'
CLOSED=False


def main(stage):
    if CLOSED:raise SystemExit('v241 execution closed; use NEXT_DECISION')
    choices={
      'test':('tests/testTetraStreamtubesV241.py',[str(BASE/'test01/report.json')],40,'offline'),
      'test2':('tests/testTetraStreamtubesV241.py',[str(BASE/'test02/report.json')],40,'offline'),
      'produce':('scripts/produceForebodyStreamtubesV241.py',[],85,'offline'),
      'audit':('scripts/auditForebodyStreamtubesV241.py',[],60,'audit'),
      'export':('scripts/exportForebodyStreamtubesV241.py',[],45,'audit'),
    }
    script,args,seconds,kind=choices[stage]
    paths=[Path(__file__).resolve(),DOC/'PLAN.md',DOC/'STREAMTUBE_PLAN.md',ROOT/script,
      ROOT/'scripts/tetraStreamtubes.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',
      ROOT/'scripts/su2PhaseAirFields.py',ROOT/'scripts/runSU2RocketPhaseV222.py',
      ROOT/'data/thermo/phase_aware_dry_air_v1.json']
    if not stage.startswith('test'):
        paths += [ROOT/'outputs/v240_validation/run01'/n for n in ('mesh.su2','state_000200.csv','original_faces.csv')]
        paths += [ROOT/'outputs/v228_validation/prepare03/region.json',ROOT/'scripts/phaseEulerRegionalVolume.py',BASE/'test02/report.json']
    if stage in ('audit','export'):
        paths+=list((BASE/'streamtubes01').glob('*'))
    if stage=='export':
        paths += [BASE/'audit01/report.json',ROOT/'web/streamtube-research-v241.js',ROOT/'web/streamtube-research-v241.html',
          ROOT/'tests/testStreamtubeWorkbenchV241.js',ROOT/'outputs/v230_validation/audit02/geometry.npz']
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    size=guard.size_bytes;guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
    # tests/附加计划/模板也按执行字节保存；原监督器只自动快照scripts和PLAN。
    snap=BASE/('inputs_'+stage);snap.mkdir(exist_ok=False)
    for p in paths:
        if p.suffix in ('.py','.md','.js','.html'):
            target=snap/(guard.sha(p)+p.suffix);target.write_bytes(p.read_bytes())
    return guard.supervise(stage,[sys.executable,str(ROOT/script)]+args,seconds,256,paths,BASE,kind)


if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
