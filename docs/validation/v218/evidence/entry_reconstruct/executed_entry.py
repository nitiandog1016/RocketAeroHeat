"""模块职责：监督完整离线参考的事前预算；当前不含原生/构建/网格/CFD入口。"""
from pathlib import Path
import json,sys,shutil
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v218_validation'

def main(stage):
    BASE.mkdir(exist_ok=True)
    import runCompleteRegionDecisionV216 as supervisor
    supervisor.BASE=BASE
    entries={'reconstruct':'inspectIndependentReferenceV218.py'}
    assert stage in entries
    source=ROOT/'scripts'/entries[stage]
    locks=[source,Path(__file__).resolve(),ROOT/'docs/validation/v218/PLAN.md',
        ROOT/'scripts/runCompleteRegionDecisionV216.py',ROOT/'outputs/v214_validation/control01/directions.csv',
        ROOT/'outputs/v217_validation/control01/stencil.csv',ROOT/'tools/runtime/eilmer-nk-geometric-v1/overlay/src/lmr/fluidfvcell.d']
    versions=BASE/('entry_'+stage);versions.mkdir(exist_ok=False);shutil.copyfile(__file__,versions/'executed_entry.py')
    used=sum(json.loads(p.read_text())['elapsedSeconds'] for p in BASE.glob('offline_*/report.json'))
    assert used<120
    r=supervisor.windows_run([sys.executable,str(source)],BASE/('offline_'+stage),min(60,120-used),locks)
    print(json.dumps({k:r[k] for k in ('passed','elapsedSeconds','peakObservedGroupRSSBytes','nativeExitCode')}));return r['passed']

if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
