"""模块职责：监督定向源码与既存数据核对；没有构建、原生查询、网格或CFD入口。"""
from pathlib import Path
import json,shutil,sys
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v219_validation'

def main(stage):
    BASE.mkdir(exist_ok=True)
    import runCompleteRegionDecisionV216 as supervisor
    supervisor.BASE=BASE
    names={'source':'inspectSU2RouteV219.py'};assert stage in names
    task=ROOT/'scripts'/names[stage]
    locks=[task,Path(__file__).resolve(),ROOT/'scripts/runCompleteRegionDecisionV216.py',ROOT/'docs/validation/v219/PLAN.md',
        ROOT/'outputs/v167_validation/build01/source_tree.json',ROOT/'tools/runtime/su2/8.6.0-rah-sst-energy-v1/runtime.json']
    entry=BASE/('entry_'+stage);entry.mkdir(exist_ok=False);shutil.copyfile(__file__,entry/'executed_caller.py')
    used=sum(json.loads(p.read_text(encoding='utf8'))['elapsedSeconds'] for p in BASE.glob('process_*/report.json'))
    assert used<120
    r=supervisor.windows_run([sys.executable,str(task)],BASE/('process_'+stage),min(60,120-used),locks)
    print(json.dumps({k:r[k] for k in ('passed','elapsedSeconds','peakObservedGroupRSSBytes','nativeExitCode')}));return r['passed']

if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
