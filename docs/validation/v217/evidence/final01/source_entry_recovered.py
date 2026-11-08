"""模块职责：监督v217离线设计和最小观察；各新增原生阶段须有独立事前合同。"""
from pathlib import Path
import json,os,sys
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v217_validation'

def main(stage):
    BASE.mkdir(exist_ok=True)
    if os.name=='nt':
        import runCompleteRegionDecisionV216 as prior
        prior.BASE=BASE
        entries={'source':'inspectRegionalClosureV217.py'}
        assert stage in entries
        source=ROOT/'scripts'/entries[stage]
        locks=[source,Path(__file__).resolve(),ROOT/'docs/validation/v217/PLAN.md',ROOT/'scripts/runCompleteRegionDecisionV216.py',
               ROOT/'outputs/v214_validation/build01/report.json',ROOT/'outputs/v216_validation/inventory.json',
               ROOT/'outputs/v216_validation/final01/decision.json']
        used=sum(json.loads(p.read_text())['elapsedSeconds'] for p in BASE.glob('offline_*/report.json'))
        assert used<120
        r=prior.windows_run([sys.executable,str(source)],BASE/('offline_'+stage),min(60,120-used),locks)
        print(json.dumps({k:r[k] for k in ['passed','elapsedSeconds','peakObservedGroupRSSBytes','nativeExitCode']}))
        return r['passed']
    raise ValueError('尚无原生阶段准入')

if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
