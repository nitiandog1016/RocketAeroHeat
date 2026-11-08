"""模块职责：监督完整离线参考的事前预算；当前不含原生/构建/网格/CFD入口。"""
from pathlib import Path
import json,sys,shutil,os,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v218_validation'

def main(stage):
    BASE.mkdir(exist_ok=True)
    used=sum(json.loads(p.read_text())['elapsedSeconds'] for p in BASE.glob('offline_*/report.json'))
    if stage=='source':
        assert os.name=='posix' and used<90
        from runEilmerProcess import run
        source=ROOT/'scripts/extractIndependentReferenceSourcesV218.py'
        locks=[source,Path(__file__).resolve(),ROOT/'docs/validation/v218/PLAN.md',
               ROOT/'outputs/v214_validation/build01/report.json']
        from extractIndependentReferenceSourcesV218 import NAMES
        locks += [Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source')/p for p in NAMES]
        entry=BASE/'entry_source';entry.mkdir(exist_ok=False);shutil.copyfile(__file__,entry/'executed_entry.py')
        shutil.copyfile(source,entry/'executed_task.py')
        work=BASE/'source_work01';work.mkdir(exist_ok=False);start=time.monotonic()
        r=run(['taskset','-c',str(min(os.sched_getaffinity(0))),sys.executable,str(source)],work,
              BASE/'offline_source',min(30,120-used),16,locks,1024,1024)
        summary=dict(passed=r['passed'],elapsedSeconds=r['elapsedSeconds'],entrySeconds=time.monotonic()-start,
                     nativeExitCode=r['nativeExitCode'],peakObservedGroupRSSBytes=r['peakObservedGroupRSSBytes'],
                     cpuLimitScope='taskset one CPU; shared supervisor upper bound two CPUs',
                     outputLimitScope='source work directory; extracted source bytes separately inventoried')
        (entry/'report.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary));return r['passed']
    import runCompleteRegionDecisionV216 as supervisor
    supervisor.BASE=BASE
    entries={'reconstruct':'inspectIndependentReferenceV218.py','assess':'assessIndependentReferenceV218.py'}
    assert stage in entries
    source=ROOT/'scripts'/entries[stage]
    locks=[source,Path(__file__).resolve(),ROOT/'docs/validation/v218/PLAN.md',
        ROOT/'scripts/runCompleteRegionDecisionV216.py',ROOT/'outputs/v214_validation/control01/directions.csv',
        ROOT/'outputs/v217_validation/control01/stencil.csv',ROOT/'tools/runtime/eilmer-nk-geometric-v1/overlay/src/lmr/fluidfvcell.d']
    if stage=='assess':
        locks += [ROOT/'outputs/v214_validation/control01/nk_samples.csv',
                  ROOT/'outputs/v214_validation/control01/lmrsim/config',ROOT/'outputs/v214_validation/control01/lmrsim/nk-config',
                  BASE/'reconstruct01/report.json',BASE/'source01/report.json']
        locks += [p for p in (BASE/'source01/source').rglob('*') if p.is_file()]
    versions=BASE/('entry_'+stage);versions.mkdir(exist_ok=False);shutil.copyfile(__file__,versions/'executed_entry.py')
    assert used<120
    r=supervisor.windows_run([sys.executable,str(source)],BASE/('offline_'+stage),min(60,120-used),locks)
    print(json.dumps({k:r[k] for k in ('passed','elapsedSeconds','peakObservedGroupRSSBytes','nativeExitCode')}));return r['passed']

if __name__ == '__main__':
    raise ValueError('v218评估已完成、冻结生产分支关闭；请按NEXT_DECISION登记新任务。')
