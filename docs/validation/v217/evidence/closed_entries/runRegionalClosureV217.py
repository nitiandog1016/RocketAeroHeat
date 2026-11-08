"""模块职责：监督v217离线设计和最小观察；各新增原生阶段须有独立事前合同。"""
from pathlib import Path
import json,os,sys,shutil,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v217_validation'

def main(stage):
    BASE.mkdir(exist_ok=True)
    if os.name=='nt':
        import runCompleteRegionDecisionV216 as prior
        prior.BASE=BASE
        entries={'source':'inspectRegionalClosureV217.py','prepare':'prepareRegionalObservationV217.py','audit':'auditRegionalStencilV217.py'}
        assert stage in entries
        source=ROOT/'scripts'/entries[stage]
        locks=[source,Path(__file__).resolve(),ROOT/'docs/validation/v217/PLAN.md',ROOT/'scripts/runCompleteRegionDecisionV216.py',
               ROOT/'outputs/v214_validation/build01/report.json',ROOT/'outputs/v216_validation/inventory.json',
               ROOT/'outputs/v216_validation/final01/decision.json']
        if stage=='prepare':
            locks += [ROOT/'docs/validation/v217/OBSERVATION_PLAN.md',ROOT/'outputs/v214_validation/prepared01/nk_probe.d']
            locks += [p for p in (ROOT/'outputs/v214_validation/control01').rglob('*') if p.is_file()]
        if stage=='audit':
            locks += [BASE/'control01/stencil.csv',ROOT/'outputs/v214_validation/control01/nk_samples.csv',
                      ROOT/'outputs/v214_validation/control01/directions.csv',ROOT/'data/thermo/phase_aware_dry_air_v1.json']
        used=sum(json.loads(p.read_text())['elapsedSeconds'] for p in BASE.glob('offline_*/report.json'))
        assert used<120
        r=prior.windows_run([sys.executable,str(source)],BASE/('offline_'+stage),min(60,120-used),locks)
        print(json.dumps({k:r[k] for k in ['passed','elapsedSeconds','peakObservedGroupRSSBytes','nativeExitCode']}))
        return r['passed']
    from runEilmerProcess import run,sha
    package=ROOT/'tools/runtime/eilmer-stencil-observation-v1'
    common=[ROOT/'docs/validation/v217/PLAN.md',ROOT/'docs/validation/v217/OBSERVATION_PLAN.md',
            Path(__file__).resolve(),ROOT/'scripts/prepareRegionalObservationV217.py']
    archive=BASE/('entry_'+stage);archive.mkdir(exist_ok=False)
    shutil.copyfile(__file__,archive/'executed_entry.py')
    shutil.copyfile(common[-1],archive/'executed_preparation.py')
    start=time.monotonic()
    if stage=='build':
        assert json.loads((BASE/'offline_prepare/report.json').read_text())['passed']
        package.mkdir(exist_ok=False)
        old=json.loads((ROOT/'outputs/v214_validation/build01/report.json').read_text())
        oldpkg=ROOT/'tools/runtime/eilmer-frozen-nk-probe-v1'
        shutil.copyfile(oldpkg/'buildinfo.json',package/'buildinfo.json')
        oldmain=ROOT/'outputs/v214_validation/prepared01/nk_probe.d';new=BASE/'prepared01/stencil_probe.d'
        cmd=[str(new) if p==str(oldmain) else p.replace(str(oldpkg),str(package)) for p in old['argv']]
        locks=[]
        for p,h in old['inputSHA256'].items():
            if Path(p).suffix not in ('.d','.di','.a','.json'):continue
            assert sha(p)==h;locks.append(p)
        locks += common+[new,package/'buildinfo.json']
        r=run(cmd,package,BASE/'build01',150,512,locks,6144,6144)
        (BASE/'build_result.json').write_text(json.dumps(dict(passed=r['passed'],binarySHA256=sha(package/'probe') if r['passed'] else None))+'\n')
    elif stage=='query':
        b=json.loads((BASE/'build_result.json').read_text());assert b['passed']
        binary=package/'probe';assert sha(binary)==b['binarySHA256'];case=BASE/'control01'
        locks=common+[binary]+[p for p in case.rglob('*') if p.is_file()]
        r=run(['taskset','-c',str(min(os.sched_getaffinity(0))),str(binary)],case,BASE/'query01',30,32,locks,2048,2048)
    else:raise ValueError('stage未登记')
    summary=dict(passed=r['passed'],elapsedSeconds=r['elapsedSeconds'],entrySeconds=time.monotonic()-start,
                 nativeExitCode=r['nativeExitCode'],peakObservedGroupRSSBytes=r['peakObservedGroupRSSBytes'])
    (archive/'report.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary));return r['passed']

if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
