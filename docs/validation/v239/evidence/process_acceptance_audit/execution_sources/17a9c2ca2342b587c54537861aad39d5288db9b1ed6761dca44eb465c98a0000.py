"""模块职责：在暂停/恢复后如实记录丢失的监督终态与现存资产；不重建原回执或重新启动CFD。"""
from pathlib import Path
import csv,datetime,json,os,sys
import psutil
import runMUSCLRampV239 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;CASE=run.CASE;DOC=run.DOC


def active_original_processes():
    engine=os.path.normcase(str(run.BUILD/'driver.exe'));found=[]
    for p in psutil.process_iter(['pid','name','exe','cmdline']):
        try:
            info=p.info;args=info['cmdline'] or []
            native=os.path.normcase(info['exe'] or '')==engine
            supervisor=len(args)>2 and Path(args[1]).name=='runMUSCLRampV239.py' and args[2]=='run'
            if native or supervisor:found.append({k:info[k] for k in ('pid','name','exe')})
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    return found


def status():
    """后处理只接受原回执或明确的当前不活动观察；二者的成功身份不混用。"""
    original=BASE/'process_rocket/report.json'
    if original.exists():
        r=read(original);assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup']
        return r
    r=read(BASE/'interruption01/report.json')
    assert r['passed'] and r['originalTerminalReceiptMissing'] and r['currentOriginalProcesses']==[]
    assert not active_original_processes()
    for n,h in r['originalOutputSHA256'].items():assert sha(CASE/n)==h,n
    return dict(passed=None,status='original-terminal-unknown-process-now-absent',
      inputsStable=None,inputsMatchAtObservation=r['allAdmittedInputsMatchNow'],liveGroupAfterCleanup=None,
      currentOriginalProcesses=[],originalExitCode=None,originalElapsedSeconds=None,originalPeakRSSBytes=None)


def main():
    out=BASE/'interruption01';out.mkdir(exist_ok=False)
    assert not (BASE/'process_rocket/report.json').exists()
    live=active_original_processes();assert not live,live
    admitted=read(BASE/'process_rocket/admission.json');changed=[]
    for p,h in admitted['inputSHA256'].items():
        if sha(p)!=h:changed.append(p)
    assert not changed,changed
    files=[p for p in CASE.iterdir() if p.is_file()];before={p.name:sha(p) for p in files}
    with (CASE/'steps.csv').open() as f:steps=list(csv.DictReader(f))
    with (CASE/'stages.csv').open() as f:stages=list(csv.DictReader(f))
    n=len(steps);assert [int(r['update']) for r in steps]==list(range(1,n+1))
    partial=[r for r in stages if int(r['update'])>n]
    assert not (CASE/'native_completion.json').exists() and not list(CASE.glob('restart_flow*.csv'))
    assert all(sha(p)==before[p.name] for p in files) and not active_original_processes()
    report=dict(passed=True,observedUTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
      originalTerminalReceiptMissing=True,currentOriginalProcesses=[],allAdmittedInputsMatchNow=True,
      originalInputStabilityForWholeExecutionUnknown=True,originalExitCode=None,originalElapsedSeconds=None,originalPeakRSSBytes=None,
      originalCFDBudgetConservativelyChargedSeconds=admitted['limitSeconds'],
      recordedCompletedUpdates=n,recordedPartialUpdateRows=len(partial),
      recordedPartialUpdatePhases=[dict(update=int(r['update']),system=int(r['system']),stage=r['stage']) for r in partial],
      originalOutputSHA256=before,availableFullFrames=sorted(int(p.stem.split('_')[-1]) for p in CASE.glob('state_*.bin')),
      availableWallFrames=sorted(int(p.stem.split('_')[-1]) for p in CASE.glob('wall_[0-9]*.csv')),
      exactTerminationCauseUnknown=True,userPauseResumeContext=True,numericalFailureProven=False,normalCompletionProven=False,
      fullFinalStateMissing=True,plannedBenefitCheckReached=False,rerun=False,physicalQualification=False)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ('passed','recordedCompletedUpdates','recordedPartialUpdateRows','originalTerminalReceiptMissing','availableFullFrames','availableWallFrames')}))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runMUSCLRampV239.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
           BASE/'process_rocket/admission.json',BASE/'process_rocket/live.json']+[p for p in CASE.iterdir() if p.is_file()]
        raise SystemExit(0 if run.guard.supervise('interruption',[sys.executable,str(Path(__file__).resolve()),'task'],40,1024,locks,BASE,'offline') else 1)
