"""验证真正的求解子进程内存可见、组内存越限被停止且无存活子进程。"""
from pathlib import Path
import argparse,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from runEilmerProcess import run,group_members


def test(out):
    start=time.monotonic();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    child='import time; x=bytearray(40*1024**2); time.sleep(1.5)'
    parent='import subprocess,sys; p=subprocess.Popen([sys.executable,"-c",'+repr(child)+']); p.wait()'
    normal=run([sys.executable,'-c',parent],out,out/'normal',10,128,group_rss_mib=128)
    assert normal['passed'] and not normal['liveGroupAfterCleanup']
    assert len(normal['observedProcessPeakRSSBytes'])>=2
    assert normal['peakObservedGroupRSSBytes']>normal['peakObservedParentRSSBytes']+30*1024**2
    child='import time; x=bytearray(48*1024**2); time.sleep(20)'
    parent='import subprocess,sys,time; ps=[subprocess.Popen([sys.executable,"-c",'+repr(child)+']) for i in range(2)]; time.sleep(20)'
    limit=run([sys.executable,'-c',parent],out,out/'limit',10,128,group_rss_mib=80)
    assert not limit['passed'] and limit['stopReason']=='process-group-rss-limit'
    assert len(limit['observedProcessPeakRSSBytes'])==3 and not limit['liveGroupAfterCleanup']
    assert not group_members(limit['processId'])
    result=dict(passed=True,childMemoryObserved=True,groupLimitStopped=True,noLiveChildren=True,
        normalPeakGroupBytes=normal['peakObservedGroupRSSBytes'],limitedPeakGroupBytes=limit['peakObservedGroupRSSBytes'],
        elapsedSeconds=time.monotonic()-start,nativeCFDRuns=0)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);print(json.dumps(test(p.parse_args().out)))
