"""监督器默认分支及新保存场拒绝真实进程组清点；两个极小系统过程，无CFD。"""
from pathlib import Path
import argparse
import json
import shutil
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from runEilmerProcess import run,sha


def verify(out):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    shutil.copyfile(__file__,out/'executed_test.py')
    ok=run([sys.executable,'-c','print("default-no-guard")'],out,out/'default',5,16)
    assert ok['passed'] and 'checkpointGuardStop' not in ok
    count=[0]
    def reject():
        count[0]+=1
        return dict(snapshot=1,invalidCells=1,fixture=True) if count[0]>=2 else None
    bad=run([sys.executable,'-c','import subprocess,time; subprocess.Popen(["sleep","10"]); time.sleep(10)'],
        out,out/'guard',5,16,checkpoint_guard=reject)
    assert not bad['passed'] and bad['stopReason']=='checkpoint-guard-rejected' and not bad['liveGroupAfterCleanup']
    assert bad['checkpointGuardStop']['fixture'] and bad['nativeExitCode']!=0
    r=dict(passed=True,defaultBranchPassed=True,guardProcessGroupReaped=True,nativeCFDRuns=0,
           elapsedSeconds=time.monotonic()-start,runnerSHA256=sha(Path(__file__).parents[1]/'scripts/runEilmerProcess.py'))
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);verify(p.parse_args().out)
