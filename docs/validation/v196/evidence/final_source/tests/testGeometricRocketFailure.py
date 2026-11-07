"""只读验证真实新失败被拒绝，并在任何原生启动前关闭本候选。"""
from pathlib import Path
import json
import shutil
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from geometricRocketEuler import execute,BASE
from hybridRocketEuler import SavedFieldGuard,validate_terminal,read_wall

def main():
    start=time.monotonic();out=BASE/'closure01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py')
    case=BASE/'rocket01';prep=BASE/'preparation01/report.json';rp=BASE/'run01/report.json'
    actions=[execute,lambda:validate_terminal(case,prep,rp),lambda:read_wall(case,1,prep,rp)]
    for action in actions:
        try:action()
        except ValueError:pass
        else:raise AssertionError('关闭候选/异常研究导出未拒绝')
    bad=SavedFieldGuard(case,64982,out/'guard.json',maximum_dt=.001)()
    assert bad['snapshot']==1 and bad['invalidCells']==1351 and bad['temperatureRangeK'][0]<0
    process=json.loads((rp.parent/'native/report.json').read_text())
    assert process['nativeExitCode']==-15 and process['stopReason']=='checkpoint-guard-rejected'
    assert process['inputsStable'] and not process['liveGroupAfterCleanup']
    report=dict(passed=True,rejections=4,nativeRuns=0,actualGuardStopConfirmed=True,
                elapsedSeconds=time.monotonic()-start,physicalQualification=False)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
