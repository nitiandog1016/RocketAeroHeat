"""回放唯一真实首步保存场，验证新增伪时间拒绝；不重复CFD。"""
from pathlib import Path
import argparse
import json
import shutil
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from hybridRocketEuler import SavedFieldGuard,execute,sha,save

def verify(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_verify.py');root=Path(__file__).resolve().parents[1]
    case=root/'outputs/v195_validation/rocket01'
    guard=SavedFieldGuard(case,64982,out/'guard.json',maximum_dt=.001);rejection=guard()
    assert rejection['kind']=='invalid-saved-pseudo-time' and rejection['snapshot']==1 and rejection['step']==1
    assert rejection['dt']>1e90 and len(guard.rows)==2 and all(row['invalidCells']==0 for row in guard.rows)
    try:execute(root/'outputs/v195_validation/admission.json',out/'must-not-run')
    except ValueError as e:assert '3e98' in str(e)
    else:raise AssertionError('关闭候选仍被放行')
    assert not (out/'must-not-run').exists()
    result=dict(passed=True,actualSavedFrameTimeRejected=True,validThermalFieldsStillReported=True,
        closedCandidateRejectedBeforeLaunch=True,nativeCFDRuns=0,
        scope='post-run saved-frame guard added after manual stop; no retroactive runtime claim',
        elapsedSeconds=time.monotonic()-start,sourceSHA256={str(p):sha(p) for p in [Path(__file__),root/'scripts/hybridRocketEuler.py']})
    save(out/'report.json',result);print(json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);verify(p.parse_args().out)
