"""新运行包只准入已证明的唯一候选；旧异常不能通过新身份读取。"""
from pathlib import Path
import copy
import json
import shutil
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from geometricRocketEuler import check_admission,BASE
from hybridRocketEuler import execute,validate_terminal,SavedFieldGuard

def main():
    start=time.monotonic();out=BASE/'admission_test01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py')
    a=json.loads((BASE/'admission.json').read_text());check_admission(a)
    mutations=[('schema','bad'),('profile','full-wall-euler-400-v1'),('budget',{}),
        ('temperatureRangeK',[60.,6000.]),('runtimeSHA256','0'*64),('evidenceSHA256',{}),
        ('case',str(BASE/'wrong')),('preparation',str(BASE/'wrong')),('record',str(BASE/'wrong'))]
    for key,value in mutations:
        wrong=copy.deepcopy(a);wrong[key]=value
        try:check_admission(wrong)
        except ValueError:pass
        else:raise AssertionError('错误合同未拒绝: '+key)
    try:execute(BASE/'admission.json',out/'must-not-run')
    except ValueError:pass
    else:raise AssertionError('旧入口重新开放')
    old=ROOT/'outputs/v195_validation'
    try:validate_terminal(old/'rocket01',old/'preparation01/report.json',old/'run01/report.json')
    except ValueError:pass
    else:raise AssertionError('旧异常结果被接受')
    rejected=SavedFieldGuard(old/'rocket01',64982,out/'guard.json',maximum_dt=.001)()
    assert rejected['kind']=='invalid-saved-pseudo-time' and rejected['snapshot']==1
    report=dict(passed=True,rejections=12,nativeRuns=0,elapsedSeconds=time.monotonic()-start,
                oldFailureKept=True,physicalQualification=False)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
