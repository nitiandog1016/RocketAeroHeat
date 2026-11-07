"""真实exit0内存失败须拒绝；正常阶段小例与旧封存回执只读兼容。"""
from pathlib import Path
import argparse,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from auditEilmerCompletion import audit_completion,completion_log_for_receipt
from readEilmerNative import export_loads


def test(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    small=ROOT/'outputs/v180_validation/small01';real=ROOT/'outputs/v180_validation/real01'
    good=audit_completion(small,small.parent/'small01_run/process00/stdout.txt');assert good['finalStep']==8
    # 精确官方停止名称的合成日志检查；不声称小例真的达到残差门槛。
    stop_names=[]
    original_log=(small.parent/'small01_run/process00/stdout.txt').read_text()
    for name in ['absolute-global-residual-target','relative-global-residual-target']:
        log=out/(name+'.txt');log.write_text(original_log.replace('STOP-REASON: maximum-steps','STOP-REASON: '+name))
        assert audit_completion(small,log)['stopReason']==name;stop_names.append(name)
    rejected=[]
    for name,case,log in [('real-exit0-oom',real,real.parent/'real01_run/process00/stdout.txt')]:
        try:audit_completion(case,log)
        except ValueError as e:rejected.append(dict(name=name,error=str(e)))
        else:raise AssertionError(name)
    broken=out/'missing-final.txt';broken.write_text('STOP-REASON: maximum-steps\n')
    try:audit_completion(small,broken)
    except ValueError as e:rejected.append(dict(name='missing-final',error=str(e)))
    else:raise AssertionError('missing-final')
    old=ROOT/'outputs/v178_validation/native03';receipt=ROOT/'outputs/v178_validation/interface01/observed_run.json'
    legacy=json.loads(receipt.read_text());log=completion_log_for_receipt(old,receipt,legacy)
    oldcheck=audit_completion(old,log);assert oldcheck['finalStep']==4
    q=export_loads(old,old/'lmrsim/loads/0000/blk-0000-bndry-0-wall.dat',out/'legacy-wall.csv',receipt)
    assert q['wallFaces']==7 and q['newQualifiedAreaM2']==0
    # 失败回执仍保持原始passed=true，当前导出必须在读取载荷前拒绝它。
    try:export_loads(real,real/'lmrsim/loads/0000/blk-0000-bndry-0-wall.dat',out/'must-not-exist.csv',real.parent/'real01_run/report.json')
    except ValueError as e:rejected.append(dict(name='reject-real-export',error=str(e)))
    else:raise AssertionError('real-export')
    assert not (out/'must-not-exist.csv').exists()
    r=dict(passed=True,normalNativeStep=8,legacyNativeStep=4,rejected=rejected,syntheticAcceptedStopNames=stop_names,oldFailedReceiptPreserved=True,
        newNativeRuns=0,physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');return r


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);print(json.dumps(test(p.parse_args().out)))
