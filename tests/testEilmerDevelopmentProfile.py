"""具名发展配置与旧入口兼容检查；错误阶段/保存合同不得进入原生长任务。"""
from pathlib import Path
import argparse,copy,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareEilmerRocketCase import validate_condition,make_job


def test(out):
    start=time.monotonic();out=Path(out);assert not out.exists()
    source=ROOT/'outputs/v178_validation/native03';m=json.loads((source/'case_manifest.json').read_text())
    c=validate_condition(m['condition']);assert make_job(m['mesh']['nativeSourcePath'],c)==(source/'job.lua').read_text()
    c['maxNewtonSteps']=8;c['development']=dict(profile='rah-nk-development-v1',initialOrderSteps=2,loadsEvery=2,snapshotsEvery=4,snapshotCount=2)
    validate_condition(c);rejected=[]
    for field,value in [('profile','unknown'),('initialOrderSteps',8),('initialOrderSteps',True),('loadsEvery',0),('snapshotCount',4)]:
        bad=copy.deepcopy(c);bad['development'][field]=value
        try:validate_condition(bad)
        except ValueError:rejected.append(field+':'+str(value))
        else:raise AssertionError('未拒绝'+field)
    oldprefix=(source/'job.lua').read_text().split('NewtonKrylovGlobalConfig')[0]
    assert make_job(m['mesh']['nativeSourcePath'],c).split('NewtonKrylovGlobalConfig')[0]==oldprefix
    out.parent.mkdir(parents=True,exist_ok=True)
    (out.parent/'small-condition.json').write_text(json.dumps(c,indent=2)+'\n')
    result=dict(passed=True,legacyJobBitwiseTextIdentical=True,physicalLuaPrefixIdentical=True,rejected=rejected,
        elapsedSeconds=time.monotonic()-start,nativeCFDRuns=0)
    out.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);print(json.dumps(test(p.parse_args().out)))
