"""资源必要条件的独立分组复算及真实启动前拒绝，不调用CFD或构建。"""
from pathlib import Path
import argparse
import copy
import json
import sys
import time
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'scripts'))
from eilmerResourceBudget import lower_bound, assess_budget, PACKAGE_SHA, sha
from readEilmerNative import audit_run_scope
from runPreparedEilmerCase import execute
from prepareEilmerDevelopmentCase import prepare


def verify(out):
    start = time.monotonic(); out = Path(out); out.mkdir(parents=True,exist_ok=False)
    result = dict(passed=False,nativeRuns=0,builds=0,physicalQualification=False)
    try:
        case = ROOT/'outputs/v187_validation/full/real01'
        m = json.loads((case/'case_manifest.json').read_text())
        bound = lower_bound()
        # 按所有权重新归并：不读取被测函数的行系数，分别计算单元/面/顶点和关联关系。
        # 每个C含6064字节，F含1930，V含576；N另250。连接CV/CF/FV分别32/92/24。
        # 完整推导与逐项去重解释见本轮RESOURCE_ACCOUNTING.md。
        oracle = (731814*6064 + 1741600*1930 + 335726*576 + 705490*250
                  + (70574*4+634916*6)*32 + (70574*4+634916*5)*92 + 1741600*3*24)
        assert oracle == 8743123908 == bound['minimumAllocationLowerBoundBytes']
        assert len({r['name'] for r in bound['rows']}) == len(bound['rows'])
        assert bound['fullRunPeakBytes'] is None and not bound['fitProven']
        original = assess_budget(m['mesh'],PACKAGE_SHA,8192)
        assert original['status']=='impossible' and original['runAuthorized'] is False
        # 下界与MiB边界交叉测试；更大预算只是不被下界排除，依然不准入。
        floor_mib = oracle//1024**2
        assert assess_budget(m['mesh'],PACKAGE_SHA,floor_mib)['status']=='impossible'
        for limit in (floor_mib+1,16384):
            r=assess_budget(m['mesh'],PACKAGE_SHA,limit)
            assert r['status']=='not-excluded-unproven' and not r['runAuthorized']
        rejected=[]
        def reject(label,fn):
            try: fn()
            except ValueError as error: rejected.append(dict(label=label,reason=str(error)))
            else: raise AssertionError('未拒绝: '+label)
        reject('wrong-package',lambda:assess_budget(m['mesh'],'0'*64,8192))
        for key,value in [('sourceSHA256','0'*64),('cells',705489),('nodes',335725),
                          ('wallFaces',19751),('farfieldFaces',6571),('wallAreaM2',1.9),
                          ('cellTypeCounts',{'10':705490,'12':0,'13':0,'14':0})]:
            bad=copy.deepcopy(m['mesh']);bad[key]=value
            reject('wrong-'+key,lambda bad=bad:assess_budget(bad,PACKAGE_SHA,8192))
        for invalid in (0,-1,8192.0,True,None):
            reject('invalid-budget-'+str(invalid),lambda invalid=invalid:assess_budget(m['mesh'],PACKAGE_SHA,invalid))
        # 使用真实入口，在任何原生任务/输出目录创建前拒绝，意外启动则直接让测试失败。
        with patch('subprocess.Popen',side_effect=AssertionError('resource test must never launch a process')):
            reject('real-run-rejected-before-launch',lambda:execute(case,'run',out/'must-not-run',660,1024))
            # 准备入口先检查旧大输入。本检查复用已审哈希作替身，不再读取大网格/初场。
            prior=ROOT/'outputs/v178_validation/real01_prepare/report.json'
            stored=json.loads(prior.read_text())['nativeInputSHA256']
            def existing_hash(path):
                return stored[str(path)] if str(path) in stored else sha(path)
            with patch('prepareEilmerDevelopmentCase.sha',side_effect=existing_hash):
                reject('prepare-rejected-old-hashes-substituted',lambda:prepare(
                    ROOT/'outputs/v178_validation/real01',prior,
                    ROOT/'outputs/v178_validation/real01_geometry.json',out/'must-not-case',
                    out/'must-not-prepare',m['condition'],120,1024))
            reject('closed-full-scope',lambda:audit_run_scope(m,case))
        assert not list(out.glob('must-not*'))
        for row in rejected[-3:]:
            assert str(oracle) in row['reason'] and '8589934592' in row['reason']
        # 此修改只影响完整域关闭分支；既存小例身份和原保存帧不重算。
        small=json.loads((ROOT/'outputs/v187_validation/small01/case_manifest.json').read_text())
        with patch('eilmerWorkspaceRuntime.audit_workspace_runtime') as check:
            audit_run_scope(small,ROOT/'outputs/v187_validation/small01')
            check.assert_called_once_with(small['condition'])
        result.update(passed=True,independentGroupedLowerBoundBytes=oracle,
                      originalResourceDecision=original,rejected=rejected,
                      largerBudgetStillUnauthorized=True,smallScopeDispatchPreserved=True,
                      preparePriorLargeHashesSubstituted=True,
                      sourceSHA256={n:sha(ROOT/n) for n in ['scripts/eilmerResourceBudget.py',
                          'scripts/readEilmerNative.py','tests/testEilmerResourceBudget.py']})
    except Exception as exc:
        result['error']=repr(exc);raise
    finally:
        result['elapsedSeconds']=time.monotonic()-start
        (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','independentGroupedLowerBoundBytes','elapsedSeconds')}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    verify(parser.parse_args().out)
