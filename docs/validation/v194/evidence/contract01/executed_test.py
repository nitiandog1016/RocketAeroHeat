# -*- coding: utf-8 -*-
"""真实Euler准备合同拒绝测试；只读现有输入，不启动原生任务。"""
from pathlib import Path
import argparse
import copy
import json
import shutil
import sys
import time
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepareHybridRocketEulerNative import verify_design, sha


def check(source, out):
    """模块一：正向读取实文件；仅替换manifest读取来验证拒绝规则，明确测试范围。"""
    start=time.monotonic(); source,out=map(Path,(source,out)); out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py')
    manifest=source/'manifest.json'; original=verify_design(source); original_read=Path.read_text
    changes=[('schema','wrong'),('profile','coarse-restart-400-to-2000'),('sourceCommit','wrong'),
             ('executionAuthorized',True),('nativeConfig.max_step',400),('condition.Mach',3.),
             ('condition.alphaDeg',0.),('condition.temperatureK',300.),('geometry.cells',705490),
             ('geometry.wallTriangles',19751),('geometry.wallAreaM2',1.)]
    rejected=[]
    for dotted,value in changes:
        modified=copy.deepcopy(original); cursor=modified
        keys=dotted.split('.')
        for key in keys[:-1]:cursor=cursor[key]
        cursor[keys[-1]]=value
        def substituted(path,*args,**kwargs):
            return json.dumps(modified) if path.resolve()==manifest.resolve() else original_read(path,*args,**kwargs)
        with patch.object(Path,'read_text',substituted):
            try:verify_design(source)
            except ValueError:rejected.append(dotted)
            else:raise AssertionError('应拒绝但通过: '+dotted)
    # 来源SHA和设计SHA走真实校验路径，仅制造错误预期摘要；文件本身完全不改。
    for key in ['filesSHA256','sourceSHA256']:
        modified=copy.deepcopy(original); first=next(iter(modified[key])); modified[key][first]='0'*64
        def substituted(path,*args,**kwargs):
            return json.dumps(modified) if path.resolve()==manifest.resolve() else original_read(path,*args,**kwargs)
        with patch.object(Path,'read_text',substituted):
            try:verify_design(source)
            except ValueError:rejected.append(key)
            else:raise AssertionError('错误SHA未拒绝')
    result=dict(passed=True,positiveRealDesign=True,rejected=rejected,
        scope='manifest substitution for contract rejection; actual file SHA reads; no native rerun or full reader fault injection',
        sourceSHA256={str(p.resolve()):sha(p) for p in [manifest,Path(__file__),Path(__file__).parents[1]/'scripts/prepareHybridRocketEulerNative.py']},
        physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(passed=True,rejections=len(rejected),elapsedSeconds=result['elapsedSeconds'])))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();check(a.source,a.out)
