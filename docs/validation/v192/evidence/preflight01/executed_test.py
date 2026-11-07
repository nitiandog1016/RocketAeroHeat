# -*- coding: utf-8 -*-
"""同场发展准备来源拒绝与精确初始场检查；不启动CFD。"""
from pathlib import Path
import argparse
import json
import sys
import shutil
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from hybridEulerEvolution import audit_evolution_inputs, local_path
from readHybridEulerNative import read_wall_frame
from readEilmerNative import sha


def test(case,prep,out):
    case,prep,out=map(Path,(case,prep,out));out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py');r=json.loads(prep.read_text())
    audit,grid,mapping=audit_evolution_inputs(case,prep)
    source=local_path(r['baselineConfig']).parents[1]
    assert sha(case/'lmrsim/snapshots/0002/fluid-0000')==sha(source/'lmrsim/snapshots/0002/fluid-0000')
    assert sha(case/'lmrsim/transient-residuals.dat')==r['historyPrefixSHA256']
    rejected=[]
    changes=[('profile',lambda d:d.update(profile='original-full-domain')),
             ('case',lambda d:d.update(case=d['case']+'-wrong')),
             ('failed',lambda d:d.update(passed=False)),
             ('source',lambda d:d['sourceSHA256'].update({next(iter(d['sourceSHA256'])):'0'*64})),
             ('input',lambda d:d['nativeInputSHA256'].update({next(iter(d['nativeInputSHA256'])):'0'*64})),
             ('budget-control',lambda d:d['expectedControl'].update(max_step=2001))]
    for name,change in changes:
        d=json.loads(prep.read_text());change(d);p=out/(name+'.json');p.write_text(json.dumps(d))
        try:audit_evolution_inputs(case,p)
        except ValueError:rejected.append(name)
        else:raise AssertionError('未拒绝 '+name)
    # 新分派不破坏旧v191严格来源读取；原生计算不重复。
    old=read_wall_frame(source,2,local_path(r['sourcePreparation']),local_path(r['sourceRun']))[0]
    assert old['algebraicInterfacePassed'] and not old['physicalQualification']
    result=dict(passed=True,sourceSnapshotBitwiseEqual=True,historyPrefixEqual=True,rejections=rejected,
                previousV191ReadingPassed=True,counts={k:audit[k] for k in ['cells','vertices','faces','ghostCells']})
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['case','preparation','out']:p.add_argument('--'+k,required=True)
    a=p.parse_args();test(a.case,a.preparation,a.out)
