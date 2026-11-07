# -*- coding: utf-8 -*-
"""新Euler合同的原生准备入口；不复用Wilcox准入，不启动发展计算。"""
from pathlib import Path
import argparse
import json
import shutil
import time
from runEilmerProcess import run, sha, COMMIT


def prepare(source, case, record):
    """模块一：逐字绑定设计输入和既有准确运行包，禁止覆盖旧算例。"""
    started = time.monotonic()
    source, case, record = map(lambda p: Path(p).resolve(), (source, case, record))
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((source/'manifest.json').read_text())
    assert manifest['schema'] == 'RAH.HybridEulerPreparation.v1'
    assert manifest['sourceCommit'] == COMMIT and manifest['geometry']['cells'] == 2400
    assert not case.exists() and not record.exists()
    for name, digest in manifest['filesSHA256'].items():
        assert sha(source/name) == digest, '设计输入已变化: '+name
    runtime_path = root/'tools/runtime/eilmer-explicit-workspace-v1/runtime.json'
    runtime = json.loads(runtime_path.read_text())
    assert runtime['sourceCommit'] == COMMIT
    for name, digest in runtime['filesSHA256'].items():
        assert sha(name) == digest, '运行包变化: '+name
    case.mkdir(parents=True); record.mkdir(parents=True)
    for name in list(manifest['filesSHA256'])+['manifest.json']:
        shutil.copyfile(source/name, case/name)
    shutil.copyfile(__file__, record/'executed_entry.py')
    locks = [case/n for n in list(manifest['filesSHA256'])+['manifest.json']]
    locks += [Path(__file__), runtime_path]+list(runtime['filesSHA256'])
    result = dict(schema='RAH.HybridEulerNativePreparation.v1', case=str(case),
                  sourceCommit=COMMIT, passed=False, physicalQualification=False,
                  runtimeSHA256=sha(runtime_path), commands=[], executionAuthorized=False)
    # 模块二：三项准备共享一个有限预算；包装退出0须结合产物审核。
    try:
        for i, command in enumerate([['prep-gas','gas.inp','rah-air.gas'],
                                     ['lmr','prep-grids','--job=grid.lua'],
                                     ['lmr','prep-sim','--job=job.lua']]):
            left = 60-(time.monotonic()-started)
            assert left > 0, '准备总预算耗尽'
            receipt = run(command, case, record/('process%02d'%i), left,
                          128, locks, memory_mib=4096, group_rss_mib=4096)
            result['commands'].append(receipt)
            assert receipt['passed'], '原生准备失败'
        for name in ['rah-air.gas','native-input-state.json','lmrsim/config','lmrsim/control',
                     'lmrsim/grid/grid-0000','lmrsim/snapshots/fluid.metadata',
                     'lmrsim/snapshots/0000/fluid-0000']:
            assert (case/name).is_file() and (case/name).stat().st_size > 0, '缺少准备产物: '+name
        result['nativeInputSHA256'] = {str(p):sha(p) for p in case.rglob('*') if p.is_file()}
        result['passed'] = True
        result['scope'] = 'native preparation only; geometry/config/EOS audit pending'
    except Exception as error:
        result['error'] = repr(error)
    result['elapsedSeconds'] = time.monotonic()-started
    (record/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True)
    p.add_argument('--case',required=True);p.add_argument('--record',required=True)
    a=p.parse_args();r=prepare(a.source,a.case,a.record)
    print(json.dumps(dict(passed=r['passed'],elapsedSeconds=r['elapsedSeconds'],error=r.get('error'))))
    raise SystemExit(0 if r['passed'] else 1)
