"""模块用途：复用已保存小场验证稀疏历史拒绝与旧入口；不重算或修改任何原任务。"""
import argparse
import copy
import csv
import io
import os
from pathlib import Path
import re
import sys
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from su2NACASpatialInterface import inventory
from su2ThermalRuntime import read,sha,resolve
from validateIdealAir import save


def check(out):
    """模块一：所有变体只改独占夹具的history和内存回执，原场逐字保持。"""
    started=time.monotonic();out=Path(out).resolve();out.relative_to(ROOT/'outputs')
    out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'outputs/v171_validation/spatial_scalar01'
    legacy=ROOT/'outputs/v169_validation/small3d03'
    receipt=read(source/'engine_process.json')
    assert receipt['status']=='terminal' and receipt['exitCode']==0
    assert all(receipt['inputsStable'].values()) and receipt['checkpoints']['normalScheduleComplete']
    frames=[source/row[name] for row in receipt['checkpoints']['complete'] for name in ('volume','surface')]
    expected_steps=[row['localIteration'] for row in receipt['checkpoints']['complete']]
    assert expected_steps==[0,2,3]
    originals=[source/name for name in ('engine_process.json','case_manifest.json','mesh.su2','history.csv')]+frames
    originals += [legacy/name for name in ('case.cfg','case_manifest.json','ideal_air_contract.json','mesh.su2','geometry.json','initial_flow.csv')]
    production=[ROOT/'scripts'/name for name in ('su2NACASpatialInterface.py','prepareNACASpatialCase.py',
        'su2NACALaminarInterface.py','su2ThermalRuntime.py','su2IdealAirInterface.py','su2IdealAirCheckpoints.py')]
    locks={str(path):sha(path) for path in originals+production+[Path(__file__).resolve()]}
    save(out/'input_manifest.json',locks)
    for path in production+[Path(__file__).resolve()]:
        (out/('executed_'+path.name)).write_bytes(path.read_bytes())
    with (source/'history.csv').open(encoding='utf8') as f:
        raw=list(csv.reader(f,skipinitialspace=True))
    header=raw[0];rows=raw[1:];step_column=[name.strip() for name in header].index('Inner_Iter')
    assert [int(float(row[step_column])) for row in rows]==[0,1,2,3]
    results=[]

    def fixture(name,steps=(0,2,3),changed_rows=None,abnormal=False):
        folder=out/name;folder.mkdir(exist_ok=False)
        for file in ('case_manifest.json','mesh.su2'):
            os.link(str(source/file),str(folder/file))
        for file in frames:
            step=int(re.search(r'_(\d+)\.',file.name)[1])
            if step in steps:os.link(str(file),str(folder/file.name))
        # history会被测试改变，必须复制并单独写，绝不建立指向原history的硬链接。
        with (folder/'history.csv').open('w',encoding='utf8',newline='') as f:
            writer=csv.writer(f);writer.writerow(header);writer.writerows(rows if changed_rows is None else changed_rows)
        local=copy.deepcopy(receipt)
        for file in frames:
            step=int(re.search(r'_(\d+)\.',file.name)[1])
            if step not in steps:local['outputSHA256'].pop(file.name,None)
        local['outputSHA256']['history.csv']=sha(folder/'history.csv')
        if abnormal:local.update(exitCode=1,timedOut=False,cancelled=False)
        save(folder/'fixture_receipt.json',local)
        return folder,local

    def status_test(name,check_result,**kwargs):
        folder,local=fixture(name,**kwargs)
        result=inventory(folder,local);check_result(result)
        save(folder/'inventory.json',result)
        results.append(dict(name=name,passed=True,normalScheduleComplete=result['normalScheduleComplete'],
            expectedLocalIterations=result['expectedLocalIterations'],missing=result['missingScheduledLocalIterations'],
            historyContiguous=result['historyContiguous'],lastFiniteHistoryLocalIteration=result['lastFiniteHistoryLocalIteration'],
            lastCompleteLocalIteration=result['lastComplete']['localIteration'] if result['lastComplete'] else None))

    def require(predicate,message):
        assert predicate,message

    status_test('normal-complete',lambda r:require(r['normalScheduleComplete'] and r['expectedLocalIterations']==[0,2,3],
        '完整正常小场被误拒绝'))
    status_test('abnormal-last-history-without-pair',lambda r:require(not r['normalTerminal'] and
        r['expectedLocalIterations']==[0,2] and not r['missingScheduledLocalIterations'] and
        r['lastFiniteHistoryLocalIteration']==3 and r['lastComplete']['localIteration']==2 and
        not r['abnormalInstantStateAvailable'],'异常瞬间不应强制完整末场'),steps=(0,2),abnormal=True)
    status_test('normal-missing-final',lambda r:require(not r['normalScheduleComplete'] and
        not r['normalFinalPairPresent'] and r['missingScheduledLocalIterations']==[3],
        '正常缺末场必须拒绝'),steps=(0,2))
    status_test('normal-missing-periodic',lambda r:require(not r['normalScheduleComplete'] and
        r['normalFinalPairPresent'] and r['missingScheduledLocalIterations']==[2],
        '有末场不能掩盖缺周期场'),steps=(0,3))
    status_test('missing-history-row',lambda r:require(not r['normalScheduleComplete'] and
        not r['historyContiguous'] and [x['localIteration'] for x in r['rejected']]==[2,3],
        '缺中间历史不能接受后续保存场'),changed_rows=[row for row in rows if int(float(row[step_column]))!=1])
    for name,wrong_step in [('negative-history-step',-1),('over-budget-history-step',4)]:
        wrong=rows[-1].copy();wrong[step_column]=str(wrong_step)
        changed=[wrong]+rows if wrong_step<0 else rows+[wrong]
        folder,local=fixture(name,changed_rows=changed)
        try:inventory(folder,local)
        except AssertionError as ex:
            record=dict(name=name,passed=True,raised='AssertionError',reason=str(ex),wrongStep=wrong_step)
            save(folder/'rejection.json',record);results.append(record)
        else:raise AssertionError('未拒绝历史步 '+str(wrong_step))

    # 模块二：直接执行当前resolve，记录旧全合同是否真的通过，不给旧来源拒绝改名。
    text=(legacy/'case.cfg').read_text(encoding='utf8')
    old=dict(case=str(legacy),fullHistoricalCaseRequalified=False)
    try:
        _,_,flavor,profile,legacy_locks=resolve(legacy,text)
    except (AssertionError,KeyError,OSError,ValueError,RuntimeError) as ex:
        old.update(resolveAccepted=False,rejection=type(ex).__name__+': '+str(ex),
            traceback=traceback.format_exc(),reasonScope='current-source/runtime/input contract rejection; not a physical result')
    else:
        old.update(resolveAccepted=True,runtimeFlavor=flavor,thermalProfileReturned=profile is not None,
            validatedLocks=len(legacy_locks),fullInputContractReadOnlyValidated=True,
            fullHistoricalCaseRequalified=False)
        assert flavor=='rah-ideal-air-observed-v1' and profile is None
    altered=text.replace('CONV_RESIDUAL_MINVAL= -8','CONV_RESIDUAL_MINVAL= -7')
    assert altered!=text
    try:resolve(legacy,altered)
    except (AssertionError,KeyError,OSError,ValueError,RuntimeError) as ex:
        old.update(relaxedResidualRejected=True,relaxedResidualReason=type(ex).__name__+': '+str(ex),
                   relaxedResidualTraceback=traceback.format_exc(),configFileModified=False)
    else:raise AssertionError('旧入口放行了改变原−8的配置文本')
    save(out/'legacy_contract.json',old)
    stable={path:sha(Path(path))==h for path,h in locks.items()};assert all(stable.values())
    report=dict(passed=True,tests=results,legacy=old,inputSHA256=locks,inputsStable=stable,
        sourceTasksModified=False,newCFDTasks=0,newPhysicalQualification=False,elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',report)
    print(dict(passed=True,historyCases=len(results),legacyResolveAccepted=old['resolveAccepted'],
               relaxedResidualRejected=old['relaxedResidualRejected'],elapsedSeconds=report['elapsedSeconds']))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    check(parser.parse_args().output)
