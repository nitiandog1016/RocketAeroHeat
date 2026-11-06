"""模块用途：有限只读核对v173异常终态的小日志、4条事件和准确终止处理源；不读取大场。"""
import csv
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'outputs/v173_validation/real01'
BUILD=ROOT/'outputs/v173_validation/build01'
PARENT=ROOT/'outputs/v171_validation/build01'
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source'
PACKAGE=ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-spatial-step-v1'
OUT=ROOT/'outputs/v173_validation/termination_review01'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text('utf8'))


def save(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')


def main():
    """模块一：绑定真实进程回执和源码身份；模块二：扫描全部日志字节，不从退出码推断原因。"""
    started=time.monotonic();OUT.mkdir(parents=True,exist_ok=False)
    (OUT/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    locks={}
    def bind(path,expected=None):
        digest=sha(path)
        if expected is not None:assert digest==expected,str(path)
        locks[str(path.resolve())]=digest
        return digest
    receipt_path=CASE/'engine_process.json';runtime_path=PACKAGE/'runtime.json'
    receipt=read(receipt_path);runtime=read(runtime_path)
    for path in (Path(__file__),receipt_path,runtime_path,BUILD/'source_tree.json',PARENT/'source_tree.json',
                 BUILD/'objects_after.json',PARENT/'objects_after.json',
                 ROOT/'outputs/v173_validation/native_real01_process.json'):
        bind(path)
    assert receipt['status']=='terminal' and receipt['exitCode']==134
    assert not any(receipt.get(k,False) for k in ('timedOut','cancelled','supervisorError','outputBudgetExceeded'))
    assert all(receipt['inputsStable'].values()) and len(receipt['inputsStable'])==85
    assert receipt['runtimeFlavor']==runtime['flavor']=='rah-ideal-air-spatial-step-v1'
    bind(runtime_path,receipt['inputSHA256'][str(runtime_path)])
    executable=PACKAGE/'bin/SU2_CFD.exe'
    bind(executable,next(row['sha256'] for row in runtime['files'] if row['path']=='bin/SU2_CFD.exe'))
    assert locks[str(executable)]==receipt['inputSHA256'][str(executable)]
    tree=read(BUILD/'source_tree.json');old_tree=read(PARENT/'source_tree.json')
    assert bind(BUILD/'source_tree.json')==runtime['sourceTreeSHA256']
    source_records={};texts={}
    names=('SU2_CFD/src/SU2_CFD.cpp','SU2_CFD/include/NativeTermination.hpp',
           'SU2_CFD/include/variables/ConservativeStepTrace.hpp','SU2_CFD/include/fluid/PropertyQueryContext.hpp')
    for relative in names:
        path=BUILD/'executed_sources'/relative;digest=bind(path,tree[str(SOURCE/relative)])
        inherited=old_tree.get(str(SOURCE/relative))==digest
        source_records[relative]=dict(path=str(path),sha256=digest,unchangedFromSpatialWallParent=inherited)
        texts[relative]=path.read_text('utf8')
    handler=texts[names[1]];main_source=texts[names[0]]
    assert source_records[names[0]]['unchangedFromSpatialWallParent']
    assert source_records[names[1]]['unchangedFromSpatialWallParent']
    assert 'rah_termination::install();' in main_source
    assert 'previous_handler = std::set_terminate(record_and_terminate);' in handler
    assert handler.index('std::_Exit(134)')<handler.index('RAH_NATIVE_TERMINATION_BEGIN')
    assert handler.count('std::fflush(stderr)')==1
    assert handler.index('RAH_NATIVE_TERMINATION_END')<handler.index('std::fflush(stderr)')
    assert 'std::fflush(trace.file)' in texts[names[2]]
    objects=read(BUILD/'objects_after.json');old_objects=read(PARENT/'objects_after.json')
    main_objects=[name for name in objects if name.replace('\\','/').endswith('/SU2_CFD.cpp.obj')]
    assert len(main_objects)==1 and objects[main_objects[0]]==old_objects[main_objects[0]]
    runner=ROOT/'scripts/runEmbeddedSU2.py'
    bind(runner,receipt['inputSHA256'][str(runner)])
    runner_text=runner.read_text('utf8')
    assert "with (folder / 'engine_stdout.log').open('wb') as stream:" in runner_text
    assert 'stdout=stream, stderr=subprocess.STDOUT' in runner_text

    # 全字节扫描包含NUL和非文本尾部；仅限制文件≤1MiB，不读取保存体场/网格。
    candidates=[p for p in CASE.iterdir() if p.is_file() and
        (p.suffix.lower() in ('.log','.txt') or any(s in p.name.lower() for s in ('stderr','stdout','termination')))]
    assert [p.name for p in candidates]==['engine_stdout.log']
    logs=[]
    for path in candidates+[ROOT/'outputs/v173_validation/native_real01.log']:
        assert path.stat().st_size<1024**2
        bind(path,receipt['outputSHA256'][path.name] if path.parent==CASE else None)
        raw=path.read_bytes()
        tokens={word:raw.count(word.encode('ascii')) for word in ('RAH_NATIVE_TERMINATION','RAH_TERMINATION',
            'RAH_PROPERTY_QUERY','exception_type=','domain_error','inverse target outside temperature range')}
        assert not any(tokens.values()),tokens
        logs.append(dict(path=str(path),bytes=len(raw),nulBytes=raw.count(b'\0'),
            nonASCIIBytes=sum(byte>=128 for byte in raw),markerCounts=tokens,
            last256BytesHex=raw[-256:].hex(),last256BytesText=raw[-256:].decode('utf8',errors='replace'),
            endsWithNewline=raw.endswith(b'\n')))
    assert logs[0]['bytes']==12791
    event_path=CASE/'conservative_steps.csv';bind(event_path,receipt['outputSHA256'][event_path.name])
    with event_path.open(encoding='utf8',newline='') as stream:events=list(csv.DictReader(stream))
    assert len(events)==4 and [int(row['event']) for row in events]==[1,2,3,4]
    observed=[]
    for row in events:
        assert row['old_valid']=='1' and float(row['theta'])==.5 and row['attempts']=='1'
        assert 'mode=rhoe phase=node.trial' in row['reason']
        observed.append({key:row[key] for key in ('event','local_point','global_point','inner_iter','output',
            'old_valid','theta','attempts','trial_e','old_e','accepted_e','reason')})
    stable={path:sha(Path(path))==digest for path,digest in locks.items()};assert all(stable.values())
    report=dict(passed=True,scope='small-log-and-recorded-event-termination-evidence-only',
        runtimeFlavor=runtime['flavor'],executableSHA256=locks[str(executable)],
        sourceTreeSHA256=runtime['sourceTreeSHA256'],exactSources=source_records,
        inheritedMainObject=dict(path=main_objects[0],sha256=objects[main_objects[0]],parentSHA256=old_objects[main_objects[0]]),
        logByteReviews=logs,caseDirectoryAlternativeTextLogFound=False,stderrMergedIntoStdout=True,
        terminationMessagePresent=False,terminationActualCauseKnown=False,terminalExceptionTypeKnown=False,
        terminalQueryModeKnown=False,terminalNodeKnown=False,terminalTemperatureDirectionKnown=False,
        operatingSystemEventLogsReviewed=False,recordedRecoveredNodeEvents=observed,
        recoveredEventsIdentifyFinalFailure=False,
        handlerLimitations=dict(cppTerminateOnly=True,concurrentOrRecursiveEntryCanExit134BeforeWrite=True,
            concurrentOrRecursiveEntryProven=False,stderrFlushOnlyAfterFullExceptionAndStackText=True,
            allTerminationPathsGuaranteedToFlush=False,exit134IdentifiesDomainError=False,
            nativeTraceExplicitlyFlushesEachEvent=True),
        nativeExitCode=receipt['exitCode'],nativeElapsedSeconds=receipt['elapsedSeconds'],
        generatedOutputBytes=receipt['generatedOutputBytes'],largeFieldsRead=0,CFDLaunches=0,builds=0,
        inputSHA256=locks,inputsStable=stable,elapsedSeconds=time.monotonic()-started)
    save(OUT/'report.json',report)
    print(json.dumps(dict(passed=True,elapsedSeconds=report['elapsedSeconds'],logs=len(logs),events=len(events),
        terminationCauseKnown=False,largeFieldsRead=0),ensure_ascii=False))


if __name__=='__main__':main()
