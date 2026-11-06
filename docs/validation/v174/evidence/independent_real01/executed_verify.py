"""v174实际候选终态独立核对：保存身份/停止语义/cvT反算，不重做全热审核。"""
import argparse
import csv
from datetime import datetime
import hashlib
import io
import json
from pathlib import Path
import shutil
import time

import numpy as np
import psutil

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'outputs/v170_validation/mesh02'
REFERENCE_FILES={
    'mesh.su2':'5562d7e6a3985ec6877585db9984d273a409178e54b1fb1d336065c93cd77152',
    'mesh_arrays.npz':'4e82f33ae66d0a19bf1a755d6c0d9210ca5cadd8db977ab5377da438517e3c01',
    'geometry.json':'62cb540c4a79d214434121a9f44ad9537de2258d006bcd808cd92ce814ba966b',
    'wall_temperature.csv':'d9b4e8916609fcb64ce4e95c10badd3793975da614d8b0f696bec3004d5ee660'}


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):result.update(block)
    return result.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def history(path):
    raw=Path(path).read_bytes();assert raw.endswith(b'\n')
    rows=csv.reader(io.StringIO(raw.decode('ascii')),skipinitialspace=True)
    keys=[v.strip() for v in next(rows)];values=[]
    assert len(keys)==len(set(keys)) and 'Inner_Iter' in keys
    for row in rows:
        if not row:continue
        v=np.asarray(row,float)
        assert len(v)==len(keys) and np.isfinite(v).all()
        assert v[keys.index('Inner_Iter')]==len(values)
        values.append(v)
    return keys,np.array(values)


def visible_owned_processes(folder,receipt):
    """只读查当前账户可见命令行，不发送信号、不按进程名终止。"""
    active=[];inaccessible=0
    own_folder=str(folder.resolve()).replace('\\','/').lower()
    for process in psutil.process_iter(['pid','cmdline','create_time']):
        try:
            args=process.info['cmdline'] or []
            normalized=[str(v).replace('\\','/').lower() for v in args]
            if any(v.endswith('/nacaconstantcpdomainguard.py') for v in normalized):
                if '--worker' in normalized:
                    index=normalized.index('--worker')
                    if index+1<len(normalized) and normalized[index+1].startswith(own_folder+'/'):
                        active.append(dict(pid=process.pid,command=args))
        except (psutil.NoSuchProcess,psutil.ZombieProcess):continue
        except psutil.AccessDenied:inaccessible+=1
    native=None
    try:
        process=psutil.Process(receipt['processID'])
        started=datetime.fromisoformat(receipt['startedUTC']).timestamp()
        ended=datetime.fromisoformat(receipt['endedUTC']).timestamp()
        if started-2<=process.create_time()<=ended+2:
            native=dict(pid=process.pid,command=process.cmdline(),created=process.create_time())
    except psutil.NoSuchProcess:pass
    return dict(activeGuardHelpers=active,originalNativeStillRunning=native,
        inaccessibleProcessCommandLineCount=inaccessible,
        scope='read-only current visible process command lines and recorded native PID/create-time; no signals')


def verify(case,output):
    start=time.monotonic();case=Path(case).resolve();output=Path(output).resolve()
    case.relative_to(ROOT/'outputs');output.relative_to(ROOT/'outputs')
    assert not output.exists();output.mkdir(parents=True)
    shutil.copy2(str(Path(__file__)),str(output/'executed_verify.py'))
    report=dict(passed=False,physicalQualified=False,experimentalQualified=False,newQualifiedArea_m2=0,
        newCFDRuns=0,newBuilds=0,fullThermalAuditRepeated=False,checks=[],frames=[])
    locked={};signatures={}
    def lock(path):
        path=Path(path).resolve();value=sha(path);locked[str(path)]=value
        stat=path.stat();signatures[str(path)]=[stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns]
        return value
    def check(name,condition,**evidence):
        assert condition,name
        report['checks'].append(dict(name=name,passed=True,**evidence))
        assert time.monotonic()-start<45,'独立终态核对45秒预算耗尽'
    try:
        lock(Path(__file__))
        receipt_path=case/'engine_process.json';lock(receipt_path);receipt=read(receipt_path)
        check('native-terminal-receipt',receipt['status']=='terminal')
        manifest_path=case/'case_manifest.json';lock(manifest_path);manifest=read(manifest_path)
        contract=manifest['idealGas'];guard=receipt['savedFieldDomainGuard']
        check('actual-fresh-contract',manifest['caseKind']=='naca-tn3344-constant-cp-spatial-wall-v1'
            and manifest['specification']['variant']=='validation-profile'
            and manifest['specification']['actualValidationMesh'] is True
            and manifest['nodes']==190334 and manifest['elements']==179568
            and contract['selectedFluidModel']=='IDEAL_GAS' and not contract['turbulent']
            and contract['energyReference']=='e=cv*T; h=cp*T; no-OTH-offset'
            and contract['developmentTemperatureRangeK']==[150.,360.])
        check('runtime-input-locks-reported-stable',bool(receipt['inputsStable']) and all(receipt['inputsStable'].values()))
        proof=ROOT/'docs/validation/v170/evidence_manifest.json';proof_hash=lock(proof);proof_value=read(proof)
        check('original-v170-manifest',proof_hash=='975599db91d573a26e480f82bf5a1214fcd4d0ee40a3ef46f693c38c91c9ca0e')
        identified={Path(v['path']).name:v['sha256'] for v in proof_value['localOnly']
                    if Path(v['path']).parent.as_posix()=='outputs/v170_validation/mesh02'}
        identified.update({Path(v['source']).name:v['sha256'] for v in proof_value['files']
                    if Path(v.get('source','')).parent.as_posix()=='outputs/v170_validation/mesh02'})
        for name,wanted in REFERENCE_FILES.items():
            check('original-four-file-'+name,identified[name]==wanted==lock(REFERENCE/name)==lock(case/name)
                ==manifest['actualValidationMesh']['fileSHA256'][name])
        init=manifest['initialization'];initial_hash=lock(case/'initial_flow.csv')
        check('fresh-prior-zero-no-old-state',init['priorSuccessfulUpdates']==0
              and init['reusesPreviousConservativeState'] is False
              and not manifest['actualValidationMesh']['previousFailedConservativeStateReused']
              and init['initialStateSHA256']==initial_hash
              and receipt['inputSHA256'][str(case/'initial_flow.csv')]==initial_hash,
              independentlyRegeneratedInitialField=False)
        check('guard-explicitly-enabled',guard['contract']==contract['savedFieldDomainGuard']
            and guard['contract']['temperatureRangeK']==[150.,360.]
            and guard['contract']['historyStrictlyAfterForLive'] is True
            and guard['contract']['helperMaximumSeconds']==15.
            and guard['contract']['totalHelperMaximumSeconds']==60.
            and guard['observesUnsavedInstantaneousState'] is False)
        check('guard-owned-helper-finished',not guard['active'] and not guard['cleanupErrors'])
        check('guard-total-budget',0<=guard['helperElapsedSeconds']<=60.)
        lock(case/'history.csv');keys,h=history(case/'history.csv');last=len(h)-1 if len(h) else None
        inventory=receipt['checkpoints'];complete=inventory['complete']
        check('independent-history-accounting',last==inventory['lastFiniteHistoryLocalIteration']
              and len(h)==inventory['historyRows'] and inventory['historyContiguous'])
        frames={row['localIteration']:row for row in complete};checked={row['step']:row for row in guard['checkedFrames']}
        check('checked-frames-are-final-complete-pairs',set(checked).issubset(frames))
        # 每一体场只读取一次，SHA与数值使用同一字节；不再调用原守卫/热审算子。
        cv=(8.31451/.02896546)/(1.4-1.)
        check('independent-native-cv',contract['gasConstant']==8.31451/.02896546 and contract['gamma']==1.4)
        for step,row in sorted(frames.items()):
            vp=case/row['volume'];sp=case/row['surface'];raw=vp.read_bytes()
            volume_hash=hashlib.sha256(raw).hexdigest();surface_hash=lock(sp)
            locked[str(vp)]=volume_hash;stat=vp.stat();signatures[str(vp)]=[stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns]
            header=next(csv.reader(io.StringIO(raw.split(b'\n',1)[0].decode('ascii'))))
            names=['PointID','Density','Momentum_x','Momentum_y','Momentum_z','Energy','Temperature','x','y','z']
            columns=[header.index(name) for name in names]
            a=np.loadtxt(io.BytesIO(raw),delimiter=',',skiprows=1,usecols=columns,ndmin=2)
            del raw
            check('saved-pair-identity-%s'%step,row['volumeSHA256']==volume_hash==receipt['outputSHA256'][vp.name]
                and row['surfaceSHA256']==surface_hash==receipt['outputSHA256'][sp.name]
                and a.shape==(190334,10) and np.isfinite(a).all() and np.array_equal(a[:,0],np.arange(190334)))
            rho=a[:,1];t=a[:,6];assert min(rho)>0 and min(t)>0
            inferred=(a[:,5]/rho-.5*np.sum((a[:,2:5]/rho[:,None])**2,axis=1))/cv
            error=float(np.max(abs(inferred-t)/np.maximum(abs(t),1.)))
            temperature_range=[float(t.min()),float(t.max())]
            inferred_range=[float(inferred.min()),float(inferred.max())]
            outside=(t<150.)|(t>360.);derived_outside=(inferred<150.)|(inferred>360.)
            check('independent-cvT-%s'%step,error<2e-12
                and np.array_equal(outside,derived_outside),relativeError=error)
            frame=dict(localIteration=step,successfulUpdates=step+1,temperatureRangeK=temperature_range,
                inferredTemperatureRangeK=inferred_range,temperatureRelativeError=error,
                outOfRangeCount=int(np.count_nonzero(outside)),volumeSHA256=volume_hash,surfaceSHA256=surface_hash,
                modelDomainPassed=not bool(np.any(outside)),minimumTemperaturePointID=int(np.argmin(t)),
                maximumTemperaturePointID=int(np.argmax(t)),observedByGuard=step in checked,
                abnormalInstantState=False,physicalQualified=False)
            if step==100:
                # 只定位已保存100；从冻结的原法向列取标签，不重建近壁或猜异常瞬间。
                meta=read(case/'geometry.json')
                with np.load(str(case/'mesh_arrays.npz')) as arrays:labels=arrays['nodeIndices'].copy()
                assert labels.shape==(len(a),3)
                h0=manifest['publicBenchmark']['conditions']['totalEnthalpy_Jkg']
                total=1.4*cv*t+.5*np.sum((a[:,2:5]/rho[:,None])**2,axis=1)
                def location(index):
                    j,i,k=[int(v) for v in labels[index]]
                    return dict(pointID=int(index),xyzM=a[index,7:10].tolist(),
                        columnLabels=dict(normalLayer=j,streamwiseStation=i,azimuthIndex=k),
                        geometryOffsetM=float(meta['offsets'][j]),geometrySM=float(meta['s'][i]),
                        offsetIsOriginalGeometryCoordinateNotSolvedBoundaryLayerThickness=True)
                frame['saved100Locations']=dict(minimumTemperature=location(np.argmin(t)),
                    maximumTemperature=location(np.argmax(t)),maximumTotalEnthalpy=location(np.argmax(total)))
                frame['saved100TotalEnthalpy']=dict(freestreamH0Jkg=h0,
                    minimumJkg=float(total.min()),maximumJkg=float(total.max()),
                    maximumRatioToFreestreamH0=float(np.max(total/h0)),
                    maximumRelativeExcessToFreestreamH0=float(np.max(total/h0-1)),
                    identity='cp*T + |momentum/rho|^2/2; native no-OTH-offset caloric model',
                    wallHeatingPresent=True,uniqueFailureCauseEstablished=False)
            if step in checked:
                observed=checked[step]
                check('guard-record-matches-final-frame-%s'%step,observed['volumeSHA256']==volume_hash
                    and observed['surfaceSHA256']==surface_hash and observed['temperatureRangeK']==temperature_range
                    and observed['outOfRangeCount']==frame['outOfRangeCount']
                    and observed['exceeded']==(not frame['modelDomainPassed'])
                    and observed['noWriterHandlesAcquired'] is True)
            report['frames'].append(frame)
        successful_attempts=[]
        for attempt in guard['helperAttempts']:
            if attempt['status']!='complete':continue
            request_path=case/attempt['request'];result_path=case/attempt['result']
            request_hash=lock(request_path);result_hash=lock(result_path)
            request=read(request_path);result=read(result_path)
            check('helper-result-binds-request-%s'%attempt['step'],result['requestSHA256']==request_hash
                and result['token']==request['token'] and request['step']==attempt['step']
                and result['step']==attempt['step'] and result==checked[attempt['step']]
                and result_hash==receipt['outputSHA256'][result_path.name]
                and request_hash==receipt['outputSHA256'][request_path.name]
                and 0<=attempt['elapsedSeconds']<=15.)
            check('helper-latest-frame-safety-%s'%attempt['step'],request['nativeTerminalAtSubmission']
                or result['historyMaximumAtInspection']>attempt['step'])
            successful_attempts.append(dict(step=attempt['step'],nativeTerminalAtSubmission=request['nativeTerminalAtSubmission'],
                historyMaximumAtInspection=result['historyMaximumAtInspection'],helperSeconds=attempt['elapsedSeconds']))
        report['helperChecks']=successful_attempts
        exceeded=bool(guard['exceeded']);requested=bool(receipt['domainGuardTerminationRequested'])
        terminal_only=bool(guard['terminalOnlyDetection']);error=guard['error']
        check('domain-outcome-identity',receipt['modelDomainExceeded']==exceeded)
        if exceeded:
            first=guard['firstExceededFrame'];step=first['step']
            check('first-observed-exceeded-frame',first==checked[step] and first['outOfRangeCount']>0
                and receipt['stopReason']==guard['stopReason']=='model-domain-exceeded'
                and all(not row['exceeded'] for s,row in checked.items() if s<step))
            check('online-terminal-stop-distinction',not (terminal_only and requested))
            classification='online-stop-requested' if requested else ('terminal-only-detection' if terminal_only
                else 'live-poll-detection-native-already-exited')
        else:
            check('no-false-domain-trigger',not requested and guard['firstExceededFrame'] is None
                  and all(not row['exceeded'] for row in checked.values()))
            classification='guard-error' if error else 'no-observed-domain-exceedance'
        if error:check('guard-error-not-hidden',bool(receipt.get('supervisorError')))
        cancel_path=case/'engine_cancel.json';matching_cancel=False
        if cancel_path.is_file():
            lock(cancel_path);cancel=read(cancel_path)
            matching_cancel=bool(receipt['request'].get('jobToken')
                and cancel.get('jobToken')==receipt['request']['jobToken'])
        check('user-cancellation-distinct',not receipt['cancelled'] or matching_cancel)
        if not cancel_path.is_file():check('no-user-cancel-without-request',receipt['cancelled'] is False)
        process_state=visible_owned_processes(case,receipt)
        check('no-visible-owned-helper-or-native-residue',not process_state['activeGuardHelpers']
              and process_state['originalNativeStillRunning'] is None)
        report['processReview']=process_state
        check('qualification-remains-zero',manifest['physicalQualified'] is False
            and receipt['physicalSolutionQualified'] is False and not inventory['physicsQualified']
            and contract['newQualifiedArea_m2']==0)
        normal=receipt['exitCode']==0 and not any(receipt.get(k,False) for k in
            ('timedOut','cancelled','outputBudgetExceeded','supervisorError','modelDomainExceeded','guardError'))
        check('normal-terminal-not-inferred-from-exit-alone',inventory['normalTerminal']==normal)
        if normal:
            check('normal-run-complete-schedule-and-guard-coverage',inventory['normalScheduleComplete']
                and inventory['normalFinalPairPresent'] and set(frames)==set(checked) and not error)
        last_saved=max(frames) if frames else None
        check('abnormal-instant-state-not-invented',inventory['abnormalInstantStateAvailable'] is False)
        report['termination']=dict(classification=classification,nativeExitCode=receipt['exitCode'],
            modelDomainExceeded=exceeded,domainGuardTerminationRequested=requested,terminalOnlyDetection=terminal_only,
            cancelled=receipt['cancelled'],matchingUserCancelRequest=matching_cancel,
            guardError=error,supervisorError=receipt.get('supervisorError'),lastFiniteHistoryLocalIteration=last,
            successfulUpdatesAtLastFiniteHistory=last+1 if last is not None else 0,
            lastCompleteSavedLocalIteration=last_saved,successfulUpdatesAtLastSaved=last_saved+1 if last_saved is not None else 0,
            firstExceededSavedLocalIteration=guard['firstExceededFrame']['step'] if exceeded else None,
            abnormalInstantStateAvailable=False,unsavedIntermediateDomainKnown=False,
            helperElapsedSeconds=guard['helperElapsedSeconds'],normalTerminal=normal,
            nativeTerminationRequestDoesNotIdentifyFirstInstantaneousExcursion=True)
        report['passed']=True
    except Exception as error:
        report['error']=type(error).__name__+': '+str(error);raise
    finally:
        report['elapsedSeconds']=time.monotonic()-start
        # 避免重复遍历大场：每项原始SHA在首次完整读取取得，末尾核对文件身份未变。
        stable={}
        for name,wanted in signatures.items():
            stat=Path(name).stat();stable[name]=[stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns]==wanted
        report['inputSHA256']=locked;report['inputsStableMetadata']=stable
        report['stabilityMeaning']='single full-read SHA bound to terminal receipts; unchanged size/mtime/ctime at review end'
        report['passed']=report['passed'] and all(stable.values())
        (output/'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
        print(json.dumps(dict(passed=report['passed'],checks=len(report['checks']),frames=len(report['frames']),
            elapsedSeconds=report['elapsedSeconds'],error=report.get('error')),ensure_ascii=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();verify(args.case,args.output)
