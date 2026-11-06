"""模块用途：只读验证已终态自适应冷重启任务，分开完整保存帧、正常完成和物理资格。"""
import argparse
import csv
from pathlib import Path
import sys
import time
import traceback

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from prepareNACAAdaptiveContinuation import KIND,INITIALIZATION,FLAVOR,ADAPTIVE,INITIAL_COLUMNS,specification
from prepareSU2NACALaminar import STAGE
from nacaSpatialWallFlux import read_temperature_table,temperature_at,verify_applied
from su2IdealAirInterface import settings
from su2ThermalRuntime import read,sha
from validateIdealAir import save


def history_values(path):
    """模块一：读取真实有限历史，保留异常末尾坏行，不补造局部CFL或松弛因子。"""
    with path.open(encoding='utf8') as stream:raw=list(csv.reader(stream,skipinitialspace=True))
    assert raw,'历史文件为空'
    keys=[key.strip() for key in raw[0]];required=['Inner_Iter','Min CFL','Max CFL','Avg CFL','LinSolRes',
        'Linear_Solver_Iterations','Min DT','Max DT','rms[Rho]','rms[RhoE]']
    assert set(required).issubset(keys) and len(keys)==len(set(keys))
    records=[];rejected=[]
    for number,row in enumerate(raw[1:],2):
        try:
            values=np.asarray(row,float)
            assert len(values)==len(keys) and np.isfinite(values).all()
            item=dict(zip(keys,values));assert item['Inner_Iter']==int(item['Inner_Iter'])
            records.append(item)
        except (AssertionError,ValueError,OverflowError) as exc:
            rejected.append(dict(line=number,reason=type(exc).__name__+': '+str(exc)))
    steps=[int(row['Inner_Iter']) for row in records]
    assert steps==list(range(len(steps))),'有限历史缺段、重复或编号异常'
    return records,rejected


def verify(folder,audit_folder,out):
    """模块二：核对来源和终态；模块三：双帧/热审核；模块四：实际CFL历史。"""
    started=time.monotonic();folder=Path(folder).resolve();audit_folder=Path(audit_folder).resolve()
    out=Path(out).resolve();out.relative_to(ROOT/'outputs');out.mkdir(parents=True,exist_ok=False)
    locks={};report=dict(reviewCompleted=False,passed=False,physicalQualified=False,experimentalQualification=False,
        newQualifiedArea_m2=0,newCFDTasks=0,uninterruptedTrajectoryEquivalent=False)

    def bind(path,expected=None):
        path=Path(path).resolve();digest=locks.get(str(path))
        if digest is None:digest=sha(path);locks[str(path)]=digest
        if expected is not None:assert digest==expected,'文件身份不符: '+str(path)
        return digest

    def mapped_hash(mapping,path):
        path=Path(path).resolve();matches=[value for key,value in mapping.items() if Path(key).resolve()==path]
        assert matches and len(set(matches))==1,'缺精确路径绑定或绑定冲突: '+str(path)
        return matches[0]

    def stable_record(item):
        stable=item.get('inputsStable');inputs=item.get('inputSHA256')
        assert isinstance(stable,dict) and stable and all(value is True for value in stable.values())
        assert isinstance(inputs,dict) and set(inputs)==set(stable)
        return inputs

    try:
        sources=[Path(__file__).resolve()]+[ROOT/'scripts'/name for name in (
            'prepareNACAAdaptiveContinuation.py','su2NACAAdaptiveInterface.py','prepareSU2NACALaminar.py',
            'nacaSpatialWallFlux.py','nacaWallFlux.py','su2IdealAirInterface.py','su2ThermalRuntime.py','validateIdealAir.py')]
        for path in sources:
            bind(path);(out/('executed_'+path.name)).write_bytes(path.read_bytes())
        receipt_path=folder/'engine_process.json';manifest_path=folder/'case_manifest.json'
        bind(receipt_path);bind(manifest_path)
        receipt=read(receipt_path);manifest=read(manifest_path)
        assert receipt.get('status')=='terminal','禁止读取未终态任务'
        assert type(receipt.get('exitCode')) is int
        receipt_inputs=stable_record(receipt)
        assert not receipt.get('supervisorError'),'监督完整性错误不得作为有效保存帧证据'
        flags={key:receipt.get(key,False) for key in (
            'timedOut','cancelled','outputBudgetExceeded','memoryBudgetExceeded','supervisorError')}
        assert all(isinstance(value,bool) for value in flags.values())
        normal=receipt['exitCode']==0 and not any(flags.values())
        report.update(normalRunCompleted=normal,actualExitCode=receipt['exitCode'],terminalFlags=flags)
        assert manifest['caseKind']==KIND and manifest['requiredRuntimeFlavor']==receipt['runtimeFlavor']==FLAVOR
        spec=specification(manifest['specification']['variant']);assert spec==manifest['specification']
        assert manifest['advancement']==ADAPTIVE and not manifest['physicalQualified']
        if not spec['actualValidationMesh']:assert normal,'小型合同检查必须正常结束'
        for path,digest in receipt_inputs.items():bind(path,digest)
        bind(manifest_path,mapped_hash(receipt_inputs,manifest_path))

        init_path=folder/'adaptive_initialization.json';init=read(init_path);contract=read(folder/'ideal_air_contract.json')
        bind(init_path,contract['initializationSHA256']);bind(folder/'ideal_air_contract.json',manifest['idealAirContractSHA256'])
        assert init==manifest['initialization'] and contract==manifest['idealAir']
        assert init['kind']==INITIALIZATION and init['localIterationStartsAt']==0
        assert init['sourceLocalIteration']==spec['sourceLocalIteration']
        prior=init['priorSuccessfulUpdates'];assert type(prior) is int and prior==spec['sourceLocalIteration']+1
        assert init['adaptiveCFLReset'] and init['solverHistoryReset']
        assert init['samePhysicalModel'] and init['sameSpatialDiscretization']
        assert not init['samePhysicsAndNumerics'] and not init['uninterruptedTrajectoryEquivalent']
        assert contract['temperatureRangeK']==[60,2000] and contract['energyEntropyReference']=='OTH'
        assert not contract['turbulent'] and not contract['physicsQualified'] and contract['newQualifiedArea']==0
        for path,digest in init['sourceSHA256'].items():bind(path,digest)
        source=Path(init['sourceDirectory']).resolve();assert source==(ROOT/spec['sourceDirectory']).resolve()
        source_receipt_path=source/'engine_process.json';source_audit_path=(ROOT/spec['sourceAudit']).resolve()
        assert Path(init['sourceReceiptPath']).resolve()==source_receipt_path
        assert Path(init['sourceAuditPath']).resolve()==source_audit_path
        bind(source_receipt_path,init['sourceReceiptSHA256']);bind(source_audit_path,init['sourceAuditSHA256'])
        sr=read(source_receipt_path);sa=read(source_audit_path)
        assert sr['status']=='terminal' and sr['runtimeFlavor']==FLAVOR
        source_frame=next(row for row in sr['checkpoints']['complete'] if row['localIteration']==spec['sourceLocalIteration'])
        assert source_frame['successfulUpdates']==source_frame['seriesSuccessfulUpdates']==prior
        assert source_frame['stage']==init['sourceStage']==STAGE
        source_volume=source/source_frame['volume'];source_surface=source/source_frame['surface']
        assert Path(init['sourceStatePath']).resolve()==source_volume and Path(init['sourceSurfacePath']).resolve()==source_surface
        bind(source_volume,init['sourceStateSHA256']);bind(source_surface,init['sourceSurfaceSHA256'])
        assert init['sourceStateSHA256']==source_frame['volumeSHA256']==sr['outputSHA256'][source_volume.name]
        assert init['sourceSurfaceSHA256']==source_frame['surfaceSHA256']==sr['outputSHA256'][source_surface.name]
        assert mapped_hash(sa['inputSHA256'],source_volume)==init['sourceStateSHA256']
        source_thermal=next(row for row in sa['frames'] if row['localIteration']==spec['sourceLocalIteration'])
        assert source_thermal['interfacePassed']
        assert init['sourceLastFiniteHistoryLocalIteration']==sr['checkpoints']['lastFiniteHistoryLocalIteration']
        assert init['sourceTermination']['status']==sr['status'] and init['sourceTermination']['exitCode']==sr['exitCode']
        assert init['sourceTermination']['timedOut']==sr.get('timedOut',False)

        runtime_path=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)/'runtime.json'
        bind(runtime_path,init['sourceRuntimeManifestSHA256']);runtime=read(runtime_path)
        executable=Path(runtime['executable']).resolve();executable.relative_to(runtime_path.parent)
        bind(executable,init['sourceExecutableSHA256'])
        for path in (runtime_path,executable):
            assert bind(path)==mapped_hash(receipt_inputs,path)==mapped_hash(sr['inputSHA256'],path)
        assert runtime['flavor']==FLAVOR and init['sourceRuntimeFlavor']==FLAVOR

        # 初场只含原全局节点与五个守恒量；无需再解EOS或重复全场热审核。
        initial_path=folder/'initial_flow.csv';bind(initial_path,init['initialStateSHA256'])
        with initial_path.open(encoding='utf8') as stream:assert next(csv.reader(stream))==INITIAL_COLUMNS
        with source_volume.open(encoding='utf8') as stream:assert next(csv.reader(stream))[:9]==INITIAL_COLUMNS
        initial=np.loadtxt(str(initial_path),delimiter=',',skiprows=1)
        original=np.loadtxt(str(source_volume),delimiter=',',skiprows=1,usecols=range(9))
        assert initial.shape==original.shape==(manifest['nodes'],9) and np.isfinite(initial).all()
        assert np.array_equal(initial[:,0],np.arange(len(initial))) and np.array_equal(initial,original)
        del original
        report.update(initialNineColumnsExactlyEqual=True,conservativeColumnMaximumAbsoluteDifference=0.,
            sourceLocalIteration=spec['sourceLocalIteration'],priorSuccessfulUpdates=prior,
            sourceNormalRunCompleted=bool(sa['interfacePassed']),sourceSelectedSavedFramePassed=True,
            sourceRuntimeAndExecutableUnchanged=True)
        for name in ('mesh.su2','mesh_arrays.npz','geometry.json','wall_temperature.csv'):
            assert bind(folder/name)==bind(source/name),'物理输入改变: '+name
        text=(folder/'case.cfg').read_text(encoding='utf8');old_text=(source/'case.cfg').read_text(encoding='utf8')
        bind(folder/'case.cfg',contract['configurationSHA256']);bind(source/'case.cfg')
        cfg=settings(text);old_cfg=settings(old_text)
        changed={key:dict(old=old_cfg.get(key),new=cfg.get(key)) for key in set(cfg)|set(old_cfg) if cfg.get(key)!=old_cfg.get(key)}
        assert set(changed).issubset({'CFL_ADAPT','CFL_ADAPT_PARAM','ITER','OUTPUT_WRT_FREQ'})
        assert cfg['CFL_ADAPT']=='YES' and float(cfg['CFL_NUMBER'])==.4 and cfg['MGLEVEL']=='0'
        parameters=[float(value) for value in cfg['CFL_ADAPT_PARAM'].strip('() ').split(',')]
        assert parameters==ADAPTIVE['parameters']
        assert int(cfg['ITER'])==spec['updates'] and int(cfg['OUTPUT_WRT_FREQ'])==spec['every']
        assert cfg['LINEAR_SOLVER_ERROR']==old_cfg['LINEAR_SOLVER_ERROR'] and float(cfg['LINEAR_SOLVER_ERROR'])==1e-5
        assert cfg['LINEAR_SOLVER_ITER']==old_cfg['LINEAR_SOLVER_ITER']=='10'
        assert cfg['CONV_RESIDUAL_MINVAL']==old_cfg['CONV_RESIDUAL_MINVAL']=='-8'
        report['configurationChanges']=changed

        # 精确NPZ已与原源SHA绑定，使用其原层标识检查全部壁点，避免重建网格。
        with np.load(str(folder/'mesh_arrays.npz')) as arrays:
            assert np.array_equal(arrays['points'],initial[:,1:4])
            labels=arrays['nodeIndices'];assert labels.shape==(len(initial),3)
            wall_ids=np.flatnonzero(labels[:,0]==0)
        table=read_temperature_table(folder/'wall_temperature.csv')
        assert np.all((table[:,1]>=60)&(table[:,1]<=2000))
        given=temperature_at(initial[wall_ids,1],table)
        applied_path=folder/'rah_wall_temperature_applied.csv';bind(applied_path)
        applied=verify_applied(applied_path,initial[:,1:4],{'ids':wall_ids},given)
        report['appliedTemperatureCheck']=applied
        del initial

        history_path=folder/'history.csv';bind(history_path,receipt['outputSHA256']['history.csv'])
        history,rejected_rows=history_values(history_path);assert history,'没有有限历史可核验'
        last=int(history[-1]['Inner_Iter']);assert last<spec['updates']
        if normal:assert not rejected_rows,'正常任务存在不完整历史行'
        cfl=np.array([[row[key] for key in ('Min CFL','Max CFL','Avg CFL')] for row in history])
        # history使用有限位文本，均值并行归约允许1e−10绝对表示误差；不改变实际上下限。
        cfl_tol=1e-10
        assert np.isfinite(cfl).all() and cfl.min()>=.4-cfl_tol and cfl.max()<=100+cfl_tol
        assert np.all(cfl[:,0]<=cfl[:,1]+cfl_tol)
        assert np.all(cfl[:,2]>=cfl[:,0]-cfl_tol) and np.all(cfl[:,2]<=cfl[:,1]+cfl_tol)
        differs_from_initial=bool(np.any(abs(cfl-.4)>cfl_tol))
        varies_in_history=bool(np.any(np.ptp(cfl,axis=0)>cfl_tol))
        if not spec['actualValidationMesh']:
            assert differs_from_initial and varies_in_history,'小任务未观察到实际CFL反馈变化'
        linear=np.array([row['LinSolRes'] for row in history])
        report['actualHistory']=dict(rows=len(history),lastFiniteLocalIteration=last,
            lastFiniteSeriesSuccessfulUpdates=prior+last+1,rejectedRows=rejected_rows,
            minimumReportedCFL=float(cfl.min()),maximumReportedCFL=float(cfl.max()),
            firstCFL=dict(zip(('min','max','avg'),map(float,cfl[0]))),lastCFL=dict(zip(('min','max','avg'),map(float,cfl[-1]))),
            firstDiffersFromInitial=bool(np.any(abs(cfl[0]-.4)>cfl_tol)),anyDiffersFromInitial=differs_from_initial,
            variesAcrossHistory=varies_in_history,textAndMeanRoundoffTolerance=cfl_tol,
            actualLinearTarget=1e-5,linearTargetReachedRows=int(np.sum(linear<=-5)),
            linearLogResidualRange=[float(linear.min()),float(linear.max())],
            adaptationControlTolerance=.001,localUnderRelaxationInferred=False,unwrittenLocalCFLReplayed=False,
            reportedCFLStage='after-update value for next update; reported DT is from current completed update')

        cp=receipt['checkpoints'];pairs=cp['complete'];steps=[item['localIteration'] for item in pairs]
        assert steps==sorted(set(steps)) and all(type(step) is int and 0<=step<=last for step in steps)
        assert cp['lastFiniteHistoryLocalIteration']==last and cp['seriesUpdatesAtLastFiniteHistory']==prior+last+1
        assert cp['normalTerminal']==normal and cp['historyContiguous'] and cp['historyRows']==len(history)
        expected=list(range(0,last+1,spec['every']))
        if normal:expected=sorted(set(expected)|{last})
        missing=sorted(set(expected)-set(steps))
        assert cp['expectedLocalIterations']==expected and cp['missingScheduledLocalIterations']==missing
        if normal:
            assert steps==expected and not cp['rejected'] and cp['normalScheduleComplete'] and cp['normalFinalPairPresent']
        else:assert not cp['normalScheduleComplete'] and not cp['normalFinalPairPresent']
        for pair in pairs:
            step=pair['localIteration']
            assert pair['successfulUpdates']==step+1 and pair['seriesSuccessfulUpdates']==prior+step+1
            assert pair['stage']==STAGE and not pair['physicalQualified']
            for key,prefix,extension in (('volume','restart_flow','csv'),('surface','surface_flow','vtk')):
                expected_name=prefix+'_%06d.'%step+extension;assert pair[key]==expected_name
                assert pair[key+'SHA256']==receipt['outputSHA256'][expected_name]
                bind(folder/expected_name,pair[key+'SHA256'])

        audit_path=audit_folder/'report.json';bind(audit_path);audit=read(audit_path)
        audit_inputs=stable_record(audit)
        for path,digest in audit_inputs.items():bind(path,digest)
        assert mapped_hash(audit_inputs,receipt_path)==bind(receipt_path)
        assert mapped_hash(audit_inputs,manifest_path)==bind(manifest_path)
        assert mapped_hash(audit_inputs,applied_path)==bind(applied_path)
        assert audit['normalRunCompleted']==normal and audit['stage']==STAGE
        assert audit['savedFrameInterfacePassed'] and audit['interfacePassed']==normal
        assert audit['auditScope']==('normal-complete-saved-frames' if normal else 'complete-saved-frames-from-abnormal-task')
        assert audit['appliedTemperatureCheck']==applied
        assert audit['iterationOrigin']['priorSuccessfulUpdates']==prior
        assert not audit['experimentalQualification'] and audit['newQualifiedArea_m2']==0
        assert [row['localIteration'] for row in audit['frames']]==steps and pairs,'严格审核与完整双帧不一致'
        for pair,frame in zip(pairs,audit['frames']):
            assert mapped_hash(audit_inputs,folder/pair['volume'])==pair['volumeSHA256']
            assert frame['interfacePassed']
            for key,limit in (('energyMaximumRelativeError',2e-11),('eosMaximumRelativeError',2e-11),
                ('paperCpMaximumRelativeError',2e-8),('viscosityMaximumRelativeError',2e-11),
                ('nativeHeatMaximumScaledError',2e-10),('powerScaledDifference',2e-9)):
                value=float(frame[key]);assert np.isfinite(value) and 0<=value<limit
        report.update(reviewCompleted=True,passed=True,normalRunPassed=normal,
            savedFrameInterfacePassed=True,completeSavedLocalIterations=steps,
            completeSavedSeriesSuccessfulUpdates=[pair['seriesSuccessfulUpdates'] for pair in pairs],
            expectedLocalIterations=expected,missingScheduledLocalIterations=missing,rejectedSavedPairs=cp['rejected'],
            normalScheduleComplete=cp['normalScheduleComplete'],lastCompleteLocalIteration=steps[-1],
            lastFiniteHistoryLocalIteration=last,lastFiniteHistoryHasFullState=steps[-1]==last,
            abnormalInstantStateAvailable=False,qualificationFromCFLChange=False)
    except Exception as exc:
        report.update(passed=False,error=type(exc).__name__+': '+str(exc),traceback=traceback.format_exc())
        raise
    finally:
        stable={path:Path(path).is_file() and sha(Path(path))==value for path,value in locks.items()}
        if not all(stable.values()):report.update(passed=False,inputMutationDetected=True)
        report.update(inputSHA256=locks,inputsStable=stable,elapsedSeconds=time.monotonic()-started)
        save(out/'input_manifest.json',locks);save(out/'report.json',report)
    assert report['passed'] and all(report['inputsStable'].values())
    print(dict(passed=True,normalRunCompleted=normal,frames=steps,
        actualCFLChanged=report['actualHistory']['variesAcrossHistory'],elapsedSeconds=report['elapsedSeconds']))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('folder','audit','output'):parser.add_argument('--'+name,required=True)
    args=parser.parse_args();verify(args.folder,args.audit,args.output)
