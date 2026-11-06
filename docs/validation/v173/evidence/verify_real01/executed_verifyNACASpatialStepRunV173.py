"""模块用途：只读验证空间Tw限步冷重启的来源、真实事件、保存帧与热审核；不启动CFD。"""
import argparse
import csv
from pathlib import Path
import re
import sys
import time
import traceback

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from prepareNACASpatialStep import (KIND,INITIALIZATION,FLAVOR,ADAPTIVE,INITIAL_COLUMNS,
    NUMERICAL_CHANGE,QUERY_CONTEXT,SPATIAL,specification)
from prepareNACAAdaptiveContinuation import FLAVOR as SOURCE_FLAVOR
from prepareSU2NACALaminar import STAGE
from nacaSpatialWallFlux import read_temperature_table,temperature_at,verify_applied
from su2IdealAirInterface import settings, CONSERVATIVE_STEP, R
from su2NACASpatialStepInterface import inventory as step_inventory
from su2IdealAirCheckpoints import output_length
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


def verify_step_events(folder,receipt,points,last,normal,out,bind):
    """模块一：独立重算事件守恒方向/热内能/EOS，真实坐标仅通过已知全局编号核对。"""
    path=folder/'conservative_steps.csv'
    if not path.exists():
        assert path.name not in receipt['outputSHA256'],'监督回执声称存在的事件文件已丢失'
        return dict(traceFilePresent=False,eventsObserved=0,actualStepEventsVerified=False,
            stepEventActivationDemonstrated=False,reason='任务没有事件文件；未实测限步事件，不推断保护不会触发',
            globalConservationProven=False)
    bind(path,receipt['outputSHA256'][path.name])
    fields=('rho','rhou','rhov','rhow','rhoE')
    required=['event','local_point','dimension','old_valid','theta','attempts','k','context_known',
        'global_point','domain','inner_iter','time_iter','mesh','rk','output','x','y','z','volume','alpha']
    required += [prefix+'_'+name for prefix in ('trial','old','accepted','linear','stored_residual') for name in fields]
    required += ['trial_e','old_e','accepted_e','reason']
    records=[];rejected=[]
    with path.open(encoding='utf8',newline='') as stream:
        reader=csv.reader(stream);assert next(reader)==required,'限步事件列合同不同'
        for line,row in enumerate(reader,2):
            if len(row)!=len(required):
                rejected.append(dict(line=line,columns=len(row),expectedColumns=len(required)))
            else:
                assert not rejected,'不完整事件后又出现完整事件，不能当作仅尾行中断'
                records.append(dict(zip(required,row)))
    if normal:assert not rejected,'正常终态含不完整限步事件'
    assert len(records)<=4096
    if not records:
        return dict(traceFilePresent=True,eventsObserved=0,rejectedRows=rejected,traceComplete=not rejected,
            actualStepEventsVerified=False,stepEventActivationDemonstrated=False,
            reason='没有完整限步事件行；不能据此验证真实限步',globalConservationProven=False)
    from idealAirProperties import IdealAir
    core=IdealAir()
    for locked,digest in core.locks.items():bind(locked,digest)
    low,high=core.state('PT',101325.,np.array([60.,2000.]))['e']
    segment_max=energy_max=location_max=query_max=0.
    linear_max=0.;linear_checked=0;known=0;old_invalid=zero=positive=unchanged_positive=0
    output_counts={};phase_counts={};global_ids=set();local_ids=set();associations={}
    accepted_rho=[];accepted_e=[];old_rho=[];old_e=[];theta_values=[];summaries=[]

    def integer(row,name):
        value=float(row[name]);assert np.isfinite(value) and value==int(value),(name,row[name])
        return int(value)

    def energy(state):
        with np.errstate(divide='ignore',invalid='ignore',over='ignore'):
            return float(state[4]/state[0]-.5*np.sum((state[1:4]/state[0])**2))

    def valid(state,e):
        return bool(np.isfinite(state).all() and state[0]>0 and np.isfinite(e) and low<=e<=high)

    def scalar_error(actual,expected):
        if np.isfinite(actual) and np.isfinite(expected):return abs(actual-expected)/max(1.,abs(expected))
        assert (np.isnan(actual) and np.isnan(expected)) or actual==expected,'非有限数类别/符号不同'
        return 0.

    for index,row in enumerate(records,1):
        assert integer(row,'event')==index
        local=integer(row,'local_point');assert 0<=local<len(points)
        assert integer(row,'dimension')==3 and float(row['k'])==0.,'必须为三维层流真实事件'
        valid_old=integer(row,'old_valid');context=integer(row,'context_known')
        assert valid_old in (0,1) and context in (0,1)
        theta=float(row['theta']);attempts=integer(row,'attempts')
        assert np.isfinite(theta) and 0<=theta<=.5
        old,trial,accepted,linear=[np.array([float(row[prefix+'_'+name]) for name in fields])
            for prefix in ('old','trial','accepted','linear')]
        energies={prefix:energy(value) for prefix,value in [('old',old),('trial',trial),('accepted',accepted)]}
        for prefix,value in energies.items():
            error=scalar_error(value,float(row[prefix+'_e']))
            assert error<2e-11,(index,prefix,error)
            energy_max=max(energy_max,error)
        assert valid(old,energies['old'])==bool(valid_old),'old_valid与独立热内能/密度判定不同'
        assert not valid(trial,energies['trial']),'记录的原试探在独立物性范围内'
        reason=row['reason']
        assert reason.startswith('RAH_PROPERTY_QUERY mode=rhoe ')
        pairs=dict(re.findall(r'\b([A-Za-z][A-Za-z0-9]*)=([^\s]+)',reason))
        expected_phase='node.trial' if valid_old else 'node.old'
        assert pairs['phase']==expected_phase and pairs['localKnown']=='1' and int(pairs['local'])==local
        # CNS查询只知道局部RCM号；真实全局号来自单独求解器上下文，不能使用query的默认0。
        assert pairs['globalKnown']==pairs['markerKnown']==pairs['innerKnown']=='0'
        query_state=trial if valid_old else old
        query_e=energies['trial'] if valid_old else energies['old']
        query_error=max(scalar_error(float(pairs['input1']),query_state[0]),
                        scalar_error(float(pairs['input2']),query_e))
        assert query_error<2e-11;query_max=max(query_max,query_error)
        phase_counts[expected_phase]=phase_counts.get(expected_phase,0)+1
        if valid_old:
            assert 1<=attempts<=32 and valid(accepted,energies['accepted'])
            if theta==0:
                assert attempts==32 and np.array_equal(accepted,old);zero+=1
                expected=old
            else:
                assert theta==2.**(-attempts)
                expected=old+theta*(trial-old);positive+=1;theta_values.append(theta)
                unchanged_positive+=int(np.array_equal(accepted,old))
            segment=float(np.max(abs(accepted-expected)/np.maximum(1.,abs(accepted))))
            assert segment<2e-13;segment_max=max(segment_max,segment)
            accepted_rho.append(accepted[0]);accepted_e.append(energies['accepted'])
            old_rho.append(old[0]);old_e.append(energies['old'])
        else:
            assert not normal and theta==0 and attempts==0
            assert np.array_equal(accepted,old,equal_nan=True);old_invalid+=1
        event=dict(event=index,localRCMPoint=local,queryPhase=expected_phase,oldValid=bool(valid_old),
            theta=theta,attempts=attempts,geometryContextKnown=bool(context))
        local_ids.add(local)
        if context:
            known+=1;global_id=integer(row,'global_point');assert 0<=global_id<len(points)
            global_ids.add(global_id)
            if local in associations:assert associations[local]==global_id
            associations[local]=global_id
            xyz=np.array([float(row[axis]) for axis in 'xyz'])
            assert np.isfinite(xyz).all()
            location_error=float(np.max(abs(points[global_id]-xyz)))
            assert location_error<1e-13;location_max=max(location_max,location_error)
            inner=integer(row,'inner_iter');output=integer(row,'output')
            assert 0<=inner<=last+1 and output in (0,1)
            assert integer(row,'mesh')==0 and integer(row,'time_iter')==0
            assert integer(row,'domain') in (0,1) and integer(row,'rk')>=0
            volume=float(row['volume']);alpha=float(row['alpha'])
            assert np.isfinite(volume) and volume>0 and np.isfinite(alpha) and 0<=alpha<=1
            output_counts[str(output)]=output_counts.get(str(output),0)+1
            if np.isfinite(old).all() and np.isfinite(trial).all() and np.isfinite(linear).all():
                increment=float(np.max(abs(trial-old-alpha*linear)/np.maximum(1.,abs(trial))))
                linear_max=max(linear_max,increment);linear_checked+=1
            event.update(globalPoint=global_id,coordinateMaximumErrorM=location_error,
                innerIteration=inner,outputPhase=bool(output),meshLevel=0,domain=bool(integer(row,'domain')))
        summaries.append(event)
    eos_error=0.;accepted_temperature=None;old_temperature=None
    if accepted_rho:
        state=core.state('rhoe',np.asarray(accepted_rho),np.asarray(accepted_e))
        previous=core.state('rhoe',np.asarray(old_rho),np.asarray(old_e))
        pressure=np.asarray(accepted_rho)*R*state['T']
        eos_error=float(np.max(abs(state['p']-pressure)/np.maximum(1.,abs(pressure))))
        assert eos_error<2e-11 and np.all((state['T']>=60)&(state['T']<=2000))
        accepted_temperature=[float(state['T'].min()),float(state['T'].max())]
        old_temperature=[float(previous['T'].min()),float(previous['T'].max())]
    save(out/'step_event_summary.json',summaries)
    return dict(traceFilePresent=True,traceComplete=not rejected,rejectedRows=rejected,eventsObserved=len(records),
        actualStepEventsVerified=True,stepEventActivationDemonstrated=True,
        validOldEvents=len(accepted_rho),invalidOldEvents=old_invalid,zeroStepEvents=zero,
        positiveThetaEvents=positive,positiveThetaWithoutRepresentableChange=unchanged_positive,
        positiveEffectiveEvents=positive-unchanged_positive,
        minimumPositiveTheta=min(theta_values) if theta_values else None,
        maximumConservativeSegmentError=segment_max,maximumThermalEnergyScaledReplayError=energy_max,
        maximumQueryInputScaledError=query_max,independentEOSChecked=bool(accepted_rho),
        independentEOSMaximumError=eos_error,acceptedTemperatureRangeK=accepted_temperature,
        oldTemperatureRangeK=old_temperature,knownGeometryEvents=known,
        coordinateMaximumErrorM=location_max,globalPointIDs=sorted(global_ids),localRCMPointIDs=sorted(local_ids),
        localAndGlobalIDsAssumedEqual=False,phaseCounts=phase_counts,outputPhaseCounts=output_counts,
        linearIncrementComparedEvents=linear_checked,maximumObservedLinearIncrementError=linear_max,
        rawLinearIncrementAgreement=bool(linear_checked and linear_max<2e-13),
        linearIncrementIsIndependentObservation=True,dualVolumeIndependentlyRecomputed=False,
        globalConservationProven=False,finalEventEqualsSavedTerminalField=False,
        physicalQualificationFromEvents=False)


def verify(folder,audit_folder,out):
    """模块二：核对来源和终态；模块三：双帧/热审核；模块四：实际CFL历史。"""
    started=time.monotonic();folder=Path(folder).resolve();audit_folder=Path(audit_folder).resolve()
    out=Path(out).resolve();out.relative_to(ROOT/'outputs');out.mkdir(parents=True,exist_ok=False)
    locks={};historical_rebinding={};report=dict(reviewCompleted=False,passed=False,physicalQualified=False,experimentalQualification=False,
        newQualifiedArea_m2=0,newCFDTasks=0,uninterruptedTrajectoryEquivalent=False)

    def bind(path,expected=None,historical_source=False):
        """模块五：当前源码照常锁定；仅历史源码或本轮计划可绑定准确SHA快照，物理文件不代换。"""
        path=Path(path).resolve();digest=locks.get(str(path))
        if digest is None:digest=sha(path);locks[str(path)]=digest
        if expected is not None and digest!=expected and historical_source:
            try:relative=path.relative_to(ROOT)
            except ValueError:relative=None
            is_source=relative is not None and relative.parts[0] in ('scripts','tests','third_party') and \
                path.suffix.lower() in ('.py','.m','.cpp','.c','.hpp','.h','.inc','.js')
            is_plan=path==(ROOT/'docs/validation/v173/PLAN.md').resolve()
            if is_source or is_plan:
                candidates=([folder/'plan_at_prepare.md'] if is_plan else [
                    ROOT/'outputs/v173_validation/build01/project_sources'/relative,
                    ROOT/'outputs/v173_validation/contract01/executed_source'/relative])
                selected=next((candidate for candidate in candidates if candidate.is_file() and sha(candidate)==expected),None)
                assert selected is not None,'缺准确历史源码快照: '+str(path)
                locks[str(selected.resolve())]=expected
                historical_rebinding[(str(path),expected)]=dict(originalPath=str(path),expectedHistoricalSHA256=expected,
                    currentSourceSHA256=digest,exactSnapshotPath=str(selected.resolve()),snapshotSHA256=expected,
                    reason='prepared-plan-snapshot' if is_plan else 'historical-source-snapshot',
                    physicalInputSubstitution=False)
                return expected
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
            'prepareNACASpatialStep.py','su2NACASpatialStepInterface.py','prepareNACAAdaptiveContinuation.py',
            'prepareSU2SpatialStepPatch.py','prepareSU2NACALaminar.py','idealAirProperties.py',
            'su2NACASpatialInterface.py','su2NACALaminarInterface.py','su2IdealAirCheckpoints.py','auditSU2SSTState.py',
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
        assert manifest['numericalChange']==NUMERICAL_CHANGE
        if not spec['actualValidationMesh']:assert normal,'小型合同检查必须正常结束'
        for path,digest in receipt_inputs.items():bind(path,digest,historical_source=True)
        bind(manifest_path,mapped_hash(receipt_inputs,manifest_path))

        init_path=folder/'step_initialization.json';init=read(init_path);contract=read(folder/'ideal_air_contract.json')
        bind(init_path,contract['initializationSHA256']);bind(folder/'ideal_air_contract.json',manifest['idealAirContractSHA256'])
        assert init==manifest['initialization'] and contract==manifest['idealAir']
        assert init['kind']==INITIALIZATION and init['localIterationStartsAt']==0
        assert init['numericalChange']==NUMERICAL_CHANGE
        assert init['sourceLocalIteration']==spec['sourceLocalIteration']
        prior=init['priorSuccessfulUpdates'];assert type(prior) is int and prior==spec['sourceLocalIteration']+1
        assert init['adaptiveCFLReset'] and init['solverHistoryReset']
        assert init['samePhysicalModel'] and init['sameSpatialDiscretization']
        assert not init['samePhysicsAndNumerics'] and not init['uninterruptedTrajectoryEquivalent']
        assert contract['temperatureRangeK']==[60,2000] and contract['energyEntropyReference']=='OTH'
        assert not contract['turbulent'] and not contract['physicsQualified'] and contract['newQualifiedArea']==0
        for path,digest in init['sourceSHA256'].items():bind(path,digest,historical_source=True)
        source=Path(init['sourceDirectory']).resolve();assert source==(ROOT/spec['sourceDirectory']).resolve()
        source_receipt_path=source/'engine_process.json';source_audit_path=(ROOT/spec['sourceAudit']).resolve()
        assert Path(init['sourceReceiptPath']).resolve()==source_receipt_path
        assert Path(init['sourceAuditPath']).resolve()==source_audit_path
        bind(source_receipt_path,init['sourceReceiptSHA256']);bind(source_audit_path,init['sourceAuditSHA256'])
        sr=read(source_receipt_path);sa=read(source_audit_path)
        assert sr['status']=='terminal' and sr['runtimeFlavor']==SOURCE_FLAVOR
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

        # 来源必须是原空间Tw包；本次使用独立限步包，明确数值改变，不要求两个exe相等。
        source_runtime_path=ROOT/'tools/runtime/su2'/('8.6.0-'+SOURCE_FLAVOR)/'runtime.json'
        bind(source_runtime_path,init['sourceRuntimeManifestSHA256']);source_runtime=read(source_runtime_path)
        source_executable=Path(source_runtime['executable']).resolve()
        source_executable.relative_to(source_runtime_path.parent)
        bind(source_executable,init['sourceExecutableSHA256'])
        for path in (source_runtime_path,source_executable):
            assert bind(path)==mapped_hash(sr['inputSHA256'],path)
        assert source_runtime['flavor']==init['sourceRuntimeFlavor']==SOURCE_FLAVOR
        runtime_path=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)/'runtime.json'
        bind(runtime_path,mapped_hash(receipt_inputs,runtime_path));runtime=read(runtime_path)
        executable=Path(runtime['executable']).resolve();executable.relative_to(runtime_path.parent)
        bind(executable,mapped_hash(receipt_inputs,executable))
        assert runtime['flavor']==runtime['numericalBaselineFlavor']==FLAVOR
        assert runtime['numericalParentFlavor']==SOURCE_FLAVOR and runtime['spatialWall']==SPATIAL
        assert runtime['conservativeStep']==CONSERVATIVE_STEP and runtime['queryContext']==QUERY_CONTEXT
        assert runtime['samePhysicsAndNumerics'] is False
        assert runtime['samePhysicsAndNumericsWithoutSpatialTable'] is False
        for row in runtime['files']:bind(runtime_path.parent/row['path'],row['sha256'])
        assert not NUMERICAL_CHANGE['originalCFLFeedbackChanged']
        assert not NUMERICAL_CHANGE['effectiveBacktrackingFedToCFL']
        assert not NUMERICAL_CHANGE['wholeQueryPathAdmissibilityGuaranteed']
        report.update(sourceRuntimeFlavor=SOURCE_FLAVOR,taskRuntimeFlavor=FLAVOR,
            sourceRuntimeAndExecutableStable=True,numericalChange=NUMERICAL_CHANGE,
            taskExecutableDiffersFromSource=bind(executable)!=bind(source_executable),
            originalCFLFeedbackUnchanged=True,effectiveBacktrackingFedToCFL=False,
            wholeQueryPathAdmissibilityGuaranteed=False)

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
            selectedSourceRuntimeIdentityPreserved=True)
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
        points=initial[:,1:4].copy();del initial

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

        report['conservativeStepEvents']=verify_step_events(folder,receipt,points,last,normal,out,bind)
        del points

        # 模块六：独立只读清点；原回执缺失仍保留缺失，不向原回执补写字段。
        cp=step_inventory(folder,receipt)
        native_checkpoints='checkpoints' in receipt
        if native_checkpoints:assert cp==receipt['checkpoints'],'离线清点与原监督回执不一致'
        save(out/'read_only_checkpoint_inventory.json',cp)
        native_output_budget='generatedOutputBytes' in receipt and 'outputBudgetExceeded' in receipt
        assert ('generatedOutputBytes' in receipt)==('outputBudgetExceeded' in receipt),'容量监督元数据仅部分存在'
        offline_bytes=output_length(folder,receipt_inputs)
        assert type(offline_bytes) is int and offline_bytes>=0
        maximum_bytes=spec['maximumOutputBytes']
        if normal or not spec['actualValidationMesh']:
            assert offline_bytes<=maximum_bytes,'离线文件长度超过本次计划上限'
        if native_output_budget:
            native_bytes=receipt['generatedOutputBytes']
            assert type(native_bytes) is int and native_bytes>=0
            assert type(receipt['outputBudgetExceeded']) is bool
            if normal:assert native_bytes<=maximum_bytes and not receipt['outputBudgetExceeded']
        if spec['actualValidationMesh']:
            assert native_checkpoints and native_output_budget,'实际候选必须形成原生清点和输出容量监督元数据'
        report.update(nativeCheckpointsRecorded=native_checkpoints,readOnlyCheckpointInventory=True,
            pipelineMetadataMissing=not(native_checkpoints and native_output_budget),
            nativeOutputBudgetRecorded=native_output_budget,offlineGeneratedOutputBytes=offline_bytes,
            maximumOutputBytes=maximum_bytes,offlineOutputWithinBudget=offline_bytes<=maximum_bytes,
            nativeGeneratedOutputBytes=receipt.get('generatedOutputBytes'),
            nativeOutputCapacityEnforcementRecorded=native_output_budget,
            nativeOutputCapacityEnforcementDemonstrated=False,
            outputCapacityReviewSemantics=('原生监督字段存在；本审核只读核对终态，不声称已测试超容量停止'
                if native_output_budget else '本小例运行时输出容量监督未激活/未实测；离线尺寸合格不能称事中强制'),
            originalReceiptModified=False)
        pairs=cp['complete'];steps=[item['localIteration'] for item in pairs]
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
        for path,digest in audit_inputs.items():bind(path,digest,historical_source=True)
        assert mapped_hash(audit_inputs,receipt_path)==bind(receipt_path)
        assert mapped_hash(audit_inputs,manifest_path)==bind(manifest_path)
        assert mapped_hash(audit_inputs,applied_path)==bind(applied_path)
        # 保存合同和热审核各有既定阶段标签；两者不可误当成同一字符串。
        assert audit['normalRunCompleted']==normal and audit['stage']=='same-flow-update-laminar-no-SST'
        report.update(checkpointStage=STAGE,thermalAuditStage=audit['stage'])
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
            abnormalInstantStateAvailable=False,qualificationFromCFLChange=False,
            qualificationFromNumericalStepProtection=False)
    except Exception as exc:
        report.update(passed=False,error=type(exc).__name__+': '+str(exc),traceback=traceback.format_exc())
        raise
    finally:
        stable={path:Path(path).is_file() and sha(Path(path))==value for path,value in locks.items()}
        if not all(stable.values()):report.update(passed=False,inputMutationDetected=True)
        report.update(inputSHA256=locks,inputsStable=stable,historicalSourceRebinding=list(historical_rebinding.values()),
            historicalSourceRebindingUsed=bool(historical_rebinding),elapsedSeconds=time.monotonic()-started)
        save(out/'input_manifest.json',locks);save(out/'report.json',report)
    assert report['passed'] and all(report['inputsStable'].values())
    print(dict(passed=True,normalRunCompleted=normal,frames=steps,
        actualCFLChanged=report['actualHistory']['variesAcrossHistory'],elapsedSeconds=report['elapsedSeconds']))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('folder','audit','output'):parser.add_argument('--'+name,required=True)
    args=parser.parse_args();verify(args.folder,args.audit,args.output)
