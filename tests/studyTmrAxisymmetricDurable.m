function report=studyTmrAxisymmetricDurable(folder,stepsThisRun,warmEvidenceFile)
%STUDYTMRAXISYMMETRICDURABLE 一接受步一检查点的真实4000格SA求稳态链。
% 模块一：从已验证normal稳态逻辑复制暖初猜，重新建立4000格初始账。
% 模块二：每次advance只接受一BE步；原方程、尺度、Newton及稳态门槛不变。
% 模块三：保存所有接受端点Q和局部原result，逐步重算全链积分/历史与真实成本。
% 模块四：原子写最新检查点；重启校验SHA/Q连接，并恢复已写端点而未更新索引的孤儿步。
% 模块五：全链仅以真正接受的末三步判断稳态，绝不把重启首行或旧中断52步计入。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v67_validation','both');end
if nargin<2,stepsThisRun=5;end
if nargin<3
    warmEvidenceFile=fullfile(root,'outputs','v67_validation','normal','tmr_sensitivity_evidence.mat');
    if ~exist(warmEvidenceFile,'file'),warmEvidenceFile=fullfile(root,'validation','v67','normal','tmr_sensitivity_evidence.mat');end
end
validateattributes(stepsThisRun,{'numeric'},{'real','finite','scalar','positive','integer'});
assert((ischar(folder)&&isrow(folder))||(isstring(folder)&&isscalar(folder)), ...
    'AxisDurable:Folder','检查点目录必须是单一路径。');folder=char(folder);
if ~exist(folder,'dir'),mkdir(folder);end
stepFolder=fullfile(folder,'steps');if ~exist(stepFolder,'dir'),mkdir(stepFolder);end
target=fullfile(folder,'tmr_sensitivity_evidence.mat');bindings=sourceBindings(root);pid=feature('getpid');
identity=struct('runLabel',getenv('ROCKET_AEROHEAT_RUN_LABEL'),'engineProcessID',pid, ...
    'launcherProcessID',NaN,'processRecordPath',getenv('ROCKET_AEROHEAT_PROCESS_RECORD'));
if ~isempty(identity.processRecordPath)&&exist(identity.processRecordPath,'file')
    runRecord=jsondecode(fileread(identity.processRecordPath));identity.launcherProcessID=runRecord.processID;
end
if exist(target,'file')
    saved=load(target,'caseDefinition','mesh','gas','boundary','Q','initialInfo','implicit','durable','report');
    assert(all(isfield(saved,{'caseDefinition','mesh','gas','boundary','Q','initialInfo','implicit','durable'})), ...
        'AxisDurable:Checkpoint','已有文件不是持久化检查点，禁止覆盖旧研究结果。');
    caseDefinition=saved.caseDefinition;mesh=saved.mesh;gas=saved.gas;boundary=saved.boundary;
    Q=saved.Q;initialInfo=saved.initialInfo;implicit=saved.implicit;durable=saved.durable;
    assert(strcmp(durable.kind,'one-accepted-BE-step-checkpoint-chain-v1')&&isequaln(durable.sourceBindings,bindings)&& ...
        strcmp(fileSHA(durable.warmBinding.path),durable.warmBinding.sha256)&& ...
        numel(durable.steps)==implicit.acceptedSteps&&isequal(Q,implicit.initialQ), ...
        'AxisDurable:Checkpoint','源码、暖源或真实接受步链已改变。');
    for k=1:numel(durable.steps)
        assert(strcmp(fileSHA(fullfile(folder,durable.steps(k).relativePath)),durable.steps(k).sha256), ...
            'AxisDurable:StepSHA','已保存端点文件SHA改变。');
    end
    [replay,replayInfo]=axisymmetricRANSResidual(implicit.Q,mesh,gas,boundary);durable.diagnosticResidualEvaluations=durable.diagnosticResidualEvaluations+1;
    assert(isequaln(replayInfo,implicit.finalInfo)&&max(abs(replay)./implicit.residualScale,[],'all')==implicit.residualNorm, ...
        'AxisDurable:Replay','重启终态不能逐值重放同一完整方程。');
else
    coarseFactory=createTmrAxisymmetricSensitivityCase('both');
    warmSHA=fileSHA(warmEvidenceFile);warm=load(warmEvidenceFile,'mesh','gas','boundary','implicit','report');
    assert(warm.implicit.steadyConverged&&warm.report.steadyConverged&&isequal(size(warm.mesh.X),[41 51])&& ...
        warm.implicit.residualNorm<=2e-6&&all(warm.implicit.history.residual(end-2:end)<=2e-6)&& ...
        isequaln(warm.gas,coarseFactory.gas)&&isequaln(warm.boundary,coarseFactory.boundary)&& ...
        isequal(warm.implicit.residualScale,coarseFactory.residualScale), ...
        'AxisDurable:WarmSource','暖源须是同气体/边界/固定尺度的真实normal稳态。');
    [warmRate,warmInfo]=axisymmetricRANSResidual(warm.implicit.Q,warm.mesh,warm.gas,warm.boundary);
    assert(isequaln(warmInfo,warm.implicit.finalInfo)&&max(abs(warmRate)./warm.implicit.residualScale,[],'all')==warm.implicit.residualNorm);
    mesh=coarseFactory.mesh;gas=coarseFactory.gas;boundary=coarseFactory.boundary;scale=coarseFactory.residualScale;
    assert(isequal(mesh.X(1:2:end,:),warm.mesh.X)&&isequal(mesh.R(1:2:end,:),warm.mesh.R));
    [ii,jj]=ind2sub(mesh.shape,(1:mesh.cellCount)');parentI=floor((ii-1)/2)+1;parentJ=jj;
    parentCellID=sub2ind(warm.mesh.shape,parentI,parentJ);Q=warm.implicit.Q(parentCellID,:);
    initialIntegral=sum(Q.*mesh.volumePerRadian,1);parentIntegral=sum(warm.implicit.Q.*warm.mesh.volumePerRadian,1);
    grouped=accumarray(parentCellID,mesh.volumePerRadian,[warm.mesh.cellCount 1]);difference=grouped-warm.mesh.volumePerRadian;
    warmBinding=struct('path',warmEvidenceFile,'sha256',warmSHA,'parentCellCount',warm.mesh.cellCount, ...
        'parentNodeShape',size(warm.mesh.X),'parentOriginalResidual',warm.implicit.residualNorm, ...
        'completeSpatialReplayPassed',true,'parentCellID',parentCellID,'parentI',parentI,'parentJ',parentJ, ...
        'logicalFactorI',2,'logicalFactorJ',1,'conservativeRemapPerformed',false, ...
        'initialConservedIntegralPerRadian',initialIntegral,'parentConservedIntegralPerRadian',parentIntegral, ...
        'physicalSolution',false,'eligibleForSurfaceCoverage',false);
    caseDefinition=coarseFactory;caseDefinition.Q=Q;caseDefinition.coarseFactoryInitialization=coarseFactory.gridStudy;
    caseDefinition.coarseFactoryParentBinding=coarseFactory.parentBinding;caseDefinition.parentBinding=warmBinding;
    g=caseDefinition.gridStudy;g.initialGuessKind='logical-copy-from-verified-normal-steady-not-conservative-remap';
    g.baselineRefinementFactorI=2;g.baselineRefinementFactorJ=2;g.factorI=2;g.factorJ=1;
    g.parentCellID=parentCellID;g.parentI=parentI;g.parentJ=parentJ;
    g.initialConservedIntegralPerRadian=initialIntegral;g.parentConservedIntegralPerRadian=parentIntegral;
    g.initialIntegralDifferencePerRadian=initialIntegral-parentIntegral;
    g.initialIntegralDifferenceOverFixedVolumeScale=(initialIntegral-parentIntegral)./(sum(warm.mesh.volumePerRadian)*scale);
    g.parentVolumePerRadianM3=sum(warm.mesh.volumePerRadian);g.volumeDifferencePerRadianM3=sum(mesh.volumePerRadian)-g.parentVolumePerRadianM3;
    g.childGroupedVolumeDifferencePerParentM3=difference;g.maximumGroupedVolumeRelativeDifference=max(abs(difference)./warm.mesh.volumePerRadian);
    caseDefinition.gridStudy=g;caseDefinition.warmInitialization=warmBinding;
    [rate,initialInfo]=axisymmetricRANSResidual(Q,mesh,gas,boundary);
    opts=struct('maxSteps',1,'initialPseudoCFL',.5,'maximumPseudoCFL',1e6,'pseudoCFLGrowth',2, ...
        'maxNewtonIterations',24,'jacobianRefreshEvery',1,'differenceScheme','central','differenceStep',2e-6, ...
        'residualScale',scale,'steadyResidualTolerance',2e-6,'steadyRequiredSteps',3,'minimumSteadySteps',3, ...
        'conservationTolerance',5e-10,'nonlinearAbsoluteTolerance',1e-11,'nonlinearRelativeTolerance',1e-9);
    implicit=initialResult(Q,mesh,gas,initialInfo,rate,opts);
    durable=struct('kind','one-accepted-BE-step-checkpoint-chain-v1','warmBinding',warmBinding, ...
        'sourceBindings',bindings,'steps',struct('number',{},'relativePath',{},'sha256',{},'matlabProcessID',{}), ...
        'nextPseudoCFL',.5,'solverElapsedSeconds',0,'diagnosticResidualEvaluations',3, ...
        'recoveredOrphanSteps',0,'runCount',0,'allAcceptedEndpointsSaved',true, ...
        'interruptedCoarseInitializationRunIncluded',false,'physicalTimeAdvancedS',0);
    durable.runIdentities=struct('runLabel',{},'engineProcessID',{},'launcherProcessID',{},'processRecordPath',{});
    assert(strcmp(warmSHA,fileSHA(warmEvidenceFile)),'AxisDurable:WarmSource','初始化期间暖源文件改变。');
end
durable.runCount=durable.runCount+1;
durable.runIdentities(end+1)=identity;
% 若上一进程写完不可变端点但未替换latest，则从真实文件恢复，不能覆写或重算冒充同一步。
while exist(fullfile(stepFolder,sprintf('step_%06d.mat',implicit.acceptedSteps+1)),'file')
    name=sprintf('steps/step_%06d.mat',implicit.acceptedSteps+1);path=fullfile(folder,strrep(name,'/',filesep));
    item=load(path,'segment','stepContext');
    assert(isequal(item.segment.initialQ,implicit.Q)&&isequaln(item.stepContext.sourceBindings,bindings), ...
        'AxisDurable:Orphan','孤儿端点不能与当前链的Q/源码连接。');
    validateOrphan(item.segment,item.stepContext,implicit,durable,mesh,gas,boundary);
    durable.diagnosticResidualEvaluations=durable.diagnosticResidualEvaluations+2;
    [implicit,durable]=appendStep(implicit,durable,item.segment,item.stepContext,mesh,path,name);
    durable.diagnosticResidualEvaluations=durable.diagnosticResidualEvaluations+item.stepContext.diagnosticResidualEvaluations;
    durable.recoveredOrphanSteps=durable.recoveredOrphanSteps+1;
end
publish();
for attempt=1:stepsThisRun
    if implicit.steadyConverged,break;end
    options=implicit.options;options.maxSteps=1;options.initialPseudoCFL=durable.nextPseudoCFL;
    [beforeRate,beforeInfo]=axisymmetricRANSResidual(implicit.Q,mesh,gas,boundary);durable.diagnosticResidualEvaluations=durable.diagnosticResidualEvaluations+1;
    assert(isequaln(beforeInfo,implicit.finalInfo)&&max(abs(beforeRate)./implicit.residualScale,[],'all')==implicit.residualNorm);
    explicit0=1/max(beforeInfo.convectiveFrequency+beforeInfo.diffusiveFrequency);
    clock=tic;segment=advanceAxisymmetricSteady(implicit.Q,mesh,gas,boundary,options);elapsed=toc(clock);
    if segment.acceptedSteps~=1
        failure=struct('status','one-step-solver-rejected-not-advanced','segment',segment,'sourceBindings',bindings, ...
            'matlabProcessID',pid,'attemptedGlobalStep',implicit.acceptedSteps+1);
        path=fullfile(folder,sprintf('rejected_%06d_pid%d.mat',implicit.acceptedSteps+1,pid));save(path,'failure');
        error('AxisDurable:Rejected','单BE步被原求解器拒绝，已保存失败，不改变已有接受链。');
    end
    [rate,info]=axisymmetricRANSResidual(segment.Q,mesh,gas,boundary);durable.diagnosticResidualEvaluations=durable.diagnosticResidualEvaluations+1;
    assert(isequaln(info,segment.finalInfo)&&max(abs(rate)./segment.residualScale,[],'all')==segment.residualNorm);
    stepContext=struct('number',implicit.acceptedSteps+1,'matlabProcessID',pid,'processIdentity',identity,'sourceBindings',bindings, ...
        'explicitStepAtInitialState',explicit0,'elapsedSolverSeconds',elapsed, ...
        'nextPseudoCFL',min(options.maximumPseudoCFL,options.pseudoCFLGrowth*segment.history.dt(end)/explicit0), ...
        'diagnosticResidualEvaluations',2,'localSteadyFlagWasNotPromoted',~segment.steadyConverged, ...
        'physicalSolution',false,'eligibleForSurfaceCoverage',false);
    name=sprintf('steps/step_%06d.mat',stepContext.number);path=fullfile(folder,strrep(name,'/',filesep));
    assert(~exist(path,'file'),'AxisDurable:Overwrite','禁止覆盖已存在的真实接受端点。');
    temporary=[tempname(stepFolder) '.mat'];save(temporary,'segment','stepContext');movefile(temporary,path);
    [implicit,durable]=appendStep(implicit,durable,segment,stepContext,mesh,path,name);
    publish();
    fprintf('Durable accepted=%d, original residual=%.9g, stored endpoint=%s, steady=%d.\n', ...
        implicit.acceptedSteps,implicit.residualNorm,name,implicit.steadyConverged);
end
publish();

    function publish()
        % 模块用途：完整最新状态经临时文件替换，所有不可变step另以SHA留存。
        report=struct('status','actual-durable-one-step-SA-chain-not-experimental-acceptance','kind','both', ...
            'gridStudy',caseDefinition.gridStudy,'parentBinding',durable.warmBinding,'cellCount',mesh.cellCount, ...
            'initialOriginalResidual',implicit.history.residual(1),'finalOriginalResidual',implicit.residualNorm, ...
            'steadyConverged',implicit.steadyConverged,'stopReason',implicit.stopReason, ...
            'acceptedSteps',implicit.acceptedSteps,'rejectedAttempts',implicit.rejectedAttempts, ...
            'solverResidualEvaluations',implicit.residualEvaluations,'diagnosticResidualEvaluations',durable.diagnosticResidualEvaluations, ...
            'elapsedSeconds',durable.solverElapsedSeconds,'fiveQuantityBalanceRelative',implicit.globalBalanceRelative, ...
            'signedWallHeatRangeWm2',[min(implicit.finalInfo.wall.heatFluxIntoWall),max(implicit.finalInfo.wall.heatFluxIntoWall)], ...
            'allAcceptedEndpointsSaved',true,'warmSourceKind','independently verified normal-grid steady; fresh logical-copy initial guess', ...
            'oldInterrupted52StepsIncluded',false,'initializationOnlyNotConservativeRemap',true, ...
            'steadyCriterionAppliedToActualAcceptedEndpoints',true,'currentMatlabProcessID',pid, ...
            'currentProcessIdentity',identity, ...
            'physicalTimeAdvanced',0,'pseudoTimeIsPhysical',false,'physicalSolution',false,'experimentalValidation',false, ...
            'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
        payload=struct('report',report,'caseDefinition',caseDefinition,'mesh',mesh,'gas',gas,'boundary',boundary, ...
            'Q',Q,'initialInfo',initialInfo,'implicit',implicit,'durable',durable);
        temporary=[tempname(folder) '.mat'];save(temporary,'-struct','payload');
        movefile(temporary,target,'f');
        writeJSON(fullfile(folder,'tmr_sensitivity_summary.json'),report);
    end
end

function validateOrphan(segment,context,previous,durable,mesh,gas,boundary)
% 模块用途：孤儿步必须重新证明同一方程、原实际BE缺陷和真实边界/源积分，不能只信文件名。
assert(isequaln(segment.gas,gas)&&isequaln(segment.geometryBinding,previous.geometryBinding)&& ...
    isequal(segment.options.residualScale,previous.residualScale)&& ...
    segment.options.initialPseudoCFL==durable.nextPseudoCFL,'AxisDurable:Orphan','孤儿步气体/几何/尺度/CFL不同。');
for key={'steadyResidualTolerance','steadyRequiredSteps','minimumSteadySteps','conservationTolerance', ...
        'nonlinearAbsoluteTolerance','nonlinearRelativeTolerance','differenceStep','differenceScheme','jacobianRefreshEvery'}
    k=key{1};assert(isequal(segment.options.(k),previous.options.(k)),'AxisDurable:Orphan','孤儿步修改原数值门槛或方法。');
end
[r0,i0]=axisymmetricRANSResidual(previous.Q,mesh,gas,boundary);
[r1,i1]=axisymmetricRANSResidual(segment.Q,mesh,gas,boundary);
assert(isequaln(i0,previous.finalInfo)&&isequaln(i1,segment.finalInfo)&& ...
    max(abs(r1)./previous.residualScale,[],'all')==segment.residualNorm, ...
    'AxisDurable:Orphan','孤儿步初末完整方程不能重放。');
dt=segment.history.dt(end);assert(isfinite(dt)&&dt>0);
assert(isequal(segment.accumulatedBoundaryIntegral,dt*i1.netBoundaryFlux)&& ...
    isequal(segment.accumulatedSourceIntegral,dt*i1.sourceIntegral)&& ...
    context.explicitStepAtInitialState==1/max(i0.convectiveFrequency+i0.diffusiveFrequency)&& ...
    context.nextPseudoCFL==min(segment.options.maximumPseudoCFL,segment.options.pseudoCFLGrowth*dt/context.explicitStepAtInitialState), ...
    'AxisDurable:Orphan','孤儿步边界/源积分或实际CFL恢复数据不一致。');
U0=previous.Q./previous.residualScale;U1=segment.Q./previous.residualScale;scaledRate=r1./previous.residualScale;
defect=U1-U0-dt*scaledRate;
target=segment.options.nonlinearAbsoluteTolerance+segment.options.nonlinearRelativeTolerance*max(abs(dt*(r0./previous.residualScale)),[],'all');
roundoff=512*eps*max(abs(U0)+abs(U1)+abs(dt*scaledRate),[],'all');
assert(target==segment.history.newtonTarget(end)&&max(abs(defect),[],'all')<=target+roundoff, ...
    'AxisDurable:Orphan','孤儿步真实BE端点缺陷未通过原Newton目标及有单位舍入预算。');
end

function [result,durable]=appendStep(previous,durable,segment,context,mesh,path,name)
% 模块用途：保留局部原result，以每步实际Q及累计边界/源重算全链history和账。
assert(segment.acceptedSteps==1&&~segment.steadyConverged&&numel(segment.history.residual)==2&& ...
    isequal(segment.initialQ,previous.Q)&&isequal(segment.residualScale,previous.residualScale)&& ...
    context.number==previous.acceptedSteps+1&&segment.globalBalanceRelative<=5e-10, ...
    'AxisDurable:StepJoin','只接受原单BE步，重复首残差不能充当新的接受步。');
result=previous;n=previous.acceptedSteps+1;result.Q=segment.Q;result.finalInfo=segment.finalInfo;
result.finalConservedIntegral=sum(segment.Q.*mesh.volumePerRadian,1);
result.accumulatedBoundaryIntegral=previous.accumulatedBoundaryIntegral+segment.accumulatedBoundaryIntegral;
result.accumulatedSourceIntegral=previous.accumulatedSourceIntegral+segment.accumulatedSourceIntegral;
scale=sum(mesh.volumePerRadian)*result.residualScale;
defect=result.finalConservedIntegral-result.initialConservedIntegral+result.accumulatedBoundaryIntegral-result.accumulatedSourceIntegral;
balance=max(abs(defect)./scale);assert(balance<=result.options.conservationTolerance,'AxisDurable:Ledger','累计五量账未通过原门槛。');
endpointDefect=result.finalConservedIntegral-previous.finalConservedIntegral+segment.accumulatedBoundaryIntegral-segment.accumulatedSourceIntegral;
assert(all(abs(endpointDefect)./scale<=result.options.conservationTolerance),'AxisDurable:Ledger','真实端点的一步五量账未通过。');
for key=fieldnames(previous.history)'
    k=key{1};row=segment.history.(k);value=row(end,:);
    if isrow(row)&&size(row,2)==2&&~strcmp(k,'stepDefect'),value=row(end);end
    switch k
        case 'pseudoTime',value=previous.pseudoTime+segment.pseudoTime;
        case 'balance',value=balance;
        case 'massRelativeDrift',value=(result.finalConservedIntegral(1)-result.initialConservedIntegral(1))/result.initialConservedIntegral(1);
        case 'residualEvaluations',value=previous.residualEvaluations+segment.residualEvaluations;
        case 'stepDefect',value=endpointDefect;
    end
    if strcmp(k,'stepDefect'),result.history.(k)(end+1,:)=value;else,result.history.(k)(end+1)=value;end
end
result.pseudoTime=previous.pseudoTime+segment.pseudoTime;result.acceptedSteps=n;
for key={'rejectedAttempts','residualEvaluations','linearSolves','factorizations','jacobianBuilds'}
    k=key{1};result.(k)=previous.(k)+segment.(k);
end
for k=1:numel(segment.rejectLog)
    entry=segment.rejectLog(k);entry.step=entry.step+previous.acceptedSteps;result.rejectLog(end+1)=entry;
end
result.residualNorm=segment.residualNorm;result.globalBalanceRelative=balance;result.options.maxSteps=max(3,n);
result.steadyConverged=n>=result.options.minimumSteadySteps&& ...
    all(result.history.residual(end-2:end)<=result.options.steadyResidualTolerance);
if result.steadyConverged,result.stopReason='steady-residual-converged';else,result.stopReason='persisted-accepted-step-budget';end
durable.nextPseudoCFL=context.nextPseudoCFL;durable.solverElapsedSeconds=durable.solverElapsedSeconds+context.elapsedSolverSeconds;
durable.steps(end+1)=struct('number',n,'relativePath',name,'sha256',fileSHA(path),'matlabProcessID',context.matlabProcessID);
end

function result=initialResult(Q,mesh,gas,info,rate,options)
% 模块用途：尚未接受任何步的真实初态；不以初始化残差伪造推进成本或稳态。
initial=sum(Q.*mesh.volumePerRadian,1);residual=max(abs(rate)./options.residualScale,[],'all');
history=struct('pseudoTime',0,'dt',0,'residual',residual,'balance',0,'newtonIterations',0, ...
    'newtonDefect',0,'newtonTarget',0,'stepDefect',zeros(1,5),'massRelativeDrift',0, ...
    'lineSearchTrials',0,'residualEvaluations',0,'linearDefect',0,'augmentedLinearDefect',0,'massConstraintDefect',0);
result=struct('Q',Q,'initialQ',Q,'pseudoTime',0,'acceptedSteps',0,'rejectedAttempts',0,'steadyConverged',false, ...
    'stopReason','initial-state-not-advanced','residualNorm',residual,'residualScale',options.residualScale, ...
    'history',history,'rejectLog',struct('step',{},'dt',{},'reason',{},'nonlinearDefect',{}), ...
    'initialConservedIntegral',initial,'finalConservedIntegral',initial,'accumulatedBoundaryIntegral',zeros(1,5), ...
    'accumulatedSourceIntegral',zeros(1,5),'globalBalanceRelative',0,'finalInfo',info,'gas',gas,'options',options, ...
    'residualEvaluations',0,'linearSolves',0,'factorizations',0,'jacobianBuilds',0,'closedMassConstraintActive',false, ...
    'geometryBinding',struct('X',mesh.X,'R',mesh.R,'cells',mesh.cells,'volumePerRadian',mesh.volumePerRadian), ...
    'pseudoTimeIsPhysical',false,'physicalTimeAdvanced',0,'wallTemperatureValueKind','prescribed', ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false);
end

function out=sourceBindings(root)
% 模块用途：持久化重启只接受同一组实际执行源码，拒绝悄悄换方程或修改驱动。
paths={'tests/studyTmrAxisymmetricDurable.m','regionalFlow/axisymmetric/advanceAxisymmetricSteady.m', ...
    'regionalFlow/axisymmetric/axisymmetricRANSResidual.m','regionalFlow/axisymmetric/axisymmetricNavierStokesResidual.m', ...
    'regionalFlow/axisymmetric/axisymmetricWallThermalGradient.m','regionalFlow/axisymmetric/axisymmetricFVJacobian.m', ...
    'regionalFlow/axisymmetric/createTmrAxisymmetricSensitivityCase.m', ...
    'regionalFlow/axisymmetric/axisymmetricSADiscretization.m','regionalFlow/axisymmetric/axisymmetricSAModel.m', ...
    'regionalFlow/axisymmetric/axisymmetricSAWallDistance.m','regionalFlow/axisymmetric/axisymmetricFVState.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVConservative.m','regionalFlow/axisymmetric/axisymmetricFVGas.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVTransport.m','regionalFlow/axisymmetric/axisymmetricFVTemperatureTrace.m', ...
    'regionalFlow/axisymmetric/buildAxisymmetricFVGeometry.m','regionalFlow/axisymmetric/loadTmrMach7Grid.m', ...
    'regionalFlow/axisymmetric/createTmrAxisymmetricAirGas.m'};
out=struct('path',{},'sha256',{});for k=1:numel(paths),out(k)=struct('path',paths{k},'sha256',fileSHA(fullfile(root,paths{k})));end
end
function writeJSON(path,value)
temporary=[path '.tmp'];fid=fopen(temporary,'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(value));clear cleanup;movefile(temporary,path,'f');
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0,'AxisDurable:Source','缺少链输入/端点：%s',path);cleanup=onCleanup(@()fclose(fid));
bytes=fread(fid,Inf,'*uint8');clear cleanup;md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
