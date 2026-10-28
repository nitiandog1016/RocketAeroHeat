function report=studyTmrAxisymmetricExplicitComparison(checkpoints,folder,maxSteps)
%STUDYTMRAXISYMMETRICEXPLICITCOMPARISON 复用真实隐式分段，不重跑已完成求解。
% 模块一：逐段核对Q连接、冻结尺度/边界、SHA及完整当前空间残差重放。
% 模块二：显式从同一最初Q起步，实际计数包括被拒绝尝试的完整残差调用。
% 模块三：固定1e-3初残差仅为工作比较目标；原2e-6稳态标准保持独立。
% 模块四：未达到共同目标时不输出调用量比，历史分段耗时不作受控加速结论。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
assert(iscell(checkpoints)&&isvector(checkpoints)&&~isempty(checkpoints));
validateattributes(maxSteps,{'numeric'},{'real','finite','scalar','integer','positive'});
bindings=struct([]);records=struct([]);calls=0;residuals=[];callHistory=[];last=[];
for j=1:numel(checkpoints)
    path=char(checkpoints{j});loaded=load(path);piece=loaded.implicit;
    assert(~piece.physicalSolution&&~piece.eligibleForSurfaceCoverage&&~piece.pseudoTimeIsPhysical&& ...
        piece.physicalTimeAdvanced==0&&~isfield(piece,'time'));
    if j==1
        mesh=loaded.mesh;gas=loaded.gas;boundary=loaded.boundary;scale=piece.residualScale;Q=piece.initialQ;
        [initialRate,~]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);
        initialResidual=max(abs(initialRate)./scale,[],'all');workTarget=1e-3*initialResidual;
    else
        assert(isequaln(piece.initialQ,last.Q)&&isequaln(loaded.mesh,mesh)&& ...
            isequaln(loaded.gas,gas)&&isequaln(loaded.boundary,boundary)&&isequaln(piece.residualScale,scale), ...
            'TmrCompare:Chain','隐式分段未在同一离散合同下连续连接。');
        assert(strcmp(loaded.report.parentSHA256,bindings(j-1).sha256),'TmrCompare:Parent','父SHA不匹配。');
    end
    [rate,info]=axisymmetricNavierStokesResidual(piece.Q,mesh,gas,boundary);
    assert(isequaln(info,piece.finalInfo)&&max(abs(rate)./scale,[],'all')==piece.residualNorm, ...
        'TmrCompare:Replay','历史隐式段不能由当前空间源码逐值重放。');
    item=struct('path',path,'sha256',fileSHA(path));
    if isempty(bindings),bindings=item;else,bindings(end+1)=item;end %#ok<AGROW>
    record=struct('sourceMAT',path,'sourceMAT_SHA256',item.sha256,'Q',piece.Q,'rate',rate, ...
        'info',info,'residualScale',scale,'mesh',mesh);
    if isempty(records),records=record;else,records(end+1)=record;end %#ok<AGROW>
    first=1;if j>1,first=2;end
    residuals=[residuals,piece.history.residual(first:end)]; %#ok<AGROW>
    callHistory=[callHistory,calls+piece.history.residualEvaluations(first:end)]; %#ok<AGROW>
    calls=calls+piece.residualEvaluations;last=piece;
end
if ~exist(folder,'dir'),mkdir(folder);end
% 当前重放单独留证；这些诊断调用不混入推进器的调用成本。
save(fullfile(folder,'tmr_current_spatial_replay.mat'),'records','mesh','gas','boundary');
opts=struct('maxSteps',maxSteps,'convectiveCFL',.25,'diffusiveCFL',.25, ...
    'residualScale',scale,'steadyResidualTolerance',2e-6,'steadyRequiredSteps',3,'minimumSteadySteps',3);
clock=tic;
try
    explicit=advanceAxisymmetricNavierStokes(Q,mesh,gas,boundary,opts);
catch exception
    failure=struct('status','explicit-comparison-failed','identifier',exception.identifier, ...
        'message',exception.message,'implicitCheckpointBindings',bindings,'physicalSolution',false);
    fid=fopen(fullfile(folder,'tmr_comparison_failure.json'),'w','n','UTF-8');assert(fid>=0);
    fprintf(fid,'%s',jsonencode(failure));fclose(fid);rethrow(exception);
end
explicitSeconds=toc(clock);
indexI=find(residuals<=workTarget,1);indexE=find(explicit.history.residual<=workTarget,1);
callsI=[];callsE=[];ratio=[];
if ~isempty(indexI),callsI=callHistory(indexI);end
if ~isempty(indexE),callsE=explicit.history.residualEvaluations(indexE);end
if ~isempty(callsI)&&~isempty(callsE),ratio=callsE/callsI;end
report=struct('status','completed-budgeted-same-initial-state-numerical-comparison', ...
    'implicitCheckpointBindings',bindings,'initialOriginalResidual',initialResidual, ...
    'fixedResidualWorkTarget',workTarget,'requiredSteadyTolerance',2e-6, ...
    'implicitOriginalResidual',last.residualNorm,'implicitSteadyConverged',last.steadyConverged, ...
    'implicitTotalResidualEvaluations',calls,'implicitCallsAtWorkTarget',callsI, ...
    'explicitOriginalResidual',explicit.residualNorm,'explicitSteadyConverged',explicit.steadyConverged, ...
    'explicitAcceptedSteps',explicit.acceptedSteps,'explicitRejectedAttempts',explicit.rejectedAttempts, ...
    'explicitTotalResidualEvaluations',explicit.residualEvaluations,'explicitCallsAtWorkTarget',callsE, ...
    'sameTargetResidualEvaluationRatio',ratio,'explicitElapsedSeconds',explicitSeconds, ...
    'diagnosticReplayResidualEvaluations',numel(checkpoints)+1, ...
    'controlledWallTimeSpeedupMeasured',false,'pseudoTimeIsPhysical',false, ...
    'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'tmr_explicit_comparison.mat'),'report','explicit','mesh','gas','boundary', ...
    'bindings','residuals','callHistory');
fid=fopen(fullfile(folder,'tmr_explicit_comparison.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('TMR same-initial explicit comparison: %.8g -> %.8g; actual calls %d; common-target ratio %s.\n', ...
    initialResidual,explicit.residualNorm,explicit.residualEvaluations,mat2str(ratio));
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
