function report=testRocketCollarNewtonDiagnostics(folder,prefixFolder,stepFile)
%TESTROCKETCOLLARNEWTONDIAGNOSTICS 用真实七步内部状态重演第八步的线搜索失败。
% 模块一：严格绑定冻结工况、旧端点与新源码；是独立诊断重演，不续写旧链。
% 模块二：同Q/CFL/门槛开启记录，完整result逐值对照旧第八步，包含所有工作量。
% 模块三：逐候选重算完整五量残差或同一状态错误，区分越界与有效但不下降。
% 模块四：保留真实试探Q及表格；内部暖初猜、给定Tw和未稳态物理解分开。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v71_validation','rocket_failure_replay');end
if nargin<2,prefixFolder=fullfile(root,'outputs','v70_validation','normal_collar');end
if nargin<3,stepFile=fullfile(root,'outputs','v70_validation','normal_collar_continued','steps','step_000008.mat');end
assert(~exist(folder,'dir'),'RocketNewtonDiagnostic:Overwrite','不能覆盖实际诊断证据。');
prefixFile=fullfile(prefixFolder,'rocket_axisymmetric_evidence.mat');
old=load(prefixFile,'caseDefinition','implicit','durable');
step=load(stepFile,'segment','stepContext');c=old.caseDefinition;baseline=step.segment;
assert(old.implicit.acceptedSteps==7&&step.stepContext.number==8&&isequal(old.implicit.Q,baseline.initialQ)&& ...
    baseline.options.initialPseudoCFL==old.durable.nextPseudoCFL&&baseline.options.differenceStep==2e-6&& ...
    baseline.rejectedAttempts==1&&strcmp(baseline.rejectLog.reason,'line-search-failed'));
fresh=createRocketAxisymmetricCollarCase(c.options);assert(isequaln(c,fresh));
bindings=old.durable.sourceBindings;currentBindings=bindings;
for k=1:numel(bindings)
    currentBindings(k).sha256=fileSHA(fullfile(root,bindings(k).path));
    assert(strcmp(currentBindings(k).sha256,bindings(k).sha256), ...
        'RocketNewtonDiagnostic:Source','冻结原计算来源改变：%s。',bindings(k).path);
end
inputs=struct('prefixFile',prefixFile,'prefixSHA256',fileSHA(prefixFile), ...
    'stepFile',stepFile,'stepSHA256',fileSHA(stepFile),'frozenSourceBindings',bindings, ...
    'currentSourceBindings',currentBindings,'initialGuessKind','internal-frozen-accepted-state-not-boundary-or-physical-solution', ...
    'observationSource',struct('path','regionalFlow/axisymmetric/advanceAxisymmetricSteadyObserved.m', ...
        'sha256',fileSHA(fullfile(root,'regionalFlow','axisymmetric','advanceAxisymmetricSteadyObserved.m'))), ...
    'oldChainModified',false,'newAcceptedChainClaimed',false);
clock=tic;
[result,diagnostics]=advanceAxisymmetricSteadyObserved(old.implicit.Q,c.mesh,c.gas,c.boundary, ...
    baseline.options,struct('enabled',true,'maxRecords',512));
seconds=toc(clock);mkdir(folder);
save(fullfile(folder,'newton_diagnostic_raw.mat'),'result','diagnostics','inputs','seconds','c','baseline','-v7');
assert(isequaln(result,baseline),'RocketNewtonDiagnostic:Changed','完整result或调用计数与冻结旧第八步不同。');
assert(diagnostics.complete&&diagnostics.seenRecords==numel(diagnostics.records)&&diagnostics.droppedRecords==0);
records=diagnostics.records;scale=result.residualScale;U0=result.initialQ./scale;
assert(result.acceptedSteps==1&&result.rejectedAttempts==1&&~result.steadyConverged&& ...
    ~result.physicalSolution&&~result.eligibleForSurfaceCoverage);
replayError=0;rows=struct([]);
for k=1:numel(records)
    item=records(k);errorID='';candidateNorm=NaN;
    try
        [rate,~]=axisymmetricRANSResidual(item.candidateQ,c.mesh,c.gas,c.boundary);
    catch exception
        errorID=exception.identifier;
    end
    if strcmp(item.status,'state-failure')
        assert(~isempty(errorID)&&strcmp(errorID,item.errorIdentifier)&&isnan(item.candidateDefect), ...
            'RocketNewtonDiagnostic:ErrorReplay','候选错误码不能从保存真实Q重演。');
    else
        assert(isempty(errorID)&&isempty(item.errorIdentifier));
        candidateNorm=max(abs(item.candidateQ./scale-U0-item.dt*(rate./scale)),[],'all');
        replayError=max(replayError,abs(candidateNorm-item.candidateDefect));
        budget=512*eps*max(abs(item.candidateQ./scale)+abs(U0)+abs(item.dt*(rate./scale)),[],'all');
        assert(abs(candidateNorm-item.candidateDefect)<=budget);
        accepted=item.candidateDefect<=item.armijoLimit||item.candidateDefect<=item.target;
        assert(accepted==strcmp(item.status,'accepted'));
    end
    assert(item.alpha==2^(-item.searchIndex)&&item.step==1&&item.jacobianRefreshedThisIteration);
    assert(max(abs(item.candidateQ./scale-item.iterateQ./scale-item.alpha*item.incrementScaled),[],'all')<128*eps* ...
        max(1,max(abs(item.candidateQ./scale)+abs(item.iterateQ./scale)+abs(item.alpha*item.incrementScaled),[],'all')));
    row=rmfield(item,{'iterateQ','candidateQ','incrementScaled','errorMessage'});
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
bad=records([records.attempt]==0);lastIteration=max([bad.newtonIteration]);
terminal=bad([bad.newtonIteration]==lastIteration);
assert(numel(terminal)==baseline.options.maxLineSearch+1&&~any(strcmp({terminal.status},'accepted')));
assert(all([terminal.dt]==result.rejectLog.dt)&&all([terminal.baseDefect]==result.rejectLog.nonlinearDefect));
stateIDs={records(strcmp({records.status},'state-failure')).errorIdentifier};
report=struct('status','passed-real-frozen-Q7-new-source-Newton-observation-replay', ...
    'matlabEngineProcessID',feature('getpid'),'cellCount',c.mesh.cellCount, ...
    'completeResultAndWorkCountsEqualFrozenV70StepEight',true,'allCandidateStatesReplayed',true, ...
    'recordCount',numel(records),'stateFailureCount',sum(strcmp({records.status},'state-failure')), ...
    'validNotDecreasingCount',sum(strcmp({records.status},'valid-not-decreasing')), ...
    'acceptedLineSearchCount',sum(strcmp({records.status},'accepted')), ...
    'stateErrorIdentifiers',{unique(stateIDs)},'terminalRejectedIteration',lastIteration, ...
    'terminalCandidateStatuses',{{terminal.status}},'terminalCandidateDefects',[terminal.candidateDefect], ...
    'terminalCandidateAlphas',[terminal.alpha],'terminalBaseDefect',terminal(1).baseDefect, ...
    'terminalTarget',terminal(1).target,'candidateReplayMaximumNormDifference',replayError, ...
    'acceptedPseudoStepS',result.history.dt(end),'rejectedPseudoStepS',result.rejectLog.dt, ...
    'residualEvaluations',result.residualEvaluations,'elapsedSecondsWithObservation',seconds, ...
    'fiveQuantityBalanceRelative',result.globalBalanceRelative, ...
    'inputBinding',inputs,'diagnosticCoversJacobianPerturbations',false, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'wholeGoalComplete',false);
save(fullfile(folder,'newton_diagnostic_evidence.mat'),'report','result','diagnostics','inputs','c','baseline','rows','-v7');
writetable(struct2table(rows),fullfile(folder,'candidate_observations.csv'));
fid=fopen(fullfile(folder,'newton_diagnostic_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;
fprintf('Real Newton diagnostic passed: %d records, %d state failures, %d valid non-decreasing; terminal Newton=%d.\n', ...
    report.recordCount,report.stateFailureCount,report.validNotDecreasingCount,lastIteration);
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
