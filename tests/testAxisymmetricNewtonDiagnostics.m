function report=testAxisymmetricNewtonDiagnostics(folder)
%TESTAXISYMMETRICNEWTONDIAGNOSTICS 小型闭域真实BE验证旁路诊断不改变任何求解结果。
% 模块一：非均匀九格环隙激活四量方程、墙热和闭域质量约束。
% 模块二：关闭/开启/截断记录逐字段比较完整result及调用计数。
% 模块三：独立重算候选BE缺陷，检查输入错误继续抛出、容量损失明确报告。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v71_validation','diagnostics_unit');end
assert(~exist(folder,'dir'),'AxisDiagnosticTest:Overwrite','测试不能覆盖既有证据。');
gas=struct('gamma',1.4,'R',1,'mu',.5,'Pr',.72,'temperatureRange',[200 500]);
[X,R]=ndgrid([0 .02 .05 .08],[.03 .035 .045 .06]);
mesh=buildAxisymmetricFVGeometry(X,R,struct('periodicI',true));
x=mesh.volumeCentroid(:,1);r=mesh.volumeCentroid(:,2);phase=2*pi*x/.08;xi=log(r/.03)/log(2);
P=[1+.03*cos(phase)+.02*sin(pi*xi),-.4+.08*sin(phase).*(.5+xi), ...
    .12*sin(pi*xi).*(1+.15*cos(phase)),300+.2*cos(phase)+.1*sin(pi*xi)];
Q=axisymmetricFVConservative(P,gas);periodic=struct('type','periodic');
boundary=struct('iMin',periodic,'iMax',periodic, ...
    'jMin',struct('type','wall','temperature',300,'velocity',[-3 0]), ...
    'jMax',struct('type','wall','temperature',302,'velocity',[0 0]));
options=struct('maxSteps',1,'maximumPseudoStep',1e-7,'residualScale',[1 20 20 800], ...
    'steadyResidualTolerance',0,'jacobianRefreshEvery',1);
plain=advanceAxisymmetricSteady(Q,mesh,gas,boundary,options);
[defaultOff,off]=advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options);
[enabled,diag]=advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options,struct('enabled',true));
[disabled,explicitOff]=advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options,struct('enabled',false));
[truncated,cap]=advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options,struct('enabled',true,'maxRecords',1));
[emptyCap,zero]=advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options,struct('enabled',true,'maxRecords',0));
assert(isequaln(plain,enabled)&&isequaln(plain,disabled)&&isequaln(plain,truncated)&&isequaln(plain,emptyCap), ...
    'AxisDiagnosticTest:Observation','诊断改变了完整结果或工作量计数。');
assert(isequaln(plain,defaultOff)&&~off.enabled&&~off.complete&&off.seenRecords==0&&isempty(off.records)&&isequaln(off,explicitOff));
assert(diag.enabled&&diag.complete&&diag.droppedRecords==0&&diag.seenRecords>1&& ...
    diag.seenRecords==numel(diag.records)&&isequal(diag.residualScale,plain.residualScale));
assert(cap.seenRecords==diag.seenRecords&&numel(cap.records)==1&&cap.droppedRecords==diag.seenRecords-1&& ...
    ~cap.complete&&isequaln(cap.records,diag.records(1)));
assert(isempty(zero.records)&&zero.droppedRecords==diag.seenRecords&&~zero.complete);
scale=diag.residualScale;maxError=0;
for k=1:numel(diag.records)
    item=diag.records(k);assert(strcmp(item.status,'accepted')&&isempty(item.errorIdentifier));
    [rate,~]=axisymmetricNavierStokesResidual(item.candidateQ,mesh,gas,boundary);
    actual=max(abs(item.candidateQ./scale-Q./scale-item.dt*(rate./scale)),[],'all');
    maxError=max(maxError,abs(actual-item.candidateDefect));
    assert(abs(actual-item.candidateDefect)<64*eps(max(1,actual))&& ...
        (actual<=item.armijoLimit||actual<=item.target));
    assert(max(abs(item.candidateQ./scale-item.iterateQ./scale-item.alpha*item.incrementScaled),[],'all')<64*eps&& ...
        item.step==1&&item.attempt==0&&item.searchIndex==0&&item.jacobianRefreshedThisIteration);
end
refusals=struct('label',{},'identifier',{});
reject('unknown-diagnostic-option',@()advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options,struct('sink',1)));
reject('nonlogical-enabled',@()advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options,struct('enabled',1)));
reject('negative-capacity',@()advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options,struct('maxRecords',-1)));
reject('noninteger-capacity',@()advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options,struct('maxRecords',1.2)));
reject('nan-capacity',@()advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options,struct('maxRecords',NaN)));
bad=boundary;bad.jMin.temperature=NaN;
reject('input-error-not-swallowed',@()advanceAxisymmetricSteadyObserved(Q,mesh,gas,bad,options,struct('enabled',true)));
report=struct('status','passed-closed-four-quantity-observational-diagnostics', ...
    'matlabEngineProcessID',feature('getpid'),'completeResultAndCountsIdentical',true, ...
    'recordCount',diag.seenRecords,'candidateDefectReplayMaximumError',maxError, ...
    'truncationAndZeroCapacityExplicit',true,'refusals',refusals, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'newton_diagnostics_evidence.mat'),'report','Q','mesh','gas','boundary','options','plain','diag','cap','zero');
fid=fopen(fullfile(folder,'newton_diagnostics_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
    function reject(label,action)
        try,action();catch exception
            assert(startsWith(exception.identifier,{'AxisSteady:Diagnostics','AxisFV:', 'MATLAB:expected', ...
                'MATLAB:validators:','MATLAB:invalidType'}),'AxisDiagnosticTest:Unexpected','非合同错误%s。',exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return;
        end
        error('AxisDiagnosticTest:NoRefusal','非法输入%s没有拒绝。',label);
    end
end
