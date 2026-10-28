function report=testRocketPressureCandidates(folder,evidence)
%TESTROCKETPRESSURECANDIDATES 只读重演已保存253个真实Q11后续候选，不重复长Newton试验。
% 模块一：核对真实当前源码、内部初猜来源与完整工况，不修改任何接受链。
% 模块二：逐条调用完整五量RANS，比较真实缺陷或同一状态错误identifier。
% 模块三：结果只证明已保存候选可重演，不能宣称HLLE已求解或全面积物理通过。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v72_validation','candidate_replay');end
if nargin<2,evidence=fullfile(root,'outputs','v72_validation','continuation_trial','continuation_evidence.mat');end
assert(~exist(folder,'dir'),'PressureCandidateTest:Overwrite','不可覆盖已有测试证据。');
d=load(evidence,'report','c','inputs','trial','diagnostics');
assert(d.report.complete&&d.report.candidateCount==253&&d.trial.acceptedSteps==1&&d.trial.rejectedAttempts==3);
for k=1:numel(d.inputs.sourceBindings)
    b=d.inputs.sourceBindings(k);assert(strcmp(fileSHA(fullfile(root,b.path)),b.sha256));
end
assert(strcmp(d.inputs.observerSHA256,fileSHA(fullfile(root,'regionalFlow','axisymmetric','advanceAxisymmetricSteadyObserved.m'))));
fresh=createRocketAxisymmetricPressureCase(d.c.options);assert(isequaln(fresh,d.c));
c=d.c;scale=d.trial.residualScale;U0=d.trial.initialQ./scale;maxError=0;failures=0;
for k=1:numel(d.diagnostics.records)
    r=d.diagnostics.records(k);identifier='';
    try,rate=axisymmetricRANSResidual(r.candidateQ,c.mesh,c.gas,c.boundary);
    catch exception,identifier=exception.identifier;end
    if strcmp(r.status,'state-failure')
        assert(strcmp(identifier,r.errorIdentifier)&&strcmp(identifier,'AxisFV:HLLC'));failures=failures+1;
    else
        assert(isempty(identifier));G=r.candidateQ./scale-U0-r.dt*(rate./scale);actual=max(abs(G),[],'all');
        budget=512*eps*max(abs(r.candidateQ./scale)+abs(U0)+abs(r.dt*(rate./scale)),[],'all');
        maxError=max(maxError,abs(actual-r.candidateDefect));assert(abs(actual-r.candidateDefect)<=budget);
    end
end
assert(failures==226);
report=struct('status','passed-all-253-saved-Q11-pressure-source-continuation-candidates', ...
    'matlabEngineProcessID',feature('getpid'),'sourceEvidenceSHA256',fileSHA(evidence), ...
    'candidateCount',253,'actualHLLCFailuresReplayed',failures,'maximumNormReplayDifference',maxError, ...
    'originalChainModified',false,'HLLEActuallyAdvanced',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);fid=fopen(fullfile(folder,'candidate_replay_summary.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
