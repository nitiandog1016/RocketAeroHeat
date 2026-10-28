function report=testRocketKrylovHistoryWork(folder,sourceFolder)
%TESTROCKETKRYLOVHISTORYWORK 用真实失败尾部计数校核下一段偏移，不运行新的NS更新。
% 模块一：读取真实接受历史与实际总调用数，要求真实末尾存在未接受工作。
% 模块二：合成有接受/零接受接续，检查尾部工作既不丢失，也不伪装为新接受状态。
% 模块三：无失败尾部时须逐字段复现原合并表达式；损坏合同明确拒绝。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<2,sourceFolder=fullfile(root,'outputs','v83_validation','krylov05');end
assert(~exist(folder,'dir'),'HistoryWorkTest:Overwrite','不得覆盖计数验证。');mkdir(folder);
s=load(fullfile(sourceFolder,'wall_velocity_krylov_evidence.mat'),'report','displayHistory');
previous=s.displayHistory;priorTotal=s.report.totalNewEquationResidualEvaluations;
assert(priorTotal>previous.residualEvaluations(end),'HistoryWorkTest:Source','本测试需要真实未接受尾部。');
local=struct;
for name=fieldnames(previous)',key=name{1};local.(key)=[previous.(key)(end),previous.(key)(end)];end
local.iteration=[0 1];local.residualEvaluations=[1 4];
[merged,work]=mergeRocketKrylovHistory(previous,priorTotal,local,7);
assert(merged.iteration(end)==previous.iteration(end)+1&&merged.residualEvaluations(end)==priorTotal+4);
assert(work.totalActualResidualEvaluations==priorTotal+7&&work.residualEvaluationsAfterLastAcceptedHistory==3);
for name=fieldnames(previous)',key=name{1};assert(isequaln(merged.(key)(1:end-1),previous.(key)));end
empty=struct;for name=fieldnames(local)',key=name{1};empty.(key)=local.(key)(1);end
[unchanged,zeroWork]=mergeRocketKrylovHistory(previous,priorTotal,empty,3);
assert(isequaln(unchanged,previous)&&zeroWork.totalActualResidualEvaluations==priorTotal+3&& ...
    zeroWork.residualEvaluationsAfterLastAcceptedHistory==priorTotal+3-previous.residualEvaluations(end));
legacyTotal=previous.residualEvaluations(end);legacy=previous;
for name=fieldnames(local)'
    key=name{1};values=local.(key)(2:end);
    if strcmp(key,'iteration'),values=values+previous.iteration(end);end
    if strcmp(key,'residualEvaluations'),values=values+legacyTotal;end
    legacy.(key)=[legacy.(key),values];
end
actual=mergeRocketKrylovHistory(previous,legacyTotal,local,4);assert(isequaln(actual,legacy));
bad=local;bad.iteration=[1 2];rejections=0;
for item={{previous,previous.residualEvaluations(end)-1,local,7,'KrylovHistory:Work'}, ...
          {previous,priorTotal,local,3,'KrylovHistory:Work'}, ...
          {previous,priorTotal,bad,7,'KrylovHistory:Iteration'}}
    x=item{1};caught=false;try,mergeRocketKrylovHistory(x{1},x{2},x{3},x{4});catch ex,caught=strcmp(ex.identifier,x{5});end
    assert(caught);rejections=rejections+1;
end
report=struct('status','passed-real-unaccepted-tail-work-offset-and-zero-acceptance-history-contract', ...
    'sourcePriorAcceptedIterations',previous.iteration(end),'sourceActualTotalCalls',priorTotal, ...
    'sourceLastAcceptedCalls',previous.residualEvaluations(end),'sourceUnacceptedTailCalls',priorTotal-previous.residualEvaluations(end), ...
    'syntheticAcceptedContinuationTotalCalls',work.totalActualResidualEvaluations, ...
    'syntheticZeroAcceptedContinuationTotalCalls',zeroWork.totalActualResidualEvaluations, ...
    'legacyCompleteHistoryExactlyReproduced',true,'invalidInputRejections',rejections, ...
    'newNSEvaluations',0,'syntheticHistoryIsPhysicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'history_work_evidence.mat'),'report','previous','priorTotal','local','merged','work','empty','zeroWork','-v7');
fid=fopen(fullfile(folder,'history_work_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
