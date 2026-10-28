function report=testRocketMomentumSummary(folder)
%TESTROCKETMOMENTUMSUMMARY 检查真实已保存迭代报告的数值身份与摘要体积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'MomentumSummaryTest:Overwrite','不得覆盖摘要测试。');mkdir(folder);
source=fullfile(root,'outputs','v83_validation','momentumlineunit04','momentum_lines_unit_evidence.mat');s=load(source,'cases');bytes=zeros(2,1);
for j=1:2
    original=s.cases{j}.detail;summary=compactRocketMomentumKrylovReport(original);raw=jsonencode(summary);bytes(j)=numel(raw);
    assert(summary.relativeAResidual==original.relativeAResidual&&summary.AApplications==original.AApplications&&summary.lineInverseApplications==original.lineInverseApplications);
    assert(isequal(summary.krylov.residualHistory,original.krylov.residualHistory)&&summary.krylov.flag==original.krylov.flag);
    assert(~contains(raw,'"V":')&&~contains(raw,'"Z":')&&~contains(raw,'"candidate":')&&bytes(j)<128*1024);
end
report=struct('status','passed-compact-momentum-summary-values-counts-history-and-no-basis-arrays','summaryBytes',bytes,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'momentum_summary_unit_evidence.mat'),'report','source','-v7');fprintf('%s\n',jsonencode(report));
end
