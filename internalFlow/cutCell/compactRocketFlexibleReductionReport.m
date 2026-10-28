function summary=compactRocketFlexibleReductionReport(detail)
%COMPACTROCKETFLEXIBLEREDUCTIONREPORT 可变内层仅导出标量账，真实V/Z及完整内解只留MAT。
summary=rmfield(detail,{'krylov','preconditionerReports'});k=detail.krylov;
if isfield(k,'cycles')
    for j=1:numel(k.cycles)
        c=k.cycles{j};remove=intersect(fieldnames(c),{'V','Z','retainedU','retainedC','xStart','rStart','weights','candidate','actualRemaining','nextU','nextC'});
        if ~isempty(remove),c=rmfield(c,remove);end;k.cycles{j}=c;
    end
end
k.basisArraysOmittedFromJSON=true;summary.krylov=k;summary.preconditionerReports=cell(size(detail.preconditionerReports));
for j=1:numel(detail.preconditionerReports)
    p=detail.preconditionerReports{j};q=rmfield(p,{'firstFReport','secondFReport'});
    for key={'firstFReport','secondFReport'}
        r=p.(key{1});q.(key{1})=struct('relativeFullFResidual',r.relativeFullFResidual,'fullFInternalTargetMet',r.fullFInternalTargetMet, ...
            'nativeIterations',r.nativeIterations,'nativeFlag',r.nativeFlag,'operatorApplications',r.operatorApplications, ...
            'adjointIdentityRelativeDifference',r.adjointIdentityRelativeDifference);
    end
    summary.preconditionerReports{j}=q;
end
assert(numel(jsonencode(summary))<128*1024,'FlexibleReduction:Summary','可变空间摘要超128KiB。');
end
