function summary=compactRocketCoupledKrylovReport(detail)
%COMPACTROCKETCOUPLEDKRYLOVREPORT 五量内部报告仅导出标量与历史，完整基留在MAT。
summary=struct;keys={'method','options','relativeCResidual','CApplications','localBlockApplications','elapsedSeconds', ...
    'baseLocalMapFixedLinear','innerSolveFixedLinearMapClaimed','requiresFlexibleOuterMethod','originalNSAcceptanceClaimed','physicalSolution','newActualNSEvaluations','newProductionSolvedAreaM2'};
for key=keys,summary.(key{1})=detail.(key{1});end
k=detail.krylov;small=struct;keys={'method','flag','stopReason','options','arnoldiSteps','modelApplications','preconditionerApplications', ...
    'basisStorageUpperBoundBytes','maximumBasisBytes','maximumBasisOrthogonalityError','absoluteResidual','relativeResidual','residualHistory'};
for key=keys,small.(key{1})=k.(key{1});end
small.cycles=cell(size(k.cycles));for j=1:numel(k.cycles)
    entry=struct;keys={'cycle','arnoldiSteps','startResidualNorm','actualResidualNorm','arnoldiRelationRelativeDifference','basisOrthogonalityError','modelApplications','preconditionerApplications'};
    for key=keys,entry.(key{1})=k.cycles{j}.(key{1});end;small.cycles{j}=entry;
end
small.basisArraysOmittedFromJSON=true;summary.krylov=small;assert(numel(jsonencode(summary))<128*1024,'CoupledSummary:Size','五量内部摘要超限。');
end
