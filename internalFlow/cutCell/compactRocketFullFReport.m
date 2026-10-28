function summary=compactRocketFullFReport(detail)
%COMPACTROCKETFULLFREPORT 仅保存完整F求解的计数、标量与历史，不在JSON复制V/Z。
summary=struct;keys={'method','options','relativeFullFResidual','FApplications','localFApplications','elapsedSeconds', ...
    'fullFOperatorIncludesCrossPatchBlocks','baseLocalFMapFixedLinear','innerSolveFixedLinearMapClaimed', ...
    'oldApproximateSchurIdentityReplaced','originalCOrNSConvergenceClaimed','physicalSolution','newActualNSEvaluations','newProductionSolvedAreaM2'};
for key=keys,summary.(key{1})=detail.(key{1});end
k=detail.krylov;small=struct;keys={'method','flag','stopReason','options','arnoldiSteps','modelApplications','preconditionerApplications', ...
    'basisStorageUpperBoundBytes','maximumBasisBytes','maximumBasisOrthogonalityError','absoluteResidual','relativeResidual','residualHistory'};
for key=keys,small.(key{1})=k.(key{1});end
small.cycles=cell(size(k.cycles));for j=1:numel(k.cycles)
    entry=struct;keys={'cycle','arnoldiSteps','startResidualNorm','actualResidualNorm','arnoldiRelationRelativeDifference','basisOrthogonalityError','modelApplications','preconditionerApplications'};
    for key=keys,entry.(key{1})=k.cycles{j}.(key{1});end;small.cycles{j}=entry;
end
small.basisArraysOmittedFromJSON=true;summary.krylov=small;assert(numel(jsonencode(summary))<128*1024,'FullFKrylov:Summary','完整F摘要超过128KiB。');
end
