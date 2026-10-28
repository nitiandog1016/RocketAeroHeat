function summary=compactRocketSpatialReductionReport(detail)
%COMPACTROCKETSPATIALREDUCTIONREPORT 空间内部摘要仅含标量、计数和历史，V/Z只留MAT。
summary=struct;keys={'method','options','relativeCResidual','CApplications','spatialReductionApplications','FLocalApplications', ...
    'coarseLocalApplications','rectangularCouplingApplications','elapsedSeconds','baseSpatialMapFixedLinear', ...
    'innerSolveFixedLinearMapClaimed','originalNSAcceptanceClaimed','physicalSolution','newActualNSEvaluations','newProductionSolvedAreaM2'};
for key=keys,summary.(key{1})=detail.(key{1});end
k=detail.krylov;small=struct;keys={'method','flag','stopReason','options','arnoldiSteps','modelApplications','preconditionerApplications', ...
    'basisStorageUpperBoundBytes','maximumBasisBytes','maximumBasisOrthogonalityError','absoluteResidual','relativeResidual','residualHistory'};
for key=keys,small.(key{1})=k.(key{1});end
small.cycles=cell(size(k.cycles));for j=1:numel(k.cycles)
    entry=struct;keys={'cycle','arnoldiSteps','startResidualNorm','actualResidualNorm','arnoldiRelationRelativeDifference','basisOrthogonalityError','modelApplications','preconditionerApplications'};
    for key=keys,entry.(key{1})=k.cycles{j}.(key{1});end;small.cycles{j}=entry;
end
small.basisArraysOmittedFromJSON=true;summary.krylov=small;assert(numel(jsonencode(summary))<128*1024,'SpatialReduction:Summary','摘要超过128KiB。');
end
