function summary=compactRocketPatchHierarchyReport(detail)
%COMPACTROCKETPATCHHIERARCHYREPORT 只导出层次迭代的标量、计数和历史，禁止V/Z进入JSON。
summary=struct;keys={'method','options','relativeCResidual','solverCApplications','VcycleApplications', ...
    'preconditionerCApplicationsByLevel','preconditionerLocalApplicationsByLevel','preconditionerTerminalApplicationsByLevel', ...
    'totalRootCApplications','elapsedSeconds','baseHierarchyFixedLinear','innerSolveFixedLinearMapClaimed', ...
    'originalNSAcceptanceClaimed','physicalSolution','newActualNSEvaluations','newProductionSolvedAreaM2'};
for key=keys,summary.(key{1})=detail.(key{1});end
k=detail.krylov;small=struct;keys={'method','flag','stopReason','options','arnoldiSteps','modelApplications','preconditionerApplications', ...
    'basisStorageUpperBoundBytes','maximumBasisBytes','maximumBasisOrthogonalityError','absoluteResidual','relativeResidual','residualHistory'};
for key=keys,small.(key{1})=k.(key{1});end
small.cycles=cell(size(k.cycles));
for j=1:numel(k.cycles)
    entry=struct;keys={'cycle','arnoldiSteps','startResidualNorm','actualResidualNorm','arnoldiRelationRelativeDifference','basisOrthogonalityError','modelApplications','preconditionerApplications'};
    for key=keys,entry.(key{1})=k.cycles{j}.(key{1});end;small.cycles{j}=entry;
end
small.basisArraysOmittedFromJSON=true;summary.krylov=small;assert(numel(jsonencode(summary))<128*1024,'PatchHierarchy:Summary','层次摘要超128KiB。');
end
