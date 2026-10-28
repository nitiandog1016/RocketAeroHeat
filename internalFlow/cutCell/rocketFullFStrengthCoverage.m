function coverage=rocketFullFStrengthCoverage(graph,factor,oldFineGroupID)
%ROCKETFULLFSTRENGTHCOVERAGE 统计新F预条件纳入的旧跨组连接，显式统一向量形状。
% 模块一：旧组编号无论行/列输入均转列，避免隐式扩展成平方规模数组。
% 模块二：只统计预条件块；原完整F矩阵及物理面积不因统计改变。
% 模块三：CV对角替代组只计原点对角块，不把被拒绝整块继续算作纳入。
nf=graph.cellCount;oldFineGroupID=double(oldFineGroupID(:));
assert(numel(oldFineGroupID)==nf&&all(isfinite(oldFineGroupID))&&all(oldFineGroupID>=1&oldFineGroupID==floor(oldFineGroupID))&&factor.cellCount==nf,'FullFCoverage:Input','旧F组编号或新因子尺寸不符。');
partition=factor.partition;owners=repelem((1:nf)',diff(double(graph.rowPointer)));rowCV=graph.order(owners);colCV=graph.order(graph.column);
rowCV=rowCV(:);colCV=colCV(:);inside=partition.groupID(rowCV)==partition.groupID(colCV);
fallback=false(partition.patchCount,1);
for k=1:numel(factor.buckets),bucket=factor.buckets{k};fallback(bucket.patchIDs)=bucket.fallback;end
effective=inside&(~fallback(partition.groupID(rowCV))|rowCV==colCV);oldCross=oldFineGroupID(rowCV)~=oldFineGroupID(colCV);
assert(isequal(size(effective),[numel(graph.column) 1])&&isequal(size(oldCross),size(effective)),'FullFCoverage:Shape','块标记必须为同尺寸列向量。');
coverage=struct('allOriginalFFBlocks',numel(graph.column),'oldWithinOriginalPatchFFBlocks',nnz(~oldCross),'newWithinStrengthPatchBlocks',nnz(inside), ...
    'newEffectiveLocalBlocksAfterFallback',nnz(effective),'oldCrossPatchBlocks',nnz(oldCross),'newlyIncludedOldCrossPatchBlocks',nnz(effective&oldCross), ...
    'countsAreAlgebraicBlocksNotPhysicalArea',true,'oldApproximateSchurIdentityClaimed',false);
end
