function reordered=reorderRocketBlockGraphByFlow(packed,centroid,flowDirection)
%REORDERROCKETBLOCKGRAPHBYFLOW 按给定来流方向投影重排完整五量块C，原矩阵值保持。
% 模块一：控制体质心沿给定U∞从上游到下游排序，同投影按原CV编号确定顺序。
% 模块二：仅更改代数顺序，不删CV/方程/连接，也不把来流投影当局部流线。
% 模块三：复制值数组前检查预算；既有分解与原NS验收规则保持。
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&& ...
    ~isfield(packed,'enrichmentReport')&&~isfield(packed,'flowReorderingReport'), ...
    'BlockFlowOrder:Kind','需要原完整未填充/未因子化的块图，不允许重复流向排序。');
n=packed.cellCount;nb=numel(packed.column);
assert(isnumeric(centroid)&&isreal(centroid)&&isequal(size(centroid),[n 3])&&all(isfinite(centroid),'all'), ...
    'BlockFlowOrder:Geometry','必须提供与全部原CV一一对应的有限SI质心。');
assert(isnumeric(flowDirection)&&isreal(flowDirection)&&numel(flowDirection)==3&&all(isfinite(flowDirection))&&norm(flowDirection)>0, ...
    'BlockFlowOrder:Direction','给定来流方向必须是非零有限三维向量。');
constructionBytes=3*nb*25*8+nb*96+n*128;
persistentBytes=packed.packingReport.persistentFactorStorageUpperBoundBytes+4*nb;
assert(persistentBytes<=packed.options.maximumPersistentBytes&&constructionBytes<=packed.options.maximumConstructionBytes, ...
    'BlockFlowOrder:Budget','重排值数组与因子/构造估算超预算，不裁原连接。');
timer=tic;direction=flowDirection(:)/norm(flowDirection);projection=centroid*direction;
[~,order]=sortrows([projection (1:n)'],[1 2]);inverse=zeros(n,1);inverse(order)=1:n;
oldOwner=repelem((1:n)',diff(double(packed.rowPointer)));
rows=inverse(packed.order(oldOwner));columns=inverse(packed.order(double(packed.column)));
[pairs,newToOld]=sortrows([rows columns],[1 2]);clear rows columns oldOwner;
rowCounts=accumarray(pairs(:,1),1,[n 1]);rowPointer=uint32([1;1+cumsum(rowCounts)]);
diagonalPosition=uint32(find(pairs(:,1)==pairs(:,2)));assert(numel(diagonalPosition)==n);
values=zeros(5,5,nb);
for first=1:4096:nb
    ids=first:min(first+4095,nb);values(:,:,ids)=packed.values(:,:,newToOld(ids));
end
reordered=packed;reordered.values=values;reordered.rowPointer=rowPointer;reordered.column=uint32(pairs(:,2));
reordered.diagonalPosition=diagonalPosition;reordered.order=order;reordered.inverseOrder=inverse;reordered.options.ordering='prescribed-freestream-projection';
reordered.flowReorderingReport=struct('method','upstream-to-downstream-prescribed-freestream-centroid-projection', ...
    'givenUnitFlowDirection',direction','minimumProjectionM',min(projection),'maximumProjectionM',max(projection), ...
    'projectionTieBreak','original-physical-cell-ID','originalAlgebraicOrdering',packed.options.ordering, ...
    'newToOldBlockPosition',uint32(newToOld),'uniqueStoredBlocks',nb,'allOriginalBlockValuesKept',true, ...
    'allOriginalCellsAndFiveEquationsKept',true,'localSolvedStreamlineClaimed',false,'physicalEquationChanged',false, ...
    'persistentFactorStorageUpperBoundBytes',persistentBytes,'algorithmArrayConstructionEstimateBytes',constructionBytes, ...
    'sourceOperatorAndLibraryWorkspaceExcluded',true,'processPeakMemoryClaimed',false,'reorderingSeconds',toc(timer), ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
fprintf('BLOCK_FLOW_ORDER cells=%d blocks=%d persistentBound=%.0f constructionEstimate=%.0f\n',n,nb,persistentBytes,constructionBytes);
end
