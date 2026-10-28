function enriched=enrichRocketBlockGraphLevel1(packed)
%ENRICHROCKETBLOCKGRAPHLEVEL1 对原完整五量块图作一级符号填充，不改变原C。
% 模块一：沿同一代数排序用P+strictLower(P)*strictUpper(P)构造一级消元模式。
% 模块二：所有原块逐值保留，新增块初值为零；不改变控制体、方程或边界。
% 模块三：先核对符号数组，再核对精确数值存储预算；估算不代表进程峰值。
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&& ...
    ~isfield(packed,'enrichmentReport'),'BlockFill:Kind','需要原始未因子化块C，不能二次填充。');
n=packed.cellCount;oldCount=numel(packed.column);timer=tic;
oldOwners=repelem((1:n)',diff(double(packed.rowPointer)));
P=spones(sparse(oldOwners,double(packed.column),1,n,n));lower=tril(P,-1);upper=triu(P,1);
productBound=full(sum(lower,1)) * full(sum(upper,2));upperCount=oldCount+productBound;
symbolicEstimate=64*(upperCount+2*oldCount)+128*n;
assert(symbolicEstimate<=packed.options.maximumConstructionBytes,'BlockFill:Budget','符号构造估算超预算。');
pattern=spones(P+lower*upper);clear P lower upper;
blockCount=nnz(pattern);assert(blockCount<2^32);
persistentBytes=blockCount*(25*8+4)+(n+1)*4+n*4+n*25*8+2*n*8;
% 同时包括输入原值、输出及因子写入副本、符号索引和批次临时量；库工作区另计。
constructionBytes=oldCount*25*8+2*blockCount*25*8+48*(blockCount+oldCount)+n*200+4096*25*8*3;
assert(persistentBytes<=packed.options.maximumPersistentBytes&&constructionBytes<=packed.options.maximumConstructionBytes, ...
    'BlockFill:Budget','一级填充数值/因子或构造数组估算超预算，不裁填充块。');
[columns,owners]=find(pattern');clear pattern;
rowCounts=accumarray(owners,1,[n 1]);rowPointer=uint32([1;1+cumsum(rowCounts)]);
position=sparse(columns,owners,(1:blockCount)',n,n);diagonalPosition=uint32(full(diag(position)));
oldToNew=uint32(full(position(sub2ind([n n],double(packed.column),oldOwners))));clear position oldOwners owners;
assert(all(oldToNew>0)&&numel(unique(oldToNew))==oldCount,'BlockFill:Original','一级填充丢失原块。');
values=zeros(5,5,blockCount);
for first=1:4096:oldCount
    ids=first:min(first+4095,oldCount);values(:,:,oldToNew(ids))=packed.values(:,:,ids);
end
enriched=packed;enriched.values=values;enriched.column=uint32(columns);enriched.rowPointer=rowPointer;enriched.diagonalPosition=diagonalPosition;
enriched.enrichmentReport=struct('method','symbolic-level-one-fill-P-plus-strictLowerP-times-strictUpperP', ...
    'originalStoredBlocks',oldCount,'enrichedStoredBlocks',blockCount,'addedInitiallyZeroBlocks',blockCount-oldCount, ...
    'originalToEnrichedPosition',oldToNew,'originalOrderingKept',true,'allOriginalBlockValuesKept',true, ...
    'minimumBlocksPerRow',min(rowCounts),'maximumBlocksPerRow',max(rowCounts),'meanBlocksPerRow',mean(rowCounts), ...
    'symbolicProductCountUpperBound',productBound,'symbolicArrayConstructionEstimateBytes',symbolicEstimate, ...
    'persistentFactorStorageUpperBoundBytes',persistentBytes,'algorithmArrayConstructionEstimateBytes',constructionBytes, ...
    'sourceOperatorAndLibraryWorkspaceExcluded',true,'processPeakMemoryClaimed',false,'enrichmentSeconds',toc(timer), ...
    'allOriginalCellsAndFiveEquationsKept',true,'physicalEquationChanged',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
fprintf('BLOCK_FILL1 cells=%d original=%d enriched=%d persistentBound=%.0f constructionEstimate=%.0f\n',n,oldCount,blockCount,persistentBytes,constructionBytes);
end
