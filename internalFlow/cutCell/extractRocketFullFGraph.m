function [graph,report]=extractRocketFullFGraph(packed,split,options)
%EXTRACTROCKETFULLFGRAPH 从完整原C逐系数提取全部五量A_FF，保留跨原组F-F连接。
% 模块一：F子系统的控制体顺序严格沿用已验split，C值通过右端处理。
% 模块二：没有删除物理CV，原完整C仍是外层验收方程；此处只取精确子块。
% 模块三：保存所有原F-F块及原CV映射，预算不足明确拒绝，不裁小系数。
if nargin<3,options=struct;end
defaults=struct('maximumGraphBytes',384*1024^2,'maximumConstructionBytes',1536*1024^2);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'FullFGraph:Options','未知完整F图选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&& ...
    strcmp(split.kind,'all-five-field-original-CV-spatial-CF-split')&&packed.cellCount==split.cellCount&&split.fineCells>0,'FullFGraph:Input','需要原C及非空完整F映射。');
n=packed.cellCount;nf=split.fineCells;nb=numel(packed.column);
owners=repelem((1:n)',diff(double(packed.rowPointer)));rowCV=packed.order(owners);colCV=packed.order(packed.column);
rowF=split.fineIndex(rowCV);colF=split.fineIndex(colCV);selected=find(rowF>0&colF>0);m=numel(selected);
bound=208*m+32*nf+65536;construction=200*nb+3*208*m+96*nb+96*n+65536;
assert(bound<=options.maximumGraphBytes&&construction<=options.maximumConstructionBytes,'FullFGraph:Budget','全部F-F子块超预算，不删跨组连接。');
[pairs,sorting]=sortrows([rowF(selected) colF(selected)]);selected=selected(sorting);
assert(all(pairs(:,1)>=1&pairs(:,1)<=nf&pairs(:,2)>=1&pairs(:,2)<=nf),'FullFGraph:Mapping','原F编号超出范围。');
diagonal=find(pairs(:,1)==pairs(:,2));assert(numel(diagonal)==nf,'FullFGraph:Diagonal','完整F缺少原点对角位置。');
same=split.originalGroupID(rowCV(selected))==split.originalGroupID(colCV(selected));
report=struct('method','exact-all-original-five-field-spatial-F-submatrix-extraction','originalCells',n,'FCells',nf,'originalBlocks',nb, ...
    'allOriginalFFBlocks',m,'withinOriginalPatchFFBlocks',nnz(same),'crossPatchFFBlocks',nnz(~same), ...
    'graphStorageUpperBoundBytes',bound,'constructionArrayEstimateBytes',construction,'processPeakMemoryClaimed',false, ...
    'allFFCoefficientsKept',true,'originalPhysicalGridChanged',false,'originalFullCReplaced',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
graph=struct('kind','complete-five-conservation-block-CSR-unshifted-steady-graph','cellCount',nf, ...
    'rowPointer',uint32([1;1+cumsum(accumarray(double(pairs(:,1)),1,[nf 1]))]),'column',uint32(pairs(:,2)), ...
    'diagonalPosition',uint32(diagonal),'order',(1:nf)','inverseOrder',(1:nf)','values',packed.values(:,:,selected), ...
    'systemScope','exact-spatial-F-submatrix-A_FF-not-original-full-C','originalPhysicalCVs',split.fineOriginalCVs,'sourceOriginalCellCount',n);
storage=whos('graph');assert(storage.bytes<=bound&&storage.bytes<=options.maximumGraphBytes,'FullFGraph:Budget','实际全部F图超过预算或上界。');
fprintf('FULL_F_GRAPH cells=%d blocks=%d cross=%d bytes=%.0f\n',nf,m,nnz(~same),storage.bytes);
end
