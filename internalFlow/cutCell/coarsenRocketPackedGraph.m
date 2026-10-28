function [coarse,report]=coarsenRocketPackedGraph(packed,partition,maximumConstructionBytes)
%COARSENROCKETPACKEDGRAPH 直接汇总全部原C系数形成C_c=R*C*P，不丢非零。
% 模块一：P为组常数延拓，R为组大小倒数的等权代数限制；R*P=I。
% 模块二：逐个五量系数累加所有原连接，不靠模型采样或旧粗像投影构造。
% 模块三：粗层仅用于预条件；完整细层C和物理控制体保持原样。
if nargin<3,maximumConstructionBytes=2048*1024^2;end
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph')&&partition.cellCount==packed.cellCount,'PatchCoarse:Input','完整C与分组不符。');
n=packed.cellCount;ng=partition.patchCount;nb=numel(packed.column);
validateattributes(maximumConstructionBytes,{'numeric'},{'scalar','real','finite','positive'});
estimate=nb*(200+80)+n*64+ng*128+65536;
assert(estimate<=maximumConstructionBytes,'PatchCoarse:Budget','原系数聚合构造数组超预算，不能丢连接或CV。');
owners=repelem((1:n)',diff(double(packed.rowPointer)));rowCV=packed.order(owners);columnCV=packed.order(packed.column);
row=double(partition.groupID(rowCV));column=double(partition.groupID(columnCV));
[pairs,~,which]=unique([row column],'rows','sorted');nc=size(pairs,1);sizes=double(partition.sizes(:));
values=zeros(5,5,nc);
for a=1:5,for b=1:5
    coefficients=reshape(packed.values(a,b,:),[],1)./sizes(row);
    values(a,b,:)=reshape(accumarray(which,coefficients,[nc 1]),1,1,[]);
end,end
assert(all(isfinite(values),'all'),'PatchCoarse:Nonfinite','原单位聚合系数非有限。');
diagonal=find(pairs(:,1)==pairs(:,2));assert(numel(diagonal)==ng,'PatchCoarse:Diagonal','粗图缺少原组对角位置。');
coarse=struct('kind','complete-five-conservation-block-CSR-unshifted-steady-graph','cellCount',ng, ...
    'rowPointer',uint32([1;1+cumsum(accumarray(pairs(:,1),1,[ng 1]))]),'column',uint32(pairs(:,2)), ...
    'diagonalPosition',uint32(diagonal),'order',(1:ng)','inverseOrder',(1:ng)','values',values);
report=struct('method','all-original-coefficients-unsmoothed-equal-average-R-C-P', ...
    'fineCells',n,'coarseCells',ng,'fineBlocks',nb,'coarseBlocks',nc,'restriction','equal-member-average-algebraic-only', ...
    'prolongation','five-conservative-field-group-constant','zeroOrSmallCoefficientsDropped',false, ...
    'physicalCellsMerged',false,'originalFineEquationChanged',false,'constructionArrayEstimateBytes',estimate, ...
    'processPeakMemoryClaimed',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
