function partition=buildRocketStrongPatches(packed,maximumPatchCells)
%BUILDROCKETSTRONGPATCHES 按完整五量矩阵块的相对强连接构造有界连通局部组。
% 模块一：权重为||Cij||F/||Cii||F的双向最大值，仅用于分组，不改原系数。
% 模块二：按强度和原CV编号确定种子/增长顺序；所有CV恰好进入一个组。
% 模块三：局部组只用于预条件，不合并物理控制体或代替已有粗空间。
assert(strcmp(packed.kind,'complete-five-conservation-block-CSR-unshifted-steady-graph'),'StrongPatches:Kind','需要原完整五量C。');
validateattributes(maximumPatchCells,{'numeric'},{'scalar','real','finite','integer','positive','<=',16});n=packed.cellCount;nb=numel(packed.column);
owners=repelem((1:n)',diff(double(packed.rowPointer)));magnitude=sqrt(reshape(sum(packed.values.^2,[1 2]),[],1));diagonal=magnitude(packed.diagonalPosition);
assert(all(isfinite(magnitude))&&all(diagonal>0),'StrongPatches:Weight','原块范数非有限或对角零，不能制造连接。');
weight=magnitude./diagonal(owners);weight(packed.diagonalPosition)=0;G=sparse(owners,double(packed.column),weight,n,n);G=max(G,G');
% 图转回原CV顺序，编号用于可追溯确定性，不使用随机种子。
G=G(packed.inverseOrder,packed.inverseOrder);[column,row,weights]=find(G');clear G magnitude weight owners;
ptr=[1;1+cumsum(accumarray(row,1,[n 1]))];score=accumarray(row,weights,[n 1],@max,0);[~,seeds]=sortrows([-score (1:n)']);clear row score;
groupID=zeros(n,1,'uint32');count=0;
for seed=seeds'
    if groupID(seed)~=0,continue;end
    count=count+1;groupID(seed)=count;member=seed;pool=[];poolWeight=[];
    for step=2:maximumPatchCells
        ids=ptr(member):ptr(member+1)-1;pool=[pool;column(ids)];poolWeight=[poolWeight;weights(ids)]; %#ok<AGROW>
        keep=groupID(pool)==0;pool=pool(keep);poolWeight=poolWeight(keep);if isempty(pool),break;end
        [candidate,~,index]=unique(pool);strength=accumarray(index,poolWeight,[],@max);[~,ranking]=sortrows([-strength candidate]);member=candidate(ranking(1));groupID(member)=count;
    end
end
assert(all(groupID>0),'StrongPatches:Coverage','有CV遗漏。');sizes=accumarray(double(groupID),1,[count 1]);assert(all(sizes<=maximumPatchCells));
[~,members]=sortrows([double(groupID) (1:n)']);groupPointer=uint32([1;1+cumsum(sizes)]);
partition=struct('kind','bounded-connected-five-field-strength-patches','cellCount',n,'patchCount',count,'maximumPatchCells',maximumPatchCells, ...
    'groupID',groupID,'groupPointer',groupPointer,'members',uint32(members),'sizes',uint8(sizes), ...
    'constructionReport',struct('weight','symmetric-max-of-full-five-field-block-Frobenius-to-row-diagonal-Frobenius', ...
    'minimumPatchCells',min(sizes),'maximumPatchCells',max(sizes),'meanPatchCells',mean(sizes),'allOriginalCellsKept',true,'physicalGridChanged',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0));
end
