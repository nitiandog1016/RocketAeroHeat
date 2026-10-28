function split=buildRocketSpatialSplit(packed,partition,coarseOriginalCVs)
%BUILDROCKETSPATIALSPLIT 每个既有局部组保留一个原CV全部五量，其余CV为F。
% 模块一：C/F只改变预条件排列，完整物理网格及所有五量方程不变。
% 模块二：默认保留相对非对角耦合最强的CV；外部给定代表须每组恰一个。
% 模块三：所有原CV恰好属于C或F，F按原组及原编号排列，支持单点组。
n=packed.cellCount;ng=partition.patchCount;
assert(partition.cellCount==n&&isequal(sort(double(partition.members)),(1:n)'),'SpatialSplit:Input','原C与完整局部组不符。');
if nargin<3
    owners=repelem((1:n)',diff(double(packed.rowPointer)));magnitude=sqrt(reshape(sum(packed.values.^2,[1 2]),[],1));
    diagonal=magnitude(packed.diagonalPosition);assert(all(diagonal>0)&&all(isfinite(magnitude)),'SpatialSplit:Scale','原对角或块范数非法。');
    strength=magnitude./diagonal(owners);strength(packed.diagonalPosition)=0;ordered=accumarray(owners,strength,[n 1],@max,0);score=zeros(n,1);score(packed.order)=ordered;
    coarseOriginalCVs=zeros(ng,1);
    for g=1:ng
        ids=double(partition.members(double(partition.groupPointer(g)):double(partition.groupPointer(g+1))-1));
        [~,pick]=max(score(ids));coarseOriginalCVs(g)=ids(pick);
    end
    policy='maximum-relative-offdiagonal-block-coupling-with-lowest-original-ID-tie';
else
    coarseOriginalCVs=double(coarseOriginalCVs(:));policy='fixed-original-CV-representatives-provided-by-source-bound-preflight';
end
assert(numel(coarseOriginalCVs)==ng&&all(isfinite(coarseOriginalCVs))&&all(coarseOriginalCVs==floor(coarseOriginalCVs))&& ...
    all(coarseOriginalCVs>=1&coarseOriginalCVs<=n)&&numel(unique(coarseOriginalCVs))==ng&& ...
    isequal(double(partition.groupID(coarseOriginalCVs)),(1:ng)'),'SpatialSplit:Representatives','C代表必须按原组顺序每组恰一个。');
coarseIndex=zeros(n,1,'uint32');coarseIndex(coarseOriginalCVs)=uint32(1:ng);nf=n-ng;
fineOriginalCVs=partition.members(coarseIndex(partition.members)==0);fineIndex=zeros(n,1,'uint32');fineIndex(fineOriginalCVs)=uint32(1:nf);
counts=double(partition.sizes(:))-1;pointer=[1;1+cumsum(counts)];fineGroupID=partition.groupID(fineOriginalCVs);
start=repelem(pointer(1:end-1),counts);finePosition=uint8((1:nf)'-start(:)+1);
assert(all((coarseIndex>0)+(fineIndex>0)==1)&&isequal(sort([double(fineOriginalCVs);coarseOriginalCVs]),(1:n)'),'SpatialSplit:Coverage','C/F遗漏或重复原CV。');
split=struct('kind','all-five-field-original-CV-spatial-CF-split','cellCount',n,'coarseCells',ng,'fineCells',nf, ...
    'coarseOriginalCVs',uint32(coarseOriginalCVs),'fineOriginalCVs',uint32(fineOriginalCVs),'coarseIndex',coarseIndex,'fineIndex',fineIndex, ...
    'originalGroupID',partition.groupID,'fineGroupID',fineGroupID,'finePosition',finePosition,'groupFCounts',uint8(counts),'groupFPointer',uint32(pointer), ...
    'representativePolicy',policy,'physicalGridChanged',false,'allOriginalCVsAndFiveFieldsKept',true);
end
