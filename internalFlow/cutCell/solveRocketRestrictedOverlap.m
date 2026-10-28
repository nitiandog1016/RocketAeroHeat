function value=solveRocketRestrictedOverlap(factor,rhs)
%SOLVEROCKETRESTRICTEDOVERLAP 扩展组参与求解，仅散回原核心五量，全部F点恰一次输出。
% 模块一：核心是扩展成员前缀，邻点解不散回其他核心，不平均物理状态。
% 模块二：所有缓存固定，无右端选择/早停，形成明确的固定线性预条件映射。
assert(strcmp(factor.kind,'fixed-five-field-overlap-inverse-return-original-core-rows-only')&& ...
    isequal(size(rhs),[factor.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'OverlapF:Input','重叠因子或全部五量右端不符。');
value=zeros(size(rhs));
for id=1:numel(factor.buckets)
    bucket=factor.buckets{id};count=numel(bucket.groupIDs);if count==0,continue;end
    k=bucket.overlapSize;m=bucket.coreSize;b=rhs(bucket.members(:),:);
    pages=reshape(permute(reshape(b,k,count,5),[3 1 2]),5*k,1,count);solved=pagemtimes(bucket.inverse,pages);
    core=solved(1:5*m,:,:);value(bucket.coreMembers(:),:)=reshape(permute(reshape(core,5,m,count),[2 3 1]),m*count,5);
end
assert(all(isfinite(value),'all'),'OverlapF:Nonfinite','重叠核心回代非有限，拒绝使用。');
end
