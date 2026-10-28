function value=solveRocketSpatialFBlocksTranspose(factor,rhs)
%SOLVEROCKETSPATIALFBLOCKSTRANSPOSE 作用已认证固定局部F映射的转置，不重新求逆。
% 模块一：使用缓存中的实际逆，包括已有明确CV对角替代的实际映射。
% 模块二：交换每页逆的行列，固定收集/散回全部五量，构成M_F的代数伴随。
assert(strcmp(factor.kind,'fixed-within-original-patch-five-field-F-inverse')&& ...
    isequal(size(rhs),[factor.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'SpatialF:Input','固定F逆转置或五量右端不符。');
value=zeros(size(rhs));
for k=1:numel(factor.buckets)
    bucket=factor.buckets{k};count=numel(bucket.groupIDs);if count==0,continue;end
    b=rhs(bucket.members(:),:);pages=reshape(permute(reshape(b,k,count,5),[3 1 2]),5*k,1,count);
    x=pagemtimes(permute(bucket.inverse,[2 1 3]),pages);
    value(bucket.members(:),:)=reshape(permute(reshape(x,5,k,count),[2 3 1]),k*count,5);
end
assert(all(isfinite(value),'all'),'SpatialF:Nonfinite','F局部逆转置作用非有限。');
end
