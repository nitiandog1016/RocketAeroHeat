function value=solveRocketSpatialFBlocks(factor,rhs)
%SOLVEROCKETSPATIALFBLOCKS 批量作用认证的组内F逆，输入输出均为固定F点顺序。
% 模块一：所有原F控制体和五量保留，不按右端选择局部逆。
% 模块二：局部F组之间的耦合由完整外层方程处理，此处不称精确全F解。
assert(strcmp(factor.kind,'fixed-within-original-patch-five-field-F-inverse')&& ...
    isequal(size(rhs),[factor.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'SpatialF:Input','固定F逆或五量右端不符。');
value=zeros(size(rhs));
for k=1:numel(factor.buckets)
    bucket=factor.buckets{k};count=numel(bucket.groupIDs);if count==0,continue;end
    b=rhs(bucket.members(:),:);pages=reshape(permute(reshape(b,k,count,5),[3 1 2]),5*k,1,count);x=pagemtimes(bucket.inverse,pages);
    value(bucket.members(:),:)=reshape(permute(reshape(x,5,k,count),[2 3 1]),k*count,5);
end
assert(all(isfinite(value),'all'),'SpatialF:Nonfinite','F局部逆作用非有限。');
end
