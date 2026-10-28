function value=solveRocketCoupledPatches(factor,rhs)
%SOLVEROCKETCOUPLEDPATCHES 全部CV的五量局部块批量作用，保留原物理顺序。
% 模块一：不作跨组长三角传播；全部跨组耦合由完整原C作用处理。
% 模块二：缓存构成固定线性映射，不按右端选择区域、截零或改局部逆。
assert(strcmp(factor.kind,'fixed-five-conservation-strong-connection-local-block-inverse'),'CoupledPatches:Kind','需要已验五量局部块因子。');
n=factor.cellCount;assert(isequal(size(rhs),[n 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'CoupledPatches:RHS','五量右端非法。');value=zeros(n,5);
for k=1:numel(factor.buckets)
    bucket=factor.buckets{k};count=numel(bucket.patchIDs);if count==0,continue;end
    b=rhs(bucket.members(:),:);pages=reshape(permute(reshape(b,k,count,5),[3 1 2]),5*k,1,count);solved=pagemtimes(bucket.inverse,pages);
    value(bucket.members(:),:)=reshape(permute(reshape(solved,5,k,count),[2 3 1]),k*count,5);
end
assert(all(isfinite(value),'all'),'CoupledPatches:Nonfinite','五量局部作用非有限，拒绝使用。');
end
