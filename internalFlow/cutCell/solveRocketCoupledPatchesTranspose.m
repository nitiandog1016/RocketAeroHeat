function value=solveRocketCoupledPatchesTranspose(factor,rhs)
%SOLVEROCKETCOUPLEDPATCHESTRANSPOSE 作用原全部五量局部块逆的转置，不重新求逆。
% 模块一：成员和原CV顺序沿用既有认证因子，包含其实际对角替代身份。
% 模块二：每页逆转置后按同一原CV顺序散回，仅作固定代数预条件伴随。
assert(strcmp(factor.kind,'fixed-five-conservation-strong-connection-local-block-inverse'),'CoupledPatches:Kind','需要原全部五量局部块因子。');
n=factor.cellCount;assert(isequal(size(rhs),[n 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'CoupledPatches:RHS','全部五量转置右端非法。');
value=zeros(n,5);
for k=1:numel(factor.buckets)
    b=factor.buckets{k};count=numel(b.patchIDs);if count==0,continue;end
    collected=rhs(b.members(:),:);pages=reshape(permute(reshape(collected,k,count,5),[3 1 2]),5*k,1,count);
    solved=pagemtimes(permute(b.inverse,[2 1 3]),pages);
    value(b.members(:),:)=reshape(permute(reshape(solved,5,k,count),[2 3 1]),k*count,5);
end
assert(all(isfinite(value),'all'),'CoupledPatches:Nonfinite','原五量局部逆转置非有限。');
end
