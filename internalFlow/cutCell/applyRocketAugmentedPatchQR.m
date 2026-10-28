function value=applyRocketAugmentedPatchQR(factor,rhs,flag)
%APPLYROCKETAUGMENTEDPATCHQR 按原CV排列批量作用固定QR右映射及其精确转置。
% 本模块只实施代数映射；不按右端选区，不重新拟合，不修改流场。
if nargin<3,flag='notransp';end
assert(strcmp(factor.kind,'fixed-complete-augmented-image-patch-QR-right-map')&& ...
    isreal(rhs)&&isequal(size(rhs),[factor.cellCount 5])&&all(isfinite(rhs),'all'), ...
    'AugmentedQR:Input','QR映射类型或五量输入非法。');
assert(ismember(flag,{'notransp','transp'}),'AugmentedQR:Flag','未知QR映射方向。');
value=zeros(size(rhs));
for k=1:numel(factor.buckets)
    b=factor.buckets{k};count=numel(b.patchIDs);if count==0,continue;end
    collected=rhs(b.members(:),:);pages=reshape(permute(reshape(collected,k,count,5),[3 1 2]),5*k,1,count);
    mapping=b.map;if strcmp(flag,'transp'),mapping=permute(mapping,[2 1 3]);end
    solved=pagemtimes(mapping,pages);
    value(b.members(:),:)=reshape(permute(reshape(solved,5,k,count),[2 3 1]),k*count,5);
end
assert(all(isfinite(value),'all'),'AugmentedQR:Nonfinite','QR映射作用非有限。');
end
