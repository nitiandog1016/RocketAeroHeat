function coarse=restrictRocketPatchResidual(partition,fine)
%RESTRICTROCKETPATCHRESIDUAL 将全部原CV五量按局部组作等权代数平均。
% 模块一：R是固定代数限制，不声称真体积积分或修改物理守恒方程。
% 模块二：与组常数延拓P满足R*P=I，任何控制体都不能被遗漏。
n=partition.cellCount;ng=partition.patchCount;
assert(isequal(size(fine),[n 5])&&isreal(fine)&&all(isfinite(fine),'all'),'PatchTransfer:Input','需要全部CV的有限五量。');
coarse=zeros(ng,5);group=double(partition.groupID);sizes=double(partition.sizes(:));
for equation=1:5,coarse(:,equation)=accumarray(group,fine(:,equation),[ng 1])./sizes;end
assert(all(isfinite(coarse),'all'),'PatchTransfer:Nonfinite','限制作用非有限。');
end
