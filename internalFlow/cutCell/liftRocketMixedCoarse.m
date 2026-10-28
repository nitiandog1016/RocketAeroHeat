function direction=liftRocketMixedCoarse(coarse,coefficients)
%LIFTROCKETMIXEDCOARSE 按逐列身份提升保守体积方向及筛选的原始变量耦合方向。
% 模块一：保守部分直接延拓原缩放五量；附加部分按同Q的rho/u/T切线提升。
% 模块二：两个部分相加是固定线性映射，不改变物理网格或产生非线性候选。
nc=5*coarse.space.groupCount;indices=coarse.selectedPrimitiveColumns(:);d=nc+numel(indices);
assert(isvector(coefficients)&&numel(coefficients)==d&&all(isfinite(coefficients)),'MixedCoarse:Coefficients','混合粗系数非法。');
values=coefficients(:).*coarse.columnLiftScale;
direction=coarse.space.prolongation*reshape(values(1:nc),coarse.space.groupCount,5);
primitive=zeros(coarse.primitiveSpace.groupCount,5);primitive(indices)=values(nc+1:end);
direction=direction+liftRocketPrimitiveCoarse(coarse.baseQ,coarse.gas,coarse.scale,coarse.primitiveSpace,primitive);
end
