function direction=liftRocketLocalCoarse(coarse,coefficients)
%LIFTROCKETLOCALCOARSE 将保留混合列及原CV单自由度方向作固定线性提升。
% 模块一：继承混合列仍按原保守量/rho-u-T公式提升，未改变旧缓存的含义。
% 模块二：附加列只扰动指定原CV的一个缩放保守量，不删除CV或修改物理网格。
nr=numel(coarse.retainedMixedColumns);nl=numel(coarse.localCellIDs);n=size(coarse.baseQ,1);
assert(isnumeric(coefficients)&&isreal(coefficients)&&isvector(coefficients)&&numel(coefficients)==nr+nl&& ...
    all(isfinite(coefficients)),'LocalCoarse:Coefficients','局部粗系数非法。');
coefficients=coefficients(:);parent=zeros(numel(coarse.mixedLift.columnLiftScale),1);
parent(coarse.retainedMixedColumns)=coefficients(1:nr);direction=liftRocketMixedCoarse(coarse.mixedLift,parent);
indices=coarse.localCellIDs+n*(coarse.localEquationIDs-1);
assert(numel(unique(indices))==nl,'LocalCoarse:Lift','局部保守列身份不能重复。');
direction(indices)=direction(indices)+coefficients(nr+1:end).*coarse.localColumnLiftScale;
end
