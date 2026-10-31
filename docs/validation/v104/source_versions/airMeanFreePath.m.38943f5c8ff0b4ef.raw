function lambda = airMeanFreePath(pressure,temperature,dynamicViscosity,gasConstant)
%AIRMEANFREEPATH 根据空气黏度估算分子平均自由程。
% 模块作用：为连续介质适用性筛选提供长度尺度，而不是求解稀薄流。
% 采用理想气体动理论常用黏度关系
% lambda=(mu/p)*sqrt(pi*R*T/2)。结果依赖输入空气模型和黏度关联式；
% 不能代替激波层内的非平衡分子碰撞或 DSMC 计算。

%% 1. 检查气体常数，并允许标量与逐面元数组组合
validateattributes(gasConstant,{'numeric'}, ...
    {'real','finite','scalar','positive'});
pressure=double(pressure);
temperature=double(temperature);
dynamicViscosity=double(dynamicViscosity);
try
    candidate=dynamicViscosity./pressure .* ...
        sqrt(pi*gasConstant*temperature/2);
catch
    error('RocketAeroHeat:MeanFreePath:Size', ...
        '压力、温度与动力黏度的数组尺寸不兼容。');
end

%% 2. 非物理或缺失状态保持 NaN，不将其伪装成连续介质
valid=isfinite(pressure) & pressure>0 & ...
    isfinite(temperature) & temperature>0 & ...
    isfinite(dynamicViscosity) & dynamicViscosity>0;
lambda=nan(size(candidate));
lambda(valid)=candidate(valid);
end
