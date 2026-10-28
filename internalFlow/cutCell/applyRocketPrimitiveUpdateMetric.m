function out=applyRocketPrimitiveUpdateMetric(metric,value,flag)
%APPLYROCKETPRIMITIVEUPDATEMETRIC 应用状态更新度量或其精确转置，不构造全局稠密矩阵。
% 模块一：正向输出[相对密度，三向速度/当地声速，相对温度]的一阶更新。
% 模块二：转置完整保留密度—动量—能量耦合，供增广最小二乘使用。
if nargin<3,flag='notransp';end
assert(strcmp(metric.kind,'current-state-linear-primitive-relative-update-metric')&& ...
    isequal(size(value),[metric.cellCount 5])&&all(isfinite(value),'all'), ...
    'PrimitiveMetric:Input','度量或五量方向非法。');
out=metric.diagonal.*value;
if strcmp(flag,'notransp')
    out(:,2:4)=out(:,2:4)+metric.densityToVelocity.*value(:,1);
    out(:,5)=out(:,5)+metric.densityToTemperature.*value(:,1)+sum(metric.momentumToTemperature.*value(:,2:4),2);
elseif strcmp(flag,'transp')
    out(:,1)=out(:,1)+sum(metric.densityToVelocity.*value(:,2:4),2)+metric.densityToTemperature.*value(:,5);
    out(:,2:4)=out(:,2:4)+metric.momentumToTemperature.*value(:,5);
else,error('PrimitiveMetric:TransposeFlag','未知状态度量作用方向。');end
end
