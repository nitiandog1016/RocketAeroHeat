function [net,derivative,convective,radiative]=evaluateBackWallHeatFlux(temperature,boundary)
%EVALUATEBACKWALLHEATFLUX 内侧散热 [W/m²] 及对背面温度的导数。
% 输入 boundary 由 wallBackBoundary 归一化；正值表示材料向内环境散热，
% 负值表示内环境加热材料。NaN 保留为未解面元，不补零。
validateattributes(temperature,{'numeric'},{'real'});
if any(isinf(temperature),'all') || any(temperature(isfinite(temperature))<=0)
    error('RocketAeroHeat:BackBoundary:Temperature','背面温度应为正有限值或 NaN。');
end
convective=zeros(size(temperature)); radiative=convective; derivative=convective;
if boundary.active
    convective=boundary.heatTransferCoefficient*(temperature-boundary.fluidTemperature);
    factor=boundary.emissivity*5.670374419e-8;
    % 分解四次方差，降低接近平衡温度时两个大数相减的精度损失。
    radiative=factor*(temperature-boundary.radiationTemperature).* ...
        (temperature+boundary.radiationTemperature).* ...
        (temperature.^2+boundary.radiationTemperature^2);
    derivative=boundary.heatTransferCoefficient+4*factor*temperature.^3;
end
net=convective+radiative;
forInvalid=isnan(temperature);
net(forInvalid)=NaN; derivative(forInvalid)=NaN;
convective(forInvalid)=NaN; radiative(forInvalid)=NaN;
if any(~isfinite(net(~forInvalid))) || any(~isfinite(derivative(~forInvalid)))
    error('RocketAeroHeat:BackBoundary:Numerical','内侧热流溢出，请检查温度和换热参数量级。');
end
end
