function [CpMax,state] = calcCpMax(freestream)
%CALCCPMAX 计算量热完全气体正激波后滞止压力及最大压力系数。
% 输入：freestream.Mach、gamma、pressure（Pa）。要求 M>=1、gamma>1。
% 输出：CpMax 为标量；state 为正激波后状态比和参考滞止压力。
% 用途：给修正牛顿模型提供压力系数上限，不代表每个面元经历正激波。
% 假设：定比热、无化学反应；高温真实气体效应不在本模块中处理。

%% 1. 检查来流物理输入
for name = {'Mach','gamma','pressure'}
    key = name{1};
    if ~isfield(freestream,key)
        error('RocketAeroHeat:CpMax:Input','缺少 freestream.%s。',key);
    end
    validateattributes(freestream.(key),{'numeric'},{'real','finite','scalar','positive'});
end
M = freestream.Mach; gamma = freestream.gamma; pinf = freestream.pressure;
if M < 1 || gamma <= 1
    error('RocketAeroHeat:CpMax:Domain','正激波参考模型要求 Mach >= 1 且 gamma > 1。');
end

%% 2. 正激波跃迁：压力比、密度比和下游马赫数
M2 = M*M;
pressureRatio = (2*gamma*M2-(gamma-1))/(gamma+1);
densityRatio = (gamma+1)*M2/((gamma-1)*M2+2);
downstreamMachSquared = (1+(gamma-1)*M2/2)/(gamma*M2-(gamma-1)/2);

%% 3. 激波后等熵减速至静止，得到参考滞止压力 p02
logP02Ratio = log(pressureRatio)+gamma/(gamma-1)* ...
    log1p((gamma-1)*downstreamMachSquared/2);
p02Ratio = exp(logP02Ratio);
qinf = 0.5*gamma*pinf*M2;
CpMax = expm1(logP02Ratio)/(0.5*gamma*M2);
state = struct('CpMax',CpMax,'stagnationPressure',pinf*p02Ratio, ...
    'stagnationPressureRatio',p02Ratio,'postShockPressureRatio',pressureRatio, ...
    'postShockDensityRatio',densityRatio, ...
    'postShockTemperatureRatio',pressureRatio/densityRatio, ...
    'postShockMach',sqrt(downstreamMachSquared),'dynamicPressure',qinf);
if any(~isfinite([CpMax,state.stagnationPressure,qinf])) || CpMax <= 0
    error('RocketAeroHeat:CpMax:NumericalRange','输入超出有限精度计算范围。');
end
end
