function Taw = airRecoveryTemperature(Te,U,recoveryFactor)
%AIRRECOVERYTEMPERATURE 按空气焓守恒求局部绝热恢复温度。
% 模块作用：求解 h(Taw)=h(Te)+r*U^2/2；使用二分法保证单调焓
% 方程的稳定性。超过 80~2500 K 的热物性适用域时返回 NaN。

%% 1. 组装待求解的焓；无效外流面元保持 NaN
Te=double(Te); U=double(U);
validateattributes(recoveryFactor,{'numeric'}, ...
    {'real','finite','scalar','>=',0,'<=',1});
if ~isequal(size(Te),size(U))
    error('RocketAeroHeat:AirRecovery:Size','外缘温度和速度尺寸不一致。');
end
Taw=nan(size(Te));
valid=isfinite(Te) & isfinite(U) & Te>=80 & Te<=2500 & U>=0;
if ~any(valid(:)), return; end
[~,hEdge]=airThermalProperties(Te(valid));
[~,hMaximum]=airThermalProperties(2500);
target=hEdge+recoveryFactor*U(valid).^2/2;
solvable=isfinite(target) & target<=hMaximum;
if ~any(solvable), return; end

%% 2. 单调二分反解温度；不在适用域外延拓
low=Te(valid); low=low(solvable);
high=2500*ones(size(low));
selectedTarget=target(solvable);
for iteration=1:40
    midpoint=(low+high)/2;
    [~,hMid]=airThermalProperties(midpoint);
    below=hMid<selectedTarget;
    low(below)=midpoint(below);
    high(~below)=midpoint(~below);
end
resolved=(low+high)/2;
ids=find(valid);
Taw(ids(solvable))=resolved;
end
