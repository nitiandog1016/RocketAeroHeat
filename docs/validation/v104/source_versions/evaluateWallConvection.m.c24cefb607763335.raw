function [h,dh]=evaluateWallConvection(hBase,Taw,Twall,law)
%EVALUATEWALLCONVECTION 壁温变化时更新工程换热系数及其导数。
% 外流状态、恢复温度和转捩类别保持当前飞行工况；仅以 Eckert 参考
% 温度下的运动黏度修正 h。层流指数 1/2、湍流指数 1/5；若提供
% 干空气焓模型，还按当前焓差等效比热更新 h 及 Newton 导数。
% 这是局部参考性质反馈，不会重算三维边界层或外部激波。

hBase=hBase(:); Taw=Taw(:); Twall=Twall(:);
if numel(hBase)~=numel(Twall) || numel(Taw)~=numel(Twall)
    error('RocketAeroHeat:WallConvection:Size','壁温、h 和恢复温度尺寸不一致。');
end
if nargin<4 || isempty(law)
    h=hBase; dh=zeros(size(hBase)); return;
end
Te=law.edgeTemperature(:); p=law.edgePressure(:);
exponent=law.exponent(:); nu0=law.referenceViscosity0(:);
if any([numel(Te),numel(p),numel(exponent),numel(nu0)]~=numel(Twall))
    error('RocketAeroHeat:WallConvection:Law','参考性质数据尺寸不一致。');
end
Tref=Te+0.5*(Twall-Te)+0.22*(Taw-Te);
if any(~isfinite(Tref) | Tref<=0 | ~isfinite(p) | p<=0)
    error('RocketAeroHeat:WallConvection:Temperature', ...
        '壁温使参考温度或外缘压力无效。');
end
mu=1.458e-6*Tref.^(3/2)./(Tref+110.4);
nu=mu.*law.gasConstant.*Tref./p;
h=hBase.*(nu./nu0).^exponent;
% ln(nu)=常数+2.5 ln(Tref)-ln(Tref+110.4)。
dh=h.*exponent.*0.5.*(2.5./Tref-1./(Tref+110.4));
if isfield(law,'thermalPropertyMode') && ...
        strcmp(law.thermalPropertyMode,'thermallyPerfect')
    [cpCurrent,cpDerivative]=airEnthalpySecantCp(Taw,Twall);
    cpInitial=law.referenceCp0(:);
    if numel(cpInitial)~=numel(Twall) || ...
            any(~isfinite(cpCurrent) | cpCurrent<=0 | ...
            ~isfinite(cpInitial) | cpInitial<=0)
        error('RocketAeroHeat:WallConvection:AirRange', ...
            '壁温或恢复温度超出干空气焓模型 80~2500 K 适用域。');
    end
    factor=cpCurrent./cpInitial;
    dh=dh.*factor+h.*cpDerivative./cpInitial;
    h=h.*factor;
end
end
