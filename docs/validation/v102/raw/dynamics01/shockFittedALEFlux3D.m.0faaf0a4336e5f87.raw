function [flux,audit]=shockFittedALEFlux3D(QL,QR,normal,gridNormalSpeed,gas)
%SHOCKFITTEDALEFLUX3D 用真实移动界面的精确常gamma Riemann求全部五量ALE通量。
% 模块一：以D*n作Galilean变换，在界面静止系精确采样，不将F-DQ套在错误采样态上。
% 模块二：质量/动量/能量通量共同转回实验室系，保留压力的移动壁面功。
% 模块三：当前精确Riemann只支持显式常gamma；NASA7明确拒绝，不静默回退。
assert(strcmp(gas.thermalPropertyMode,'caloricallyPerfect'),'ShockFit:FluxGas', ...
    '当前ALE精确Riemann只支持常gamma；不能用于NASA7或静默切换物性。');
N=size(QL,1);assert(isequal(size(QL),[N 5])&&isequal(size(QR),[N 5])&& ...
    isequal(size(normal),[N 3])&&isequal(size(gridNormalSpeed),[N 1])&&isreal(gridNormalSpeed)&& ...
    all(isfinite(gridNormalSpeed))&&all(abs(sum(normal.^2,2)-1)<1e-12),'ShockFit:FluxGeometry','五量、法向和网格速度尺寸不符。');
w=gridNormalSpeed.*normal;
left=boost(QL,w);right=boost(QR,w);
constant=struct('gamma',gas.gamma,'R',gas.R,'temperatureRange',gas.temperatureRangeK);
[relative,riemann]=cutCellExactEulerFlux(left,right,normal,constant);
flux=[relative(:,1),relative(:,2:4)+w.*relative(:,1), ...
    relative(:,5)+sum(w.*relative(:,2:4),2)+.5*sum(w.^2,2).*relative(:,1)];
audit=struct('method','移动法向系精确常gamma Godunov；Galilean变换全部五量ALE通量', ...
    'riemann',riemann,'gridNormalSpeedMS',gridNormalSpeed,'physicalQualification',false);
end

function shifted=boost(Q,w)
shifted=[Q(:,1),Q(:,2:4)-Q(:,1).*w, ...
    Q(:,5)-sum(Q(:,2:4).*w,2)+.5*Q(:,1).*sum(w.^2,2)];
end
