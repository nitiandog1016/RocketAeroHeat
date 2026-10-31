function [flux,audit]=shockFittedSlipWallFlux3D(Q,normal,wallNormalSpeed,gas)
%SHOCKFITTEDSLIPWALLFLUX3D 从内部状态与移动滑移壁面的镜像Riemann求壁压。
% 模块一：镜像相对法向速度，压力由两侧波系相容求出，不把输入壁压当求解。
% 模块二：壁面质量通量为0、动量为p*n、能量为p*D；零值来自壁面恒等式。
% 模块三：仅无黏滑移条件，不计算黏性热流、壁温或真实边界层外缘。
state=shockFittedState3D(Q,gas,'conservative');
assert(isequal(size(normal),size(state.velocity))&&isequal(size(wallNormalSpeed),size(state.density)), ...
    'ShockFit:WallGeometry','壁面状态、法向和移动速度尺寸不符。');
relative=sum(state.velocity.*normal,2)-wallNormalSpeed;
mirrorVelocity=state.velocity-2*relative.*normal;
mirror=shockFittedState3D([state.density mirrorVelocity state.temperature],gas,'primitive');
[~,solver]=shockFittedALEFlux3D(Q,mirror.Q,normal,wallNormalSpeed,gas);
p=solver.riemann.pressureStarPa;
assert(all(p>0&isfinite(p)),'ShockFit:WallPressure','壁面Riemann没有得到正压力。');
flux=[zeros(size(p)),p.*normal,p.*wallNormalSpeed];
audit=struct('pressurePa',p,'interiorPressurePa',state.pressure,'inputRelativeNormalVelocityMS',relative, ...
    'relativeStarNormalVelocityMS',solver.riemann.normalVelocityStarMs, ...
    'method','镜像相对法向速度/精确常gamma壁面Riemann；能量通量为求得壁压乘壁面速度', ...
    'pressureGivenAsBoundaryCondition',false,'heatFluxSolved',false,'wallTemperatureSolved',false);
end
