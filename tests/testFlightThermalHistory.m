function testFlightThermalHistory()
% 变工况壁温：常工况须回归旧求解器，变工况须满足离散能量平衡。
addpath(fileparts(fileparts(mfilename('fullpath'))));
cfg.flight.Twall=300;
mesh.F=[1 2 3];
mesh.localGeometry.area=2;
t=[0 1 2];
material=struct('density',1000,'specificHeat',1000, ...
    'thickness',0.01,'thermalConductivity',100, ...
    'emissivity',0,'radiationTemperature',300);
profile=struct('time',t,'initialWallTemperature',300, ...
    'heatTransferCoefficient',100*ones(1,3), ...
    'recoveryTemperature',500*ones(1,3), ...
    'eligibleFace',true(1,3));
history=calcWallThermalHistory(mesh,profile,material,cfg);
single=mesh;
single.heating=struct('heatTransferCoefficient',100, ...
    'recoveryTemperature',500,'eligibleFace',true, ...
    'wallTemperature',300);
single=calcWallThermalResponse(single,t,material,cfg);
assert(max(abs(history.wallThermal.temperature- ...
    single.wallThermal.temperature),[],'all')<1e-9);
profile.heatTransferCoefficient=[100 200 300];
variable=calcWallThermalHistory(mesh,profile,material,cfg);
variable=calcAerothermalExposure(variable);
T=variable.wallThermal.temperature;
assert(T(end)>history.wallThermal.temperature(end));
C=material.density*material.specificHeat*material.thickness;
residual=C*(T(3)-T(2))-(t(3)-t(2))* ...
    variable.wallThermal.netHeatFlux(1,3);
assert(abs(residual)<1e-8*max(1,abs(C*(T(3)-T(2)))));
assert(abs(variable.aerothermalExposure.energyBalanceResidualPerArea) ...
    <1e-8*max(1,variable.aerothermalExposure.cumulativeNetHeat(1,end)));
recalculated=calcWallThermalHistory(variable,profile,material,cfg);
assert(~isfield(recalculated,'aerothermalExposure'));
profile.eligibleFace(2)=false;
unknown=calcWallThermalHistory(mesh,profile,material,cfg);
assert(~unknown.wallThermal.validFace && ...
    all(isnan(unknown.wallThermal.temperature)));

% 厚壁也与固定工况一维有限体积解对照。
material.thickness=0.05; material.thermalConductivity=0.2;
profile.eligibleFace(:)=true;
profile.heatTransferCoefficient(:)=100;
thick=calcWallThermalHistory(mesh,profile,material,cfg);
thick=calcAerothermalExposure(thick);
single=calcWallThermalResponse(single,t,material,cfg);
assert(thick.wallThermal.oneDimensionalFace);
assert(max(abs(thick.wallThermal.temperature- ...
    single.wallThermal.temperature),[],'all')<1e-8);
assert(abs(thick.aerothermalExposure.energyBalanceResidualPerArea) ...
    <1e-8*max(1,thick.aerothermalExposure.cumulativeNetHeat(1,end)));
% 对常值 h、Taw 的解析一阶热容解，较密时间站应降低隐式步进误差。
thinMaterial=material;
thinMaterial.thickness=0.01;
thinMaterial.thermalConductivity=100;
coarseProfile=struct('time',[0 10], ...
    'initialWallTemperature',300, ...
    'heatTransferCoefficient',[100 100], ...
    'recoveryTemperature',[500 500], ...
    'eligibleFace',[true true]);
denseProfile=coarseProfile;
denseProfile.time=0:10;
denseProfile.heatTransferCoefficient=100*ones(1,11);
denseProfile.recoveryTemperature=500*ones(1,11);
denseProfile.eligibleFace=true(1,11);
coarseWall=calcWallThermalHistory(mesh,coarseProfile,thinMaterial,cfg);
denseWall=calcWallThermalHistory(mesh,denseProfile,thinMaterial,cfg);
capacity=thinMaterial.density*thinMaterial.specificHeat* ...
    thinMaterial.thickness;
exact=500-(500-300)*exp(-100*10/capacity);
coarseError=abs(coarseWall.wallThermal.temperature(end)-exact);
denseError=abs(denseWall.wallThermal.temperature(end)-exact);
assert(denseError<coarseError/2);
fprintf('Flight wall-history constant/variable and thick-wall checks passed.\n');
end
