function testWallThermalResponse()
% 无辐射隐式解析递推、辐射降温与 Biot 适用性。
addpath(fileparts(fileparts(mfilename('fullpath'))));
m.F=[1 2 3];
m.localGeometry.area=2;
m.heating.heatTransferCoefficient=100;
m.heating.recoveryTemperature=1000;
m.heating.eligibleFace=true;
cfg.flight.Twall=300;
mat=struct('density',1000,'specificHeat',1000,'thickness',0.01, ...
    'thermalConductivity',10,'emissivity',0, ...
    'radiationTemperature',300);
time=[0 10 20];
out=calcWallThermalResponse(m,time,mat,cfg);
out=calcAerothermalExposure(out);
C=mat.density*mat.specificHeat*mat.thickness;
expected=(C/10*300+100*1000)/(C/10+100);
assert(out.wallThermal.validFace);
assert(abs(out.wallThermal.temperature(1,2)-expected)<1e-10);
assert(abs(out.wallThermal.temperature(1,3)- ...
    (C/10*expected+100*1000)/(C/10+100))<1e-10);
assert(all(out.wallThermal.radiativeHeatFlux==0));
expectedDose=sum(diff(time).*out.wallThermal.convectiveHeatFlux(1,2:end));
assert(abs(out.aerothermalExposure.cumulativeNetHeat(1,end)- ...
    expectedDose)<1e-9*expectedDose);
assert(abs(out.aerothermalExposure.summary.totalNetEnergyResolved- ...
    2*expectedDose)<1e-9*expectedDose);
assert(out.aerothermalExposure.summary.resolvedAreaFraction==1);
recalculated=calcWallThermalResponse(out,time,mat,cfg);
assert(~isfield(recalculated,'aerothermalExposure'));
mat.emissivity=0.9;
radiating=calcWallThermalResponse(m,time,mat,cfg);
radiating=calcAerothermalExposure(radiating);
assert(radiating.wallThermal.temperature(1,end)<out.wallThermal.temperature(1,end));
assert(radiating.aerothermalExposure.cumulativeRadiativeHeat(1,end)>0);
assert(abs(radiating.aerothermalExposure.cumulativeNetHeat(1,end)- ...
    radiating.aerothermalExposure.cumulativeConvectiveHeat(1,end)+ ...
    radiating.aerothermalExposure.cumulativeRadiativeHeat(1,end))<1e-8);
assert(abs(radiating.aerothermalExposure.energyBalanceResidualPerArea)<1e-7);
mat.thermalConductivity=0.1;
thick=calcWallThermalResponse(m,time,mat,cfg);
assert(thick.wallThermal.validFace && thick.wallThermal.oneDimensionalFace);
assert(isfinite(thick.wallThermal.temperature(1,end)));
assert(thick.wallThermal.temperature(1,end)>out.wallThermal.temperature(1,end));
assert(isequal(size(thick.wallThermal.oneDimensional.finalDepthTemperature),[1 12]));
% 绝热背面且无辐射时，材料储热应等于隐式时刻对流输入积分。
mat.emissivity=0;
energy=calcWallThermalResponse(m,time,mat,cfg);
energy=calcAerothermalExposure(energy);
profile=energy.wallThermal.oneDimensional.finalDepthTemperature;
dx=mat.thickness/(size(profile,2)-1);
stored=mat.density*mat.specificHeat*dx* ...
    (sum(profile(2:end-1))+0.5*(profile(1)+profile(end)) ...
    -300*(size(profile,2)-1));
inputEnergy=sum(diff(time).*energy.wallThermal.convectiveHeatFlux(1,2:end));
assert(abs(stored-inputEnergy)<1e-8*max(1,inputEnergy));
assert(abs(energy.aerothermalExposure.energyBalanceResidualPerArea) ...
    <1e-8*max(1,inputEnergy));
fprintf('Wall thermal response passed: Tw(10 s)=%.3f K, Bi=%.3f.\n', ...
    expected,out.wallThermal.Biot);
end
